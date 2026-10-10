// SPDX-License-Identifier: GPL-3.0-or-later
// Experimental Cartesian spherical mechanical verification solver, M2A-01.
#include "fvCFD.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "MechanicalBackend.H"
#include "SphereSurface.H"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <chrono>

using namespace Foam;
static vector asFoamVector(const double* p) {return vector(p[0],p[1],p[2]);}
static vector analytical(const vector& x,const vector& centre,scalar radius,
                         const vector& speed,const vector& omega,bool translation)
{
    const vector r=x-centre;
    const scalar distance=mag(r);
    if(distance<=radius) return translation ? vector::zero : (omega^r);
    if(!translation) return pow3(radius/distance)*(omega^r);
    const vector n=r/distance;
    const scalar a=radius/distance;
    return (1-0.75*a-0.25*pow3(a))*speed +(-0.75*a+0.75*pow3(a))*(speed&n)*n;
}

int main(int argc,char* argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    int initialized=0;MPI_Initialized(&initialized);
    bool ownMPI=!initialized;
    if(ownMPI) MPI_Init(&argc,&argv);
    int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try
    {
        if(ranks!=Pstream::nProcs()) throw std::runtime_error("MPI/OpenFOAM rank mismatch");
        IOdictionary config(IOobject("mechanicalProperties",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
        const word mode(config.lookup("mode"));
        if(mode!="fixed" && mode!="rotate" && mode!="translate" && mode!="free")
            throw std::runtime_error("Unknown mechanical mode");
        const scalar rho=readScalar(config.lookup("rho")), viscosity=readScalar(config.lookup("nu"));
        const scalar demDt=readScalar(config.lookup("demDeltaT"));
        const scalar dt=runTime.deltaTValue(), penalty=readScalar(config.lookup("penalty"));
        const label q=readLabel(config.lookup("quadrature")), correctors=readLabel(config.lookup("correctors"));
        const label maxCorrectors=config.getOrDefault<label>("maxCorrectors",correctors);
        const scalar momentumTolerance=config.getOrDefault<scalar>("momentumImpulseTolerance",-1);
        const word constraintScheme=config.getOrDefault<word>("constraintScheme","volumePenalty");
        const bool surfaceExtension=constraintScheme=="surfaceExtension";
        if(!surfaceExtension&&constraintScheme!="volumePenalty") throw std::runtime_error("Unknown constraintScheme");
        if(surfaceExtension&&mode!="fixed"&&mode!="rotate")
            throw std::runtime_error("Surface extension candidate is stationary single-sphere only");
        const label maxSurfaceCorrectors=config.getOrDefault<label>("maxSurfaceCorrectors",32);
        const scalar targetTolerance=config.getOrDefault<scalar>("surfaceTargetTolerance",1e-7);
        const word reconstruction=config.getOrDefault<word>("surfaceReconstruction","linear");
        const scalar surfaceRelaxation=config.getOrDefault<scalar>("surfaceRelaxation",1);
        const scalar continuityTolerance=config.getOrDefault<scalar>("continuityTolerance",-1);
        if(reconstruction!="linear"&&reconstruction!="quadratic") throw std::runtime_error("Unknown surface reconstruction");
        if(surfaceRelaxation<=0||surfaceRelaxation>1) throw std::runtime_error("Invalid surface relaxation");
        if(surfaceExtension&&(maxSurfaceCorrectors<1||targetTolerance<=0))
            throw std::runtime_error("Invalid surface target controls");
        if(rank==0) Info<<"M2A_CONSTRAINT scheme="<<constraintScheme<<endl;
        if(surfaceExtension&&rank==0) Info<<"M2A_RECONSTRUCTION order="<<reconstruction<<endl;
        const int substeps=static_cast<int>(std::llround(dt/demDt));
        if(rho<=0||viscosity<=0||penalty<=0||q<2||q>8||correctors<1||maxCorrectors<correctors
            || substeps<1||std::abs(substeps*demDt-dt)>dt*1e-10)
            throw std::runtime_error("Invalid mechanical controls or nonintegral DEM substeps");
        const dimensionedScalar nu("nu",dimViscosity,viscosity);
        volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
        volScalarField p(IOobject("p",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
        if(U.dimensions()!=dimVelocity || p.dimensions()!=sqr(dimVelocity))
            throw std::runtime_error("M2A uses velocity and kinematic pressure, not Pa");
        surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),fvc::flux(U));
        mesh.setFluxRequired(p.name());
        volScalarField solid(IOobject("solidFraction",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar("zero",dimless,0));
        volScalarField lambda(IOobject("constraintRate",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("zero",dimless/dimTime,0));
        volVectorField rigid(IOobject("rigidVelocity",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedVector("zero",dimVelocity,vector::zero));
        lpbfM2::MechanicalBackend backend("input.dem");
        auto particles=backend.state();
        // First integrated package isolates one sphere; API does not require consecutive IDs.
        if(particles.size()!=1) throw std::runtime_error("M2A-01 requires one noncontact sphere");
        const label np=particles.size();
        List<vector> previousP(np,vector::zero),previousL(np,vector::zero),previousCentre(np,vector::zero);
        labelList ids(np);
        for(label a=0;a<np;++a) {ids[a]=particles[a].id;previousCentre[a]=asFoamVector(particles[a].x);}
        const scalar h=std::cbrt(gMin(mesh.V().field()));
        if(gMax(mesh.V().field())>pow3(h)*(1+1e-10)) throw std::runtime_error("Uniform cubic mesh required");
        forAll(mesh.cells(),cell)
        {
            const boundBox box(mesh.points(),mesh.cellPoints()[cell],false);
            const vector span=box.span();
            if(mag(span-vector(h,h,h))>h*1e-8) throw std::runtime_error("Axis-aligned cubic cells required");
        }
        std::unique_ptr<lpbfM2::SphereSurface> surface;
        if(surfaceExtension) surface.reset(new lpbfM2::SphereSurface(mesh,h,reconstruction=="quadratic"));
        auto mapGeometry=[&]()
        {
            if(surface) surface->gather(U,p);
            solid=dimensionedScalar("zero",dimless,0);
            rigid=dimensionedVector("zero",dimVelocity,vector::zero);
            forAll(mesh.C(),cell)
            {
                const vector centre=mesh.C()[cell], rp=asFoamVector(particles[0].x);
                if(mag(centre-rp)>particles[0].radius+0.866026*h) continue;
                label hits=0;vector average=vector::zero;
                for(label i=0;i<q;++i) for(label j=0;j<q;++j) for(label k=0;k<q;++k)
                {
                    const vector x=centre+h*vector((i+0.5)/q-0.5,(j+0.5)/q-0.5,(k+0.5)/q-0.5);
                    if(mag(x-rp)<=particles[0].radius)
                    {++hits;average+=asFoamVector(particles[0].v)+(asFoamVector(particles[0].omega)^(x-rp));}
                }
                solid[cell]=scalar(hits)/(q*q*q);
                if(hits) rigid[cell]=average/hits;
            }
            solid.correctBoundaryConditions();rigid.correctBoundaryConditions();
            lambda=dimensionedScalar("rate",dimless/dimTime,penalty/dt)*solid;
            if(surface)
                forAll(mesh.C(),cell)
                {
                    const vector rp=asFoamVector(particles[0].x);
                    const bool inside=mag(mesh.C()[cell]-rp)<=particles[0].radius;
                    lambda[cell]=inside ? penalty/dt:0;
                    if(inside) rigid[cell]=surface->target(mesh.C()[cell],rp,particles[0].radius,
                        asFoamVector(particles[0].v),asFoamVector(particles[0].omega));
                }
            rigid.correctBoundaryConditions();
        };
        auto interior=[&](vector& momentum,vector& angular)
        {
            momentum=vector::zero;angular=vector::zero;
            forAll(mesh.C(),cell)
            {
                const vector m=rho*solid[cell]*mesh.V()[cell]*U[cell];
                momentum+=m;angular+=(mesh.C()[cell]-asFoamVector(particles[0].x))^m;
            }
            reduce(momentum,sumOp<vector>());reduce(angular,sumOp<vector>());
        };
        if(runTime.value()==0)
        {
            forAll(mesh.C(),cell)
            {
                const vector r=mesh.C()[cell]-asFoamVector(particles[0].x);
                const scalar distance=mag(r), radius=particles[0].radius;
                if(mode=="fixed")
                {
                    U[cell]=analytical(mesh.C()[cell],asFoamVector(particles[0].x),radius,vector(0.01,0,0),vector::zero,true);
                    p[cell]=distance>radius ? -1.5*viscosity*radius*(vector(0.01,0,0)&r)/pow3(distance) : 0;
                }
                else if(mode=="rotate") U[cell]=analytical(mesh.C()[cell],asFoamVector(particles[0].x),radius,vector::zero,asFoamVector(particles[0].omega),false);
                else
                {
                    U[cell]=asFoamVector(particles[0].v)-analytical(mesh.C()[cell],asFoamVector(particles[0].x),radius,asFoamVector(particles[0].v),vector::zero,true);
                }
            }
        }
        // Exact Dirichlet velocity for stationary Stokes references removes
        // the usual finite-box velocity-boundary error; closed walls otherwise.
        forAll(U.boundaryField(),patch)
            if(!mesh.boundary()[patch].coupled())
                forAll(U.boundaryField()[patch],face)
                    U.boundaryFieldRef()[patch][face]=(mode=="fixed"||mode=="rotate")
                        ? analytical(mesh.Cf().boundaryField()[patch][face],asFoamVector(particles[0].x),particles[0].radius,
                                     vector(0.01,0,0),asFoamVector(particles[0].omega),mode=="fixed") : vector::zero;
        U.correctBoundaryConditions();
        if(runTime.value()==0) phi=fvc::flux(U);
        mapGeometry();
        IOdictionary restore(IOobject("couplingState",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE,false));
        if(runTime.value()>0)
        {
            if(restore.getOrDefault<word>("constraintScheme","volumePenalty")!=constraintScheme)
                throw std::runtime_error("Restart constraint scheme differs from checkpoint");
            if(surfaceExtension&&restore.getOrDefault<word>("surfaceReconstruction","linear")!=reconstruction)
                throw std::runtime_error("Restart surface reconstruction differs from checkpoint");
            if(!restore.found("ids")) throw std::runtime_error("Missing common coupling checkpoint");
            labelList saved(restore.lookup("ids"));
            if(saved!=ids) throw std::runtime_error("Restart ID mismatch");
            restore.lookup("interiorMomentum")>>previousP;
            restore.lookup("interiorAngular")>>previousL;
            restore.lookup("mappedCentres")>>previousCentre;
            if(previousP.size()!=np||previousL.size()!=np||previousCentre.size()!=np)
                throw std::runtime_error("Restart state length mismatch");
        }
        else interior(previousP[0],previousL[0]);
        if(std::abs(backend.time()-runTime.value())>1e-12) throw std::runtime_error("Restart clocks differ");
        std::ofstream history;
        if(rank==0)
        {
            history.open("mechanical-history.csv",runTime.value()>0 ? std::ios::app : std::ios::out);
            if(!history) throw std::runtime_error("Cannot write mechanical history");
            history<<std::setprecision(17);
            if(runTime.value()==0)
                history<<"time,id,owner,x,y,z,vx,vy,vz,wx,wy,wz,fx,fy,fz,tx,ty,tz,volume_error,slip_rms,div_max,momentum_residual,angular_impulse_residual,clock_error,covered_ranks,force_ratio,torque_ratio,fluid_px,fluid_py,fluid_pz,old_fluid_px,old_fluid_py,old_fluid_pz,boundary_fx,boundary_fy,boundary_fz,constraint_fx,constraint_fy,constraint_fz,inertia_fx,inertia_fy,inertia_fz,constraint_tx,constraint_ty,constraint_tz,inertia_tx,inertia_ty,inertia_tz,pressure_correctors,momentum_equation_impulse_L1,support_volume_ratio,surface_correctors,surface_target_defect,pressure_correctors_total,bulk_target_slip_rms,stress_pressure_fx,stress_pressure_fy,stress_pressure_fz,stress_viscous_fx,stress_viscous_fy,stress_viscous_fz,stress_pressure_tx,stress_pressure_ty,stress_pressure_tz,stress_viscous_tx,stress_viscous_ty,stress_viscous_tz,window_wall_seconds\n";
        }
        const scalar refPressure=p[0];
        while(runTime.run())
        {
            const auto wallStart=std::chrono::steady_clock::now();
            if(mag(runTime.deltaTValue()-dt)>dt*1e-12) throw std::runtime_error("Adaptive dt unsupported in first package");
            const auto before=particles;
            const vector oldPhysical=gSum(U.primitiveField()*mesh.V().field())*rho-previousP[0];
            ++runTime;
            mapGeometry();
            surfaceScalarField advecting(IOobject("advecting",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),phi);
            label usedCorrectors=0,surfaceCorrectors=0,totalCorrectors=0;
            scalar equationImpulse=GREAT,targetDefect=0,continuityDefect=GREAT;
            for(label outer=0;outer<(surface ? maxSurfaceCorrectors:1);++outer)
            {
                fvVectorMatrix equation(fvm::ddt(U)+fvm::div(advecting,U)-fvm::laplacian(nu,U)+fvm::Sp(lambda,U)==lambda*rigid);
                // One momentum predictor per outer target update. Each inner PISO
                // loop corrects H without re-solving the predictor. All outer solves
                // retain the SAME old-time field and frozen advecting flux/window.
                solve(equation==-fvc::grad(p));
                volScalarField rAU("rAU",1.0/equation.A());
                for(label c=0;c<maxCorrectors;++c)
                {
                    volVectorField HbyA(constrainHbyA(rAU*equation.H(),U,p));
                    surfaceScalarField predicted("predicted",fvc::flux(HbyA));
                    // Standard transient collocated flux consistency, using the
                    // same inverse diagonal that includes the solid constraint.
                    predicted+=fvc::interpolate(rAU)*fvc::ddtCorr(U,phi);
                    if(p.needReference()) adjustPhi(predicted,U,p);
                    constrainPressure(p,U,predicted,rAU);
                    fvScalarMatrix pressure(fvm::laplacian(rAU,p)==fvc::div(predicted));
                    if(p.needReference()) pressure.setReference(rank==0 ? 0:-1,refPressure);
                    pressure.solve();
                    phi=predicted-pressure.flux();
                    U=HbyA-rAU*fvc::grad(p);U.correctBoundaryConditions();
                    const scalar maxU=gMax(mag(U)().primitiveField());
                    const scalar maxP=gMax(mag(p)().primitiveField());
                    continuityDefect=gMax(mag(fvc::div(phi))().primitiveField());
                    // residual() is source - A*U, integrated over each cell.
                    // Include the pressure source and sum magnitudes, so opposite
                    // local defects cannot cancel in a global momentum ledger.
                    vectorField defect(equation.residual());
                    // In this v2512 residual() path, addBoundarySource includes
                    // coupled neighbours and lduMatrix::residual includes their
                    // interface contribution again. Remove one explicit copy;
                    // the momentum matrix and pressure solve remain unchanged.
                    forAll(U.boundaryField(),patch)
                        if(U.boundaryField()[patch].coupled())
                        {
                            const vectorField neighbour(U.boundaryField()[patch].patchNeighbourField());
                            const labelUList& cells=mesh.boundary()[patch].faceCells();
                            forAll(cells,face)
                                defect[cells[face]]-=cmptMultiply(equation.boundaryCoeffs()[patch][face],neighbour[face]);
                        }
                    defect-=mesh.V().field()*fvc::grad(p)().primitiveField();
                    equationImpulse=rho*dt*gSum(mag(defect));
                    usedCorrectors=c+1;
                    ++totalCorrectors;
                    if(rank==0) Info<<"M2A_CORRECTION time="<<runTime.value()
                        <<" iteration="<<c+1<<" maxU="<<maxU<<" maxP="<<maxP
                        <<" equationImpulseL1="<<equationImpulse<<" continuityMax="<<continuityDefect<<endl;
                    if(!std::isfinite(maxU)||!std::isfinite(maxP)||!std::isfinite(equationImpulse)||!std::isfinite(continuityDefect))
                        throw std::runtime_error("Nonfinite CFD state before DEM feedback");
                    if(usedCorrectors>=correctors && (momentumTolerance<0||equationImpulse<=momentumTolerance)
                        && (continuityTolerance<0||continuityDefect<=continuityTolerance)) break;
                }
                if(momentumTolerance>=0 && equationImpulse>momentumTolerance)
                    throw std::runtime_error("CFD momentum equation did not converge before DEM feedback");
                if(continuityTolerance>=0&&continuityDefect>continuityTolerance)
                    throw std::runtime_error("CFD continuity did not converge before DEM feedback");
                surfaceCorrectors=outer+1;
                if(!surface) break;
                surface->gather(U,p);
                vectorField updated(rigid.primitiveField());targetDefect=0;
                forAll(mesh.C(),cell) if(lambda[cell]>0)
                {
                    updated[cell]=surface->target(mesh.C()[cell],asFoamVector(before[0].x),before[0].radius,
                        asFoamVector(before[0].v),asFoamVector(before[0].omega));
                    targetDefect=Foam::max(targetDefect,mag(updated[cell]-rigid[cell]));
                }
                reduce(targetDefect,maxOp<scalar>());
                if(rank==0) Info<<"M2A_SURFACE time="<<runTime.value()<<" iteration="<<surfaceCorrectors
                    <<" targetDefect="<<targetDefect<<endl;
                if(!std::isfinite(targetDefect)) throw std::runtime_error("Nonfinite surface target");
                // Retain exactly the target used by the last solved momentum
                // equation for its force ledger. Do not silently substitute a
                // newly reconstructed source after declaring convergence.
                if(targetDefect<=targetTolerance) break;
                if(surfaceCorrectors==maxSurfaceCorrectors)
                    throw std::runtime_error("Surface target did not converge before DEM feedback");
                rigid.primitiveFieldRef()+=surfaceRelaxation*(updated-rigid.primitiveField());
                rigid.correctBoundaryConditions();
            }
            lpbfM2::SphereSurface::Diagnostics stress;
            if(surface) stress=surface->diagnose(asFoamVector(before[0].x),before[0].radius,
                asFoamVector(before[0].v),asFoamVector(before[0].omega),rho,viscosity);
            vector constraint=vector::zero, constraintTorque=vector::zero;
            scalar volume=0,slip=0,supportVolume=0;label covered=0;
            forAll(mesh.C(),cell)
            {
                const vector f=rho*lambda[cell]*(rigid[cell]-U[cell])*mesh.V()[cell];
                constraint+=f;constraintTorque+=(mesh.C()[cell]-asFoamVector(before[0].x))^f;
                volume+=solid[cell]*mesh.V()[cell];slip+=solid[cell]*mesh.V()[cell]*magSqr(U[cell]-rigid[cell]);
                if(lambda[cell]>0) supportVolume+=mesh.V()[cell];
            }
            covered=volume>0 ? 1:0;reduce(covered,sumOp<label>());
            reduce(constraint,sumOp<vector>());reduce(constraintTorque,sumOp<vector>());
            reduce(volume,sumOp<scalar>());reduce(slip,sumOp<scalar>());
            reduce(supportVolume,sumOp<scalar>());
            if(volume<=0) throw std::runtime_error("Sphere outside resolved mesh");
            vector newP,newL;interior(newP,newL);
            const vector inertia=(newP-previousP[0])/dt;
            const vector rotationalInertia=(newL-previousL[0]+((asFoamVector(before[0].x)-previousCentre[0])^previousP[0]))/dt;
            const vector force=-constraint+inertia, torque=-constraintTorque+rotationalInertia;
            const vector newPhysical=gSum(U.primitiveField()*mesh.V().field())*rho-newP;
            // Boundary momentum ledger uses exactly the selected Gauss-linear
            // convection and orthogonal Laplacian operators, not an extra force.
            surfaceVectorField faceU(fvc::interpolate(U)), diffusion(nu*fvc::snGrad(U)*mesh.magSf());
            surfaceScalarField faceP(fvc::interpolate(p));
            vector boundary=vector::zero;
            forAll(mesh.boundary(),patch)
                if(!mesh.boundary()[patch].coupled())
                    boundary+=rho*(sum(diffusion.boundaryField()[patch])
                        -sum(faceP.boundaryField()[patch]*mesh.Sf().boundaryField()[patch])
                        -sum(advecting.boundaryField()[patch]*faceU.boundaryField()[patch]));
            reduce(boundary,sumOp<vector>());
            const scalar radius=before[0].radius, inertiaSphere=0.4*before[0].mass*radius*radius;
            const scalar exactVolume=4.0/3.0*constant::mathematical::pi*pow3(radius);
            std::vector<std::array<double,6>> loads(1);
            for(int k=0;k<3;++k) {loads[0][k]=mode=="free" ? force[k]:0;loads[0][k+3]=mode=="free" ? torque[k]:0;}
            backend.advance(before,loads,substeps,demDt);
            particles=backend.state();
            if(particles.size()!=before.size()||particles[0].id!=before[0].id) throw std::runtime_error("Particle identity changed");
            const vector particleImpulse=before[0].mass*(asFoamVector(particles[0].v)-asFoamVector(before[0].v));
            // Fixed/prescribed particle fixtures have an external support force.
            const vector support=mode=="free" ? vector::zero : -force;
            const scalar momentumResidual=mag(newPhysical-oldPhysical+particleImpulse-(boundary+support)*dt);
            const scalar angularResidual=mode=="free" ? mag(inertiaSphere*(asFoamVector(particles[0].omega)-asFoamVector(before[0].omega))-torque*dt) : 0;
            const scalar clock=std::abs(backend.time()-runTime.value());
            const scalar divMax=gMax(mag(fvc::div(phi))().primitiveField());
            const scalar forceRef=6*constant::mathematical::pi*rho*viscosity*radius*0.01;
            const scalar torqueRef=8*constant::mathematical::pi*rho*viscosity*pow3(radius)*0.1;
            if(rank==0)
            {
                history<<runTime.value()<<','<<particles[0].id<<','<<particles[0].owner;
                for(double v:particles[0].x) history<<','<<v;
                for(double v:particles[0].v) history<<','<<v;
                for(double v:particles[0].omega) history<<','<<v;
                for(int k=0;k<3;++k) history<<','<<force[k];
                for(int k=0;k<3;++k) history<<','<<torque[k];
                history<<','<<std::abs(volume/exactVolume-1)<<','<<(surface ? stress.slip:std::sqrt(slip/volume))
                    <<','<<divMax<<','<<momentumResidual<<','<<angularResidual<<','<<clock<<','<<covered
                    <<','<<force.x()/forceRef<<','<<-torque.z()/torqueRef;
                for(const vector& data:{newPhysical,oldPhysical,boundary,constraint,inertia,constraintTorque,rotationalInertia})
                    for(int k=0;k<3;++k) history<<','<<data[k];
                history<<','<<usedCorrectors<<','<<equationImpulse<<','<<supportVolume/exactVolume
                    <<','<<surfaceCorrectors<<','<<targetDefect<<','<<totalCorrectors<<','<<std::sqrt(slip/volume);
                for(const vector& data:{stress.pressureForce,stress.viscousForce,stress.pressureTorque,stress.viscousTorque})
                    for(int k=0;k<3;++k) history<<','<<data[k];
                history<<','<<std::chrono::duration<double>(std::chrono::steady_clock::now()-wallStart).count()<<'\n';
                history.flush();
            }
            previousP[0]=newP;previousL[0]=newL;previousCentre[0]=asFoamVector(before[0].x);
            runTime.write();
            if(runTime.writeTime())
            {
                IOdictionary checkpoint(IOobject("couplingState",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false));
                checkpoint.add("ids",ids);checkpoint.add("interiorMomentum",previousP);
                checkpoint.add("constraintScheme",constraintScheme);
                checkpoint.add("surfaceReconstruction",reconstruction);
                checkpoint.add("interiorAngular",previousL);checkpoint.add("mappedCentres",previousCentre);
                checkpoint.regIOobject::write();
                if(rank==0) std::filesystem::create_directories(runTime.timeName().c_str());
                MPI_Barrier(MPI_COMM_WORLD);
                backend.writeRestart(std::string(runTime.timeName().c_str())+"/dem.restart");
            }
        }
        if(rank==0) Info<<"M2A_EXECUTION_COMPLETE"<<endl;
    }
    catch(const std::exception& e)
    {
        std::cerr<<"M2A_FAIL rank="<<rank<<": "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);return 1;
    }
    if(ownMPI) MPI_Finalize();
    return 0;
}
