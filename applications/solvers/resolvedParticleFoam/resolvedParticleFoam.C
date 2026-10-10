// SPDX-License-Identifier: GPL-3.0-or-later
// Experimental Cartesian spherical mechanical verification solver, M2A-01.
#include "fvCFD.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "MechanicalBackend.H"
#include <fstream>
#include <filesystem>
#include <iomanip>

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
        auto mapGeometry=[&]()
        {
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
                history<<"time,id,owner,x,y,z,vx,vy,vz,wx,wy,wz,fx,fy,fz,tx,ty,tz,volume_error,slip_rms,div_max,momentum_residual,angular_impulse_residual,clock_error,covered_ranks,force_ratio,torque_ratio,fluid_px,fluid_py,fluid_pz,old_fluid_px,old_fluid_py,old_fluid_pz,boundary_fx,boundary_fy,boundary_fz,constraint_fx,constraint_fy,constraint_fz,inertia_fx,inertia_fy,inertia_fz,constraint_tx,constraint_ty,constraint_tz,inertia_tx,inertia_ty,inertia_tz,pressure_correctors,momentum_equation_impulse_L1,support_volume_ratio\n";
        }
        const scalar refPressure=p[0];
        while(runTime.run())
        {
            if(mag(runTime.deltaTValue()-dt)>dt*1e-12) throw std::runtime_error("Adaptive dt unsupported in first package");
            const auto before=particles;
            const vector oldPhysical=gSum(U.primitiveField()*mesh.V().field())*rho-previousP[0];
            ++runTime;
            mapGeometry();
            surfaceScalarField advecting(IOobject("advecting",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),phi);
            fvVectorMatrix equation(fvm::ddt(U)+fvm::div(advecting,U)-fvm::laplacian(nu,U)+fvm::Sp(lambda,U)==lambda*rigid);
            // One momentum predictor per window; inner PISO corrections update
            // H from corrected U without repeatedly solving the predictor.
            solve(equation==-fvc::grad(p));
            volScalarField rAU("rAU",1.0/equation.A());
            label usedCorrectors=0;
            scalar equationImpulse=GREAT;
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
                // residual() is source - A*U, integrated over each cell.
                // Include the pressure source and sum magnitudes, so opposite
                // local defects cannot cancel in a global momentum ledger.
                vectorField defect(equation.residual());
                defect-=mesh.V().field()*fvc::grad(p)().primitiveField();
                equationImpulse=rho*dt*gSum(mag(defect));
                usedCorrectors=c+1;
                if(rank==0) Info<<"M2A_CORRECTION time="<<runTime.value()
                    <<" iteration="<<c+1<<" maxU="<<maxU<<" maxP="<<maxP
                    <<" equationImpulseL1="<<equationImpulse<<endl;
                if(!std::isfinite(maxU)||!std::isfinite(maxP)||!std::isfinite(equationImpulse))
                    throw std::runtime_error("Nonfinite CFD state before DEM feedback");
                if(usedCorrectors>=correctors && (momentumTolerance<0||equationImpulse<=momentumTolerance)) break;
            }
            if(momentumTolerance>=0 && equationImpulse>momentumTolerance)
                throw std::runtime_error("CFD momentum equation did not converge before DEM feedback");
            vector constraint=vector::zero, constraintTorque=vector::zero;
            scalar volume=0,slip=0,supportVolume=0;label covered=0;
            forAll(mesh.C(),cell)
            {
                const vector f=rho*lambda[cell]*(rigid[cell]-U[cell])*mesh.V()[cell];
                constraint+=f;constraintTorque+=(mesh.C()[cell]-asFoamVector(before[0].x))^f;
                volume+=solid[cell]*mesh.V()[cell];slip+=solid[cell]*mesh.V()[cell]*magSqr(U[cell]-rigid[cell]);
                if(solid[cell]>0) supportVolume+=mesh.V()[cell];
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
                history<<','<<std::abs(volume/exactVolume-1)<<','<<std::sqrt(slip/volume)
                    <<','<<divMax<<','<<momentumResidual<<','<<angularResidual<<','<<clock<<','<<covered
                    <<','<<force.x()/forceRef<<','<<-torque.z()/torqueRef;
                for(const vector& data:{newPhysical,oldPhysical,boundary,constraint,inertia,constraintTorque,rotationalInertia})
                    for(int k=0;k<3;++k) history<<','<<data[k];
                history<<','<<usedCorrectors<<','<<equationImpulse<<','<<supportVolume/exactVolume<<'\n';
                history.flush();
            }
            previousP[0]=newP;previousL[0]=newL;previousCentre[0]=asFoamVector(before[0].x);
            runTime.write();
            if(runTime.writeTime())
            {
                IOdictionary checkpoint(IOobject("couplingState",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false));
                checkpoint.add("ids",ids);checkpoint.add("interiorMomentum",previousP);
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
