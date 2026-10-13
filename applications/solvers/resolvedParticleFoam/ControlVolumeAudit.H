// SPDX-License-Identifier: GPL-3.0-or-later
// Static Cartesian audit. Never feeds a second load to the particle backend.
#ifndef LPBF_CONTROL_VOLUME_AUDIT_H
#define LPBF_CONTROL_VOLUME_AUDIT_H
#include "fvCFD.H"
#include <fstream>
#include <iomanip>

namespace lpbfM2
{
inline void controlVolumeAudit
(
    const Foam::fvMesh& mesh, const Foam::volVectorField& U,
    const Foam::volScalarField& p, const Foam::surfaceScalarField& advecting,
    const Foam::surfaceVectorField& faceU, const Foam::surfaceScalarField& faceP,
    const Foam::surfaceVectorField& diffusion, const Foam::vectorField& source,
    const Foam::vectorField& matrixDefect, const Foam::vector& centre,
    Foam::scalar radius, Foam::scalar h, Foam::scalar rho, Foam::scalar viscosity,
    Foam::scalar dt, const Foam::vector& fullSource, const Foam::vector& fullTorque,
    std::ofstream& history, bool master
)
{
    using namespace Foam;
    // Reconstruct the explicit form of the ACTUAL Euler/Gauss-linear equation.
    // This is checked against the retained final coupled matrix residual.
    const volVectorField convection(fvc::div(advecting,U));
    const volVectorField viscous(fvc::laplacian(dimensionedScalar("auditNu",dimViscosity,viscosity),U));
    const volVectorField gradient(fvc::grad(p));
    const surfaceVectorField net(diffusion-faceP*mesh.Sf()-advecting*faceU);
    // A second, diagnostic traction uses the symmetric viscous stress. Its
    // difference from the solved vector-Laplacian flux is NOT a second source.
    const surfaceVectorField symmetricDiffusion
    (
        2.0*viscosity*(fvc::interpolate(symm(fvc::grad(U)))&mesh.Sf())
    );
    const surfaceVectorField symmetricExtra(symmetricDiffusion-diffusion);
    vectorField storage(mesh.nCells()), cellNet(mesh.nCells()), imbalance(mesh.nCells());
    scalar match=0;
    forAll(storage,cell)
    {
        storage[cell]=rho*mesh.V()[cell]*(U[cell]-U.oldTime()[cell])/dt;
        cellNet[cell]=rho*mesh.V()[cell]*(viscous[cell]-gradient[cell]-convection[cell]);
        imbalance[cell]=storage[cell]-source[cell]-cellNet[cell];
        match+=mag(imbalance[cell]+rho*matrixDefect[cell])*dt;
    }
    reduce(match,sumOp<scalar>());
    if(!std::isfinite(match)||match>1e-12)
        throw std::runtime_error("Control volume explicit operators differ from final matrix residual");
    for(const int padding:{6,8})
    {
        const scalar half=radius+padding*h;
        boolList inside(mesh.nCells(),false);
        label cells=0;
        vector inventory=vector::zero, inventoryTorque=vector::zero;
        vector integratedSource=vector::zero, sourceTorque=vector::zero;
        vector operatorNet=vector::zero, operatorTorque=vector::zero;
        vector residual=vector::zero, residualTorque=vector::zero;
        scalar impulseL1=0;
        forAll(inside,cell)
        {
            const vector arm=mesh.C()[cell]-centre;
            inside[cell]=mag(arm.x())<half&&mag(arm.y())<half&&mag(arm.z())<half;
            if(!inside[cell]) continue;
            ++cells; inventory+=storage[cell]; inventoryTorque+=arm^storage[cell];
            integratedSource+=source[cell]; sourceTorque+=arm^source[cell];
            operatorNet+=cellNet[cell]; operatorTorque+=arm^cellNet[cell];
            residual+=imbalance[cell]; residualTorque+=arm^imbalance[cell];
            impulseL1+=mag(imbalance[cell])*dt;
        }
        label exteriorCells=0;
        forAll(mesh.boundary(),patch) if(!mesh.boundary()[patch].coupled())
            for(const label cell:mesh.boundary()[patch].faceCells()) if(inside[cell]) ++exteriorCells;
        reduce(exteriorCells,sumOp<label>());
        if(exteriorCells) throw std::runtime_error("Control volume touches physical outer boundary");
        vector faceNet=vector::zero, faceTorque=vector::zero, armCorrection=vector::zero;
        vector convective=vector::zero, pressure=vector::zero, diffusive=vector::zero;
        vector symmetricExtraForce=vector::zero, symmetricExtraTorque=vector::zero;
        auto side=[&](label cell,const vector& face,const vector& flux,
                      const vector& conv,const vector& press,const vector& diff,
                      const vector& symmetricCorrection)
        {
            if(!inside[cell]) return;
            faceNet+=rho*flux; faceTorque+=rho*((face-centre)^flux);
            armCorrection+=rho*((mesh.C()[cell]-face)^flux);
            convective+=rho*conv; pressure+=rho*press; diffusive+=rho*diff;
            symmetricExtraForce+=rho*symmetricCorrection;
            symmetricExtraTorque+=rho*((face-centre)^symmetricCorrection);
        };
        // Cell-side enumeration: internal faces cancel if both cells selected.
        // Processor faces are visited on each local cell side; opposing outward
        // fluxes cancel after reduction. No processor face is counted as an
        // extra exterior traction or ownership-independent particle force.
        const auto& owner=mesh.owner(); const auto& neighbour=mesh.neighbour();
        forAll(neighbour,face)
        {
            const vector conv=-advecting[face]*faceU[face];
            const vector press=-faceP[face]*mesh.Sf()[face];
            side(owner[face],mesh.Cf()[face],net[face],conv,press,diffusion[face],symmetricExtra[face]);
            side(neighbour[face],mesh.Cf()[face],-net[face],-conv,-press,-diffusion[face],-symmetricExtra[face]);
        }
        forAll(mesh.boundary(),patch)
        {
            const labelUList& patchCells=mesh.boundary()[patch].faceCells();
            forAll(patchCells,face)
                side(patchCells[face],mesh.Cf().boundaryField()[patch][face],
                     net.boundaryField()[patch][face],
                     -advecting.boundaryField()[patch][face]*faceU.boundaryField()[patch][face],
                     -faceP.boundaryField()[patch][face]*mesh.Sf().boundaryField()[patch][face],
                     diffusion.boundaryField()[patch][face],symmetricExtra.boundaryField()[patch][face]);
        }
        reduce(cells,sumOp<label>());reduce(impulseL1,sumOp<scalar>());
        for(vector* value:{&inventory,&inventoryTorque,&integratedSource,&sourceTorque,
            &operatorNet,&operatorTorque,&residual,&residualTorque,&faceNet,&faceTorque,
            &armCorrection,&convective,&pressure,&diffusive,&symmetricExtraForce,
            &symmetricExtraTorque}) reduce(*value,sumOp<vector>());
        const vector symmetricFaceNet(faceNet+symmetricExtraForce);
        const vector symmetricFaceTorque(faceTorque+symmetricExtraTorque);
        const scalar faceError=mag(faceNet-operatorNet)*dt;
        const scalar angularError=mag(faceTorque+armCorrection-operatorTorque)*dt;
        const scalar sourceError=mag(integratedSource-fullSource);
        const scalar torqueError=mag(sourceTorque-fullTorque);
        if(master)
        {
            history<<mesh.time().value()<<','<<padding<<','<<cells<<','<<match<<','<<impulseL1
                <<','<<faceError<<','<<angularError<<','<<sourceError<<','<<torqueError;
            for(const vector& value:{inventory,inventoryTorque,integratedSource,sourceTorque,
                convective,pressure,diffusive,faceNet,faceTorque,armCorrection,residual,residualTorque,
                symmetricExtraForce,symmetricExtraTorque,symmetricFaceNet,symmetricFaceTorque})
                for(int a=0;a<3;++a) history<<','<<value[a];
            history<<'\n';history.flush();
        }
        if(cells==0||!std::isfinite(faceError)||!std::isfinite(angularError)
            ||!std::isfinite(sourceError)||!std::isfinite(torqueError)||!std::isfinite(impulseL1)
            ||!std::isfinite(mag(symmetricFaceNet))||!std::isfinite(mag(symmetricFaceTorque))
            ||faceError>1e-12||angularError>radius*1e-12||sourceError>1e-12||torqueError>1e-14)
            throw std::runtime_error("Control volume face/operator/source identity failed");
    }
}
}
#endif
