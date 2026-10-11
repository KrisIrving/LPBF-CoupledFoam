// SPDX-License-Identifier: GPL-3.0-or-later
// Static verification-only J / volume-weighted transpose S. Replicated markers.
#ifndef LPBF_JOINT_SPHERE_H
#define LPBF_JOINT_SPHERE_H
#include "SphereSurface.H"
#include "WallGMRES.H"

namespace lpbfM2
{
class JointSphere
{
    const Foam::fvMesh& mesh_;
    SphereSurface& grid_;
    std::vector<Foam::vector> points_;
    std::vector<SphereSurface::Stencil> stencils_;
    std::vector<std::array<int,8>> localStencils_;
    bool sparse_;
    double minPivot_=1;
public:
    JointSphere(const Foam::fvMesh& mesh,SphereSurface& grid,const Foam::vector& centre,double radius,double h,bool sparse=false)
    :mesh_(mesh),grid_(grid),sparse_(sparse)
    {
        const double pi=Foam::constant::mathematical::pi;
        const int count=static_cast<int>(std::ceil(4*pi*radius*radius/(2.25*h*h)));
        if(count<8||count>512) throw std::runtime_error("Joint static marker budget outside [8,512]");
        for(int i=0;i<count;++i)
        {
            const double z=1-2*(i+.5)/count,phi=pi*(3-std::sqrt(5.0))*i;
            const Foam::vector n(std::sqrt(1-z*z)*std::cos(phi),std::sqrt(1-z*z)*std::sin(phi),z);
            points_.push_back(centre+radius*n);stencils_.push_back(grid_.stencil(points_.back()));
        }
        // Owner lookup is frozen with the static geometry/partition. The grid
        // constructor already checks exactly one global owner for every cell.
        std::vector<int> owner(grid_.globalCells(),-1);
        forAll(mesh_.C(),cell) owner[grid_.globalIndex(cell)]=cell;
        for(const auto& stencil:stencils_)
        {
            std::array<int,8> local;
            for(int k=0;k<8;++k) local[k]=owner[stencil.cells[k]];
            localStencils_.push_back(local);
        }
        // Normalised J J^T Cholesky without diagonal regularisation. This checks
        // interpolation-row independence only, not pressure-projected Schur rank.
        std::vector<double> diagonal(count),factor(count*count,0);
        for(int i=0;i<count;++i)
            for(double w:stencils_[i].weights) diagonal[i]+=w*w;
        for(int i=0;i<count;++i) for(int j=0;j<=i;++j)
        {
            double value=0;
            for(int a=0;a<8;++a) for(int b=0;b<8;++b)
                if(stencils_[i].cells[a]==stencils_[j].cells[b])
                    value+=stencils_[i].weights[a]*stencils_[j].weights[b];
            value/=std::sqrt(diagonal[i]*diagonal[j]);
            for(int k=0;k<j;++k) value-=factor[i*count+k]*factor[j*count+k];
            if(i==j)
            {
                minPivot_=std::min(minPivot_,value);
                if(!std::isfinite(value)||value<=1e-8)
                    throw std::runtime_error("Dependent or ill-conditioned wall interpolation rows");
                factor[i*count+j]=std::sqrt(value);
            }
            else factor[i*count+j]=value/factor[j*count+j];
        }
    }
    int count() const {return points_.size();}
    std::vector<double> normalProbe(const Foam::vector& centre,double amplitude) const
    {
        std::vector<double> result(3*count());
        for(int i=0;i<count();++i)
        {
            const Foam::vector n=(points_[i]-centre)/Foam::mag(points_[i]-centre);
            for(int a=0;a<3;++a) result[3*i+a]=amplitude*n[a];
        }
        return result;
    }
    double minPivot() const {return minPivot_;}
    void loads(const std::vector<double>& forces,const Foam::vector& centre,
               Foam::vector& total,Foam::vector& torque) const
    {
        total=Foam::vector::zero;torque=Foam::vector::zero;
        for(int i=0;i<count();++i)
        {
            const Foam::vector f(forces[3*i],forces[3*i+1],forces[3*i+2]);
            total+=f;torque+=(points_[i]-centre)^f;
        }
    }
    std::vector<double> interpolateDense(const Foam::volVectorField& U,const Foam::volScalarField& p)
    {
        grid_.gather(U,p);std::vector<double> result(3*count());
        for(int i=0;i<count();++i)
        {
            const Foam::vector u=grid_.interpolateVelocity(points_[i]);
            for(int a=0;a<3;++a) result[3*i+a]=u[a];
        }
        return result;
    }
    std::vector<double> interpolateSparse(const Foam::volVectorField& U) const
    {
        std::vector<double> result(3*count(),0);
        for(int i=0;i<count();++i) for(int k=0;k<8;++k)
        {
            const int cell=localStencils_[i][k];if(cell<0) continue;
            for(int a=0;a<3;++a) result[3*i+a]+=stencils_[i].weights[k]*U[cell][a];
        }
        MPI_Allreduce(MPI_IN_PLACE,result.data(),static_cast<int>(result.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        return result;
    }
    std::vector<double> interpolate(const Foam::volVectorField& U,const Foam::volScalarField& p)
    {return sparse_ ? interpolateSparse(U):interpolateDense(U,p);}
    std::vector<double> desired(const Foam::vector& centre,const Foam::vector& v,const Foam::vector& omega) const
    {
        std::vector<double> result(3*count());
        for(int i=0;i<count();++i)
        {
            const Foam::vector ub=v+(omega^(points_[i]-centre));
            for(int a=0;a<3;++a) result[3*i+a]=ub[a];
        }
        return result;
    }
    Foam::vectorField accelerationDense(const std::vector<double>& integratedForces,double rho) const
    {
        if(integratedForces.size()!=3*points_.size()) throw std::runtime_error("Wall load dimensions differ");
        std::vector<double> global(3*grid_.globalCells(),0);
        for(int i=0;i<count();++i) for(int k=0;k<8;++k) for(int a=0;a<3;++a)
            global[3*stencils_[i].cells[k]+a]+=stencils_[i].weights[k]*integratedForces[3*i+a];
        Foam::vectorField result(mesh_.nCells(),Foam::vector::zero);
        forAll(result,cell) for(int a=0;a<3;++a)
            result[cell][a]=global[3*grid_.globalIndex(cell)+a]/(rho*mesh_.V()[cell]);
        return result;
    }
    Foam::vectorField accelerationSparse(const std::vector<double>& integratedForces,double rho) const
    {
        if(integratedForces.size()!=3*points_.size()) throw std::runtime_error("Wall load dimensions differ");
        Foam::vectorField result(mesh_.nCells(),Foam::vector::zero);
        for(int i=0;i<count();++i) for(int k=0;k<8;++k)
        {
            const int cell=localStencils_[i][k];if(cell<0) continue;
            for(int a=0;a<3;++a) result[cell][a]+=stencils_[i].weights[k]*integratedForces[3*i+a];
        }
        forAll(result,cell) result[cell]/=rho*mesh_.V()[cell];
        return result;
    }
    Foam::vectorField acceleration(const std::vector<double>& integratedForces,double rho) const
    {return sparse_ ? accelerationSparse(integratedForces,rho):accelerationDense(integratedForces,rho);}
    std::vector<double> inverseDiagonalDense(const Foam::volScalarField& rAU,double rho) const
    {
        // J rAU S acceleration response, ignoring pressure, for right scaling.
        std::vector<double> global(grid_.globalCells(),0),result(count(),0);
        forAll(mesh_.C(),cell) global[grid_.globalIndex(cell)]=rAU[cell]/(rho*mesh_.V()[cell]);
        MPI_Allreduce(MPI_IN_PLACE,global.data(),grid_.globalCells(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        for(int i=0;i<count();++i) for(int k=0;k<8;++k)
            result[i]+=stencils_[i].weights[k]*stencils_[i].weights[k]*global[stencils_[i].cells[k]];
        for(double value:result) if(!std::isfinite(value)||value<=0)
            throw std::runtime_error("Invalid joint wall scaling");
        return result;
    }
    std::vector<double> inverseDiagonalSparse(const Foam::volScalarField& rAU,double rho) const
    {
        std::vector<double> result(count(),0);
        for(int i=0;i<count();++i) for(int k=0;k<8;++k)
        {
            const int cell=localStencils_[i][k];if(cell<0) continue;
            const double w=stencils_[i].weights[k];
            result[i]+=w*w*rAU[cell]/(rho*mesh_.V()[cell]);
        }
        MPI_Allreduce(MPI_IN_PLACE,result.data(),static_cast<int>(result.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        for(double value:result) if(!std::isfinite(value)||value<=0)
            throw std::runtime_error("Invalid joint wall scaling");
        return result;
    }
    std::vector<double> inverseDiagonal(const Foam::volScalarField& rAU,double rho) const
    {return sparse_ ? inverseDiagonalSparse(rAU,rho):inverseDiagonalDense(rAU,rho);}
    std::array<double,3> auditTransfer(const Foam::volVectorField& U,const Foam::volScalarField& p,
                                     const Foam::volScalarField& rAU,double rho)
    {
        auto relative=[](const std::vector<double>& a,const std::vector<double>& b)
        {double error=0,scale=1e-30;for(std::size_t i=0;i<a.size();++i){error=std::max(error,std::abs(a[i]-b[i]));scale=std::max(scale,std::abs(b[i]));}return error/scale;};
        const double interpolation=relative(interpolateSparse(U),interpolateDense(U,p));
        std::vector<double> force(3*count());
        for(std::size_t i=0;i<force.size();++i) force[i]=1e-6*std::sin(.73*(i+1));
        const Foam::vectorField dense(accelerationDense(force,rho)),sparse(accelerationSparse(force,rho));
        const double spreading=Foam::gMax(Foam::mag(sparse-dense))/Foam::max(Foam::gMax(Foam::mag(dense)),Foam::scalar(1e-30));
        const double scaling=relative(inverseDiagonalSparse(rAU,rho),inverseDiagonalDense(rAU,rho));
        return {interpolation,spreading,scaling};
    }
};

}
#endif
