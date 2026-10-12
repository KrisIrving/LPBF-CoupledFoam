// SPDX-License-Identifier: GPL-3.0-or-later
// Static verification-only J / volume-weighted transpose S. Replicated markers.
#ifndef LPBF_JOINT_SPHERE_H
#define LPBF_JOINT_SPHERE_H
#include "SphereSurface.H"
#include "WallGMRES.H"
#include "CompleteWallBasis.H"
#include <set>
#include <cstring>
#include <string>

namespace lpbfM2
{
class JointSphere
{
    const Foam::fvMesh& mesh_;
    SphereSurface& grid_;
    std::vector<Foam::vector> points_;
    std::vector<SphereSurface::Stencil> stencils_;
    std::vector<std::vector<int>> localStencils_;
    bool sparse_;
    bool complete_;
    std::vector<std::vector<double>> wallCoordinates_;
    std::vector<double> basisSums_;
    std::vector<Foam::vector> basisFirst_;
    double radius_,coverage_=0,orthogonality_=0;
    double minPivot_=1;
public:
    JointSphere(const Foam::fvMesh& mesh,SphereSurface& grid,const Foam::vector& centre,double radius,double h,bool sparse=false,bool fluidCubic=false,bool complete=false)
    :mesh_(mesh),grid_(grid),sparse_(sparse),complete_(complete),radius_(radius)
    {
        const double pi=Foam::constant::mathematical::pi;
        int count=static_cast<int>(std::ceil(4*pi*radius*radius/(2.25*h*h)));
        if(count<8||count>512) throw std::runtime_error("Joint static marker budget outside [8,512]");
        if(complete_)
        {
            points_=grid_.wallPoints(centre,radius);std::vector<SphereSurface::Stencil> full;std::set<int> support;
            for(const auto& point:points_)
            {full.push_back(grid_.stencil(point));for(int cell:full.back().cells) support.insert(cell);}
            std::vector<int> cells(support.begin(),support.end());std::vector<std::vector<double>> original;
            for(const auto& row:full)
            {std::vector<double> values(cells.size(),0);for(std::size_t k=0;k<row.cells.size();++k)
             values[std::lower_bound(cells.begin(),cells.end(),row.cells[k])-cells.begin()]=row.weights[k];
             original.push_back(values);}
            const CompleteWallBasis basis(original);wallCoordinates_=basis.coordinates;
            coverage_=basis.coverage;orthogonality_=basis.orthogonality;
            for(const auto& row:basis.rows)
            {SphereSurface::Stencil stencil;stencil.cells=cells;stencil.weights=row;stencils_.push_back(stencil);
             double sum=0;Foam::vector first=Foam::vector::zero;
             for(std::size_t k=0;k<cells.size();++k){sum+=row[k];first+=row[k]*grid_.cellCentre(cells[k]);}
             basisSums_.push_back(sum);basisFirst_.push_back(first);}
            count=stencils_.size();
            if(count<8||count>512) throw std::runtime_error("Complete wall modal budget outside [8,512]");
        }
        else for(int i=0;i<count;++i)
        {
            const double z=1-2*(i+.5)/count,phi=pi*(3-std::sqrt(5.0))*i;
            const Foam::vector n(std::sqrt(1-z*z)*std::cos(phi),std::sqrt(1-z*z)*std::sin(phi),z);
            points_.push_back(centre+radius*n);
            stencils_.push_back(fluidCubic ? grid_.fluidStencil(points_.back(),centre,radius):grid_.stencil(points_.back()));
        }
        // Owner lookup is frozen with the static geometry/partition. The grid
        // constructor already checks exactly one global owner for every cell.
        std::vector<int> owner(grid_.globalCells(),-1);
        forAll(mesh_.C(),cell) owner[grid_.globalIndex(cell)]=cell;
        for(const auto& stencil:stencils_)
        {
            std::vector<int> local(stencil.cells.size());
            for(std::size_t k=0;k<stencil.cells.size();++k) local[k]=owner[stencil.cells[k]];
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
            std::size_t a=0,b=0;
            while(a<stencils_[i].cells.size()&&b<stencils_[j].cells.size())
            {
                const int ca=stencils_[i].cells[a],cb=stencils_[j].cells[b];
                if(ca==cb){value+=stencils_[i].weights[a]*stencils_[j].weights[b];++a;++b;}
                else if(ca<cb) ++a;else ++b;
            }
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
    std::array<double,5> reconstructionStats() const
    {
        double moment=0,l1=0,minimum=1e30,maximum=0,expanded=0;
        for(const auto& stencil:stencils_)
        {moment=std::max(moment,stencil.moment);l1=std::max(l1,stencil.l1);
         minimum=std::min(minimum,double(stencil.cells.size()));maximum=std::max(maximum,double(stencil.cells.size()));
         if(stencil.expanded) ++expanded;}
        return {moment,l1,minimum,maximum,expanded};
    }
    int count() const {return stencils_.size();}
    bool completeWall() const {return complete_;}
    Foam::word basisSignature() const
    {
        unsigned long long hash=14695981039346656037ULL;
        auto add=[&](unsigned long long bits){for(int i=0;i<8;++i){hash^=(bits>>(8*i))&255ULL;hash*=1099511628211ULL;}};
        for(const auto& row:stencils_) for(std::size_t i=0;i<row.cells.size();++i)
        {add(row.cells[i]);unsigned long long bits=0;static_assert(sizeof(bits)==sizeof(double),"64-bit basis signature");
         std::memcpy(&bits,&row.weights[i],sizeof(bits));add(bits);}
        return Foam::word("b"+std::to_string(hash));
    }
    std::array<double,4> basisStats() const {return {double(points_.size()),double(count()),coverage_,orthogonality_};}
    std::array<double,2> rigidResidual(const std::vector<double>& values,
        const Foam::vector& centre,const Foam::vector& v,const Foam::vector& omega) const
    {
        double sum=0,maximum=0;
        const int samples=complete_ ? points_.size():count();
        for(int i=0;i<samples;++i)
        {
            const Foam::vector target=v+(omega^(points_[i]-centre));
            Foam::vector value=Foam::vector::zero;
            if(complete_) for(int j=0;j<count();++j) for(int a=0;a<3;++a) value[a]+=wallCoordinates_[i][j]*values[3*j+a];
            else value=Foam::vector(values[3*i],values[3*i+1],values[3*i+2]);
            const double error=Foam::mag(value-target);
            sum+=error*error;maximum=std::max(maximum,error);
        }
        return {std::sqrt(sum/samples),maximum};
    }
    std::vector<double> normalProbe(const Foam::vector& centre,double amplitude) const
    {
        std::vector<double> result(3*count());
        for(int i=0;i<count();++i)
        {
            const Foam::vector n=complete_ ? (basisFirst_[i]-basisSums_[i]*centre)/radius_:
                (points_[i]-centre)/Foam::mag(points_[i]-centre);
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
            if(complete_){total+=basisSums_[i]*f;torque+=(basisFirst_[i]-basisSums_[i]*centre)^f;}
            else {total+=f;torque+=(points_[i]-centre)^f;}
        }
    }
    std::vector<double> interpolateDense(const Foam::volVectorField& U,const Foam::volScalarField& p)
    {
        grid_.gather(U,p);std::vector<double> result(3*count());
        for(int i=0;i<count();++i)
        {
            const Foam::vector u=grid_.interpolateStencil(stencils_[i]);
            for(int a=0;a<3;++a) result[3*i+a]=u[a];
        }
        return result;
    }
    std::vector<double> interpolateSparse(const Foam::volVectorField& U) const
    {
        std::vector<double> result(3*count(),0);
        for(int i=0;i<count();++i) for(std::size_t k=0;k<stencils_[i].cells.size();++k)
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
            const Foam::vector ub=complete_ ? basisSums_[i]*v+(omega^(basisFirst_[i]-basisSums_[i]*centre)):
                v+(omega^(points_[i]-centre));
            for(int a=0;a<3;++a) result[3*i+a]=ub[a];
        }
        return result;
    }
    Foam::vectorField accelerationDense(const std::vector<double>& integratedForces,double rho) const
    {
        if(integratedForces.size()!=3*std::size_t(count())) throw std::runtime_error("Wall load dimensions differ");
        std::vector<double> global(3*grid_.globalCells(),0);
        for(int i=0;i<count();++i) for(std::size_t k=0;k<stencils_[i].cells.size();++k) for(int a=0;a<3;++a)
            global[3*stencils_[i].cells[k]+a]+=stencils_[i].weights[k]*integratedForces[3*i+a];
        Foam::vectorField result(mesh_.nCells(),Foam::vector::zero);
        forAll(result,cell) for(int a=0;a<3;++a)
            result[cell][a]=global[3*grid_.globalIndex(cell)+a]/(rho*mesh_.V()[cell]);
        return result;
    }
    Foam::vectorField accelerationSparse(const std::vector<double>& integratedForces,double rho) const
    {
        if(integratedForces.size()!=3*std::size_t(count())) throw std::runtime_error("Wall load dimensions differ");
        Foam::vectorField result(mesh_.nCells(),Foam::vector::zero);
        for(int i=0;i<count();++i) for(std::size_t k=0;k<stencils_[i].cells.size();++k)
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
        for(int i=0;i<count();++i) for(std::size_t k=0;k<stencils_[i].cells.size();++k)
            result[i]+=stencils_[i].weights[k]*stencils_[i].weights[k]*global[stencils_[i].cells[k]];
        for(double value:result) if(!std::isfinite(value)||value<=0)
            throw std::runtime_error("Invalid joint wall scaling");
        return result;
    }
    std::vector<double> inverseDiagonalSparse(const Foam::volScalarField& rAU,double rho) const
    {
        std::vector<double> result(count(),0);
        for(int i=0;i<count();++i) for(std::size_t k=0;k<stencils_[i].cells.size();++k)
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
