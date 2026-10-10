// SPDX-License-Identifier: GPL-3.0-or-later
// Verification-only replicated Cartesian interpolation. Not a scalable backend.
#ifndef LPBF_SPHERE_SURFACE_H
#define LPBF_SPHERE_SURFACE_H
#include <array>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace lpbfM2
{
class SphereSurface
{
    const Foam::fvMesh& mesh_;
    double h_;
    std::array<double,3> origin_;
    std::array<int,3> size_;
    std::vector<int> indices_;
    std::vector<double> velocity_, pressure_;
    int count_;

    int index(int i,int j,int k) const {return (i*size_[1]+j)*size_[2]+k;}
    std::array<double,4> sample(const Foam::vector& point) const
    {
        std::array<int,3> lower;
        std::array<double,3> fraction;
        for(int a=0;a<3;++a)
        {
            const double x=(point[a]-origin_[a])/h_;
            lower[a]=static_cast<int>(std::floor(x));fraction[a]=x-lower[a];
            if(lower[a]<0||lower[a]+1>=size_[a])
                throw std::runtime_error("Surface probe outside Cartesian interpolation grid");
        }
        std::array<double,4> value{};
        for(int i=0;i<2;++i) for(int j=0;j<2;++j) for(int k=0;k<2;++k)
        {
            const int cell=index(lower[0]+i,lower[1]+j,lower[2]+k);
            const double weight=(i?fraction[0]:1-fraction[0])
                *(j?fraction[1]:1-fraction[1])*(k?fraction[2]:1-fraction[2]);
            for(int a=0;a<3;++a) value[a]+=weight*velocity_[3*cell+a];
            value[3]+=weight*pressure_[cell];
        }
        return value;
    }
    static Foam::vector velocity(const std::array<double,4>& v)
    {return Foam::vector(v[0],v[1],v[2]);}

public:
    struct Diagnostics
    {
        double slip=0;
        Foam::vector pressureForce=Foam::vector::zero,viscousForce=Foam::vector::zero;
        Foam::vector pressureTorque=Foam::vector::zero,viscousTorque=Foam::vector::zero;
    };
    SphereSurface(const Foam::fvMesh& mesh,double h):mesh_(mesh),h_(h)
    {
        double localMin[3]={1e100,1e100,1e100},localMax[3]={-1e100,-1e100,-1e100};
        forAll(mesh_.C(),cell) for(int a=0;a<3;++a)
        {
            localMin[a]=std::min(localMin[a],double(mesh_.C()[cell][a]));
            localMax[a]=std::max(localMax[a],double(mesh_.C()[cell][a]));
        }
        double maximum[3];
        MPI_Allreduce(localMin,origin_.data(),3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(localMax,maximum,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        long long total=1;
        for(int a=0;a<3;++a) {size_[a]=std::llround((maximum[a]-origin_[a])/h_)+1;total*=size_[a];}
        if(total<1||total>500000) throw std::runtime_error("Surface reference limited to 500000 Cartesian cells");
        count_=static_cast<int>(total);
        std::vector<int> ownership(count_,0);indices_.resize(mesh_.nCells());
        forAll(mesh_.C(),cell)
        {
            int ijk[3];
            for(int a=0;a<3;++a)
            {
                const double x=(mesh_.C()[cell][a]-origin_[a])/h_;
                ijk[a]=std::llround(x);
                if(ijk[a]<0||ijk[a]>=size_[a]||std::abs(x-ijk[a])>1e-7)
                    throw std::runtime_error("Invalid Cartesian centre index");
            }
            indices_[cell]=index(ijk[0],ijk[1],ijk[2]);++ownership[indices_[cell]];
        }
        MPI_Allreduce(MPI_IN_PLACE,ownership.data(),count_,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        for(int owners:ownership) if(owners!=1) throw std::runtime_error("Missing or duplicate Cartesian cell owner");
        velocity_.resize(3*count_);pressure_.resize(count_);
    }
    void gather(const Foam::volVectorField& U,const Foam::volScalarField& p)
    {
        std::fill(velocity_.begin(),velocity_.end(),0);std::fill(pressure_.begin(),pressure_.end(),0);
        forAll(mesh_.C(),cell)
        {
            const int id=indices_[cell];pressure_[id]=p[cell];
            for(int a=0;a<3;++a) velocity_[3*id+a]=U[cell][a];
        }
        MPI_Allreduce(MPI_IN_PLACE,velocity_.data(),3*count_,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,pressure_.data(),count_,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    // Evaluate ghost target at an actual cell centre, using only fluid-side
    // interpolation for the image point R+2h (all stencil nodes lie outside R).
    Foam::vector target(const Foam::vector& point,const Foam::vector& centre,
                        double radius,const Foam::vector& v,const Foam::vector& omega) const
    {
        const Foam::vector r=point-centre;
        const double distance=Foam::mag(r);
        // Cover every interior node of a trilinear wall stencil (sqrt(3)h).
        if(distance<=radius-std::sqrt(3.0)*h_ || distance<1e-30) return v+(omega^r);
        const Foam::vector n=r/distance, wall=centre+radius*n;
        const Foam::vector wallVelocity=v+(omega^(radius*n));
        return wallVelocity+(distance-radius)/(2*h_)
            *(velocity(sample(wall+2*h_*n))-wallVelocity);
    }
    Diagnostics diagnose(const Foam::vector& centre,double radius,
                         const Foam::vector& v,const Foam::vector& omega,double rho,double nu) const
    {
        Diagnostics result;
        const double pi=Foam::constant::mathematical::pi;
        // 12-point Gauss-Legendre in cos(theta), 24 periodic azimuth points.
        for(int iz=0;iz<12;++iz)
        {
            double z=std::cos(pi*(iz+.75)/(12+.5)), derivative=0;
            bool converged=false;
            for(int iteration=0;iteration<32;++iteration)
            {
                double a=1,b=0;
                for(int degree=1;degree<=12;++degree)
                {const double c=b;b=a;a=((2*degree-1)*z*b-(degree-1)*c)/degree;}
                derivative=12*(z*a-b)/(z*z-1);
                const double change=a/derivative;z-=change;
                if(std::abs(change)<1e-15) {converged=true;break;}
            }
            if(!converged) throw std::runtime_error("Surface quadrature root did not converge");
            const double weight=2/((1-z*z)*derivative*derivative)*2*pi/24*radius*radius;
            for(int ip=0;ip<24;++ip)
            {
                const double phi=2*pi*(ip+.5)/24;
                const Foam::vector n(std::sqrt(1-z*z)*std::cos(phi),std::sqrt(1-z*z)*std::sin(phi),z);
                const Foam::vector wall=centre+radius*n, ub=v+(omega^(radius*n));
                const auto a=sample(wall+2*h_*n), b=sample(wall+3*h_*n), c=sample(wall+4*h_*n);
                const Foam::vector radial=(-13.0/12*ub+3*velocity(a)-8.0/3*velocity(b)+.75*velocity(c))/h_;
                // At a rigid spherical wall, tangential derivatives of Ub
                // give (gradU+gradU^T)n = radial+n*(radial.n)-omega x n.
                const Foam::vector viscous=rho*nu*(radial+n*(radial&n)-(omega^n));
                const double wallPressure=6*a[3]-8*b[3]+3*c[3]; // kinematic -> Pa via rho
                const Foam::vector pressure=-rho*wallPressure*n;
                result.pressureForce+=weight*pressure;result.viscousForce+=weight*viscous;
                result.pressureTorque+=weight*((radius*n)^pressure);
                result.viscousTorque+=weight*((radius*n)^viscous);
                result.slip+=weight*Foam::magSqr(velocity(sample(wall))-ub);
            }
        }
        result.slip=std::sqrt(result.slip/(4*pi*radius*radius));
        return result;
    }
};
}
#endif
