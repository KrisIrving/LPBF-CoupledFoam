// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LPBF_WALL_GMRES_H
#define LPBF_WALL_GMRES_H
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <vector>
namespace lpbfM2
{
// Unrestarted bounded GMRES. FV face projection need not make the wall Schur
// response symmetric; do not use CG based on the toy adjoint algebra.
inline std::vector<double> wallGMRES(const std::function<std::vector<double>(const std::vector<double>&)>& apply,
    const std::vector<double>& rhs,double tolerance,int budget,int& iterations)
{
    if(budget<1||!std::isfinite(tolerance)||tolerance<=0)
        throw std::runtime_error("Invalid joint Krylov controls");
    for(double value:rhs) if(!std::isfinite(value)) throw std::runtime_error("Nonfinite joint Krylov right-hand side");
    const int n=rhs.size(),limit=std::min(budget,n);
    auto dot=[](const std::vector<double>& a,const std::vector<double>& b)
    {double result=0;for(std::size_t i=0;i<a.size();++i) result+=a[i]*b[i];return result;};
    const double beta=std::sqrt(dot(rhs,rhs));iterations=0;
    if(beta<=tolerance) return std::vector<double>(n,0);
    std::vector<std::vector<double>> basis(1,rhs);
    for(double& value:basis[0]) value/=beta;
    std::vector<double> H((limit+1)*limit,0),cs(limit),sn(limit),g(limit+1,0);g[0]=beta;
    int used=0;
    for(int j=0;j<limit;++j)
    {
        auto w=apply(basis[j]);
        if(w.size()!=rhs.size()) throw std::runtime_error("Invalid wall operator dimensions");
        for(double value:w) if(!std::isfinite(value)) throw std::runtime_error("Nonfinite joint wall response");
        // Twice modified Gram-Schmidt for near-dependent marker responses.
        for(int sweep=0;sweep<2;++sweep) for(int i=0;i<=j;++i)
        {
            const double d=dot(w,basis[i]);H[i*limit+j]+=d;
            for(int k=0;k<n;++k) w[k]-=d*basis[i][k];
        }
        const double length=std::sqrt(dot(w,w));H[(j+1)*limit+j]=length;
        if(length>1e-30) {for(double& x:w) x/=length;basis.push_back(w);}
        for(int i=0;i<j;++i)
        {
            const double a=H[i*limit+j],b=H[(i+1)*limit+j];
            H[i*limit+j]=cs[i]*a+sn[i]*b;H[(i+1)*limit+j]=-sn[i]*a+cs[i]*b;
        }
        const double a=H[j*limit+j],b=H[(j+1)*limit+j],scale=std::hypot(a,b);
        if(!std::isfinite(scale)||scale<1e-30) throw std::runtime_error("Joint wall Schur breakdown");
        cs[j]=a/scale;sn[j]=b/scale;H[j*limit+j]=scale;H[(j+1)*limit+j]=0;
        g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];used=j+1;
        if(std::abs(g[j+1])<=tolerance) break;
        if(length<=1e-30) throw std::runtime_error("Joint wall Krylov breakdown before convergence");
    }
    if(std::abs(g[used])>tolerance)
    {
        std::ostringstream message;message.precision(17);
        message<<"Joint wall iteration budget exhausted: residual="<<std::abs(g[used])
               <<" tolerance="<<tolerance<<" iterations="<<used;
        throw std::runtime_error(message.str());
    }
    std::vector<double> y(used),answer(n,0);
    for(int i=used-1;i>=0;--i)
    {
        double value=g[i];for(int j=i+1;j<used;++j) value-=H[i*limit+j]*y[j];
        y[i]=value/H[i*limit+i];
        for(int k=0;k<n;++k) answer[k]+=y[i]*basis[i][k];
    }
    const auto check=apply(answer);double residual=0;
    if(check.size()!=rhs.size()) throw std::runtime_error("Invalid explicit wall operator dimensions");
    for(int k=0;k<n;++k) residual+=(check[k]-rhs[k])*(check[k]-rhs[k]);
    if(!std::isfinite(residual)||std::sqrt(residual)>2*tolerance)
    {
        std::ostringstream message;message.precision(17);
        message<<"Joint wall explicit residual exceeds tolerance: residual="<<std::sqrt(residual)
               <<" limit="<<2*tolerance<<" iterations="<<used;
        throw std::runtime_error(message.str());
    }
    iterations=used;return answer;
}
}
#endif
