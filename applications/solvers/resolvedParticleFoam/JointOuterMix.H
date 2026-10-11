// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LPBF_JOINT_OUTER_MIX_H
#define LPBF_JOINT_OUTER_MIX_H
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <stdexcept>
#include <utility>
namespace lpbfM2
{
// Minimise ||R_previous + alpha (R_current-R_previous)||_2. The
// caller supplies global dot products and must recheck the actual candidate.
inline double jointOuterFactor(double previousDotDifference,double differenceSquared)
{
    if(!std::isfinite(previousDotDifference)||!std::isfinite(differenceSquared)
       ||differenceSquared<=0) return 1;
    const double alpha=-previousDotDifference/differenceSquared;
    return std::isfinite(alpha) ? std::clamp(alpha,1.0,4.0):1;
}
// Small residual subspace, QR with twice modified Gram-Schmidt. The global dot
// callback makes coefficients identical on all MPI ranks; local vector sizes
// may differ. Return baseline if directions/coefficient growth are unsafe.
inline std::vector<double> jointOuterWeights(
    const std::vector<std::vector<double>>& history,const std::vector<double>& current,
    const std::function<double(const std::vector<double>&,const std::vector<double>&)>& dot)
{
    const int m=history.size();
    if(m<1||m>4) throw std::runtime_error("Outer residual history outside1..4");
    std::vector<double> weights(m+1,0);weights.back()=1;
    const double reference=std::sqrt(dot(current,current));
    if(!std::isfinite(reference)||reference<=0) return weights;
    std::vector<std::vector<double>> basis;
    std::vector<int> selected;
    std::vector<double> upper(m*m,0),rhs;
    // Prefer recent directions if older directions are nearly dependent.
    for(int h=m-1;h>=0;--h)
    {
        if(history[h].size()!=current.size()) throw std::runtime_error("Outer residual size mismatch");
        std::vector<double> v=current;
        for(std::size_t i=0;i<v.size();++i) v[i]=history[h][i]-current[i];
        const int j=basis.size();
        for(int k=0;k<=j;++k) upper[k*m+j]=0;
        for(int sweep=0;sweep<2;++sweep) for(int k=0;k<j;++k)
        {
            const double projection=dot(basis[k],v);upper[k*m+j]+=projection;
            for(std::size_t i=0;i<v.size();++i) v[i]-=projection*basis[k][i];
        }
        const double norm=std::sqrt(dot(v,v));
        if(!std::isfinite(norm)||norm<=reference*1e-10) continue;
        upper[j*m+j]=norm;
        for(double& value:v) value/=norm;
        rhs.push_back(-dot(v,current));basis.push_back(std::move(v));selected.push_back(h);
    }
    const int n=basis.size();std::vector<double> gamma(n,0);
    for(int j=n-1;j>=0;--j)
    {
        double value=rhs[j];for(int k=j+1;k<n;++k) value-=upper[j*m+k]*gamma[k];
        gamma[j]=value/upper[j*m+j];
        if(!std::isfinite(gamma[j])) return weights;
    }
    double sum=0;
    for(int j=0;j<n;++j){weights[selected[j]]=gamma[j];sum+=gamma[j];}
    weights.back()=1-sum;
    double norm=0;for(double value:weights) norm+=std::abs(value);
    if(!std::isfinite(norm)||norm>8)
    {std::fill(weights.begin(),weights.end(),0);weights.back()=1;}
    return weights;
}
}
#endif
