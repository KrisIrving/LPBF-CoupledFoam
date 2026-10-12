// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LPBF_POLYNOMIAL_STENCIL_H
#define LPBF_POLYNOMIAL_STENCIL_H
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace lpbfM2
{
inline std::array<double,20> cubicBasis(const std::array<double,3>& p)
{
    const double x=p[0],y=p[1],z=p[2];
    return {1,x,y,z,x*x,x*y,x*z,y*y,y*z,z*z,
        x*x*x,x*x*y,x*x*z,x*y*y,x*y*z,x*z*z,y*y*y,y*y*z,y*z*z,z*z*z};
}
struct PolynomialWeights {std::vector<double> weights;double moment=0,l1=0;};
inline PolynomialWeights cubicWeights(const std::vector<std::array<double,3>>& points)
{
    const int n=points.size();
    std::vector<std::array<double,20>> basis;std::vector<double> d;
    for(const auto& p:points)
    {basis.push_back(cubicBasis(p));d.push_back(std::exp(-(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2));}
    std::vector<std::vector<double>> q;
    double r[20][20]{};
    for(int j=0;j<20;++j)
    {
        std::vector<double> v(n);double original=0;
        for(int k=0;k<n;++k){v[k]=std::sqrt(d[k])*basis[k][j];original+=v[k]*v[k];}
        for(int pass=0;pass<2;++pass) for(int i=0;i<j;++i)
        {
            double projection=0;for(int k=0;k<n;++k) projection+=q[i][k]*v[k];
            r[i][j]+=projection;for(int k=0;k<n;++k) v[k]-=projection*q[i][k];
        }
        double norm=0;for(double x:v) norm+=x*x;r[j][j]=std::sqrt(norm);
        if(!std::isfinite(r[j][j])||r[j][j]<=std::max(std::sqrt(original),1.)*1e-10)
            throw std::runtime_error("Deficient fluid cubic stencil");
        for(double& x:v) x/=r[j][j];q.push_back(v);
    }
    double z[20]{};
    for(int i=0;i<20;++i)
    {double rhs=i==0 ? 1:0;for(int j=0;j<i;++j) rhs-=r[j][i]*z[j];z[i]=rhs/r[i][i];}
    PolynomialWeights result;result.weights.resize(n);
    for(int k=0;k<n;++k)
    {
        double value=0;for(int j=0;j<20;++j) value+=q[j][k]*z[j];
        result.weights[k]=std::sqrt(d[k])*value;result.l1+=std::abs(result.weights[k]);
    }
    for(int j=0;j<20;++j)
    {double value=0;for(int k=0;k<n;++k) value+=result.weights[k]*basis[k][j];
     result.moment=std::max(result.moment,std::abs(value-(j==0 ? 1:0)));}
    if(!std::isfinite(result.l1)||result.l1>8||!std::isfinite(result.moment)||result.moment>1e-10)
        throw std::runtime_error("Fluid cubic moment/coefficient gate failed");
    return result;
}
}
#endif
