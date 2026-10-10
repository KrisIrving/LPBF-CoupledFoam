// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LPBF_JOINT_RESPONSE_AUDIT_H
#define LPBF_JOINT_RESPONSE_AUDIT_H
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>
namespace lpbfM2
{
struct JointResponseSample
{
    std::vector<double> wall,velocity,pressure,flux;
};
struct JointResponseAudit
{
    double zero=0,repeat=0,linear=0,scale=0;
};
inline double responseDifference(const JointResponseSample& actual,
    const JointResponseSample& a,double ca,const JointResponseSample& b,double cb)
{
    double error=0;
    const std::vector<double>* x[]={&actual.wall,&actual.velocity,&actual.pressure,&actual.flux};
    const std::vector<double>* y[]={&a.wall,&a.velocity,&a.pressure,&a.flux};
    const std::vector<double>* z[]={&b.wall,&b.velocity,&b.pressure,&b.flux};
    // Separate relative scales: pressure must not hide a velocity/flux error.
    for(int field=0;field<4;++field)
    {
        if(x[field]->size()!=y[field]->size()||x[field]->size()!=z[field]->size())
            throw std::runtime_error("Joint response audit dimensions changed");
        double difference=0,denominator=1e-30;
        for(std::size_t i=0;i<x[field]->size();++i)
        {
            const double value=(*x[field])[i],ya=ca*(*y[field])[i],zb=cb*(*z[field])[i];
            if(!std::isfinite(value)||!std::isfinite(ya)||!std::isfinite(zb))
                throw std::runtime_error("Nonfinite joint operator audit sample");
            difference=std::max(difference,std::abs(value-ya-zb));
            denominator=std::max(denominator,std::max(std::abs(value),std::abs(ya)+std::abs(zb)));
        }
        error=std::max(error,difference/denominator);
    }
    return error;
}
inline JointResponseAudit auditJointResponse(
    const std::function<JointResponseSample(const std::vector<double>&)>& apply,int size)
{
    if(size<1) throw std::runtime_error("Empty joint response audit");
    std::vector<double> zero(size,0),a(size),b(size),combination(size),half(size);
    for(int i=0;i<size;++i)
    {
        a[i]=1e-3*std::sin(.73*(i+1));b[i]=1e-3*std::cos(.41*(i+1));
        combination[i]=a[i]-.3*b[i];half[i]=.5*a[i];
    }
    const auto z=apply(zero),x=apply(a),y=apply(b);
    const auto combined=apply(combination),scaled=apply(half),again=apply(a);
    JointResponseAudit result;
    result.zero=responseDifference(z,x,0,y,0);
    result.repeat=responseDifference(again,x,1,y,0);
    result.linear=responseDifference(combined,x,1,y,-.3);
    result.scale=responseDifference(scaled,x,.5,y,0);
    return result;
}
}
#endif
