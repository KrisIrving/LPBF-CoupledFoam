// Actual shared C++ kernel regression, executed only by the user's test script.
#include "WallGMRES.H"
#include "JointResponseAudit.H"
#include <iostream>
int main()
{
    auto apply=[](const std::vector<double>& x)
    {return std::vector<double>{3*x[0]+2*x[1],x[0]+4*x[1]+x[2],2*x[1]+5*x[2]};};
    const std::vector<double> exact{1,-2,3};int iterations=0;
    const auto answer=lpbfM2::wallGMRES(apply,apply(exact),1e-12,3,iterations);
    for(int i=0;i<3;++i) if(std::abs(answer[i]-exact[i])>1e-11) return 1;
    if(iterations!=3) return 2;
    bool exhausted=false;
    try {lpbfM2::wallGMRES(apply,apply(exact),1e-12,1,iterations);}
    catch(const std::runtime_error&) {exhausted=true;}
    if(!exhausted) return 3;
    bool singular=false;
    try {lpbfM2::wallGMRES([](const std::vector<double>& x){return std::vector<double>(x.size(),0);},exact,1e-12,3,iterations);}
    catch(const std::runtime_error&) {singular=true;}
    if(!singular) return 4;
    auto linear=[](const std::vector<double>& x)
    {
        lpbfM2::JointResponseSample s;
        s.wall=x;s.velocity=x;s.pressure=x;s.flux=x;
        return s;
    };
    const auto audit=lpbfM2::auditJointResponse(linear,9);
    if(audit.zero!=0||audit.repeat!=0||audit.linear>1e-12||audit.scale>1e-12) return 5;
    auto nonlinear=[&](const std::vector<double>& x)
    {auto s=linear(x);s.velocity[0]+=100*x[0]*x[0];return s;};
    if(lpbfM2::auditJointResponse(nonlinear,9).linear<1e-4) return 6;
    int calls=0;
    auto stale=[&](const std::vector<double>& x)
    {auto s=linear(x);s.flux[0]+=1e-3*(++calls);return s;};
    const auto staleAudit=lpbfM2::auditJointResponse(stale,9);
    if(staleAudit.zero<.1||staleAudit.repeat<.1) return 7;
    std::cout<<"JOINT_RESPONSE_AUDIT_NATIVE_KERNEL_PASS\n";
    std::cout<<"WALL_GMRES_NATIVE_KERNEL_PASS\n";
}
