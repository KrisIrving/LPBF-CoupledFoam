// Actual shared C++ kernel regression, executed only by the user's test script.
#include "WallGMRES.H"
#include "JointResponseAudit.H"
#include "JointOuterMix.H"
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
    // Two feasible states of u0+u1=3 with residual A*u+G*p-S*F-b.
    // Same affine mix must carry p and the load F, not only velocity.
    const std::vector<double> previous{1,2,2,1},current{1.2,1.8,1.6,1.2};
    auto residual=[](const std::vector<double>& x)
    {return std::vector<double>{3*x[0]+x[1]+x[2]-x[3]-4,
                               x[0]+4*x[1]-x[2]-2*x[3]-3};};
    const auto rp=residual(previous),rc=residual(current);
    double dot=0,norm=0;
    for(int i=0;i<2;++i){const double d=rc[i]-rp[i];dot+=rp[i]*d;norm+=d*d;}
    const double factor=lpbfM2::jointOuterFactor(dot,norm);
    std::vector<double> mixed(4);
    for(int i=0;i<4;++i) mixed[i]=previous[i]+factor*(current[i]-previous[i]);
    const auto rm=residual(mixed);
    if(std::abs(mixed[0]+mixed[1]-3)>1e-12) return 8;
    for(int i=0;i<2;++i) if(std::abs(rm[i]-rp[i]-factor*(rc[i]-rp[i]))>1e-12) return 9;
    if(rm[0]*rm[0]+rm[1]*rm[1]>=rc[0]*rc[0]+rc[1]*rc[1]) return 11;
    if(factor<1||factor>4||lpbfM2::jointOuterFactor(-5,1)!=4
       ||lpbfM2::jointOuterFactor(1,1)!=1||lpbfM2::jointOuterFactor(1,0)!=1) return 10;
    std::cout<<"JOINT_OUTER_MIX_NATIVE_KERNEL_PASS\n";
    auto dotProduct=[](const std::vector<double>& a,const std::vector<double>& b)
    {double value=0;for(std::size_t i=0;i<a.size();++i) value+=a[i]*b[i];return value;};
    const std::vector<std::vector<double>> directions{{1,0},{1,0},{0,1}};
    const auto weights=lpbfM2::jointOuterWeights(directions,{1,1},dotProduct);
    double affine=0;for(double w:weights) affine+=w;
    std::vector<double> mixedResidual{weights.back(),weights.back()};
    for(int h=0;h<3;++h) for(int i=0;i<2;++i) mixedResidual[i]+=weights[h]*directions[h][i];
    if(std::abs(affine-1)>1e-12||std::hypot(mixedResidual[0],mixedResidual[1])>1e-12) return 12;
    const auto scaledWeights=lpbfM2::jointOuterWeights(directions,{1,1},
        [&](const auto& a,const auto& b){return 7*dotProduct(a,b);});
    for(int i=0;i<4;++i) if(std::abs(weights[i]-scaledWeights[i])>1e-12) return 13;
    const auto unsafe=lpbfM2::jointOuterWeights({{.9,0}},{1,0},dotProduct);
    if(unsafe[0]!=0||unsafe[1]!=1) return 14;
    std::cout<<"JOINT_OUTER_SUBSPACE_NATIVE_KERNEL_PASS\n";
    std::cout<<"JOINT_RESPONSE_AUDIT_NATIVE_KERNEL_PASS\n";
    std::cout<<"WALL_GMRES_NATIVE_KERNEL_PASS\n";
}
