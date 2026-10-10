// Actual shared C++ kernel regression, executed only by the user's test script.
#include "WallGMRES.H"
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
    std::cout<<"WALL_GMRES_NATIVE_KERNEL_PASS\n";
}
