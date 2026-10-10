# M2A-02C 原生压力—壁面联合校正候选

2026-10-10。已实现原生候选及集中用户测试入口，**原生编译、共享C++内核运行及CFD运行均待用户反馈**。37项开发端离线检查通过，六组独立几何算子审计通过，旧二次反馈摘要复算完全一致，Windows Git Bash仅检查脚本语法。没有调用开发端WSL或执行CFD/DEM。M2A仍未关闭。

## 本次交付与方案边界

`constraintScheme jointSurface`限静止单球的fixed/rotate，启用`compatibleGauss`外边界。新方案取消球内ghost/体积惩罚，lambda恒零；真实表面节点与流体通过三线性J交换速度，通过同一权重的转置交换**节点积分力F（N）**。它与上一契约的牵引τ等价，F=Wτ；流体加速度`SaF=JᵀF/(ρV)`，没有漏除ρ或重复乘面积。q4内部虚拟流体Pi/Li库存保留。

节点采用固定Fibonacci球面序列，`N=ceil(4πR²/(2.25h²))`，三网格分别51、90、159个，各节点间距量级1.5h。这些约束节点不同于原12×24诊断积分节点。六种网格/相位的归一化`JJᵀ`Cholesky最小主元0.9064–0.9410，大于预设1e−8；没有加入对角正则化。**这仅验证J行独立，不证明压力投影后的Schur算子良态。** 后者通过GMRES预算、breakdown及显式残差检查保留失败。

由分片单位性和坐标线性再现可得节点/网格合力与一阶矩一致；伴随关系给出功交换一致。原生每窗同时审计节点合力、力矩、功与实际网格源，门槛分别1e−12 N、1e−14 N·m、1e−14 W，超标不得推进DEM。

## 原生求解关系

同一耦合窗口冻结旧时间场及对流phi。完整动量预测后，沿用原PISO，用当前已累计源构造HbyA和实际phiH（含ddtCorr），进行压力校正。若真实节点壁速残差超过1e−7 m/s，构造如下压力投影响应：

1. 节点积分力增量δF，经`δa=JᵀδF/(ρV)`扩散。
2. `δUH=rAU δa`，外边界速度增量为零，处理器边界正常交换；`δphiH=flux(δUH)`。
3. 用**同一rAU**解`laplacian(rAU,δp)=div(δphiH)`，压力参考增量为零，fixedFluxPressure通过constrainPressure更新。显式登记jointDeltaP的fluxRequired，并使用原p线性求解控制。
4. `δphi=δphiH−pressureMatrix.flux()`，`δU=δUH−rAU grad(δp)`。以`JδU`作为壁面Schur响应。
5. 有界GMRES求`JδU=Ub−JU`；右缩放用未投影`J rAU Sa`对角。实际FV响应不假定对称，未使用CG。两次修正Gram–Schmidt，最多64次，不靠增加预算掩盖未收束。Krylov估计之后重新计算完整算子残差，内目标为1e−8 m/s的L2范数量级。
6. 同时更新U、p、phi、累计网格源和动量矩阵source；再次实际采样节点速度并验证最大矢量残差≤1e−7 m/s。

这不是一次性对完整动量矩阵求逆的精确块分解：响应采用rAU对角近似，动量非对角部分通过后续PISO/H更新迭代。因此仍要求完整动量方程冲量L1≤1e−11 kg·m/s、实际phi散度≤1e−7 /s、至少8最多96次压力校正，全部满足后才允许DEM回传。近似响应能否在原预算内收束是本次原生测试的重要待验证项。

## 同源载荷、重启与诊断

受力账本读取最后已求解的累计网格源：`C=Σρ a V`，颗粒力`Fh=−C+(Pi_new−Pi_old)/dt`，力矩同样使用源的一阶矩及内部虚拟流体惯性修正。节点源不再独立重复施力；压力/黏性表面应力依旧只作诊断。

当前每窗从零节点源重新求解，未采用上一窗力热启动；重启同样从零求源，所以不会缺失未保存的源历史。检查点核对constraintScheme、boundaryTreatment、节点数和既有ID/内部库存/时钟。共同重启仍需原生全历史对照确认，不因这项设计就宣布通过。

CSV新增joint_iterations、joint_pressure_solves、joint_wall_residual、joint_markers、joint_min_pivot以及力/力矩/功交换误差。普通pressure_correctors保持原PISO计数；额外压力响应次数单独记录，避免隐藏成本。`surface_target_defect`在本方案记实际节点最大残差；原`slip_rms`仍是**独立12×24真实球面插值RMS**。节点约束成功不等于整个球面无滑移或受力精度通过，分析器同时检查两者。

## 一次集中原生反馈

入口`scripts/m2a-joint-test.sh`重用既有OF/LIGGGHTS依赖，不重装、不重跑通信demo。先运行六组几何预检，再在用户机器用g++编译并执行**与求解器共享的WallGMRES.H**内核测试，覆盖非对称系统、预算耗尽和奇异响应。通过后wmake机械求解器，再运行原三网格×两相位的固定/转动球与固定球MPI、20+20窗共同重启，共14组。执行失败继续收集其余轨迹，原生预检/构建失败则直接保留归档。

原40窗、三网格/两相位、球面滑移、最细受力/应力10%、相位2%、细化增量0.02、线动量及串并行/共同重启门槛均保留。新增节点/交换/迭代计数和原生边界兼容性检查；不因离线应力偏差降低门槛。单次求解资源预算仍900秒，可配置JOINT_CASE_TIMEOUT，但当前不建议先延长。超时保留exit_status=124，不当作完整验收。

```bash
cd ~/LPBF-CoupledFoam
git pull --ff-only
export OPENFOAM_BASHRC="$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
source .build/demo-liggghts.env
bash scripts/m2a-joint-test.sh
```

归档前缀为m2a-joint。完整U/p/phi/网格仍保留在用户.runs，反馈保留日志、控制、CSV及共同检查点。若失败，下一步优先完整核对新算子的编译/API、压力响应残差、实际phi和动量缺陷、球面滑移、反力/应力及成本；不只挑通过组或反复单项微调。

## 当前限制

均匀轴对齐立方网格、恒ρ、等温单相、静止球心、单个无接触球；全局Cartesian场及节点复制仅用于验证，尚非生产分布式映射。额外压力响应可能昂贵，native成本未测量。没有声称完成动态GCL、完整流体角动量、多球接触、热/蒸发或LPBF物理验证。

主文件：applications/solvers/resolvedParticleFoam/JointSphere.H、WallGMRES.H、resolvedParticleFoam.C；测试入口与分析器位于scripts/prepare-m2a-joint.py、m2a-joint-test.sh、summarize-m2a-joint.py。上阶段[公式契约](m2a-pressure-wall-contract.zh-CN.md)保留为设计历史，由本文件覆盖“原生联合算子尚未实现”的旧状态。
