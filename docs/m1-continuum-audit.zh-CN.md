# M1 连续相源码审计与下一轮检查

## 状态和范围

物理源码基线为 LaserbeamFoam V3.0，提交 `f42e08a0bfc1749675beadcf7c6a90334d7aa9f8`；项目的 `6f23a32b` 尚未更改原物理方程。本文基于公开源码静态阅读。新增诊断的编译和运行等待用户反馈，M1 的物理验证尚未完成。

本地研究资料的逐项对照另存于工作区 `output`，不上传原始书稿或该对照报告。本文件只记录公开源码分析、项目设计和检查任务。

2026-10-08 用户反馈已确认 `66dc771d` 的诊断代码在 v2512 编译通过，基线串行/双进程 MPI 正常运行，各产生 21 个样本。M1.1 诊断接入通过；这不覆盖可压缩分支运行或熔化/蒸发物理验证。具体数值见开发记录。初始激光功率来自上游场初始化值，首次激光更新前不能作物理解读。

## 主干选择建议

继续使用 V3.0 作为基础，保留 `laserbeamFoam` 为已通过短时运行检查的参考分支。对显式蒸发—气体连续相主干，优先评估 `compressibleLaserbeamFoam`：它已有相变体积分数、密度和压力耦合接口，但必须先通过源项配对、能量与气体响应检查。

在这些检查之前，主干状态为**候选**，不开始固定 DEM 适配。源码中有可压缩方程不等于已经验证了实际 LPBF 蒸汽羽流。源项驱动的低 Mach 气体路线保留为后续可选模式；需统一蒸发接口及验证指标后再实现，不能只把现有两相不可压缩程序改个名称。

## 不可压缩分支的方程与源项

| 项目 | 实际实现 | 审计解释 |
|---|---|---|
| 相分数与密度 | [createFields.H](../applications/solvers/laserbeamFoam/createFields.H)，[MULES 源项](../applications/solvers/laserbeamFoam/MULES/alphaSuSp.H)，[isoAdvector 源项](../applications/solvers/laserbeamFoam/isoAdvector/alphaSuSp.H) | 两个界面输运路径的 `Su/Sp` 为零；密度由金属/气体相分数和固定相密度构成，未将 `Qv/Lv` 作为显式蒸发质量源加入 |
| 动量 | [UEqn.H](../applications/solvers/laserbeamFoam/UEqn.H) | 变量密度动量、黏性、Darcy 阻力、Marangoni 项；`pVap` 经界面梯度作用于动量预测 |
| 压力 | [pEqn.H](../applications/solvers/laserbeamFoam/pEqn.H) | 反冲项也参与压力预测通量，之后压力修正速度；这是同一个力在压力—速度算法中的一致使用，不能仅因出现两次就认定重复施力 |
| 热与相变 | [TEqn.H](../applications/solvers/laserbeamFoam/TEqn.H) | 温度方程含激光、蒸发冷却与熔化潜热；通过 `epsilon1` 迭代处理熔化/凝固 |
| 物性 | [updateKappaCp.H](../applications/solvers/laserbeamFoam/updateKappaCp.H) | 对 `cp` 和 `kappa` 进行体积分数混合，再形成 `rhoCp=rho*cp`；需要审计界面单元热容的物理定义 |
| 网格 | [isoAdvector/firstIter.H](../applications/solvers/laserbeamFoam/isoAdvector/firstIter.H)，[MULES/firstIter.H](../applications/solvers/laserbeamFoam/MULES/firstIter.H) | 主循环通过包含文件调用 `mesh.update()`；isoAdvector 路径还包含重构和映射，动态网格功能应继续单独检查 |

蒸发项可整理为：

\[
p_{sat}(T)=p_0\exp\left[\frac{L_v M_m}{R}\left(\frac1{T_v}-\frac1T\right)\right],\quad
p_{Vap}=0.54p_{sat},\quad
Q_v=0.82L_vp_{sat}\sqrt{\frac{M_m}{2\pi RT}}.
\]

`pVap` 为 Pa；`Qv` 为 W/m²，乘 `mag(gradAlpha)` 与无量纲 `thermalDamper` 后进入温度方程，成为 W/m³。它提供蒸发冷却，但没有对应的显式蒸汽相质量生成、物种输运或指定出口蒸汽动量通量接口。因此，当前参考分支无法直接作为已闭合的蒸汽—动态颗粒模型。

温度方程内部未见显式辐射或给定 `h(T-Ta)` 源项；通过气体能量方程和边界条件传热的途径另需核查，不与给定换热系数模型自动等同。

## 可压缩分支的方程与源项

实际顺序见 [主循环](../applications/solvers/compressibleLaserbeamFoam/compressibleLaserbeamFoam.C)：

1. 网格更新；调用 `mixture.solve(&mass_dot)` 更新相分数和相变项。
2. 用相分数通量形成 `mixture.rhoPhi()`，求解密度输运。
3. 更新物性、光学界面和激光沉积；求解动量、温度。
4. 压力修正并更新各相热力学密度，最终令 `rho=mixture.rho()`。

| 量 | 单位 | 去向与意义 |
|---|---|---|
| `evaprate/condrate` | s⁻¹ | 由饱和压力、当前压力、温度、密度及 `interface_thickness` 构造的相变速率系数 |
| 每相 `alphagen` | s⁻¹ | 进入相分数输运；不是直接的 kg/(m³ s) 质量源 |
| `PCR/vDot` | s⁻¹ | 由相源项求和构成，进入 [压力方程](../applications/solvers/compressibleLaserbeamFoam/pEqn.H) |
| `massgen/mass_dot` | kg K/(m³ s) | 经过潜热和各相 `Cv` 换算，进入 [温度方程](../applications/solvers/compressibleLaserbeamFoam/TEqn.H)；这个变量名容易被误解 |
| `rhoPhi` | kg/s，每个面 | 各相质量通量求和，用于密度、动量和热输运 |

温度方程还包含压力功/动能相关项。动量与压力源码中未见不可压缩分支的 `pVap` 经验反冲项，也没有指定蒸汽法向射流速度的显式源接口。压力驱动流动能否达到目标，需要验证；不应直接复制另一个分支的反冲项。

## 需要先解决的风险

### R1：蒸发系数与实际转移质量（纠正凝结判断）

M1.2 核查实际源码发现：上一版将凝结密度比转录反了。实际液相源使用 rho_v/rho_l，撤回此前的凝结质量不配对判断。

固定密度、相同限幅速率下，蒸发源为 S_l=-e*alpha_l*rho_v/rho_l、S_v=e*alpha_l；凝结源为 S_l=c*alpha_v*rho_v/rho_l、S_v=-c*alpha_v。两者均满足 rho_l*S_l+rho_v*S_v=0，体积源分别膨胀与收缩。

待审计的是蒸发系数含 1/(interface_thickness*rho_l)，而实际质量率为 mdot=rho_v*S_v。相对以液相密度解释的供体速率，多一个 rho_v/rho_l 因子。需要追溯物理定义、界面分布和限幅，暂不修改生产求解器。

新增 scripts/check-phase-pair.py 直接提取实际 alphagen/massgen 表达式，执行密度比 2/1000、蒸发/凝结、限幅开/关的 8 个合成检查。报告包含源码哈希、表达式、质量和体积源、温度源及潜热参考功率。温度源单位与潜热功率不同，不能直接比较。此检查不包含输运、PCR、EOS 或能量推进，不证明全局 CFD 守恒。

```bash
python3 scripts/check-phase-pair.py --output .runs/m1-phase-pair.json
```

此命令只需要 Python 3，无需加载 OpenFOAM 或运行仿真。

### R2：子循环的源项时间平均

`solve()` 在每个子步调用 `solveAlphas()`；后者重置 `massdotterm`。外层对 `rhoPhi` 有 `dt_sub/dt_total` 加权，但返回的 `PCR` 被每次赋值覆盖，温度项 `massdotterm` 也只保留最后一个子步的结果。

在源项随子步改变时，这三者的时间代表值不同。需要用子步累计的相变质量、体积变化和能量交换统一外层源项，并检查 `nAlphaSubCycles=1/2/4`。此时还有随时间步变化的速率限幅 `maxrate=0.5/dt`，应同时记录限幅激活程度。

### R3：界面单元的能量定义

不可压缩分支计算 `rho=(alpha*rho_l+(1-alpha)*rho_g)`，`cp=(alpha*cp_l+(1-alpha)*cp_g)`，然后相乘。这与 `alpha*rho_l*cp_l+(1-alpha)*rho_g*cp_g` 的相热容加权结果不同。差值为：

\[
\alpha(1-\alpha)(\rho_l-\rho_g)(c_{p,g}-c_{p,l}).
\]

需要结合界面混合、`thermalDamper` 和相变焓选择一致的能量定义。变比热下，`rhoCp*T` 不能直接当作完整显热，更不能当作总能量。可压缩分支的 `Cv/rCv/mass_dot` 也要从同一能量方程重新核对。

### R4：网格、光线及命名依赖

- 速率使用给定 `interface_thickness`，需要检查它与真实界面厚度、网格和界面压缩的关系。
- 激光的高斯核为 `f*P/(pi*r²)*exp(-f*R²/r²)`，默认 `Radius_Flavour=2`。比较其他模型前要统一高斯宽度定义。
- 光线初始化截断于 `1.5*r`，径向面积离散、初始化功率和截止剩余功率都需要能量账本。`integral(deposition)` 只说明沉积功率，不说明输入功率全部去向。
- 气相通过名称包含 `vapour` 或等于 `air` 判定；液气配对依赖 `liquidName+vapour`。应把这种隐式命名约定整理为显式相类型/相对声明，避免新增保护气或材料时误分类。
- [可压缩 Test1](../tutorials/compressiblelaserbeamFoam/Test1/system/controlDict) 仍有旧应用名，且为多材料配置；不能直接当作 316L 文献验证算例。

## 本轮新增的可关闭诊断

`controlDict` 中 `continuumDiagnostics true;` 启用诊断，默认关闭。两套主循环在初始时刻及每个完成的时间步输出 `M1_CONTINUUM`：

- `massKg=integral(rho dV)`，`outwardMassFluxKgPerS=sum(rhoPhi)`（排除处理器和周期内部界面）。
- 动能、吸收激光功率、温度范围、最大速度。
- 所有全局积分经 MPI 归约后由主进程打印。

诊断不改变原方程。本轮不实现蒸发源项修复，也不计算完整焓、边界热通量、界面能量交换或机械功的收支。

[摘要脚本](../scripts/summarize-continuum.py) 生成 `m1-summary.json`：记录初末库存、库存差、右端点估计的质量残差和吸收能量积分，并列出串行/并行最终总质量差。该残差不是任意时间格式、孔隙率、网格变化或外部质量源下的一般离散残差，不能自动宣称物理验证通过。积分一致也不证明场级并行一致。

## M1 后续测试顺序和验收

| 顺序 | 测试 | 需要证明的内容 |
|---|---|---|
| M1.1 | 原小算例加诊断，串行/2 进程 | 新诊断可编译运行，输出归约正确，21 个初始/步末观测可解析；不是物理验证 |
| M1.2 | 无激光、封闭单液气对的蒸发/凝结检查 | 相间质量转移配对、体积变化符号和潜热对应关系；检查密度比和限幅 |
| M1.3 | 同条件子循环 1/2/4 | 累计交换量与外层时间推进一致，误差随步长收敛 |
| M1.4 | 恒定/变比热导热、固液相变、激光平面吸收 | 定义相焓，核对显热、潜热、边界热通量和激光输入/沉积/逸出 |
| M1.5 | 无颗粒的界面蒸发与气体响应、裸板单道 | 检查压力、法向质量/动量/能量通量、气体速度和熔池几何 |

完成这些证据后才能批准连续相主干进入 M2。此次提交只完成静态审计与 M1.1 诊断准备。

下一轮沿用 `bash scripts/remote-test.sh`，报告包会新增 `m1-summary.json`。用户运行前需更新并重新编译，不能复用上一轮二进制来验证新诊断。Python 摘要脚本使用标准库，需 `python3`。
