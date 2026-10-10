# M2A 原生联合算子诊断包

2026-10-11。按文献/方法复审实施，生成器revision4；原生编译、共享C++诊断内核和实际CFD诊断均待用户运行。开发端没有执行WSL或仿真。45项M2A离线检查通过；首轮joint及revision2失败反馈摘要复算完全一致。

## 实际调用路径

jointOperatorAudit=true时，在每次求解器启动后的首次壁面校正中，使用本次窗口同一个rAU、网格、齐次增量边界和压力参考，直接调用生产response(input,false)。每次探针清除lastInput，绕过同输入结果缓存，实际执行压力Poisson求解，而非用几何或小型矩阵替代OpenFOAM算子。审查不推进DEM、不改变节点力种子；完整矩阵作用检查暂时扰动当前U并恢复，不修改旧时间、压力、源或几何。

共享JointResponseAudit.H定义6次探针：零输入、确定性a/b、a−0.3b、0.5a及隔次重新计算a。快照包含节点壁速、内部U、内部p、全部内部/边界phi。重复/叠加/缩放以各字段分别归一化的最大差衡量，并对所有MPI rank取最大值，避免压力尺度掩盖速度/通量问题。

随后增加一次径向缩放输入探针，记录normalGain；再用实际response(a)生成右端，并通过原WallGMRES求解，显式检查结果响应。制造解检查使用原64次预算及1e−8m/s求解容差，显式残差上限2e−8m/s。径向探针的缩放输入不是均匀物理压力牵引；它只能提供方向性响应信息。一个制造右端可解不证明完整Schur全秩、唯一解或良好条件数，不会自动削除法向模式/添加正则项。

每次独立启动审查一次，因此共同重启的第一段和重启段都必须提供各自的审查证据。原14组每组有自己的真实网格/相位/rank响应，包含静止与规定转动的实际窗口状态。

## 记录和门槛

| 字段 | 定义/门槛 |
|---|---|
| zero/repeat/linear/scale | 四字段最大相对差，上限1e−8；零输入期望为严格零，任何非零响应会触发该检查 |
| div | 各探针实际div(deltaPhi)最大值，原连续性上限1e−7/s |
| gauge | rank0参考单元增量压力，除max(全局压力幅值,1)，上限1e−8；固定参考值为0 |
| diagonal | deltaU+rAU·grad(deltaP)−rAU·deltaAcceleration的相对最大缺陷，上限1e−8 |
| fullImpulse | 通过实际动量矩阵残差差分获得A·deltaU+V·grad(deltaP)−V·deltaAcceleration，再计算rho·dt·L1；仅诊断，不要求对角近似同时满足完整A逆 |
| normalGain | 缩放径向探针的壁面响应增益，仅诊断，不作秩门槛 |
| manufactured | 制造右端的实际壁面响应L2残差，上限2e−8m/s，保留迭代数 |
| pressureSolves/seconds | 审查额外压力求解数/墙钟秒数，独立记录 |

首次真实response(...,true)另记M2A_OPERATOR_COMMIT：检查速度、压力、内部phi、加速度及equation.source增量是否匹配已显式核验的deltaU/deltaP/deltaPhi/deltaAcceleration。误差分别相对于该字段基态与增量的最大幅值归一化，均不超过1e−8，避免极小增量减大基态的舍入噪声被放大；这不是极小增量的相对精度证明。它核验提交代码一致性，不代表完整动量方程收敛；后者仍由原有动量门槛独立检查。提交的边界phi由整体场加法更新，重复/叠加探针包含边界phi。

M2A_OPERATOR_PROBE在制造求解之前输出，基本探针失败立即停止该工况；M2A_OPERATOR_AUDIT为完整诊断结果。原生摘要要求全部字段有限、门槛成立以及提交证据，缺少诊断不能当作通过。失败工况仍不进入物理/比较通过，并在摘要operator_logs保留探针、完整诊断、提交、最后校正和错误。

审查额外压力求解计数从历史CSV joint_pressure_solves中分开，避免改变物理校正预算；窗口wall_seconds包含审查成本，应结合独立seconds解释启动开销。没有修改压力容差、40窗、物理受力/应力门槛或900秒单次预算。

## 明显发散保护

jointFullPredictor=true在当前对角响应模式下被求解器直接拒绝，防止再次启用已失稳组合。若校正后maxU超过参考壁速/外流速度的1e6倍，则提前失败；固定球参考0.01m/s，转动球参考max(|Vp|+R|Omega|,1e−12)。该保护只缩短明显异常轨迹，仍算失败，不替代收敛或物理门槛。低于保护阈值的失稳仍受原校正次数和900秒限制。GMRES失败信息新增真实残差、容差/上限和迭代数，不改变算法。

## 用户入口与结果解释

```bash
cd ~/LPBF-CoupledFoam
git pull --ff-only
export OPENFOAM_BASHRC="$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
source .build/demo-liggghts.env
BUILD_JOBS=8 JOINT_CASE_JOBS=4 bash scripts/m2a-joint-test.sh
```

现有入口先编译运行共享WallGMRES及JointResponseAudit内核；后者含线性通过、非线性拒绝、带状态/偏移响应拒绝。然后编译实际求解器，运行含新诊断的原14组。4组并发策略和资源上限保留，不重装LIGGGHTS，不重复通信demo。返回新的m2a-joint归档即可。

若探针失败，先根据字段定位边界/缓存/压力响应，不继续扩大迭代；若探针通过而完整动量缺陷不收束，说明对角联合响应/外层迭代仍不足，应重建一致压力—节点力块及预条件，而非将探针通过称为CFD验收。若轨迹稳定，再审查最细受力/应力和细化、MPI/共同重启。所有通过均限定于本静止、等温、无接触单球包，不扩展为自由运动、GCL、热、蒸发或LPBF验证。
