# 平台选择与解析 CFD–DEM 开发计划（待讨论确认）

## 用户确认后的执行路线（覆盖下文待定建议）

2026-10-09 用户已同意保留 OpenCFD v2512，定向移植 CFDEM 必要算法，先做审查与 OpenFOAM–LIGGGHTS 通信 demo。不再进行 OF10 平台比较。P0 平台方向已定，P1 改为 v2512 最小嵌入通信包，具体范围、依赖、固定门槛及用户运行操作见 [demo 文档](dem-communication-demo.zh-CN.md)。

已完成当前公开 C API/字段/共享库需求审查并交付 demo；真实编译及运行待反馈。解析 IB 的完整离散、低 Mach 闭合和部分熔化仍按后续模型包核实，不因本次同意而声称所有模型细节已经确定。下文 OF10 优先评估及暂停实现是此前讨论稿，保留决策历史；当前以本节和 demo 文档为准。

2026-10-09。用户要求先详细讨论并敲定计划，再动手。本页是建议方案，未切换平台、未移植 CFDEM，也未继续修改求解器。现有 v2512 工作线冻结在已测试提交 `98befc42ae4187d9217fd4793ef427d3fcee58b6`，历史测试继续保留。

## 结论及证据边界

如果准备较多复用 CFDEM-PUBLIC 的解析浸入边界实现，Foundation OF10 是优先评估候选；若仅借用通信思想、自己实现边界约束，则当前已工作的 OpenCFD v2512 更有连续性优势。OF10 的优势是可能减少同一 Foundation 系列内的 API 移植差异，这是工程推断，尚未通过编译量化。不能由版本名称推出完整耦合兼容。

| 组合 | 核实的事实 | 对选择的含义 |
|---|---|---|
| 本项目 V3.0 / v2512 | 用户已编译并完成多包连续相测试，仍有物理和守恒缺口 | 保留为工程基线，测试并非通用物理验证 |
| LaserbeamFoam / OF10 | 上游 README 声明 OF10 支持；实际分支 `OpenFoam_Org_main`，查询到提交 `38dd5fb967f1fd661dd914b45f6915d95c8a32e5` | 候选分支，未在本项目编译；目录含 laserbeamFoam、multiComponentlaserbeamFoam，不能认为与当前 compressibleLaserbeamFoam 文件一一对应 |
| CFDEMcoupling-PUBLIC / LIGGGHTS-PUBLIC | 官方 README 声明 OF6，并称近期不再更新；versionInfo 声明 CFDEM 3.8.1 / LIGGGHTS 3.8.0 | PUBLIC → OF10 仍需移植，不能直接承诺可用 |
| 扩展 CFDEM / OF10 | Aspherix 耦合安装文档要求 OF10，同时要求 Aspherix 安装及许可 | 不是本项目的开源 LIGGGHTS-PUBLIC 依赖方案 |

来源：[LaserbeamFoam](https://github.com/laserbeamfoam/LaserbeamFoam)、[OF10 求解器目录](https://github.com/laserbeamfoam/LaserbeamFoam/tree/OpenFoam_Org_main/applications/solvers)、[CFDEM-PUBLIC](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC)、[versionInfo](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/src/lagrangian/cfdemParticle/cfdTools/versionInfo.H)、[扩展版安装说明](https://doc.aspherix-dem.com/coupling/Section_installation.html)。查询日期为本页日期，实施前固定各依赖提交。

## 借用哪些算法，哪些必须重新推导

CFDEM-PUBLIC 不只有未解析经验阻力模型。[cfdemSolverIB](https://www.cfdem.com/media/CFDEM/docu/cfdemSolverIB.html) 明确面向颗粒大于网格的解析浸入边界/虚拟域求解；可借用颗粒定位、覆盖单元、几何映射、刚体速度约束、流体力/力矩、嵌入 LIGGGHTS 及 MPI 交换结构。它是算法起点，不是现成 LPBF 求解器。

1. **连续相方程需统一。** 原 IB 求解器为单相不可压缩形式；我们的蒸发体积源、热膨胀/组分变密度和移动固体排流必须进入同一压力约束。不得将原 IB 压力修正和现有蒸发压力修正直接前后叠加，否则后一投影可能破坏前一约束。
2. **受力可以用等价体积分。** ShirgaonkarIB 源码累加覆盖单元内力密度乘体积，并可累加力矩；不应误称为直接几何表面积分，也不应因此当成经验阻力。必须核对虚拟内部流体惯性、压力单位、变密度、黏性应力和反作用账本。其 forceSubModel 不可压缩分支对运动学压力的梯度乘密度；当前求解器的 Pa 压力不能照搬该表达式。不要额外叠加一次经验拖曳。
3. **几何不能混拼。** 书稿中的开口体积/面面积方法与虚拟域 IB 不是同一种离散。先选定一种，定义固体占据率和液/汽 VOF 的控制体含义，再推导几何守恒；不能同时乘两套孔隙率。IBVoidFraction 的 alphaMin、scaleUpVol 是算法参数，不能未经检验用人工放大半径代表真实粉末。
4. **热与熔化需开发。** 单个 DEM 温度不足以替代颗粒内部导热。激光吸收、固体内温度/焓、接触传热、固液界面与颗粒退出规则须具有唯一质量和能量归属。整颗粒熔化后转换可作为首个受限验证案例，但不能代表部分熔化/润湿行为已经实现。
5. **通信与物理分开。** 优先封装 particleBackend，借用 twoWayMPI 嵌入库交换；固定稳定 ID、位置/速度/角速度、力/力矩、热量和共同时间窗。MPI 数据可交换不等于模型已守恒，也不等于 48 核效率足够。

源码依据：[IB 求解器](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/applications/solvers/cfdemSolverIB/cfdemSolverIB.C)、[解析受力](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/src/lagrangian/cfdemParticle/subModels/forceModel/ShirgaonkarIB/ShirgaonkarIB.C)、[力密度](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/src/lagrangian/cfdemParticle/subModels/forceModel/forceSubModels/forceSubModel/forceSubModel.C)、[体积映射](https://www.cfdem.com/media/CFDEM/docu/voidFractionModel_IBVoidFraction.html)、[twoWayMPI](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/src/lagrangian/cfdemParticle/subModels/dataExchangeModel/twoWayMPI/twoWayMPI.C)。

## 有限开发包与退出条件

以下均待用户认可计划后实施。每包一次集中反馈，失败以原因定位补丁收束，不无期限追加均匀微测试；不同时维护两套生产平台。

| 包 | 交付及测试范围 | 退出条件 |
|---|---|---|
| P0 平台与方程方案（当前讨论） | 固定 OF10/v2512 候选、CFDEM/LIGGGHTS 来源与提交；比较 VOF、热物性、激光、压力和构建接口；给出单一移植清单 | 用户确认平台评估方向、IB 或开口面表示及第一版物理范围；本阶段只记录和审计 |
| P1 平台可行性 | 若选 OF10：在独立评估分支保留原始上游，先验证 LaserbeamFoam 基例与最小单相 IB–LIGGGHTS 交换/固定球；若保留 v2512：只移植必要模型部件 | 用户侧集中构建与代表运行完成，明确必改接口、依赖、MPI ABI 和保留功能；评估不成立则回到 v2512，不开展两平台全面开发 |
| P2 连续相核心（M1 收束） | 修正已确认的摘要时间对齐问题；统一界面质量、相库存、压力及材料能量账本；随后保护气/蒸汽组分及动量闭合 | 同一集成包验证质量/能量、边界收支、空间/时间敏感性、代表性 MPI/重启；不因源项两侧相消而宣称总量守恒 |
| P3 解析颗粒机械/热模块（M2） | 固定/平移/转动球、几何排流、力/力矩与反作用；内部导热、跨 rank 迁移及共同重启 | 与可追溯单球阻力/热传导参考比较，检查几何与总动量/能量；网格小于粒径仅为必要起点，需分辨率敏感性证据 |
| P4 固液颗粒交互（M3） | 明确部分熔化表示、接触/润湿、激光遮挡与吸收；从受限整粒转换到实际目标范围 | 金属质量、动量、焓不重复；液体/固体所有权及接触退出明确；不能在粉床完整演算中掩盖转换错误 |
| P5 文献集成（M4） | 单材料参数核对 → 无粉床单道 → 固定粉床 → 动态粉床；用同一账本和版本复现 | 参数缺口公开登记；熔池尺寸与飞溅统计分别评估，跨求解器比较与实验验证分开 |
| P6 可选加速（M5） | 先量测光线、流体、几何定位、DEM 和通信耗时，再选择固体远场只求热、耦合窗口/区域策略 | 同物理模型的全耦合参考下报告误差、守恒和速度；保护气远场保留必要流动，不默认整片远场只求热 |

P1 的孤立单相机械可行性检查与 P2 连续相收束可分开准备，但完整 LPBF 双向耦合必须等二者通过。气体优先评估带源、变密度低 Mach 模式；当前可压缩解保留作比较。不可压缩速度约束必须含蒸发排体积，不能采用无源 div(U)=0。正式选择要根据目标工况 Mach 数、压力响应和质量/能量验证，OF10 本身不决定气体模型。

## 讨论时需要固定的选择

- 平台：是否接受“OF10 优先做有限可行性评估，失败回退 v2512”，或直接留 v2512 定向移植？当前建议前者，前提是较多复用 IB。
- 几何：首版是否采用 CFDEM 虚拟域 IB，并按 LPBF 方程重构约束，还是优先复现书稿的开口体积/面方法？建议 IB 作为首版工程候选，不能称为书稿算法的等同复现。
- 范围：首版验证球形刚体、内部导热和单组分金属/保护气；部分熔化和液体润湿必须在 P4 明确模型，不能跳过后直接宣称动态粉床验证。
- 依赖：若已有适配 OF10 的 CFDEM-PUBLIC fork，提供仓库/提交即可改变平台工作量判断；现有信息只有 LIGGGHTS 可执行文件，没有已证实的耦合库或 OF10 PUBLIC 安装。

上述建议尚未成为已批准的平台迁移决策；下一步仍是讨论、澄清及更新本页。
