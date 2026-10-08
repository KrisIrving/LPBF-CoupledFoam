# V3.0 源码入口与第一次守恒审计清单

范围：固定上游提交的静态阅读，不是物理正确性认证。以下相对链接均指本仓库文件；后续重构时应同步更新。

| 模块 | 代码入口 | 目前观察与开发含义 |
|---|---|---|
| 不可压缩主循环 | [laserbeamFoam.C](../applications/solvers/laserbeamFoam/laserbeamFoam.C) | 在 PIMPLE 循环内推进相分数、物性、激光沉积、动量、温度和压力；颗粒交换的时序必须与这些更新协调 |
| 激光 | [laserHeatSource.C](../src/laserHeatSource/laserHeatSource.C) | 共用激光热源类；主循环通过 `laser.updateDeposition` 更新沉积场。动态颗粒的遮挡/吸收需要定义与 VOF 金属的表示关系 |
| 温度与蒸发冷却 | [TEqn.H](../applications/solvers/laserbeamFoam/TEqn.H) | `laser.deposition()` 加热、`Qv*mag(gradAlpha)*thermalDamper` 蒸发冷却、`TRHS` 熔化潜热。存在蒸发冷却不代表已建立显式蒸汽质量守恒 |
| 动量与反冲 | [UEqn.H](../applications/solvers/laserbeamFoam/UEqn.H) | 含 Darcy 阻力、Marangoni 项与 `pVap` 反冲作用。显式蒸汽源加入后必须重新核对压力跳跃与反冲是否重复 |
| 可压缩主循环 | [compressibleLaserbeamFoam.C](../applications/solvers/compressibleLaserbeamFoam/compressibleLaserbeamFoam.C) | `mixture.solve(&mass_dot)` 返回 `vDot`；另有密度方程。需把相间交换、压力和热源作为整体审计 |
| 蒸发/凝结交换 | [multiphaseMixtureThermo.C](../applications/solvers/compressibleLaserbeamFoam/multiphaseMixtureThermo/multiphaseMixtureThermo.C) | 有饱和压力、温度、压力差、界面厚度与限幅相关实现；应检查系数、量纲、适用压力范围和网格依赖 |
| 压力闭合 | [pEqn.H](../applications/solvers/compressibleLaserbeamFoam/pEqn.H) | `vDot` 进入压力方程，不能把蒸发模型仅当作温度方程的一项移植 |
| 粉床初始化 | [setSolidFraction.C](../applications/utilities/setSolidFraction/setSolidFraction.C) | 从粒子位置/半径文件初始化金属相分数；这个工具本身不证明存在运行时双向 CFD–DEM 耦合 |

## 不能按变量名直接接耦合接口

[createFields.H](../applications/solvers/compressibleLaserbeamFoam/createFields.H) 中：

- `vDot` 的实际量纲为 `[0 0 -1 0 0]`，即 s⁻¹；附近注释写了质量源单位，不能据此解释。
- `mass_dot` 的实际量纲为 `[1 -3 -1 1 0]`，含温度维度，并出现在[温度方程](../applications/solvers/compressibleLaserbeamFoam/TEqn.H)中；不能直接视为 kg/(m³ s) 的颗粒质量源。

首次审计应对每个源项登记：物理定义、量纲、正负号、单位面积/单位体积转换、密度或比热换算、在哪个方程加入、在哪次外迭代更新。现有注释不一致只是审计线索，不足以单独判定求解器存在量纲错误。

## 下一步可验收交付

1. 分别给两套求解器画出质量、动量和能量源项流向，并标明未闭合的界面交换。
2. 单独验证激光功率收支、相变焓、静态界面蒸发与气体响应；通过后再选择连续相主干。
3. 定义颗粒交换数据结构，明确气体—固体颗粒、熔池—颗粒和颗粒熔化三类作用的所有权，避免同一颗粒同时作为 VOF 固体和 DEM 固体重复贡献质量与热容。
4. 加入总量诊断，再实现双向交换；加速模式复用相同诊断。

尚未完成：逐项方程推导、原始模型常数核验、并行守恒测试及可压缩算例运行。此文件不替代这些工作。
