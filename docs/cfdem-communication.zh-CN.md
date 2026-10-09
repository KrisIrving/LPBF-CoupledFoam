# LaserbeamFoam—LIGGGHTS 通信与耦合讨论

## 2026-10-09 解析算法与平台说明

CFDEM-PUBLIC 已有 cfdemSolverIB 解析浸入边界候选，并非仅未解析经验拖曳。它的单相不可压缩压力修正不能直接与蒸发/变密度压力修正叠加。物理表面应力可用直接面积求积或经验证的等价 IB 体积分取得；ShirgaonkarIB 实际采用覆盖单元力密度体积分及可选力矩，需审计运动学/Pa 压力、虚拟内部流体惯性与作用反作用，不能再叠加一次经验拖曳。书稿开口面几何与虚拟域 IB 不能未经推导混用。

PUBLIC 官方支持 OF6，OF10 是待评估候选；用户要求先讨论并敲定计划，暂停实现。详细来源、模块契约与有限开发包见 [平台与 CFDEM 计划](platform-and-cfdem-plan.zh-CN.md)。本节补充下文“应力/边界积分”措辞，不宣称耦合已实现。

状态：设计草案，2026-10-09；未链接或运行新的 DEM 后端。以下区分已经读取的上游源码事实与本项目的拟议设计。

## 现有功能与建议路线

本仓库 V3.0 的 `tutorials/laserbeamFoam/LPBF_small/Allrun` 先运行 LIGGGHTS 铺粉，读取导出的 locations，调用 setSolidFraction 后运行 laserbeamFoam。它是初始粉床生成链；没有仿真中每步的颗粒受力、位置和热量双向交换。setSolidFraction 的单元中心占据判断也不能直接充当动态颗粒体积映射。

公开 CFDEMcoupling 的 [twoWayMPI 源码](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/src/lagrangian/cfdemParticle/subModels/dataExchangeModel/twoWayMPI/twoWayMPI.C) 在 CFD 进程中构造 LIGGGHTS/LAMMPS 对象，复制 MPI communicator，调用输入脚本/运行命令及 getData/giveData。这不是两个独立命令行程序靠每步文件交换。LIGGGHTS 的 [fix_cfd_coupling.cpp](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/master/src/fix_cfd_coupling.cpp) 有 MPI/file 数据交换及颗粒属性注册入口，MPI 模式可由外部程序控制交换时刻。

**拟议首选：由 OpenFOAM 求解器驱动、LIGGGHTS 作为库嵌入，外面封装 particleBackend。** 先在串行和相同 MPI rank 数下建立最小交换闭环，再做分区映射和扩展。是否移植整个 CFDEMcoupling 或只采用必要的通信/模型部件，在兼容性审计后决定；不能把“CFD–DEM 方法”与“必须整体安装 CFDEMcoupling-PUBLIC”混为一件事。

公开仓库的 [versionInfo.H](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/src/lagrangian/cfdemParticle/cfdTools/versionInfo.H) 当前声明 CFDEM 3.8.1、LIGGGHTS 3.8.0、特定 OpenFOAM 6 提交。它不能证明直接兼容本项目的 OpenCFD v2512。独立 DEM 可执行文件能运行也不能证明已有可链接的 shared library、相同 MPI ABI 或所需 fix。后续需确认用户安装的 fork/commit、shared library 路径、编译 MPI 和启用的 coupling/heat 模块；不要求先切换已工作的 v2512 环境。

## 通信前必须决定：颗粒怎样占据流体网格

LPBF 熔池界面常需细网格；如果网格尺寸小于粉末直径，不能默认套用把颗粒当作单元内点源的未解析 CFD–DEM 模型。用户已经选择解析颗粒；下面其他路线只保留为方法比较，不作为第一实现目标：

| 表示 | 适用思路 | 对本项目的影响 |
|---|---|---|
| 未解析、体积平均颗粒 | 流体网格/平均体积能容纳颗粒统计 | 需要孔隙率、颗粒体积核映射及适合气体/多相环境的受力模型；单纯 VOF 原方程加拖曳不足 |
| 解析颗粒或部分解析耦合 | 网格解析颗粒几何/附近流动 | 需要移动固体边界、浸入边界或一致体积分数约束；受力从应力积分取得，计算成本与网格运动更复杂 |
| 混合表示 | 近熔池与远处颗粒使用不同表示 | 需守恒转换及明确切换标准，后续研究再实现；不能在第一版静默混用 |

VOF 液相体积分数、DEM 固体占据率和流体孔隙率是不同量。它们的关系必须随选定控制体方程定义，不能把 alpha.metal 同时当成熔池液体与 DEM 颗粒体积，否则质量、阻力、导热或激光遮挡可能重复计算。

## 拟议交换契约

| 方向 | 最小数据 | 所有权与约束 |
|---|---|---|
| DEM → CFD | 稳定全局 ID、位置、半径、质量、平动/转动速度、温度/焓、固相状态 | LIGGGHTS 管接触与颗粒运动；数组排序及 rank 迁移不能改变 ID 含义 |
| CFD → 颗粒耦合层 | 壁面应力、表面热通量、颗粒内部温度/焓场及界面状态 | 解析边界积分给出力/力矩与热量；集总温度不能代替内部导热场 |
| 耦合层 → DEM | 一次定义的流体力/力矩、热交换功率及未来的质量变化 | 受力和热模型位置固定；两端不能同时再算一次拖曳或热通量 |
| 耦合层 → CFD | 移动壁面、开口体积/面积及相容边界反作用，转换时的质量/动量/能量 | 满足几何排流和累计交换收支；已通过边界施加的反作用不能再撒点添加一次 |

拖曳/界面力的冲量反馈应等量反向；热交换反馈亦然。压力、排体积及浮力等项必须与选定连续相方程一致，不是对所有粒子力盲目取负。能量要分清传热、流体对颗粒做功、接触耗散去向和质量携带焓。初版可以先固定粒径且不传质量，但接口保留未来的熔化转移；固定粒径测试不代表熔化转换实现。

激光能量必须有唯一去向：若 OpenFOAM 已按固体占据率吸收，则颗粒接收应从该项守恒分配；若直接在粒子端做光线吸收，连续相不能同时重复吸收。蒸发反冲与解析蒸汽压力/动量源同样需要避免重复。

## 时间同步、并行与重启

第一版采用明确耦合窗口：交换颗粒状态 → 定位并映射 → CFD 推进并计算粒子作用 → DEM 推进到同一窗口终点 → 汇总实际交换量和反作用。先用固定时间步/固定粒径；DEM 是否需要比 CFD 更小的子步由碰撞与接触稳定性决定，不能预设总是几十个子步。每窗终点必须同一物理时刻。

PIMPLE 外层迭代只更新同一个 CFD 时间步，不能每次外层都让 DEM 不可逆地向前再跑一遍。若需要强耦合，应有可恢复的窗初态和迭代收敛策略；显式耦合则必须记录滞后及窗口敏感性。后续自适应 CFD dt 需要独立的耦合时钟，处理不能整除 DEM dt 的末段，禁止两端时间漂移。

同 MPI rank 数不等于 CFD/DEM 空间分区一致。需稳定 ID、粒子所有者、CFD 单元所属 rank、迁移后的更新；跨多个单元或分区的映射守恒。首版可以采用已有交换机制验证正确性，但全局收集/广播不能未经性能测量就当成 48 核生产方案。

checkpoint 必须同时保存 CFD 场、DEM restart、全局 ID/颗粒状态、耦合窗口时钟及已累计交换量；只恢复 OpenFOAM 时间目录会重复或漏算交换。已有 M1 重启门槛为此提供连续相基础，尚未证明 DEM 重启可用。

## 讨论与实现顺序

1. 确定 Δx/dp、近场颗粒表示和相体积分数/孔隙率定义，这是模型层选择。
2. 审计实际 LIGGGHTS 库、fork/commit、MPI 与 coupling/heat 功能；决定移植范围。
3. 建立 ID/状态交换、时间同步、位置迁移及共同 checkpoint 的最小后端；暂不加完整 LPBF 物理。
4. 一次综合包验证单颗粒受力/升温、双向交换、两个 rank 迁移和重启，再进入动态粉床。
5. M3 实现颗粒到液相的质量/动量/焓转换及接触退出；VOF 液滴与残余 DEM 固体不能重复保留相同金属质量。

本页的首选方案为工程建议，尚未固定后端依赖或实现解析颗粒。详细讨论应先解决第一、二项，再写完整 coupling solver；通信成功只说明数据能交换，不说明模型守恒或粉末飞溅可信。

## 原文复核后的补充

用户于 2026-10-09 明确选择解析颗粒，此前几何占据率 + 经验受力第一候选撤回。第一实现目标为移动固体边界、几何守恒、压力/剪切力与力矩及内部导热；LIGGGHTS 管刚体运动和接触。气体是否低 Mach 与颗粒解析程度分别判断。详见 [气体与解析颗粒契约](resolved-gas-formulation.zh-CN.md) 和 [D04–D07 决策](model-decisions.zh-CN.md)。
