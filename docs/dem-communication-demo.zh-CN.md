# v2512–LIGGGHTS 最小通信 demo 与审查结论

2026-10-09。用户已确认保留 OpenCFD v2512、定向移植 CFDEM 必要算法，并要求先审查再完成通信 demo。不再评估 OF10。本包已实现代码、运行入口与独立摘要，开发端仅做离线检查；真实编译/运行待用户反馈。

## 审查与实现选择

1. 官方 LIGGGHTS-PUBLIC `library.h` 提供带 MPI communicator 的 lammps_open、命令执行、颗粒数组读取及变量接口；不需要先移植整套 CFDEM。实际读取了官方 library.h/library.cpp/atom.cpp，核实 id、x、v、omega、rmass、radius 和 nlocal。查询上游 HEAD 为 `3d5c00f20519e6bb6eb6756f51f1ad36564e649d`；这是审查来源，不强迫替换用户已有 fork，运行报告记录实际源提交、头文件及共享库 SHA。
2. CFDEM twoWayMPI 的嵌入库和 communicator 复制思路可复用。本 demo 使用公开 C API，自编适配层，不复制整个 CFDEM 源码，不建立对 OF6 cloud 的依赖。
3. 现有 LIGGGHTS 可执行文件不等于共享库。需要 MPI-enabled `.so` 和同一源码的 library.h；头文件/二进制不匹配、MPI ABI 不一致可能导致崩溃。预检保存符号与 ldd，拒绝缺符号、缺依赖、MPI-stub 库及明显 MPI 路径不一致；这些检查不能证明所有 ABI 完全一致。
4. 官方 gather API 要求连续 ID；首包限定两个固定 ID=1/2、无删除插入，使用全局 gather。生产版必须换为稳定 ID 的分布式交换，不能拿此 demo 证明稀疏 ID、规模扩展或 CFD/DEM 网格映射已实现。当前 PUBLIC 的本地 tag 按 int 读取，不支持 BIGBIG/tag64 或不相容 API 的 LAMMPS fork。
5. 原 CFDEM IB 的单相压力修正、经验模型或热模型没有混入本包。M1 现有未通过项保留；此处独立验证工程通信，不改变 M1 结论。

源码依据：[库接口](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/3d5c00f20519e6bb6eb6756f51f1ad36564e649d/src/library.h)、[库实现](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/3d5c00f20519e6bb6eb6756f51f1ad36564e649d/src/library.cpp)、[twoWayMPI](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC/blob/master/src/lagrangian/cfdemParticle/subModels/dataExchangeModel/twoWayMPI/twoWayMPI.C)、[共享库构建](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/3d5c00f20519e6bb6eb6756f51f1ad36564e649d/src/Makefile.shlib)。

## 本包证明什么

程序是真正以 wmake 构建的 OpenFOAM 应用，使用 OpenFOAM Time 和并行初始化；在同一进程/rank 集合里通过库调用推进 LIGGGHTS。它不求解流体方程，不需要生成 CFD 网格，不要求 CFDEM 安装。

每个 0.01 s 的 OpenFOAM 窗口：读取两颗球的状态 → 根据实际速度计算反馈力 → 将力通过 LIGGGHTS 的 group/fix addforce 命令回传 → DEM 以 0.001 s 推进十步 → 读取末状态和独立 DEM 时钟。球质量由密度 2500 kg/m³、半径 0.005 m 确定；这些是通信合成数据，不是 LPBF 材料参数。两个球 y 位置错开、不接触，沿 x 反向移动并越过双进程 DEM 分区。

反馈为 F_x=m·2·(v_target−v_x)，目标速度为 ±0.6 m/s；每窗力固定。按分段恒加速度解检查位移、速度及实际动量增量。反向冲量只记账，没有施加到真实流体；不能称为完整双向 CFD–DEM 动量守恒验证。转速只读取并检查零值，未回传力矩；无热交换、重启、解析壁面受力或碰撞验证。

固定门槛：20 窗×2 颗粒；位置误差 ≤1e−10 m，速度误差 ≤1e−10 m/s，冲量残差 ≤1e−12 kg·m/s，时钟误差 ≤1e−12 s；双进程两粒子均须出现所有者迁移。Python 摘要独立从轨迹重算反馈、运动和时钟，不信任日志自报误差；串行/双进程逐窗位置和速度差 ≤1e−10（分别按 m、m/s）。原有 LPBF Allwmake 不加入 DEM 依赖，只单独编译本应用。

## 用户侧操作

先运行以下命令。脚本优先使用显式变量，否则寻找 `~/CFDEM/LIGGGHTS/src` 及附近常见库路径。找不到库也会生成诊断归档，不会自动修改或重编你的 LIGGGHTS。

```bash
cd ~/LPBF-CoupledFoam
git pull --ff-only
export OPENFOAM_BASHRC=/home/kris/OpenFOAM/OpenFOAM-v2512/etc/bashrc
source "$OPENFOAM_BASHRC"
export LIGGGHTS_SRC="$HOME/CFDEM/LIGGGHTS/src"
# 改成实际共享库路径；若脚本能自动找到，可不设置此项。
# export LIGGGHTS_LIB="$HOME/CFDEM/LIGGGHTS/src/libliggghts.so"
bash scripts/dem-communication-test.sh
```

若仅安装可执行文件，先找库，不要盲目重编正在使用的安装：

```bash
find "$HOME/CFDEM/LIGGGHTS" -maxdepth 4 \
  \( -name 'libliggghts.so*' -o -name 'liblammps.so*' -o -name 'liblmp*.so*' \)
```

确认找到的库与 library.h 属于同一构建。若确实没有库，官方 Makefile 路线通常为在独立的、已配置 MPI 的源码副本 `src` 下执行 `make makeshlib`、`bash Make.sh models`、`make -f Makefile.shlib auto -j8`。这不是已核实适用于用户现有 CMake fork 的命令；现有 fork 应根据其实际构建选项生成 MPI 共享库，依赖错误归档可供下一步判断。不要使用 MPI STUBS，也不要为此更换 OpenFOAM。

反馈包在 `.runs/dem-communication-时间戳-PID.tar.gz`，包含实际提交、环境、库/头文件哈希、导出符号、依赖、构建日志、串行/MPI 日志与 CSV、独立 summary。依赖或构建失败同样归档；通过前不宣称路径已由运行证明。

## 后续边界

本包通过后才能判定 v2512 与该 LIGGGHTS 安装的最小嵌入通信可行。下一开发包先设计/验证解析颗粒几何与受力约束，建立真正 CFD 反作用；并行分布式映射、共同重启随后纳入集中验收。M1 守恒修复独立收束，完整蒸发–颗粒耦合需两者通过。用户仍执行所有真实编译和仿真，开发端不调用 WSL。
