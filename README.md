# LPBF-CoupledFoam

Research framework under development for coupled laser absorption, evaporation,
gas flow, melt-pool dynamics and moving powder particles in LPBF.

**Status: research framework under development.** Continuous-phase audit candidates
and a standalone OpenFOAM/LIGGGHTS communication demo are implemented. Full
resolved particle coupling and optional acceleration remain unimplemented;
the minimal demo has passed build, serial and two-rank communication checks.
An experimental single-sphere mechanical solver and eight-case integrated
package are now implemented; native build and all eight individual case gates
pass, including serial/MPI and serial common restart comparisons. The package
still fails fixed-sphere refinement. The mechanical engineering baseline is
frozen while actual-surface geometry and force discretization are audited.
An independent six-fixture [geometry reference](docs/m2a-geometry-reference.zh-CN.md)
passes offline checks; it has not changed the CFD constraint or force scheme.
An optional stationary true-surface extension candidate and a fourteen-case
[verification package](docs/m2a-surface-candidate.zh-CN.md) are now delivered.
Its first native build and all fourteen executions pass, but only three
individual physical gate sets pass. Wall slip and coarse-grid phase sensitivity
remain open. A quadratic reconstruction/continuity revision awaits native testing.
Full CFD/DEM physical verification remains pending.
This is not validation of resolved CFD forces or heat. See the
[demo audit](docs/dem-communication-demo.zh-CN.md) and
[M2A scope, fixed gates and user-side instructions](docs/m2a-resolved-coupling.zh-CN.md).

Based on [LaserbeamFoam](https://github.com/laserbeamfoam/LaserbeamFoam),
tag **V3.0**, commit `f42e08a0bfc1749675beadcf7c6a90334d7aa9f8`.
Original source, attribution and GPL-3.0 license are retained.
See the [upstream README](README.upstream.md) for solver documentation.

## Build and test

Baseline build and short serial/MPI runs have passed with **OpenCFD OpenFOAM
v2506 and v2512**. See the verification record for their scope. Install OpenFOAM,
GCC and OpenMPI before running the commands; Foundation OpenFOAM 9/10 is not
interchangeable. The scripts respect an already loaded OpenFOAM environment
unless `OPENFOAM_BASHRC` is explicitly set.

```bash
git clone --depth 1 https://github.com/KrisIrving/LPBF-CoupledFoam.git
cd LPBF-CoupledFoam
export OPENFOAM_BASHRC=/usr/lib/openfoam/openfoam2506/etc/bashrc
BUILD_JOBS=4 bash scripts/build.sh
bash scripts/smoke.sh
# Optional parallel execution check:
SMOKE_NPROCS=2 bash scripts/smoke.sh
```

The shallow clone is sufficient for testing on another machine. Use
`git fetch --unshallow origin` when full development history is needed.

For the Ubuntu v2512 machine, use the
[remote testing instructions](docs/remote-testing.zh-CN.md).
One command builds and runs serial/MPI checks, then packages logs on success or
failure. Use a separate fresh checkout for each OpenFOAM version: intermediate
objects in source directories are shared even though final binaries use
version-specific directories.

The wrapper places executables and libraries in ignored `.build/`, without
replacing other user solver installations. Each smoke run copies the upstream
Plate2D tutorial into a new ignored `.runs/` directory. It uses 3,200 cells and
runs to 20 microseconds. This is an execution check, not an LPBF validation case
or a speed benchmark. The optional MPI check depends on working local MPI.

To use the built executables in another case:

```bash
source scripts/environment.sh
```

For substantial WSL production runs, prefer the Linux filesystem for the working
copy and results. The initial smoke check was performed on a Windows-mounted
workspace. Do not infer performance from its elapsed wall-clock values.

## Research and reproducibility

- [项目状态与接续入口](docs/project-status.zh-CN.md)
- [测试台账与结论](docs/test-register.zh-CN.md) · [机器可读结果](docs/test-results.json)
- [M1-B 空间界面与激光综合测试](docs/m1b-integrated.zh-CN.md)
- [M1B-02 重启、并行与分项账本](docs/m1b-portability.zh-CN.md)
- [M1 主要验收收束清单](docs/m1-closeout.zh-CN.md)
- [论文/书稿依据与未定模型决策](docs/model-decisions.zh-CN.md)
- [LaserbeamFoam–LIGGGHTS 通信与耦合讨论](docs/cfdem-communication.zh-CN.md)
- [解析颗粒、蒸汽注入与气体方程权衡](docs/resolved-gas-formulation.zh-CN.md)
- [M1C-01 界面源与气体响应：七组集中测试](docs/m1c-interface.zh-CN.md)
- [中文研究路线、文献基准与可选加速](docs/research-plan.zh-CN.md)
- [Development and verification record](docs/development-log.md)
- [V3.0 源码入口与守恒审计清单](docs/source-map.zh-CN.md)
- [M1 连续相源码审计与诊断检查](docs/m1-continuum-audit.zh-CN.md)
- [M1.2 质量率、潜热与可压缩测试](docs/m1-phase-transfer.zh-CN.md)
- [Ubuntu v2512 测试与反馈流程](docs/remote-testing.zh-CN.md)

Initial material target: 316L. The proposed reference is Zhang et al.,
Acta Materialia 288 (2025), 120816,
[DOI](https://doi.org/10.1016/j.actamat.2025.120816).
Parameter gaps must be resolved before claiming reproduction.

Only code, original project documentation and public references belong here.
Local manuscripts, private reference documents, credentials, compiled binaries
and simulation outputs are excluded from project additions.

### 最新反馈（2026-10-10，744930ac）

二次表面包整体未通过：10/14组完成、4组细网格转动个例通过，固定球精度仍失败。
三组偏移固定球的散度平台已由独立边界通量审计重现，另有一个固定球超时。
暂停原包重复测试；下一步集中审查边界通量与壁面/压力/受力离散相容性。
数值、失败原因及接续计划见[项目状态](docs/project-status.zh-CN.md)。当前未修改求解器。

### M2A-02C审查交付

已准备可选边界通量修复并完成完整图像点插值/应力与联合约束代数审计。
原生边界修复未测试，联合壁面求解器仍待实现，暂不重复运行原14组。
下一原生开发契约见[压力—壁面联合约束](docs/m2a-pressure-wall-contract.zh-CN.md)。

### M2A-02C原生联合候选已交付

新增真实球面插值/反力与实际phi压力响应的联合校正；37项离线检查通过。
原生编译/运行、精度和成本待验证，M2A仍未关闭。集中用户入口：
`bash scripts/m2a-joint-test.sh`。完整命令与限制见[原生联合候选](docs/m2a-joint-candidate.zh-CN.md)。
旧surface包保留作历史对照，当前测试使用新joint包。
