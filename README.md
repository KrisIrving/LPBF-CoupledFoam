# LPBF-CoupledFoam

Research framework under development for coupled laser absorption, evaporation,
gas flow, melt-pool dynamics and moving powder particles in LPBF.

**Status: upstream baseline only.** New particle coupling and optional acceleration
are planned, not implemented or validated in this revision.

Based on [LaserbeamFoam](https://github.com/laserbeamfoam/LaserbeamFoam),
tag **V3.0**, commit `f42e08a0bfc1749675beadcf7c6a90334d7aa9f8`.
Original source, attribution and GPL-3.0 license are retained.
See the [upstream README](README.upstream.md) for solver documentation.

## Build and test

Target: Ubuntu 22.04 with **OpenCFD OpenFOAM v2506**, GCC and OpenMPI.
WSL2 is supported by the local baseline workflow. Install these prerequisites
before running the commands; Foundation OpenFOAM 9/10 is not interchangeable.

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

- [中文研究路线、文献基准与可选加速](docs/research-plan.zh-CN.md)
- [Development and verification record](docs/development-log.md)
- [V3.0 源码入口与守恒审计清单](docs/source-map.zh-CN.md)

Initial material target: 316L. The proposed reference is Zhang et al.,
Acta Materialia 288 (2025), 120816,
[DOI](https://doi.org/10.1016/j.actamat.2025.120816).
Parameter gaps must be resolved before claiming reproduction.

Only code, original project documentation and public references belong here.
Local manuscripts, private reference documents, credentials, compiled binaries
and simulation outputs are excluded from project additions.
