# Ubuntu v2512 测试与反馈

后续工作采用开发端修改并推送 GitHub、用户 Ubuntu 端编译测试、日志反馈的流程。开发端仅做必要的静态检查；本次脚本修改没有在本地运行 OpenFOAM。耗时的 GitHub 构建工作流改为手动触发。

## 环境选择

本仓库之前在 v2506 上通过了基线构建及短时运行。v2512 当前为待验证环境；已有 laserbeamFoam 程序不能替代对本仓库提交的编译验证。脚本使用项目独立目录中的新程序，保留已有安装。机器配置与本机安装路径不登记在公开文档中。

## 第一轮操作

在普通交互终端加载环境；`of2512` 如果是 alias/function，不能直接放在普通非交互脚本内调用。

```bash
of2512
export OPENFOAM_BASHRC="$WM_PROJECT_DIR/etc/bashrc"
printf 'Selected OpenFOAM: %s\n' "$WM_PROJECT_VERSION"

git clone --depth 1 https://github.com/KrisIrving/LPBF-CoupledFoam.git ~/LPBF-CoupledFoam-v2512
cd ~/LPBF-CoupledFoam-v2512
BUILD_JOBS=8 SMOKE_NPROCS=2 bash scripts/remote-test.sh
```

确认输出选择的是 `v2512`。如果目标目录已经存在且是本项目工作副本，只需进入该目录并 `git pull --ff-only`；已有未提交修改应先保留。每个 OpenFOAM 版本使用独立工作副本，避免复用其他版本的编译中间文件。

这轮测试只包括：项目编译、3200 单元的 Plate2D 串行检查、相同短时工况的双进程运行检查。模拟终点是 20 μs。默认不需要 LIGGGHTS、不需要 sudo，也不进行大规模 LPBF 计算。8 个编译线程和 2 个 MPI 进程用于初始兼容性排查，不是机器性能的最优配置。

## 反馈内容

脚本结束时会打印 `Feedback archive: .../.runs/report-日期-进程号.tar.gz`。包内包含：

- 提交号、工作区状态、OpenFOAM/编译器/MPI 版本和硬件信息。
- 构建日志、串行和并行求解器日志。
- 最后执行阶段及退出码；失败时也会生成包。

反馈该压缩包，并附终端最后几行即可；不需要发送整个计算目录或二进制文件。如果在加载 OpenFOAM 环境时就失败、尚未建立报告目录，请反馈该错误文字。

程序及库在 `.build/v2512/$WM_OPTIONS/`；运行结果与报告位于 `.runs/`。这些目录均被 Git 忽略。日志包不收集完整环境变量、密码或令牌，但会含本机路径和系统信息，仅供当前开发反馈。

## 后续迭代

```bash
of2512
export OPENFOAM_BASHRC="$WM_PROJECT_DIR/etc/bashrc"
cd ~/LPBF-CoupledFoam-v2512
git pull --ff-only
BUILD_JOBS=8 SMOKE_NPROCS=2 bash scripts/remote-test.sh
```

只有这轮通过后，才能登记 v2512 基线兼容。之后按具体修改提供更有针对性的测试，不会每轮都要求运行完整大算例。LIGGGHTS 后续接入时再核对版本和编译选项；alias 不会自动传入测试脚本。
