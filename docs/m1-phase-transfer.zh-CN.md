# M1.2 质量率、潜热与子循环测试

## 统一接口的推导

约定 mDot > 0 表示液体变为蒸汽，单位 kg/(m³ s)。固定相密度下：

```
S_liquid = -mDot/rho_liquid           [1/s]
S_vapour = +mDot/rho_vapour           [1/s]
volumeSource = mDot*(1/rho_vapour - 1/rho_liquid)
latentPower = -Lv*mDot               [W/m³]
```

质量加权源项和为零。蒸发膨胀、凝结收缩。这里的 latentPower 是采用共同显热参考、将相间能差单独列为潜热时的热源表示；若总能/相焓本身已经包含该能差，则不能再次叠加。质量携带的显热、动能和压力功也必须按所用能量变量处理。

V3.0 的原始相源可从 mDot=rho_vapour*S_vapour 推得上述配对，已通过八项独立代数检查。问题不在配对本身，而在蒸发系数除以 rho_liquid、实际转移量乘以 rho_vapour 的归一化解释。因此本轮不擅自改动交换率。

## 温度源为什么还不能直接替换

生产方程求解 rho*T，并以 rCv=Σ(alpha_i/Cv_i) 换算功率、压力功和动能项。它不是 1/Σ(alpha_i*Cv_i)。上游 mass_dot 已包含各相密度及 Cv 换算，单位 kg K/(m³ s)。

在独立检查的固定 Cv 条件下，记 r=rho_v/rho_l：

- 蒸发的两相 massgen 总和为 Lv*mDot*(1-r)/Cv_v。
- 凝结的总和为 Lv*mDot*(rho_l/rho_v-1)/Cv_l（mDot 为负）。

这不是可直接与 Lv*mDot 数值比较的功率。尤其凝结方向，不能只把 massgen 的符号当作“吸热/放热”结论：相质量变化也进入 rho*T 的左端。下一步需将该温度方程、密度变化、EOS 和参考能一起还原成能量账本，再决定替换方式。仅将 mass_dot 改为 rCv*Lv*mDot 并不足以证明守恒。

## 本轮的可选子循环修正

加入 alpha 求解控制项 averagePhaseChangeSources，默认 false，保留上游行为。true 时，PCR 和温度相变源与 rhoPhi 一样，按 dt_sub/dt_total 累积，再用于外层压力/温度求解。nAlphaSubCycles=1 不改变行为。

这解决三个外层量采用不同时间代表值的问题，尚不能保证整体离散守恒。速率限幅 0.5/dt 随子步变化、相分数迭代与 EOS 反馈仍需验证。当前不修改原始瞬时质量/热源公式。

## 用户侧最小运行门槛

```bash
cd ~/LPBF-CoupledFoam
git pull --ff-only
export OPENFOAM_BASHRC=/home/kris/OpenFOAM/OpenFOAM-v2512/etc/bashrc
source "$OPENFOAM_BASHRC"
BUILD_JOBS=8 bash scripts/m1-compressible-test.sh
```

脚本编译，然后复制上游多相 Test1 到新的 .runs 目录，修正应用名，保持 125 个单元、静态网格，以 1 ns 时间步运行 10 步。比较平均开/关和子循环 1/2/4 的六组串行运行。它包含原有多相、压力/EOS 及温度方程，不是 316L、裸板熔池或受控单次相变基准。保留原始材料和边界，只用于首次可压缩运行与开关检查，不从该短时间变化中推断准确度。

成功或失败均输出反馈压缩包，包含编译日志、每组字典和场、运行日志及源码代数报告。若某组失败，停止并保留失败组；开发端根据日志修复，用户无需自行调整材料参数。新增 C++ 代码尚未在开发端编译，六组 CFD 均等待用户测试。
