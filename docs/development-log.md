# Development and verification record

## 2026-10-09: M1C feedback reviewed; implementation paused for platform discussion

Reviewed user archive at tested commit 98befc42; build and all seven solver runs complete. Original summary fails all cases because diagnostic timestamp strings use different precision. Separate offline step-aligned review preserves thresholds: off, half-dt, kinetic and condensation pass their limited gates; prescribed coarse/fine and open fail inventory checks, with open also failing the initial-area threshold. Original failure retained. Detailed evidence and archive SHA are in the register and machine-readable record; no universal mass/energy or plume validation claimed.

User explicitly requests a detailed plan before further implementation. Added conditional OF10 versus v2512 platform plan and PUBLIC/extended-CFDEM distinction. PUBLIC supports OF6; OF10 LaserbeamFoam is a candidate, not proven drop-in compatibility. Corrected the resolved-force contract to permit validated equivalent IB volume integration as well as geometric traction integration. No source/script changes, WSL, compilation or CFD execution. Subsequent coding and migration await plan agreement.

## M1C-01 first interface-source package prepared: 2026-10-09

Added default-legacy interfacePhaseChangeModel with opt-in prescribed or pure-vapour kinetic flux. Evaluate geometry and flux once before all phase updates, map paired liquid/vapour volume sources and shared reference-energy source, bypass the legacy volume-distributed pair closure. Joint occupancy limiter logs requested/actual flux; existing pressure vDot uses the sum of paired volume sources. Restricted to equal constant Cv, reference source and one alpha subcycle, with explicit errors for active extra phases. No new Knudsen exit momentum, shielding-gas partial-pressure transport or low-Mach solver.

One user-side build and seven cases: off, prescribed mesh/time sensitivity, kinetic evaporation/condensation, open outlet. Flat smooth synthetic interface, no laser/gravity/surface tension. Diagnostics include area, pair mass/energy source, limiter, phase boundary flux and pressure response. Summary integrates last outer trial once and independently checks total/liquid/vapour inventory balances; thresholds fixed before runtime. Missing diagnostics are failures. Neither source cancellation nor power identity is claimed as universal physical conservation.

Ten offline tests pass, covering wrong inventory despite paired source, missing flux, outer-iteration counting, source off, condensation and open-field generation, plus existing portability tests. All seven cases generated and Python AST/metadata, shell syntax, JSON, document links and diff checks pass; no WSL, build or CFD executed. User-side compile/runtime remains pending. Plan/status/entrypoints updated; M1C-01 remains incomplete until runtime and remaining momentum/gas scope is resolved.

## M1B-02 accepted; resolved particles and gas formulation clarified: 2026-10-09

User tested a7395edddbe02701a968679e6b5744b8079c7ade. Build and all three forty-step trajectories complete; independent integration gates, restart and two-rank history comparisons pass unchanged thresholds. Restart maximum temperature difference is 9.999894e-9 K and reference-energy-change difference 5.273559e-15 J; MPI energy difference 3.209238e-16 J. The prior spurious fusion regression is absent. Q08 is closed for this fixed-mesh fully-liquid synthetic fixture only; combined checkpoint changes do not isolate each cause of the original mismatch. Earlier failures preserved, no strong-evaporation or material validation claim.

User explicitly selects resolved particles with mesh below diameter. Superseded prior geometry-aware empirical-force first candidate across canonical plan, decisions, communication and continuity instructions. Documented paired evaporation mass flux, source-driven volume constraint, variable-density low-Mach versus compressible pressure closures, moving solid geometric conservation, surface stress/torque and internal conduction contracts. Current compressible branch remains engineering baseline; low-Mach variable-density mode is prioritized for evaluation, not implemented or promised faster. Next package remains M1C-01 interface evaporation/gas response, then M1C-02 material energy and M2 resolved boundaries.

Read-only feedback archive analysis and local upstream source review; private archive/manuscripts remain outside Git. Offline reanalysis and JSON/document consistency checks only; no WSL, compilation or CFD execution. This commit changes curated records and design documentation, not solver equations.

## Checkpoint read-constructor regression; literature decisions and staged closeout: 2026-10-09

User tested `6a16c10b228043577250b0d3df95cf321245ab3c`: build and three 40-step runs complete, MPI comparison passes, restart comparison fails. New candidate causes 337.5 K cooling and about 0.000675 J false fusion uptake. Archive has saved epsilon1=1; restart warnings identify MUST_READ with default-value constructors. Direct read-only inspection of installed v2512 GeometricField.C confirms readIfPresent warns and does not read required-mode fields in these constructors. This regression was introduced by our prior candidate, not evidence that checkpoint serialization was missing.

Repaired candidate validates required phi/rho/epsilon1 headers separately and uses READ_IF_PRESENT in fallback-value constructors. Saved epsilon min/max are now logged and checked for this fully liquid synthetic fixture. Missing saved files remain fatal; original comparison thresholds unchanged. Build/CFD and the effect on the earlier flux mismatch remain user-side pending.

Reviewed supplied c03/c02 text plus critical formula/figure/material pages visually and the 2025 author paper/supplement. Recorded decisions D01–D08, proof table discrepancies, geometry-aware empirical particle route versus resolved surface stresses, conservative evaporation contract, variable-property energy scope, backend/MPI and full-melt transfer. Planned M1C split into two bounded packages; neither is implemented. No private PDF/page image or raw report added to Git. Offline archive reanalysis reproduces MPI PASS / restart FAIL; five Python log tests (including restored-fraction checks), candidate generation, Python AST, JSON parsing, shell syntax and diff checks pass. No developer-side compilation or simulation.

## M1B-02 MPI accepted, restart fails; bounded closeout and DEM design: 2026-10-09

User-tested `4a258dc5400a575fcf5cf6a6617b5146a5e87bac`. All three forty-step trajectories complete and independently pass integration gates. MPI history comparison passes. Restart mass relative difference 8.2578e-8, phase-volume/domain difference 4.1295e-8 and reference-energy history difference 6.2807e-10 J fail unchanged gates. Temperature/absorbed histories pass. Original failure retained.

Prepared default-off checkpointContinuumState candidate for compressible solver: saved phi header required at nonzero time, continuity rho written/read strictly, saved epsilon1 preserved, and startup zero-divergence CorrectPhi skipped while rAU creation retained. Existing dgdt fields already use READ_IF_PRESENT/AUTO_WRITE; no claim they were absent. Same three trajectories, thresholds unchanged. Runtime build and repair effect remain user-side pending.

User requested rapid M1 main verification: documented a bounded closeout matrix for meaningful evaporation/gas response, melting/moving laser and resolution sensitivity. Existing audit evidence reused with explicit scope; no new uniform micro-tests. M1 completion is not claimed. Added source-backed CFD–DEM discussion documenting upstream offline bed generation, embedded library/MPI proposal, unresolved-versus-resolved particle representation and conservative/time/checkpoint contracts. No DEM implementation or linked-library compatibility claim. Offline feedback reanalysis reproduces the original failure; five Python log tests, fixture generation, Python AST and shell syntax checks pass. No developer-side WSL, build or simulation.

## M1B-01 short integration accepted; M1B-02 prepared: 2026-10-09

User-tested `4438d4aada380e670974584262890309137c8f36`: all eight cases complete twenty steps and pass unchanged gates. Absorption is 2.57714867562 nJ from 10 nJ incident; advection errors are 0.10183205/0.08099039. Q06 is resolved only for the synthetic single-metal optical candidate; no multi-metal or physical validation. Laser phase exchange remains very weak; mesh peak temperatures differ substantially.

Prepared one milestone package: same fine laser-phase fixture over forty steps, continuous serial, checkpoint/restart at twenty steps, and two MPI ranks. Whole global inventory/temperature histories compared against predeclared tolerances. Raw restart segments retained; duplicated checkpoint diagnostics removed only in derived merged log. Explicit common Cv metadata enables signed temperature-equation budget diagnostics, never a new conservation gate. No C++ equation or production-default changes.

Development verification consists of offline feedback reanalysis, synthetic merge/comparison checks, case generation, Python AST and shell syntax checks. No WSL, compilation or CFD execution performed here. M1B-02 runtime results remain pending.

## M1B-01 completed runs expose transport/optics limits: 2026-10-09

User-tested commit `685261887b2a199613fae024207eecddda14f9f9`: all eight
cases ran twenty steps and ended normally. Temperature, mass and alpha-bound/
sum gates passed in all cases. Overall integration gate failed: coarse
advection error was 0.143000779994598 against 0.12 (fine 0.10183205477657);
four laser cases absorbed zero energy and had no temperature response. Both
spatial evaporation cases and fine advection passed all requested gates.
The initial CorrectPhi dictionary failure is resolved by runtime evidence.

Archived fields and source audit show electrical_resistivity was accumulated
from vapour/air in update.H, leaving pure metal essentially zero. Rays launched
with 1 W input, so no source-input failure explains the zero absorption.
Prepared default-off condensedPhaseOptics: sum positive condensed fractions
divided by conductivity, normalized by condensed fraction. Single-metal limit
is 1/sigma; only laser fixtures opt in. Added optional spatial M1_OPTICS values
and summary observations. Multi-metal mixing and actual absorption runtime
are not yet validated. Default optics behavior and phase/energy equations are
preserved. Algebra checks show zero resistivity gives zero normal-incidence
Drude absorption and the candidate's positive resistivity gives positive
absorption; this is not CFD validation.

Prepared revision 2 of the same eight-family matrix. Advection grids become
80/160 cells instead of 40/80 to resolve the observed interface diffusion;
old failures remain recorded. Other grids and every numerical threshold are
unchanged. Corrected documentation of actual ray sampling: default polar mode
uses 5 radial by 30 angular rays, not per-face subdivisions. State, plan and
registers distinguish the completed failed integration run from pending
candidate runtime. No WSL process or CFD run was invoked on the development
side; static generation/parser/algebra checks only.

## M1B-01 startup configuration failure and repair: 2026-10-09

User-tested commit `5e1f81da2667bb4a9dc5eb46b04accd4ea73c12f` built
successfully. All eight meshes and initial fields were prepared, but every
solver exited during initial CorrectPhi with a missing pcorrFinal dictionary
entry. No physical time steps were executed; spatial/laser behavior is not
evaluated by this archive. Failure provenance and per-case scope are retained
in the integrated test register, separately from M1-A history.

The generated fvSolution had omitted upstream pcorr.* when replacing its
solver section. Restored its PCG/GAMG controls without changing upstream
tolerance or any integration gate. Added a generator check for startup and
time-loop solver name coverage; the previous fixture is rejected and all eight
repaired dictionaries satisfy the offline check. Failed-case summaries now
include missing-entry details from the solver log. This is a configuration
repair to the same package, not a new physics test or equation change.

Updated state, plan and human/JSON registers. Repaired runtime remains pending
user execution of the same integrated wrapper. No WSL process or CFD run was
invoked on the development side.

## M1B-01 spatial/laser integration package prepared: 2026-10-09

Prepared one user-run integration wrapper with eight synthetic cases: two
resolutions each for matched-density slab advection, spatial phase exchange,
laser heating with phase exchange disabled, and laser heating with exchange
enabled. Meshes are 40/80 cells for one-dimensional translation and 216/1728
for the other families. Short mode runs twenty steps; one build and one archive
collect all cases, retaining failures rather than stopping after one case.

Added default-true phaseChangeEnabled to the existing liquid/vapour rate paths;
default numerical behavior is retained. Added default-off spatialDiagnostics
for phase bounds, sum, volumes, moments and exact cell-average translation
error on the specified Cartesian slab. Matched-density translation fixes a
pressure reference plane because this upstream branch lacks an all-incompressible
pressure reference path. Actual pair compression and diffusion tables are
controlled in generated fixtures, not merely the unused cAlpha setting.

Summaries apply predeclared execution, time-history, T, mass, phase-bound and
excitation gates, and report coarse/fine observations. Spatial reference energy
minus absorbed laser energy remains observational; no uniform energy gate is
reused. The existing volume-distributed phase law is not represented as a
validated interface-localized flux. Material reproduction, melting, full
boundary-energy balance, moving laser, restart and MPI checks remain pending.

Development-side checks: generated all eight fixtures, inspected dictionary
structure and periodic field conditions, checked shell/Python syntax and raw
phase-pair algebra. Synthetic parser logs exercise positive and negative gates,
including missing samples, phase bounds, source disabling, translation error
and failed-case retention. These are not CFD tests. Compilation/runtime await
user feedback; no WSL process or simulation was invoked. State, plan and pending
test register are synchronized with this package.

## Outer-iteration feedback and milestone documentation: 2026-10-09

User-tested commit `6464b2499d9ef63f2f2d1afdd5710103dee308f7`: all seven
cases completed, with temperature/mass and synthetic energy gates passing.
For condensation with subcycles 4, maximum energy changes at 5/7/9 outer
correctors were 38.6499/38.6599/38.6599 nJ; final phase inventory mass changes
were -1.25616e-14/-1.25715e-14/-1.25715e-14 kg. More outer correctors did not
meaningfully improve these errors. The hypothesis that increasing outer
iterations resolves this drift is not supported by this comparison.
No general cause is established; mass/energy-zero sensitivity remains an
explanation of the ledger's magnitude, not a mass correction.

In response to the user's request to avoid endless local tests, M1-A is now
closed as a limited foundational audit with Q01 retained. M1 remains incomplete;
the next development package is integrated spatial-interface and laser
verification (M1-B), not another outer/subcycle micro-test. No solver equations
or gate thresholds are changed in this documentation update.

Created a current-state handoff page, a human-readable test register and a
curated JSON register of all 14 locally available M1 compressible feedback
archives. Archive SHA-256 and tested commits preserve provenance; raw archives
and host details remain outside Git. Historical failed or unevaluated gates
are retained. Updated the research plan to v2512 and milestone packages with
acceptance requirements, and added project instructions to maintain these
records after feedback. The latest user result is distinguished from later
documentation-only commits. No WSL process or CFD simulation was invoked.

## Temperature budget feedback and mass/energy sensitivity: 2026-10-09

User-tested commit `6bf0a82b4f688200d392c4b7edeed7776376b615`: compilation
and all seven diagnostic cases completed; requested gates passed. Integrated
unrelaxed rho*T equation residuals were about 5.3--8.9e-12 J in condensation.
Mechanical, transport and diffusion terms were negligible in this uniform
zero-velocity fixture; laser was zero and fusion terms below 1.1e-13 J.
Thus the observed multi-subcycle ledger residual is not a large residual
of the reconstructed final thermal equation.

Summed phase inventories at 15-digit log precision reveal final mass changes
at deltaT 0.1 ns of -3.3177e-16, -8.4819e-15 and -1.2562e-14 kg for subcycles
1/2/4. Initial native specific energy is about 3.04228e6 J/kg. Its product
with these mass changes accounts for about -1.009, -25.804 and -38.216 nJ
of the original -1.452, -26.242 and -38.650 nJ endpoint energy changes.
After subtracting this common energy-zero sensitivity, endpoint changes are
about -0.442, -0.437 and -0.434 nJ. This does not correct mass or establish
conservation; it explains why an extremely small relative mass drift matters
when energy is compared on a much smaller phase-transfer scale.

Added a mass_energy_sensitivity observation based on original phase inventories.
The original energy ledger and all gate thresholds remain unchanged. The
pressure equation overwrites continuity rho with mixture.rho(); insufficient
outer coupling is a hypothesis, not a demonstrated cause. Next user script
m1-outer-energy-test.sh compares 5/7/9 outer correctors with subcycles 1/4,
fixed deltaT 0.1 ns, ten steps and one shared equilibrium reference: seven
small cases. It retains temperature budget diagnostics and does not change
solver equations or defaults. Static/parser checks only were run on the
development side; no WSL process or simulation was invoked.

## Energy refinement feedback and temperature budget diagnostics: 2026-10-09

User-tested commit `5c321a132a85f5c4a43e290c2abc0625d1d6c4dd`: all 13
cases completed with temperature, mass and synthetic energy gates passing.
Maximum closed-mass relative residual was 3.999784626294979e-12.
Condensation maximum energy residuals at deltaT 0.1 ns were 1.4734, 26.2417,
and 38.6499 nJ for subcycles 1/2/4; at 0.05 ns they were 0.7394, 13.1258,
and 19.3559 nJ. Halving deltaT approximately halves each residual, consistent
with a first-order temporal contribution, but two steps do not establish an
asymptotic convergence order. Evaporation stayed below 6.28e-11 J.

Archive-only integration of the final outer temperature source, multiplied by
the fixture's common Cv=800 J/(kg K), gives a condensation reference-inventory
versus source mismatch of about -0.453 nJ at 0.1 ns for every subcycle count.
The much larger ledger residual for subcycles 2/4 therefore needs additional
equation-term diagnostics; source averaging alone does not explain it. This
comparison does not isolate pressure work, fusion or transport effects.

Added default-off temperatureBudgetDiagnostics before mixture.correct in
TEqn.H. It records integrated storage, transport, diffusion, mechanical,
laser, fusion and phase-source terms of the existing rho*T equation, plus
current/old rho*T inventories and the reconstructed unrelaxed residual.
It does not alter the equation or solver defaults. Energy summaries integrate
only the final thermal solve per time, without counting PIMPLE iterations
as additional physical steps, and convert units only with explicit common
constant-Cv fixture metadata. Explicit reconstruction is a diagnostic for
this uniform case, not a general implicit-matrix conservation measurement.

Prepared m1-temperature-budget-test.sh: one equilibrium reference and six
condensation cases at subcycles 1/2/4 and deltaT 0.1/0.05 ns. This focuses
the next user run on the identified residual. Compilation and runtime of the
new C++ diagnostic await user feedback. Development-side checks are static
and archive/parser checks only; no WSL process or CFD run was invoked.

## Internal-energy reference feedback and refinement gate: 2026-10-09

User-tested commit `03639c9fbae1d97782aaa7b1e7255f8aea5ffa09`: all six
uniform synthetic cases completed, with temperature convergence and closed
mass checks passing. With referenceInternalEnergySource enabled, maximum
absolute reference-energy change was 1.4733814168721437e-9 J for condensation
and 5.820766091346741e-11 J for evaporation. Equilibrium gave
5.002220859751105e-11 J in both modes. Disabled-source active controls gave
1.6323913405358326e-4 J and 7.677044777665287e-9 J respectively.

Added an optional per-sample synthetic energy regression gate:
abs(delta E) <= 1e-9 J + 1e-4 * abs(Lv * delta vapour mass).
The absolute floor exceeds observed printed-energy cancellation noise; the
relative threshold is a project regression choice, not a physical accuracy
standard. It avoids normalizing by the much larger initial energy. Retrospective
checks accept the enabled cases and reject both disabled active controls.
The ledger now records the full diagnostic history, and rejects non-finite or
negative tolerances. Default observation-only behavior is retained.

Prepared m1-energy-refinement-test.sh: one equilibrium reference plus twelve
active runs, two directions, alpha subcycles 1/2/4, deltaT 0.1/0.05 ns and
end time 1 ns. All use the same consistent synthetic fixture, five outer
correctors, source averaging, old-time protection and the reference candidate.
Only one equilibrium run is retained to reduce redundant work. Runtime is
pending user feedback; this extends uniform numerical regression, not material
validation or full M1 completion. Solver defaults and equations are unchanged
in this update. Development-side work uses archive parsing and script checks;
no WSL process or CFD run is invoked.

## Initial baseline: 2026-10-08

- Base: LaserbeamFoam V3.0, `f42e08a0bfc1749675beadcf7c6a90334d7aa9f8`.
- Host: WSL2, Ubuntu 22.04.5, OpenCFD OpenFOAM v2506,
  GCC 11.4.0, OpenMPI 4.1.2; 32 logical CPUs, approximately 31 GiB RAM.
- Build: four jobs; all four expected programs produced:
  `laserbeamFoam`, `laserMeltFoam`, `compressibleLaserbeamFoam`,
  `setSolidFraction`.
- Execution: copied upstream Plate2D, 3,200 cells, end time 20 microseconds,
  initial and maximum time step 1 microsecond; upstream physical inputs retained.
- Serial and two-rank MPI runs reached `End`; MPI also reported
  `Finalising parallel run`. No fatal errors or non-finite values were detected
  in the solver logs by the smoke script.
- Local build and run logs remain in ignored `.build/` and `.runs/`.
- CI now builds and runs the same serial smoke script. Remote CI execution
  is separate from these local results; no remote pass is claimed here.

### Environment issues found

OpenFOAM's environment setup contains nonzero-returning probes; load it within
a conditional so shell `errexit` does not terminate setup prematurely.

OpenMPI startup stalled in hwloc's graphics discovery while connecting to a
local X server. A two-process `hostname` check and the solver both ran after
setting `HWLOC_COMPONENTS=-gl`. The environment wrapper applies this default
only on WSL and preserves an explicit user override. No system MPI settings
were changed.

The mounted Windows filesystem produced clock-skew warnings during compilation
and inconsistent wall-clock times in a short run. All required executables were
built and exercised, but these runs must not be used to report performance.
Prefer a WSL Linux filesystem checkout and a stable clock for timing studies.

An initial log check incorrectly matched the normal “floating point exception
trapping enabled” startup message. The check was corrected and rerun. Solver
exit status, normal completion, fatal errors and non-finite log values are checked.

### Scope of evidence

These checks establish build and short-run execution only. They do not validate
evaporation, gas compressibility, melt-pool dimensions, moving powder, spatter,
parallel numerical equivalence or long-time stability. Only `laserbeamFoam`
was run; the other programs were compiled. The 316L literature reproduction,
CFD–DEM implementation and acceleration modules are pending.

## User-side Ubuntu testing setup: 2026-10-08

The remote-testing workflow adds an OpenCFD v2512 compatibility candidate.
Scripts now preserve a selected environment, separate final
binaries by OpenFOAM version, reject recorded build-environment changes, and
package build/serial/MPI logs on success or failure. Source intermediate objects
still require separate checkouts per version.

The v2512 route is a compatibility candidate pending the user's first report.
No OpenFOAM compilation or simulation was run locally for this update. Heavy
GitHub build checks are now manually triggered; lightweight shell syntax and
diff checks are performed before pushing code.

## Reported v2512 baseline verification: 2026-10-08

Tested code: `9e33b47a4c93a9edd616664f59cd04237ee4f67b`.
The submitted feedback archive was inspected without invoking WSL or running
any solver on the development side. Original reports and machine metadata
remain outside this public repository.

## Temperature linear-solver comparison prepared: 2026-10-08

Added a two-case comparison of upstream smoothSolver and PBiCGStab/DILU
using the same absolute tolerance 1e-8, zero relative tolerance, 1000-iteration
limit and identical case physics. Only exact T/TFinal controls differ. This
is a diagnostic candidate, not a confirmed convergence fix. Source averaging
is disabled and alpha subcycles fixed to one to isolate the comparison.

The wrapper now records individual T solve residuals and enforces a separate
temperature convergence gate. The comparison explicitly allows upstream
failure as a reference but requires the candidate to attain tolerance. General
compressible tests require every case to pass the T check. Script syntax and
synthetic convergence-log checks pass; OpenFOAM execution awaits user feedback.
No WSL command, CFD run or tolerance loosening was used in development.

## Temperature candidate verified; controlled phase tests prepared: 2026-10-08

User report for `2e3a8ee6db1a51a029abe8cc655aedd0964e7cd7` confirms the
PBiCGStab/DILU T candidate reaches residuals between approximately 1.6e-17
and 3.3e-16 in one iteration for every one of ten solves. The upstream reference
still reaches 1000 iterations with residuals above 1e-8. Both terminate normally
and global observations remain nearly identical. This validates the candidate
for this small case, not universal solver behavior or a physical energy balance.

Added optional per-phase EOS mass and alpha-volume observations, plus a
controlled synthetic homogeneous pair preparation script. Initial alpha is 0.5
for metal1 and metal1vapour, zero for other phases; all nine upstream phases
remain in the thermo model. T=4101 K makes the active pair's initial Psat=P0.
Initial p=80000/120000 Pa drives evaporation/condensation respectively. Zero
laser power and gravity, no-slip velocity and insulated temperature boundaries
make a closed box. EOS, latent values and instantaneous source formula remain
upstream synthetic values, not 316L properties.

The next user wrapper runs twelve cases (two directions, averaging off/on,
subcycles 1/2/4) with PBiCGStab at unchanged tolerance. Inventory changes include
EOS effects and must not be equated to raw mass source. Compilation and CFD
execution await user feedback; preparation and parser checks are static only.

## Closed-box report: mass drift detected; coupling audit prepared: 2026-10-08

User report for `fe006dd58b8b3eea79a33797502f06f3c64de80d` contains twelve
normal completions with all T solves attaining tolerance. Closed-boundary mass
flux is zero throughout. Condensation increases mass by about `6.700058e-08`
kg (relative `2.671065e-05` for one subcycle), while evaporation decreases it
by `2.3107e-10` kg (relative `9.244475e-08`). Phase alpha changes have the
expected direction. Source averaging does not materially change these drifts.
Thus runtime and T convergence pass, but closed mass conservation fails.

Added a closed-mass gate for synthetic phase cases with a diagnostic relative
threshold 1e-9. This is a project regression threshold above printed precision,
not a universal scientific accuracy criterion. Applying it to the supplied
logs rejects all twelve cases. Zero-flux conserved synthetic logs pass.

Added optional couplingDiagnostics checkpoints after alpha, continuity,
temperature and each pressure correction. They report continuity-density and
EOS-density inventories, mean pressure, volume source and temperature source
integrals without changing equations. The suspected inconsistency is in full
alpha/pressure/EOS coupling; the report does not identify a proven root cause.
The next wrapper compares one versus three PIMPLE outer correctors for both
directions (four cases), holding alpha subcycles at one and source averaging
enabled. New C++ compilation and execution await user feedback. Reading the
archive and static checks used Windows only, without invoking WSL.

## Coupling iteration report and temporal checks: 2026-10-08

The user report for `260695052802c4aad20bd580ce65fe29705b65b8` contains all
four normal completions. Overall gate failure is expected because the two
one-outer-corrector references retain closed mass drift. With three outer
correctors, evaporation mass remains unchanged at logged precision and
condensation changes by approximately `4e-13` kg, relative `1.59465e-10`.
Both three-corrector cases satisfy the 1e-9 diagnostic gate over all logged
times, not just at the endpoint.

Checkpoints show alpha updates change EOS inventory while the continuity
density keeps its closed-domain inventory. The first pressure/EOS correction
overwrites rho with a different mixture inventory. Additional outer coupling
iterations reconcile these states to the gate tolerance in these cases.
This supports a coupling-lag explanation; it does not prove three iterations
are sufficient for arbitrary meshes, material closures or transient regimes.
The instantaneous phase-change formulas remain unchanged.

Added a twelve-case follow-up at three outer correctors, alpha subcycles
1/2/4 and time steps 1 ns/0.5 ns, both directions. Total duration stays 10 ns;
cases have ten/twenty steps. Source averaging stays enabled. The closed mass
gate now uses the maximum accumulated residual at every logged timestep,
preventing endpoint cancellation from hiding drift. Re-analysis of the four
submitted logs confirms the two expected passes/two failures; a synthetic
history with transient drift and a conserved endpoint is rejected. No WSL
command or solver was run during development; user execution remains pending.

## Temporal/subcycle report; old-time protection candidate: 2026-10-08

User report for `f93862db7c4302458e5665258d0543dafd4b92fd` completes all
twelve cases with converged T solves. At one alpha subcycle, both time steps
pass closed mass checks: condensation maximum relative residuals are
1.59465e-10 and 4.38530e-11; evaporation is unchanged at printed precision.
Condensation with two/four subcycles fails: maximum relative residual is
about 9.32e-8 at 1 ns and 7.06e-8 at 0.5 ns. Evaporation residuals stay below
1e-9, but vapour transfer also varies with subcycling. The overall failure
is therefore real and not a temperature convergence failure.

The solver only wraps the first alpha in subCycle. Other phase fields lack
the same preservation/restoration of global old-time histories across outer
iterations. OpenCFD v2512 subCycle.H was inspected through a Windows read-only
share: subCycleField preserves old/old-old values and supports updateTimeIndex.
Added optional protectAllPhaseOldTimes (default false) which wraps the remaining
phase fields with that guard. Time indices are protected once before time
advancement; old histories are restored after the subcycle time object ends.
No rate, latent or pressure formula was changed. This is a hypothesis-driven
candidate, not a confirmed fix. New C++ compilation/runtime await user feedback.

The four-case follow-up compares condensation subcycles 2/4 with protection
off/on, at three outer correctors and 1 ns. Known unprotected failures remain
in the overall mass gate; individual case summaries must be inspected even
if the wrapper exits nonzero. No WSL command or CFD run is used in development.

## Old-time protection candidate passes targeted comparison: 2026-10-08

User report for `cb915b72f894574f94e62440e8703feaeddff79b` confirms compilation
and four normal completions with converged T solves. Protected condensation
subcycles 2/4 have maximum relative mass residuals `6.77727e-11` and
`1.99331e-11`, both passing the 1e-9 gate. Unprotected references remain at
about 9.32e-8 and fail, so the overall archive status remains nonzero as intended.
Protected vapour mass changes are approximately `-2.40491505e-9` and
`-2.40491485e-9` kg, close to the earlier one-subcycle `-2.40491545e-9` kg.
This supports the history-protection candidate for this test; the default
remains false pending broader regression.

Prepared a twelve-case protected regression across both directions, 1/2/4
subcycles and two time steps. Added native per-phase thermo he inventory
observations. The current case thermos use sensibleInternalEnergy; no latent
or reference offsets are added, and no energy gate or conservation claim is
made. Physical phase-reference energies must be reconciled before a total
energy balance can be constructed. New diagnostics compilation and regression
await user execution; no WSL command or solver was run by the developer.

## Protected regression passes; reference energy ledger prepared: 2026-10-08

User report for `964fb4c5ffcb9825e6aee0243fc5db8344284efd` passes all twelve
runtime, T convergence and closed mass gates. Maximum history-relative mass
residual is 1.59465e-10. At fixed dt, protected subcycles give closely matching
vapour inventory changes. Halving dt changes one-subcycle condensation vapour
transfer by about 0.02%; this is sensitivity evidence, not a demonstrated order
of time convergence. Native sensible energy changes are about -2.09 mJ for
condensation and +6.23 microjoules for evaporation. These are phase-native
energies without matched latent reference offsets, not proof of physical loss.

Added a three-case equilibrium/evaporation/condensation energy observation
wrapper, using protected histories, three outer correctors, one subcycle and
PBiCGStab. At equilibrium T=4101 K, p=P0=100000 Pa, the active pair has no
initial saturation-pressure driving. The other two pressures remain 80000 and
120000 Pa. Diagnostic output precision increases to 15 digits to reduce energy
inventory cancellation. Original energy equations and source rates are unchanged.

The ledger computes a constant vapour energy offset from the equilibrium
initial phase energies and specific volumes, interpreting Lv as a reference
enthalpy gap. It enforces that reference relation once, then applies the same
offset to all cases without fitting their time histories. It records shifted
internal energy plus kinetic energy, not an energy pass/fail gate. Thermodynamic
EOS consistency and the off-reference phase energy gap still require audit.
Synthetic reference-calibration and equilibrium preparation checks pass; new
C++ diagnostic compilation and CFD execution await user feedback. No WSL
command or solver was invoked on the development side.

## Reference-energy report and common-latent candidate: 2026-10-08

User report for `55dbe8077f0c9f2d06eebd25c1e51db67caea736` passes all three
runtime/T/mass gates. Equilibrium reference energy changes by about -2e-11 J;
evaporation by -3.25155e-7 J and condensation by +1.24748e-4 J. The fixed
vapour offset is -921476.8488544316 J/kg. Non-equilibrium energy changes remain
after reference calibration and need equation/EOS analysis; total reported
energy is about 7.6 kJ, so normalization only by initial energy would hide
the much smaller phase-transfer scale. No physical energy pass is recorded.

Added optional commonLatentHeatSource, default false. For each liquid/vapour
pair, it forms Lv times the signed mass rate implied by the original capped
alpha sources, counts that pair once in the liquid branch, then uses the
existing rCv operator to convert power into the temperature-source units.
The temperature equation subtracts the resulting source, so evaporation cools
and condensation heats. Instantaneous alpha sources and pressure formulas
are unchanged. This isolates latent conversion; it does not supply full
phase-energy transport, EOS-consistent caloric behavior or pressure-work closure.
Pair molecular-weight/rate equivalence is only verified in the present synthetic
pair, not arbitrary multi-material settings.

Prepared six cases: equilibrium/evaporation/condensation with legacy/common
latent conversion. All use the same initial equilibrium calibration offset;
the ledger rejects differing initial reference states. Static shell checks,
sign/unit reasoning and reproduction of the submitted ledger offset pass.
Candidate C++ compilation and CFD runs await user feedback. No WSL command
or solver is invoked on the development side.

## Common-latent comparison and thermo identity audit: 2026-10-08

User report for `2e248420013fa7afe0a541d9d4cc73ee6b2cb9d9` passes all six
runtime/T/mass gates. Common latent conversion reduces condensation reference
energy change from +124.748 to +110.646 microjoules (about 11.3%), while the
evaporation value remains -0.325155 microjoules at recorded precision.
Equilibrium stays stable. The candidate does not close the energy ledger and
remains disabled by default. A successful wrapper status still means the
runtime/T/mass gates, not physical energy validation.

Inspected OpenCFD v2512 caloric/EOS source via Windows read-only sharing.
For the synthetic liquid hConst + adiabaticPerfectFluid, native sensible
internal energy is Hs-p/rho, with Hs independent of pressure. The EOS has
rho_T=0 and psi=rho/[gamma*(p+B)]. Thus the native derivative is
-1/rho+p*psi/rho^2, while the thermodynamic pressure-energy identity requires
p*psi/rho^2. The gap -1/rho cannot be repaired by a constant reference offset.
The synthetic vapour eConst + perfectFluid passes this pressure derivative
identity; that alone is not full thermodynamic validation.

Added a standard-library analytic audit that reads the actual synthetic phase
dictionaries, records their hashes and checks pressures 80/100/120 kPa.
Finite-difference checks reproduce the differentiated native liquid pressure
term. This establishes an incompatibility in the fixture's caloric/EOS pairing,
not the entire cause of solver energy error. No production EOS, source or
equation is changed. User-side script reproduction is next; no CFD rerun is
needed. No WSL command was invoked during analysis or development.

## User thermo audit verified; independent active-pair fixture: 2026-10-08

User thermo JSON matches the development report exactly, including both phase
dictionary hashes and all pressure derivative values. This confirms the legacy
liquid fixture's pressure-energy identity failure; it is not a full energy
solver diagnosis. The legacy tutorial and default physical solver remain intact.

Added an independent consistent active-pair fixture: liquid eConst/rhoConst,
rho=5000 kg/m3, Cv=800 J/(kg K); vapour eConst/perfectGas, molecular weight
50 kg/kmol, Cv=800 J/(kg K). The active pair reference vaporization enthalpy
is 2e6 J/kg at 4101 K and 100000 Pa. Both active native energies are independent
of pressure and satisfy the pressure-energy identity; the ideal gas R follows
molecular weight. The reference offset is still measured from equilibrium
thermo observations. Remaining upstream phases are initialized inactive.
Constant Lv away from the reference, phase transport and solver energy closure
are still under examination; this is not a real material model.

The user wrapper prepares three closed cases at 0.1 ns steps for a total 1 ns,
five outer correctors, one alpha subcycle, history protection and common-latent
candidate enabled. This new benchmark changes several fixture properties and
is not a single-parameter ablation of legacy behavior. OpenCFD v2512 registration
of both thermo combinations and rhoConst E=0 were checked by read-only Windows
source inspection. Preparation/analytic/shell checks pass; runtime is pending.
No WSL command or CFD run was performed during development.

## Consistent fixture report; internal-reference energy candidate: 2026-10-08

User report for `5786e340265758eb70efad1b48a483d78b5b28e0` passes all three
runtime/T/mass gates. Equilibrium reference energy changes by 5e-11 J.
Condensation energy change is +163.239 microjoules, compared with a latent
transfer scale 478.762 microjoules; evaporation is -7.677e-9 J compared with
2.267e-8 J. The roughly 34.1% ratio matches the reference p*delta-specific-volume
fraction of Lv. This is evidence for the enthalpy/internal-energy distinction
in this closed equal-Cv fixture, not a general energy-solver proof.

Added optional referenceInternalEnergySource (default false), requiring common
source conversion. It uses explicitly configured phase-reference energy
differences instead of Lv in the signed mass-exchange power. For this equal-Cv
pair, the reference difference C=Lv-p_ref*(1/rho_v-1/rho_l) is
1318067.1651452065 J/kg. It is taken from the submitted equilibrium initial
thermo calibration, not fitted to time-history errors. A modern SI constant
substitution gave a slightly different value; the implementation deliberately
uses the measured OpenCFD v2512 reference and checks it independently in the
ledger. Other OpenFOAM constant configurations require recalibration.

The consistent test case writes all phase offsets, with zero for inactive
phases. A six-case follow-up compares reference mode off/on for equilibrium,
evaporation and condensation. Rates, pressure equation and thermo properties
remain unchanged. This reference term is appropriate to the equal-Cv, common
native-energy fixture; variable or unequal Cv, species energy fluxes and general
moving-interface flows require a fuller conservative energy treatment.
Static preparation, energy-cancellation algebra and shell checks pass; C++
compilation and CFD execution remain pending. No WSL command was invoked.

- The build completed successfully and the wrapper checked all four expected
  executables. The run selected this project's version-specific executable.
- Serial and two-rank MPI Plate2D runs each used 3,200 cells and completed
  20 time steps, reaching `2e-05` seconds.
- Overall, serial and MPI exit statuses were all zero. Both solver logs ended
  normally; MPI also reported `Finalising parallel run`.
- Inspection found no fatal-error markers or NaN/Inf values in either solver
  log. No timing comparison is used to claim a speedup.

M0 build and short-run execution checks now pass on v2512. This records baseline
execution compatibility, not physical validation or field-level serial/MPI
equivalence. Only `laserbeamFoam` was exercised; other executables were built.
M1 conservation/source-term auditing and dedicated physical tests remain next.

## M1 source audit and optional observations: 2026-10-08

Added a public-source continuum audit and optional `continuumDiagnostics`
logging to both laserbeamFoam and compressibleLaserbeamFoam. The flag defaults
to false; smoke cases enable it. Logs report global mass, outward boundary mass
flux, absorbed laser power, kinetic energy, temperature range and maximum speed.
The feedback wrapper creates `m1-summary.json` with endpoint estimates; this
does not assess full energy balance, field equivalence or physical validation.

Static checks cover shell syntax, Python syntax, an analytic mass/power history
and rejection of non-finite observations. No WSL command, OpenFOAM compilation
or solver run was performed for this update. Compilation of these new additions
and user-side execution remain pending. The physics equations are unchanged.

## User-side M1.1 diagnostic verification: 2026-10-08

Tested commit: `66dc771d3dd06a50af8726b0a453084cabbcfb73` on OpenCFD v2512.
The supplied archive was read through Windows file sharing; no WSL command or
solver was executed on the development side. Build and serial/two-rank MPI
execution completed with zero exit statuses. Both logs have 21 diagnostic
samples from time zero to 20 microseconds and normal termination.

Both runs report initial/final mass `2.400356e-05` kg at printed precision,
absorbed energy estimated by right-endpoint integration `0.002154103865976` J,
and final absorbed power `110.486740967` W. Maximum endpoint mass residual
estimates are `4.301555069e-13` kg/s (serial) and `3.96838751553e-13` kg/s (MPI).
Final maximum temperatures are `380.334145290` K and `380.334144791` K;
field-level equivalence has not been evaluated.

Initial reported laser power is `-5e-09` W because the upstream constructor
initializes Deposition to `-1` W/m3 before the first laser update. It is an
initialization placeholder, not physical negative absorption; right-endpoint
integration does not use this initial value. Restart initialization likewise
must not be interpreted as a physical laser observation before update.

M1.1 diagnostic integration passes for this baseline. The run is below melting
and evaporation temperatures and does not validate phase-change sources,
compressible runtime behavior, full energy balance or gas/particle coupling.
Next is M1.2: isolated evaporation and condensation source-pair checks.

## M1.2 isolated source algebra: 2026-10-08

Corrected the earlier condensation audit: the actual liquid source uses
rho_v/rho_l, not rho_l/rho_v. Both raw phase-source pairs conserve mass under
fixed densities and identical paired rates. The earlier condensation failure
claim is withdrawn. Evaporation rate normalization still needs physical review:
the rate definition divides by liquid density but transferred mass uses vapour
density. No production source term was changed.

Added a standard-library Python audit which extracts actual alphagen/massgen
assignments, restricts evaluation to scalar arithmetic and min, and records
the source SHA256 and expressions. Eight synthetic checks cover evaporation,
condensation, density contrasts and rate limiting. These checks passed on the
development side without OpenFOAM or WSL. The feedback wrapper includes the
report. This is isolated algebra, not a full CFD or energy validation; M1.2
source integration and thermodynamic validation remain pending.

## User-side M1.2 algebra report verified: 2026-10-08

The submitted phase-pair JSON has source SHA256
`48497f9bca75c00d959098062316b5dc44648f4f46f05c405705fd78ea84942a`,
identical to the local source bytes. All extracted expressions and eight case
records match the development-side report exactly. Every isolated mass residual
is zero; evaporation has positive volume source and condensation negative.
Rate limiting halves the synthetic transfer as expected. Evaporation transfer
relative to donor-density rate interpretation is 0.5 and 0.001 for density
ratios 2 and 1000; condensation ratios are 1. These are synthetic algebra
observations, not measured material properties or CFD results.

M1.2 algebra checks are verified on both sides. The remaining tasks are to
derive the mass-rate normalization and thermal source consistently, then test
the integrated compressible solver with pressure/EOS and subcycling enabled.
No production physical source term was changed and no WSL command was invoked
to inspect this report. Raw user reports remain outside the public repository.

## M1 phase-transfer derivation and compressible runtime gate: 2026-10-08

Documented a shared signed mass-rate interface and latent-power convention,
including why the existing rho*T temperature source cannot be replaced from
units or signs alone. No mass-rate or instantaneous latent formula was changed.
Added optional alpha-control averagePhaseChangeSources (default false) to time
average PCR and temperature source with the same subcycle weights as rhoPhi.
The one-subcycle route remains unchanged.

Added a user-run wrapper for six short serial compressible runs, comparing
source averaging off/on and 1/2/4 alpha subcycles. It copies upstream Test1,
uses its 125-cell mesh without refinement and runs ten 1 ns steps. This is
a runtime gate using upstream synthetic materials, not a controlled material
validation case. Build and runtime are pending user feedback. Developer-side
checks are limited to script syntax, source algebra and weighted-source math;
no WSL command or CFD run is used.

## User-side compressible runtime gate: 2026-10-08

Tested commit `3753c9e0ecad7049f56d0120214d776671cffb12`. Build and all six
serial cases completed normally, with overall exit status zero. Archived
dictionaries confirm source averaging false/true and alpha subcycles 1/2/4.
Each case has 11 continuum observations, reaching 10 ns. No WSL command or
solver was invoked to inspect the user archive.

All cases report final mass `0.00384063832621` kg, versus initial
`0.00384063832664` kg at printed precision. Final maximum speed is
`0.0014733850156` m/s; internal temperature range is approximately
`949.999949468` to `950.000000002` K. Post-update absorbed laser power is zero.
Within each subcycle count, averaging on/off continuum histories are identical
at logged precision. This is not evidence that averaging is unnecessary:
vapour fractions remain essentially unchanged and phase-change excitation is
too weak for an informative source-averaging validation.

**Numerical qualification:** every one of the 60 T linear solves reports
1000 iterations with final residual above the configured `1e-8` tolerance.
For example, the final residual is `2.74518e-08`. Thus the runtime gate passes,
but linear convergence does not. The current wrapper checks termination and
fatal/non-finite errors, not attainment of every linear solver tolerance.
No source-accuracy, energy-conservation or subcycle convergence claim follows.
Next work must resolve the T solve and provide controlled nonzero phase change
before changing mass-rate or latent formulas. Raw fields and host information
remain outside this public repository.
