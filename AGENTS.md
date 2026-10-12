# Project continuity

Before continuing development, read `docs/project-status.zh-CN.md`,
`docs/research-plan.zh-CN.md`, and the latest entries in
`docs/test-register.zh-CN.md` and `docs/development-log.md`.

The user runs compilation and CFD tests. Do not invoke WSL or run simulations
on the development machine without new user authorization. Reading supplied
feedback archives is allowed. Develop, inspect statically, document and push
within the user's authorized project scope.

Use milestone development packages and integrated test matrices. Do not
continually add isolated uniform-case micro-tests. M1-A is closed as a limited
audit, not universal physical validation. Next work is M1-B as described in the
research plan. Q01 remains unresolved; more outer correctors did not fix it.

For each feedback archive, record its tested commit, archive SHA-256, test
scope, outcome, numerical evidence, limitations and changed hypotheses in
`docs/test-results.json` and `docs/test-register.zh-CN.md`. Update project
status and development log; update the plan when scope or decisions change.
Do not overwrite historical failures with later pass results. Missing or
unevaluated checks are not passes. Do not claim unexecuted runtime validation.

Keep raw reports, machine details, credentials and private manuscripts out of
the public repository. Curated result values and project-authored documentation
may be committed. Do not alter numerical acceptance thresholds after observing
results merely to make a test pass.

The user fixed resolved particles on 2026-10-09: mesh is smaller than particle
diameter. Target moving solid boundaries, surface stress/torque, geometric
conservation and internal conduction; the previous empirical-force first
candidate is superseded. See docs/resolved-gas-formulation.zh-CN.md. M1B-02
passes its limited fixture; proceed to M1C-01, not more uniform micro-tests.

Latest user instruction (2026-10-09): v2512 selective CFDEM algorithm port
is approved; no OF10 assessment. Audit then implement an OpenFOAM/LIGGGHTS
communication demo. This supersedes the earlier implementation pause for
the approved route. Read docs/dem-communication-demo.zh-CN.md and the latest
status. The user runs native builds/CFD/DEM; no developer-side WSL. Keep
communication-only evidence separate from resolved-force, fluid reaction,
thermal, restart and production distributed-mapping validation. PUBLIC
officially supports OF6; equivalent IB force volume integration still needs
LPBF equation and unit audits. Preserve unresolved M1C failures.

Communication package passed at user-tested 026247c9 on 2026-10-09.
Do not ask to repeat it without backend changes. Next bounded scope is
docs/m2a-resolved-coupling.zh-CN.md: actual CFD geometry/constraints/forces
and representative MPI/common restart, not more communication micro-tests.
Communication pass does not close M1C failures or validate resolved physics.

M2A-01 implementation and eight-case user-side entry are prepared; native
build/runtime remain untested. Read the fixed scope/limits in the M2A document
before changing it. Keep this experimental isothermal noncontact single-sphere
solver separate from production LPBF. No full angular-fluid ledger, moving-mask
GCL, contact validation or coupled thermal model is claimed. Evaluate the full
returned package before requesting additional runs; preserve failed gates.

M2A build passed at user-tested 1bcf1413, but all eight cases became unstable
in one to three windows. Pressure/PISO structure revision is pending native
feedback, not an established physical fix. Preserve failed partial trajectories
and mark absent restart/comparison evidence not_evaluated. No thermal/GCL
expansion before assessing this revised concentrated package.

At user-tested f4598d84 all eight cases complete and restart/MPI comparisons
pass. Startup momentum and fixed-sphere refinement gates fail. Residual-controlled
PISO candidate is pending native feedback. Do not erase startup windows, loosen
thresholds or call nonmonotone two-grid errors convergence. Keep refinement
failure open while investigating geometry/support and force discretization.

At user-tested 7354a5ab six cases pass. Translation reaches budget64; MPI
residual double-counts coupled neighbours in the audited v2512 call path.
Diagnostic correction and budget96 await native validation; target unchanged.
Keep fixed refinement failure open; read the geometry-force audit before
new surface/geometry work. New geometry module is not implemented yet.

At user-tested f6f31a4f all eight cases pass individually, as do serial/MPI
and serial common restart comparisons. Overall fixed refinement still fails.
Freeze the current engineering baseline; no more pressure iteration micro-runs.
Read docs/m2a-geometry-reference.zh-CN.md: independent Python geometry oracle
and six offline fixtures pass, but CFD constraints/forces remain unchanged.
Accurate intersections may enlarge cell support; replacing chi alone is not
a surface-constraint fix. Next concentrated task is actual-surface/pressure/
force compatibility, then extended static-grid/phase and mechanical regression
package. No rerun of unchanged eight cases, no developer-side WSL/native CFD.

M2A-02B static surfaceExtension candidate and fourteen-case package delivered.
Read docs/m2a-surface-candidate.zh-CN.md before edits. Native build/runtime
pending. True-radius fluid-side ghost target, bounded outer updates, q4
fictitious inventory and last-solved-source force ledger; stress is diagnostic,
never a duplicate applied force. Keep default volumePenalty and old failures.
No moving/free-sphere mode in this candidate until static evidence passes.
Predeclared fine64/stress10% and phase2% gates must not be relaxed after feedback.
User runs scripts/m2a-surface-test.sh; no WSL/developer-side CFD/DEM. All case
failures retained. Full global Cartesian field replication is validation-only,
not production distributed mapping or acceleration.

User-tested cd6ddd85 surface package: all14 execute, 3 individual pass, overall
failed. MPI/restart and momentum pass; wall slip, coarse phase sensitivity,
one divergence and one fine stress mismatch fail. Read latest records.
Quadratic extension with relaxed target updates and explicit continuity guard
is a pending native candidate. Preserve linear mode, raw target norm and old
failures. Same14 cases/physics thresholds; no arbitrary radius/drag fit or
transfer of static evidence to moving/GCL/contact/heat. Analytic geometry audit
omits CFD and Cartesian image interpolation; never label it CFD validation.

User-tested744930ac quadratic surface package fails:10/14 complete,4 individual
fine-rotation passes,3 fixed-offset continuity aborts,1 fixed-finer timeout after14
windows. MPI/restart comparisons pass. Fixed drag/slip/stress worsen versus linear.
Offline boundary face-centre flux/cell volume reproduces3 divergence plateaus;
this explains boundary compatibility, not wall accuracy. Freeze both candidates;
do not request unchanged14-case reruns or blind budget/tolerance tuning. Read latest
status/plan: M2A-02C must audit boundary flux and pressure/wall/force compatibility
before selecting/implementing next operator. No solver fix/native pass claimed.

M2A-02C audit/contract delivered: read docs/m2a-pressure-wall-contract.zh-CN.md.
Opt-in compatibleGauss static analytic outer boundary candidate now exists;
native compilation/runtime pending. Default and old surface generator retained.
Full Cartesian image interpolation and stress oracle plus small joint-constraint
algebra pass offline32 M2A tests, not CFD evidence. Next native wall direction is
actual FV phi-pressure coupled with paired J/S surface constraints; not implemented.
Audit node independence, collocated flux compatibility and last-solved source
before delivering one concentrated native package. No unchanged14-case rerun,
no empirical force fits, no gate relaxation or thermal/GCL expansion.

M2A-02C native jointSurface candidate delivered; read docs/m2a-joint-candidate.zh-CN.md.
True-surface J and integrated-force transpose spread, compatible boundary, actual
FV phi pressure response and bounded nonsymmetric GMRES; full momentum residual
via PISO because rAU response is diagonal approximation. No ghost penalty in joint
mode. Same solved source for force/inertia ledger, native exchange audits/cost and
marker checkpoint.37 offline tests and6 geometric rank fixtures pass, historical
quadratic summary exact; native build/kernel/CFD still pending. User entry is
scripts/m2a-joint-test.sh (14 cases and actual shared kernel), not old surface
runner. Preserve all physical gates, budgets, failures; no more review-only turns
before processing this full feedback. No developer-side WSL/native simulations.

User-tested0e4cc6cc joint: native build/kernel pass,6/14 complete and individual
pass,6 timeout124,2 fixed64 momentum-budget96 failures. MPI/restart and exchange
pass; coarse load errors13.6/26% do not close precision, fine comparisons missing.
Read latest docs. Revision2 groups GAMG unchanged pressure tolerance, checkpointed
integrated-force warm start, full momentum predictor per joint correction and
reuse of explicitly checked response.39 offline checks pass; native revision
pending. Same scripts/m2a-joint-test.sh14-case entry and900s/96/64/physical gates.
Never relax gates or label timeout as convergence. No developer-side WSL/native.

User requested parallel test execution. Joint runner now supports bounded independent-case concurrency (JOINT_CASE_JOBS, default1, cap4 with runtime Linux resource limits) and BUILD_JOBS. Preserve serial/MPI2/restart definitions. All workers must complete before archive; concurrent timing is not directly comparable to historical isolated cost. Three Windows Git Bash queue tests pass; native concurrency pending. No developer-side WSL.

2026-10-11: user-tested4aa89072 revision2 fails first window in all14 cases (11 explicit residual failures,3 timeouts); build/kernel pass, comparisons not_evaluated. Revision3 disables repeated full momentum predictor while retaining diagonal Schur response, GAMG/warm seed/response reuse/concurrency.42 offline M2A tests and independent archive summary pass; native pending. Read latest docs, preserve budgets and no developer-side WSL.

2026-10-11 user requested renewed literature/method review after repeated failures. Read docs/m2a-literature-method-review.zh-CN.md before further solver changes. Current IB is not the book open-volume/area algorithm; actual FV response linearity/repeatability/pressure gauge/full residual/projected-Schur audit is missing. Revision3 remains pending, not a validated fix. Prepare operator diagnostics with the original concentrated matrix before another run request; no immediate unchanged14-case rerun, no gate/budget relaxation, no developer-side WSL.

2026-10-11: revision4 actual FV joint response/commit diagnostics delivered. Read docs/m2a-native-operator-audit.zh-CN.md. Same14-case user entry now requires native zero/repeat/linearity/scaling/gauge/divergence and manufactured-solve evidence; full impulse/radial gain are diagnostic only, no rank proof.45 offline M2A checks and two old joint summaries pass; native build/kernel/CFD pending. No WSL/developer simulation; original gates/budgets retained.

2026-10-11: user-tested2aa86d9d revision4 build/kernel and14 first-window operator/commit audits pass;5 complete individual passes,9 timeout124.1211/10320 pressure solves hit1000 iterations above original1e-14 tolerance; historical summary did not gate these. Read docs/m2a-joint-revision5.zh-CN.md. Revision5 uses GAMG80 then PCG/DIC within original1000 total budget, same tolerance and per-solve gate;46 offline tests pass, native pending. No gate/budget increase or WSL/developer simulation. Finest precision and MPI/restart comparisons remain open.

2026-10-11: user-testedc17d0040 revision5 native build/kernel pass,15924 pressure solves meet1e-14,10 full trajectories/individual gates pass plus MPI/common restart. Two fixed64 budget96 failures~4e-11 vs1e-11; two rotate64 timeouts20/24windows. Read docs/m2a-joint-revision6.zh-CN.md. Revision6 safeguarded full-state affine outer mixing within same frozen window updates U/p/all phi/F/acceleration/source; accept only actual residual decrease, original wall/div and affine residual check.48 offline checks pass; new native kernel/CFD pending. Original gates/budgets retained; no WSL/developer simulations; finest spatial gates remain open.

2026-10-11: user-tested522e62d5 revision6 native build/kernels pass,10 complete individual passes plus MPI/restart. Same two fixed64 budget96 failures and two rotate64 timeout22/22;1669s vs1686s, no material speedup. Read docs/m2a-joint-revision7.zh-CN.md. Revision7 bounded four-state global residual QR mixes coherent U/p/all phi/F/acceleration/source; same actual safeguards/budgets/gates,49 offline checks, native pending. If ineffective stop coefficient/history micro-tuning and design complete coupled operator; if finest precision fails audit spatial force/geometry. No WSL/developer simulations or thermal/GCL expansion.

2026-10-11: user-testedf3c5785b revision7 improves iterations:fixed64 first33 correctors and6/7 windows,rotate64 34/34;10 full individual passes plus MPI/restart, four finest timeout. Keep effective depth4, no coefficient micro-tuning. Read docs/m2a-joint-revision8.zh-CN.md. Revision8 sparse local-owner J/S/scaling, one scaling per frozen matrix, actual dense reference<=1e-12 and rank0 phase profiling/checkpoint evidence.51 offline checks pass, native pending. Original budgets/gates retained; partial load errors diagnostic only. Global Cartesian/stress mapping still replicated; no production/GCL/heat claim or WSL/developer simulations.

2026-10-12: user-tested0fd9e979 sparse/dense15 audits and14319 pressure solves pass;10 full individual/MPI/restart pass, four64^3 timeout6/8/35/35. Projection99.33..99.56% of recorded response;1829s vs1837s no material speedup. Read docs/m2a-joint-revision9.zh-CN.md. Revision9 retains solver revision8; opt-in64^3 MPI1/2/4/8 with decomposition/native ranks/owner/coverage checks, resource rank reservation and separate diagnostic progress-cost.json retaining timeout windows without passes.53 offline M2A and3 mock queue tests; native finer MPI pending. Preserve budgets/gates, no empirical force or WSL/developer simulation. Next if timeout split pressure assembly/solve and audit cache lifecycle, if full trajectory precision fails audit spatial geometry/force.
