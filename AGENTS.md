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
