# Development and verification record

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
