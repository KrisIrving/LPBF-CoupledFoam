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
