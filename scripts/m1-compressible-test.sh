#!/usr/bin/env bash
# Small runtime gate for the upstream multiphase case; not material validation.
set -eo pipefail
source "$(dirname "$0")/environment.sh"
cd "$COUPLED_ROOT"
report_dir="$COUPLED_ROOT/.runs/m1-compressible-$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$report_dir"
stage=build
finish() {
    status=$?
    printf 'last_stage=%s\nexit_status=%s\n' "$stage" "$status" > "$report_dir/result.txt"
    tar -czf "$report_dir.tar.gz" -C "$(dirname "$report_dir")" "$(basename "$report_dir")" || status=1
    echo "Feedback archive: $report_dir.tar.gz"
    exit "$status"
}
trap finish EXIT
git rev-parse HEAD > "$report_dir/commit.txt"
git status --short >> "$report_dir/commit.txt"
printf 'version=%s\noptions=%s\n' "$WM_PROJECT_VERSION" "$WM_OPTIONS" > "$report_dir/environment.txt"
BUILD_JOBS="${BUILD_JOBS:-8}" bash scripts/build.sh > "$report_dir/log.build" 2>&1
test -x "$FOAM_USER_APPBIN/compressibleLaserbeamFoam"
python3 scripts/check-phase-pair.py --output "$report_dir/m1-phase-pair.json"
for temperature_solver in ${M1_T_SOLVERS:-upstream}; do
for mode in ${M1_AVERAGING_MODES:-false true}; do
    for cycles in ${M1_SUBCYCLES:-1 2 4}; do
        stage="T-$temperature_solver-average-$mode-subcycles-$cycles"
        case_dir="$report_dir/$stage"
        mkdir -p "$case_dir"
        cp -a tutorials/compressiblelaserbeamFoam/Test1/{initial,constant,system} "$case_dir/"
        (
            cd "$case_dir"
            cp -a initial 0
            foamDictionary system/controlDict -entry application -set compressibleLaserbeamFoam >/dev/null
            foamDictionary system/controlDict -entry startFrom -set startTime >/dev/null
            foamDictionary system/controlDict -entry startTime -set 0 >/dev/null
            foamDictionary system/controlDict -entry endTime -set 1e-8 >/dev/null
            foamDictionary system/controlDict -entry deltaT -set 1e-9 >/dev/null
            foamDictionary system/controlDict -entry maxDeltaT -set 1e-9 >/dev/null
            foamDictionary system/controlDict -entry adjustTimeStep -set false >/dev/null
            foamDictionary system/controlDict -entry writeControl -set timeStep >/dev/null
            foamDictionary system/controlDict -entry writeInterval -set 10 >/dev/null
            foamDictionary system/controlDict -entry continuumDiagnostics -set true >/dev/null
            python3 - "$cycles" "$mode" "$temperature_solver" <<'PY'
import re
import sys
from pathlib import Path
p = Path('system/fvSolution')
text, count = re.subn(r'nAlphaSubCycles\s+\d+\s*;',
                     f'nAlphaSubCycles {int(sys.argv[1])};\n'
                     f'        averagePhaseChangeSources {sys.argv[2]};',
                     p.read_text())
if count != 1:
    raise SystemExit('Expected exactly one alpha subcycle control')
if sys.argv[3] == 'PBiCGStab':
    # Exact T/TFinal entries override upstream regex entries, retaining tolerance.
    text = text.replace('solvers\n{', '''solvers
{
    T
    {
        solver PBiCGStab;
        preconditioner DILU;
        tolerance 1e-8;
        relTol 0;
        maxIter 1000;
    }
    TFinal
    {
        $T;
    }
''', 1)
    if 'solver PBiCGStab;' not in text:
        raise SystemExit('Failed to insert exact T solver control')
elif sys.argv[3] != 'upstream':
    raise SystemExit('Unsupported temperature solver')
p.write_text(text)
PY
            foamDictionary constant/dynamicMeshDict -entry dynamicFvMesh -set staticFvMesh >/dev/null
            blockMesh > log.blockMesh 2>&1
            setFields > log.setFields 2>&1
            compressibleLaserbeamFoam > log.compressibleLaserbeamFoam 2>&1
            grep -q '^End$' log.compressibleLaserbeamFoam
            if grep -Ei 'FOAM FATAL|SIGFPE signal|(^|[^a-z])(nan|inf)([^a-z]|$)' log.compressibleLaserbeamFoam; then
                exit 1
            fi
        ) > "$case_dir/log.wrapper" 2>&1
        echo "Completed: $stage"
    done
done
done
stage=linear-convergence
python3 scripts/summarize-m1-temperature.py "$report_dir" ${M1_REFERENCE_FLAG:-}
stage=complete
echo 'Compressible runtime and requested temperature convergence checks passed; physical validation remains pending.'
