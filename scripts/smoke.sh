#!/usr/bin/env bash
set -eo pipefail
source "$(dirname "$0")/environment.sh"
test -x "$FOAM_USER_APPBIN/laserbeamFoam" || {
    echo 'Build this project first; the pre-existing user solver is not used.'; exit 1;
}
case_dir="$COUPLED_ROOT/.runs/plate2d-$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$case_dir"
collect_smoke_logs() {
    smoke_status=$?
    if [ -n "${SMOKE_REPORT_DIR:-}" ]; then
        mkdir -p "$SMOKE_REPORT_DIR"
        printf 'case=%s\nexit_status=%s\n' "$case_dir" "$smoke_status" > "$SMOKE_REPORT_DIR/result.txt"
        find "$case_dir" -maxdepth 1 -type f -name 'log.*' -exec cp {} "$SMOKE_REPORT_DIR/" \;
    fi
    exit "$smoke_status"
}
trap collect_smoke_logs EXIT
cp -a "$COUPLED_ROOT/tutorials/laserbeamFoam/Plate2D/." "$case_dir/"
cd "$case_dir"
cp -a initial 0
foamDictionary system/controlDict -entry endTime -set 0.00002 >/dev/null
foamDictionary system/controlDict -entry deltaT -set 0.000001 >/dev/null
foamDictionary system/controlDict -entry maxDeltaT -set 0.000001 >/dev/null
foamDictionary system/controlDict -entry writeInterval -set 0.00002 >/dev/null
foamDictionary system/controlDict -entry continuumDiagnostics -set true >/dev/null
blockMesh >log.blockMesh 2>&1
setFields >log.setFields 2>&1
if [ "${SMOKE_NPROCS:-1}" -gt 1 ]; then
    foamDictionary system/decomposeParDict -entry numberOfSubdomains -set "$SMOKE_NPROCS" >/dev/null
    foamDictionary system/decomposeParDict -entry method -set scotch >/dev/null
    decomposePar >log.decomposePar 2>&1
    mpirun -np "$SMOKE_NPROCS" laserbeamFoam -parallel >log.laserbeamFoam 2>&1
else
    laserbeamFoam >log.laserbeamFoam 2>&1
fi
grep -q '^End$' log.laserbeamFoam
if grep -Ei 'FOAM FATAL|SIGFPE signal|(^|[^a-z])(nan|inf)([^a-z]|$)' log.laserbeamFoam; then
    echo 'Smoke test failed: non-finite value or fatal error'; exit 1
fi
echo "Smoke test completed: $case_dir (execution check, not physical validation)"
