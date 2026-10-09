#!/usr/bin/env bash
# One build and one feedback archive. Invoked by the user on Linux.
set -eo pipefail
source "$(dirname "$0")/environment.sh"
cd "$COUPLED_ROOT"
report_dir="$COUPLED_ROOT/.runs/m1c-interface-$(date +%Y%m%d-%H%M%S)-$$"
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
echo 'Building interface candidate; output is retained in the feedback archive.'
BUILD_JOBS="${BUILD_JOBS:-8}" bash scripts/build.sh > "$report_dir/log.build" 2>&1
for selection in off prescribed-coarse prescribed-fine prescribed-half-dt kinetic condensation open; do
    stage="$selection"
    case_dir="$report_dir/$selection"
    if (
        set -e
        python3 scripts/prepare-m1c-interface.py "$case_dir" "$selection" || exit $?
        cd "$case_dir" || exit $?
        blockMesh > log.blockMesh 2>&1 || exit $?
        setFields > log.setFields 2>&1 || exit $?
        compressibleLaserbeamFoam > log.compressibleLaserbeamFoam 2>&1 || exit $?
        grep -q '^End$' log.compressibleLaserbeamFoam || exit $?
    ) > "$report_dir/log.$selection.wrapper" 2>&1; then
        echo 'exit_status=0' > "$case_dir/case-result.txt"
        echo "Completed: $selection"
    else
        status=$?
        mkdir -p "$case_dir"
        printf 'exit_status=%s\n' "$status" > "$case_dir/case-result.txt"
        echo "Failed: $selection; continuing to retain other case feedback."
    fi
done
stage=interface-summary
python3 scripts/summarize-m1c-interface.py "$report_dir"
stage=complete
echo 'Interface source gates passed; further momentum/material validation remains pending.'
