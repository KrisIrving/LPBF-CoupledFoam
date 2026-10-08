#!/usr/bin/env bash
# One build, eight spatial cases, one feedback archive; CFD runs only by the user.
set -eo pipefail
source "$(dirname "$0")/environment.sh"
cd "$COUPLED_ROOT"
profile="${M1B_PROFILE:-short}"
case "$profile" in short|verification) ;; *) echo 'M1B_PROFILE must be short or verification'; exit 2 ;; esac
report_dir="$COUPLED_ROOT/.runs/m1b-integrated-$(date +%Y%m%d-%H%M%S)-$$"
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
for family in advection spatial-evaporation laser-heat laser-phase; do
    for level in coarse fine; do
        stage="$family-$level"
        case_dir="$report_dir/$stage"
        # Retain failures and continue collecting the other case families.
        if (
            set -e
            python3 scripts/prepare-m1b-case.py "$case_dir" "$family" "$level" --profile "$profile" || exit $?
            cd "$case_dir" || exit $?
            blockMesh > log.blockMesh 2>&1 || exit $?
            setFields > log.setFields 2>&1 || exit $?
            compressibleLaserbeamFoam > log.compressibleLaserbeamFoam 2>&1 || exit $?
            grep -q '^End$' log.compressibleLaserbeamFoam || exit $?
        ) > "$report_dir/log.$stage.wrapper" 2>&1; then
            echo 'exit_status=0' > "$case_dir/case-result.txt"
            echo "Completed: $stage"
        else
            status=$?
            mkdir -p "$case_dir"
            printf 'exit_status=%s\n' "$status" > "$case_dir/case-result.txt"
            echo "Failed: $stage; continuing to collect remaining cases."
        fi
    done
done
stage=integrated-summary
python3 scripts/summarize-m1b.py "$report_dir"
stage=complete
echo 'Spatial integration gates passed; material/physical energy validation remains pending.'
