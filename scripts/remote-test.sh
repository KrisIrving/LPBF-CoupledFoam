#!/usr/bin/env bash
# Run on the user's Ubuntu machine after `of2512`.
set -eo pipefail
source "$(dirname "$0")/environment.sh"
cd "$COUPLED_ROOT"
report_dir="$COUPLED_ROOT/.runs/report-$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$report_dir"
stage=environment
finish_report() {
    report_status=$?
    printf 'last_stage=%s\nexit_status=%s\n' "$stage" "$report_status" > "$report_dir/result.txt"
    if tar -czf "$report_dir.tar.gz" -C "$(dirname "$report_dir")" "$(basename "$report_dir")"; then
        echo "Feedback archive: $report_dir.tar.gz"
    else
        echo "Could not archive report; logs remain at: $report_dir"
        report_status=1
    fi
    exit "$report_status"
}
trap finish_report EXIT
{
    date -Is
    git rev-parse HEAD
    git status --short
    cat /etc/os-release
    uname -m
    lscpu
    free -h
    df -h "$COUPLED_ROOT"
    printf 'WM_PROJECT=%s\nWM_PROJECT_VERSION=%s\nWM_PROJECT_DIR=%s\nWM_OPTIONS=%s\n' \
        "$WM_PROJECT" "$WM_PROJECT_VERSION" "$WM_PROJECT_DIR" "$WM_OPTIONS"
    printf 'BUILD_JOBS=%s\nSMOKE_NPROCS=%s\n' "${BUILD_JOBS:-8}" "${SMOKE_NPROCS:-2}"
    command -v wmake blockMesh foamDictionary mpirun
    g++ --version | head -n 1
    mpirun --version | head -n 3
} > "$report_dir/environment.txt"
stage=build
echo 'Building project; logs are collected in the feedback archive.'
BUILD_JOBS="${BUILD_JOBS:-8}" bash scripts/build.sh > "$report_dir/log.build" 2>&1
printf 'laserbeamFoam=%s\n' "$(command -v laserbeamFoam)" >> "$report_dir/environment.txt"
stage=serial
echo 'Running serial smoke check.'
SMOKE_NPROCS=1 SMOKE_REPORT_DIR="$report_dir/serial" bash scripts/smoke.sh > "$report_dir/log.serial-wrapper" 2>&1
stage=parallel
echo 'Running MPI smoke check.'
SMOKE_NPROCS="${SMOKE_NPROCS:-2}" SMOKE_REPORT_DIR="$report_dir/parallel" bash scripts/smoke.sh > "$report_dir/log.parallel-wrapper" 2>&1
stage=diagnostics
python3 scripts/check-phase-pair.py --output "$report_dir/m1-phase-pair.json"
python3 scripts/summarize-continuum.py "$report_dir/serial/log.laserbeamFoam" \
    "$report_dir/parallel/log.laserbeamFoam" --output "$report_dir/m1-summary.json"
stage=complete
echo 'Build, serial and MPI execution checks passed. Physical validation is still pending.'
