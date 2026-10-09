#!/usr/bin/env bash
# User-run milestone: continuous, checkpoint/restart and two MPI ranks.
set -eo pipefail
source "$(dirname "$0")/environment.sh"
cd "$COUPLED_ROOT"
report_dir="$COUPLED_ROOT/.runs/m1b-portability-$(date +%Y%m%d-%H%M%S)-$$"
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
for mode in continuous restart mpi2; do
    stage="$mode"
    case_dir="$report_dir/$mode"
    if (
        python3 scripts/prepare-m1b-case.py "$case_dir" laser-phase fine --profile verification || exit $?
        cd "$case_dir" || exit $?
        blockMesh > log.blockMesh 2>&1 || exit $?
        setFields > log.setFields 2>&1 || exit $?
        if [ "$mode" = restart ]; then
            foamDictionary system/controlDict -entry endTime -set 1e-8 > log.checkpoint-config 2>&1 || exit $?
            foamDictionary system/controlDict -entry writeInterval -set 20 >> log.checkpoint-config 2>&1 || exit $?
            compressibleLaserbeamFoam > log.first 2>&1 || exit $?
            grep -q '^End$' log.first || exit $?
            foamDictionary system/controlDict -entry startFrom -set latestTime > log.restart-config 2>&1 || exit $?
            foamDictionary system/controlDict -entry endTime -set 2e-8 >> log.restart-config 2>&1 || exit $?
            compressibleLaserbeamFoam > log.restart 2>&1 || exit $?
            grep -q '^End$' log.restart || exit $?
            python3 "$COUPLED_ROOT/scripts/summarize-m1b-portability.py" "$case_dir" --merge || exit $?
        elif [ "$mode" = mpi2 ]; then
            cat > system/decomposeParDict <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
EOF
            decomposePar > log.decomposePar 2>&1 || exit $?
            mpirun -np 2 compressibleLaserbeamFoam -parallel > log.compressibleLaserbeamFoam 2>&1 || exit $?
        else
            compressibleLaserbeamFoam > log.compressibleLaserbeamFoam 2>&1 || exit $?
        fi
        grep -q '^End$' log.compressibleLaserbeamFoam || exit $?
    ) > "$report_dir/log.$mode.wrapper" 2>&1; then
        echo 'exit_status=0' > "$case_dir/case-result.txt"
        echo "Completed: $mode"
    else
        status=$?
        mkdir -p "$case_dir"
        printf 'exit_status=%s\n' "$status" > "$case_dir/case-result.txt"
        echo "Failed: $mode; collecting remaining cases."
    fi
done
stage=portability-summary
python3 scripts/summarize-m1b-portability.py "$report_dir"
stage=complete
echo 'Synthetic restart/MPI inventory gates passed; physical validation remains pending.'
