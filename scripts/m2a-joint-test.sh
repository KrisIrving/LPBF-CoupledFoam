#!/usr/bin/env bash
# User-side integrated package; never rebuilds DEM or runs LPBF thermophysics.
set -eo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
report="$repo/.runs/m2a-joint-$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$report"
stage=environment
finish() {
    rc=$?
    trap - EXIT
    printf 'last_stage=%s\nexit_status=%s\n' "$stage" "$rc" > "$report/result.txt"
    # Full meshes/fields remain in .runs; feedback retains logs, controls,
    # histories and common-checkpoint metadata/binaries, avoiding a huge upload.
    tar --exclude='*/constant/polyMesh' --exclude='*/U' --exclude='*/p' \
        --exclude='*/phi' --exclude='*/solidFraction' --exclude='*/core' \
        --exclude='*/core.*' -czf "$report.tar.gz" \
        -C "$(dirname "$report")" "$(basename "$report")" || rc=1
    echo "Feedback archive: $report.tar.gz"
    exit "$rc"
}
trap finish EXIT
cd "$repo"
git rev-parse HEAD > "$report/commit.txt"
git status --short > "$report/working-tree.txt"
source scripts/environment.sh > "$report/log.environment" 2>&1
if [ "$WM_PROJECT_VERSION" != v2512 ]; then
    echo 'M2A-02C currently requires OpenCFD v2512.' | tee -a "$report/log.environment"
    exit 1
fi
if [ -z "${LIGGGHTS_SRC:-}" ] && [ -z "${LIGGGHTS_LIB:-}" ] && [ -f .build/demo-liggghts.env ]; then
    source .build/demo-liggghts.env
fi
stage=dependency-audit
test -f "${LIGGGHTS_SRC:-}/library.h" && test -f "${LIGGGHTS_LIB:-}" || {
    echo 'Load .build/demo-liggghts.env or set both LIGGGHTS_SRC and LIGGGHTS_LIB.' | tee "$report/dependency-error.txt"
    exit 1
}
export LIGGGHTS_SRC="$(realpath "$LIGGGHTS_SRC")"
LIGGGHTS_LIB="$(realpath "$LIGGGHTS_LIB")"
if [[ "$LIGGGHTS_SRC$LIGGGHTS_LIB" =~ [[:space:]] ]]; then exit 1; fi
export LPBF_MPI_INC="$(mpicc --showme:compile)"
export LPBF_MPI_LIBS="$(mpicc --showme:link)"
test -n "$LPBF_MPI_INC" && test -n "$LPBF_MPI_LIBS"
{
    printf 'WM_PROJECT_VERSION=%s\nWM_OPTIONS=%s\nWM_MPLIB=%s\n' "$WM_PROJECT_VERSION" "$WM_OPTIONS" "${WM_MPLIB:-}"
    printf 'LPBF_MPI_INC=%s\nLPBF_MPI_LIBS=%s\n' "$LPBF_MPI_INC" "$LPBF_MPI_LIBS"
    printf 'LIGGGHTS_SRC=%s\nLIGGGHTS_LIB=%s\n' "$LIGGGHTS_SRC" "$LIGGGHTS_LIB"
    mpirun --version | head -n 3
    sha256sum "$LIGGGHTS_SRC/library.h" "$LIGGGHTS_LIB"
    git -C "$LIGGGHTS_SRC" rev-parse HEAD
} > "$report/environment.txt" 2>&1
ldd "$LIGGGHTS_LIB" > "$report/ldd.dem.txt"
if grep -q 'not found' "$report/ldd.dem.txt"; then exit 1; fi
nm -D --defined-only "$LIGGGHTS_LIB" > "$report/symbols.txt"
for symbol in lammps_open lammps_close lammps_file lammps_command lammps_free \
    lammps_extract_atom lammps_extract_global lammps_extract_variable; do
    grep -Eq "[[:space:]]${symbol}$" "$report/symbols.txt" || exit 1
done
linkdir="$repo/.build/$WM_PROJECT_VERSION/$WM_OPTIONS/dem-link"
mkdir -p "$linkdir"
ln -sfn "$LIGGGHTS_LIB" "$linkdir/libliggghts.so"
export LIGGGHTS_LIBDIR="$linkdir"
export LD_LIBRARY_PATH="$linkdir:$(dirname "$LIGGGHTS_LIB"):${LD_LIBRARY_PATH:-}"
stage=build
signature="$WM_PROJECT_DIR|$WM_OPTIONS|$LPBF_MPI_INC|$LPBF_MPI_LIBS|$LIGGGHTS_SRC|$LIGGGHTS_LIB|$(sha256sum "$LIGGGHTS_SRC/library.h" | awk '{print $1}')"
signature_file="$linkdir/m2a-environment.signature"
if [ -f "$signature_file" ] && [ "$(cat "$signature_file")" != "$signature" ]; then
    (cd applications/solvers/resolvedParticleFoam && wclean) > "$report/log.clean" 2>&1
fi
printf '%s\n' "$signature" > "$signature_file"
stage=operator-preflight
python3 scripts/audit-m2a-joint-stencils.py --output "$report/stencil-audit.json" > "$report/log.stencil-audit" 2>&1
g++ -std=c++17 -O2 -Iapplications/solvers/resolvedParticleFoam tests/native/wall-gmres.C -o "$report/wall-gmres-kernel" > "$report/log.kernel-build" 2>&1
"$report/wall-gmres-kernel" > "$report/log.kernel-test" 2>&1
stage=build
echo 'Building the experimental joint wall/pressure solver.'
(cd applications/solvers/resolvedParticleFoam && wmake) > "$report/log.build" 2>&1
app="$FOAM_USER_APPBIN/resolvedParticleFoam"
test -x "$app"
ldd "$app" > "$report/ldd.solver.txt"
if grep -q 'not found' "$report/ldd.solver.txt"; then exit 1; fi
dem_mpi="$(awk '$1 ~ /^libmpi\.so/ {print $3}' "$report/ldd.dem.txt")"
app_mpi="$(awk '$1 ~ /^libmpi\.so/ {print $3}' "$report/ldd.solver.txt")"
test -n "$dem_mpi" && test -n "$app_mpi"
if [ "$(realpath "$dem_mpi")" != "$(realpath "$app_mpi")" ]; then
    echo 'MPI library mismatch.' | tee "$report/dependency-error.txt"
    exit 1
fi
run_case() {
    local case_dir="$1" selection="$2"
    cd "$case_dir" || return
    blockMesh > log.blockMesh 2>&1 || return
    if [ "$selection" = fixed-mpi2 ]; then
        decomposePar > log.decomposePar 2>&1 || return
        timeout --signal=TERM --kill-after=10 "$case_timeout" mpirun -np 2 "$app" -parallel > log.solver 2>&1 || return
    elif [ "$selection" = fixed-restart ]; then
        python3 "$repo/scripts/prepare-m2a-joint.py" "$case_dir" --half || return
        timeout --signal=TERM --kill-after=10 "$case_timeout" "$app" > log.first 2>&1 || return
        test -s 0.002/couplingState && test -s 0.002/dem.restart || return
        test -s 0.002/U && test -s 0.002/p && test -s 0.002/phi || return
        cp input.dem input.first.dem
        cp log.dem log.first.dem
        cp system/controlDict system/controlDict.first
        python3 "$repo/scripts/prepare-m2a-joint.py" "$case_dir" --restart || return
        timeout --signal=TERM --kill-after=10 "$case_timeout" "$app" > log.restart 2>&1 || return
    else
        timeout --signal=TERM --kill-after=10 "$case_timeout" "$app" > log.solver 2>&1 || return
    fi
}
case_timeout="${JOINT_CASE_TIMEOUT:-900}"
if ! [[ "$case_timeout" =~ ^[1-9][0-9]*$ ]]; then echo "Invalid JOINT_CASE_TIMEOUT"; exit 1; fi
command -v timeout >/dev/null
failures=0
for selection in fixed-coarse fixed-coarse-offset fixed-fine fixed-fine-offset fixed-finer fixed-finer-offset rotate-coarse rotate-coarse-offset rotate-fine rotate-fine-offset rotate-finer rotate-finer-offset fixed-restart fixed-mpi2; do
    stage="$selection"
    case_dir="$report/$selection"
    python3 "$repo/scripts/prepare-m2a-joint.py" "$case_dir" --selection "$selection"
    echo "Running $selection (per solver invocation budget: ${case_timeout}s)"
    if (run_case "$case_dir" "$selection"); then
        printf 'exit_status=0\n' > "$case_dir/result.txt"
    else
        rc=$?
        printf 'exit_status=%s\n' "$rc" > "$case_dir/result.txt"
        failures=$((failures+1))
        echo "Failed: $selection; collecting remaining cases."
    fi
done
stage=summary
python3 "$repo/scripts/summarize-m2a-joint.py" "$report" --output "$report/summary.json"
test "$failures" -eq 0
stage=complete
echo 'M2A-02C static joint wall/pressure candidate gates passed; full physical validation remains pending.'
