#!/usr/bin/env bash
# User-side only. Does not rebuild LIGGGHTS or run the LPBF solvers.
set -eo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
report="$repo/.runs/dem-communication-$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$report"
stage=environment
finish() {
    rc=$?
    trap - EXIT
    printf 'last_stage=%s\nexit_status=%s\n' "$stage" "$rc" > "$report/result.txt"
    tar -czf "$report.tar.gz" -C "$(dirname "$report")" "$(basename "$report")" || rc=1
    echo "Feedback archive: $report.tar.gz"
    exit "$rc"
}
trap finish EXIT
cd "$repo"
git rev-parse HEAD > "$report/commit.txt"
git status --short > "$report/working-tree.txt"
source "$repo/scripts/environment.sh" > "$report/log.environment" 2>&1
if [ "$WM_PROJECT_VERSION" != v2512 ]; then
    echo 'This demo currently targets OpenCFD v2512.' | tee -a "$report/log.environment"
    exit 1
fi
stage=dependency-audit
{
    printf 'WM_PROJECT_VERSION=%s\nWM_OPTIONS=%s\nWM_MPLIB=%s\n' \
        "$WM_PROJECT_VERSION" "$WM_OPTIONS" "${WM_MPLIB:-}"
    command -v wmake mpirun mpicc python3 nm ldd
    mpirun --version | head -n 3
    mpicc --showme
} > "$report/environment.txt" 2>&1
if [ -z "${LIGGGHTS_SRC:-}" ]; then
    for candidate in "$HOME/CFDEM/LIGGGHTS/src" "$HOME/LIGGGHTS-PUBLIC/src"; do
        if [ -f "$candidate/library.h" ]; then LIGGGHTS_SRC="$candidate"; break; fi
    done
fi
if [ ! -f "${LIGGGHTS_SRC:-}/library.h" ]; then
    echo 'Set LIGGGHTS_SRC to the source directory containing library.h.' | tee "$report/dependency-error.txt"
    exit 1
fi
LIGGGHTS_SRC="$(realpath "$LIGGGHTS_SRC")"
if [ -z "${LIGGGHTS_LIB:-}" ]; then
    for candidate in "$LIGGGHTS_SRC/libliggghts.so" \
        "$(dirname "$LIGGGHTS_SRC")/src-build/libliggghts.so" \
        "$(dirname "$LIGGGHTS_SRC")/src-build/liblammps.so" \
        "$(dirname "$LIGGGHTS_SRC")/build/libliggghts.so"; do
        if [ -f "$candidate" ]; then LIGGGHTS_LIB="$candidate"; break; fi
    done
fi
if [ ! -f "${LIGGGHTS_LIB:-}" ]; then
    echo 'Shared library missing. Set LIGGGHTS_LIB to its full .so path; the executable alone is insufficient. See docs/dem-communication-demo.zh-CN.md.' | tee "$report/dependency-error.txt"
    exit 1
fi
LIGGGHTS_LIB="$(realpath "$LIGGGHTS_LIB")"
# wmake whitespace in external paths is deliberately unsupported in this demo.
if [[ "$LIGGGHTS_SRC$LIGGGHTS_LIB" =~ [[:space:]] ]]; then
    echo 'Use dependency paths without whitespace.' | tee "$report/dependency-error.txt"
    exit 1
fi
{
    printf 'LIGGGHTS_SRC=%s\nLIGGGHTS_LIB=%s\n' "$LIGGGHTS_SRC" "$LIGGGHTS_LIB"
    sha256sum "$LIGGGHTS_SRC/library.h" "$LIGGGHTS_LIB"
    git -C "$LIGGGHTS_SRC" rev-parse HEAD || true
    git -C "$LIGGGHTS_SRC" status --short || true
} > "$report/dependency.txt" 2>&1
ldd "$LIGGGHTS_LIB" > "$report/ldd.dem.txt" 2>&1
if grep -q 'not found' "$report/ldd.dem.txt"; then
    echo 'Unresolved shared-library dependencies; inspect ldd.dem.txt.' | tee "$report/dependency-error.txt"
    exit 1
fi
nm -D --defined-only "$LIGGGHTS_LIB" > "$report/symbols.txt"
for symbol in lammps_open lammps_close lammps_file lammps_command lammps_free \
    lammps_extract_atom lammps_extract_global lammps_extract_variable \
    lammps_gather_atoms lammps_get_natoms; do
    if ! grep -Eq "[[:space:]]${symbol}$" "$report/symbols.txt"; then
        echo "Missing library symbol: $symbol" | tee "$report/dependency-error.txt"
        exit 1
    fi
done
if ! grep -Eq 'libmpi\.so' "$report/ldd.dem.txt"; then
    echo 'An MPI-enabled shared library is required; serial/MPI-stub libraries are unsupported.' | tee "$report/dependency-error.txt"
    exit 1
fi
linkdir="$repo/.build/$WM_PROJECT_VERSION/$WM_OPTIONS/dem-link"
mkdir -p "$linkdir"
ln -sfn "$LIGGGHTS_LIB" "$linkdir/libliggghts.so"
export LIGGGHTS_SRC LIGGGHTS_LIBDIR="$linkdir"
export LD_LIBRARY_PATH="$linkdir:$(dirname "$LIGGGHTS_LIB"):$LD_LIBRARY_PATH"
stage=build
signature="$WM_PROJECT_DIR|$WM_OPTIONS|$LIGGGHTS_SRC|$LIGGGHTS_LIB|$(sha256sum "$LIGGGHTS_SRC/library.h" | awk '{print $1}')"
if [ -f "$linkdir/environment.signature" ] && [ "$(cat "$linkdir/environment.signature")" != "$signature" ]; then
    (cd applications/utilities/foamDemCommunicationDemo && wclean) > "$report/log.clean-demo" 2>&1
fi
printf '%s\n' "$signature" > "$linkdir/environment.signature"
echo 'Building only the OpenFOAM communication demo.'
(cd applications/utilities/foamDemCommunicationDemo && wmake) > "$report/log.build" 2>&1
app="$FOAM_USER_APPBIN/foamDemCommunicationDemo"
test -x "$app"
ldd "$app" > "$report/ldd.demo.txt"
if grep -q 'not found' "$report/ldd.demo.txt"; then exit 1; fi
# Reject an obvious OpenFOAM/LIGGGHTS MPI library-path mismatch.
dem_mpi="$(awk '$1 ~ /^libmpi\.so/ {print $3}' "$report/ldd.dem.txt")"
app_mpi="$(awk '$1 ~ /^libmpi\.so/ {print $3}' "$report/ldd.demo.txt")"
if [ -z "$app_mpi" ] || [ "$(realpath "$dem_mpi")" != "$(realpath "$app_mpi")" ]; then
    echo 'MPI library mismatch between DEM and demo; use one MPI toolchain.' | tee "$report/dependency-error.txt"
    exit 1
fi
failures=0
for mode in serial mpi; do
    stage="$mode"
    case_dir="$report/$mode"
    mkdir -p "$case_dir"
    cp -R "$repo/tutorials/communicationDemo/." "$case_dir/"
    if [ "$mode" = mpi ]; then
        mkdir -p "$case_dir/processor0" "$case_dir/processor1"
        if (cd "$case_dir" && mpirun -np 2 "$app" -parallel) > "$case_dir/log.demo" 2>&1; then
            printf 'exit_status=0\n' > "$case_dir/result.txt"
        else
            printf 'exit_status=1\n' > "$case_dir/result.txt"; failures=$((failures+1))
        fi
    else
        if (cd "$case_dir" && "$app") > "$case_dir/log.demo" 2>&1; then
            printf 'exit_status=0\n' > "$case_dir/result.txt"
        else
            printf 'exit_status=1\n' > "$case_dir/result.txt"; failures=$((failures+1))
        fi
    fi
done
stage=summary
python3 "$repo/scripts/summarize-dem-communication.py" "$report" --output "$report/summary.json"
test "$failures" -eq 0
stage=complete
echo 'Communication gates passed. Resolved CFD, heat, torque feedback and restart remain unvalidated.'
