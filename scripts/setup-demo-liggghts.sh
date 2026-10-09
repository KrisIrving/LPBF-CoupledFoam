#!/usr/bin/env bash
# User-side dependency build. No CFDEM install, sudo or system-wide changes.
set -eo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
report="$repo/.runs/dem-setup-$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$report"
stage=environment
finish() {
    rc=$?
    trap - EXIT
    printf 'last_stage=%s\nexit_status=%s\n' "$stage" "$rc" > "$report/result.txt"
    if [ -f "$report/log.communication-wrapper" ]; then
        child="$(sed -n 's/^Feedback archive: //p' "$report/log.communication-wrapper" | tail -n 1)"
        case "$child" in
            "$repo"/.runs/dem-communication-*.tar.gz)
                if [ -f "$child" ]; then cp "$child" "$report/" || rc=1; fi ;;
        esac
    fi
    tar -czf "$report.tar.gz" -C "$(dirname "$report")" "$(basename "$report")" || rc=1
    echo "Setup feedback archive: $report.tar.gz"
    exit "$rc"
}
trap finish EXIT
cd "$repo"
git rev-parse HEAD > "$report/commit.txt"
source scripts/environment.sh > "$report/log.environment" 2>&1
if [ "$WM_PROJECT_VERSION" != v2512 ] || [ "${WM_MPLIB:-}" != SYSTEMOPENMPI ]; then
    echo 'Load OpenCFD v2512 with SYSTEMOPENMPI before setup.'
    exit 1
fi
stage=prerequisites
{
    for tool in git make g++ mpicxx mpirun python3 nm ldd; do command -v "$tool"; done
    mpicxx --showme
    mpirun --version | head -n 3
    g++ --version | head -n 1
    printf 'DEPS_JOBS=%s\n' "${DEPS_JOBS:-8}"
} > "$report/environment.txt" 2>&1
if ! [[ "${DEPS_JOBS:-8}" =~ ^[1-9][0-9]*$ ]]; then echo 'DEPS_JOBS must be positive.'; exit 1; fi
pin=3d5c00f20519e6bb6eb6756f51f1ad36564e649d
dep_root="${LPBF_DEPS_ROOT:-$HOME/LPBF-dependencies}"
mkdir -p "$dep_root"
dep_root="$(realpath "$dep_root")"
if [[ "$dep_root" =~ [[:space:]] ]]; then echo 'Use a dependency path without whitespace.'; exit 1; fi
tree="$dep_root/LIGGGHTS-PUBLIC-${pin:0:12}"
stage=checkout
if [ ! -e "$tree" ]; then
    echo 'Fetching pinned LIGGGHTS-PUBLIC into an independent user directory.'
    git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git "$tree" > "$report/log.checkout" 2>&1
    git -C "$tree" checkout --detach "$pin" >> "$report/log.checkout" 2>&1
fi
if [ ! -f "$tree/src/library.h" ] || [ "$(git -C "$tree" rev-parse HEAD)" != "$pin" ]; then
    echo 'Existing target is not the expected pinned source; use a new LPBF_DEPS_ROOT.'
    exit 1
fi
# Never reset, clean or overwrite an unrelated checkout or local source edits.
if [ -n "$(git -C "$tree" diff --name-only)" ] || [ -n "$(git -C "$tree" diff --cached --name-only)" ]; then
    echo 'Pinned dependency has tracked source edits; inspect them or use a new LPBF_DEPS_ROOT.'
    exit 1
fi
stage=configuration
cd "$tree/src"
# Official postfix option isolates our settings and object directory.
cp MAKE/Makefile.user_default MAKE/Makefile.user_lpbf
cat >> MAKE/Makefile.user_lpbf <<'EOF'

# Minimal spherical MPI communication fixture; no optional visualisation libraries.
USE_MPI = "ON"
USE_FPIC = "ON"
USE_VTK = "OFF"
USE_CATALYST = "OFF"
USE_SUPERQUADRICS = "OFF"
USE_JPG = "OFF"
USE_GZIP = "OFF"
USE_CONVEX = "OFF"
USE_MFEM = "OFF"
AUTOINSTALL_VTK = "OFF"
AUTOINSTALL_CONVEX = "OFF"
AUTOINSTALL_MFEM = "OFF"
EOF
printf 'MPICXX_USR=%s\n' "$(command -v mpicxx)" >> MAKE/Makefile.user_lpbf
cp MAKE/Makefile.user_lpbf "$report/Makefile.user_lpbf.txt"
signature="$(command -v mpicxx)|$(mpicxx --showme)|$(g++ --version | head -n 1)"
if [ -f .lpbf-toolchain.signature ] && [ "$(cat .lpbf-toolchain.signature)" != "$signature" ]; then
    echo 'Toolchain changed. Use a new LPBF_DEPS_ROOT rather than mixing objects.'
    exit 1
fi
printf '%s\n' "$signature" > .lpbf-toolchain.signature
make makeshlib > "$report/log.configure" 2>&1
bash Make.sh models >> "$report/log.configure" 2>&1
stage=shared-library-build
echo "Building MPI shared library with ${DEPS_JOBS:-8} jobs; logs are retained."
make -f Makefile.shlib auto postfix=lpbf -j "${DEPS_JOBS:-8}" > "$report/log.build-liggghts" 2>&1
export LIGGGHTS_SRC="$tree/src"
export LIGGGHTS_LIB="$(realpath "$tree/src/libliggghts.so")"
test -f "$LIGGGHTS_LIB"
{
    printf 'source_commit=%s\nLIGGGHTS_SRC=%s\nLIGGGHTS_LIB=%s\n' "$pin" "$LIGGGHTS_SRC" "$LIGGGHTS_LIB"
    sha256sum "$LIGGGHTS_SRC/library.h" "$LIGGGHTS_LIB"
    ldd "$LIGGGHTS_LIB"
} > "$report/dependency.txt" 2>&1
cd "$repo"
mkdir -p .build
printf 'export LIGGGHTS_SRC=%q\nexport LIGGGHTS_LIB=%q\n' \
    "$LIGGGHTS_SRC" "$LIGGGHTS_LIB" > .build/demo-liggghts.env
stage=communication-demo
echo 'Shared library built; running the existing serial/MPI communication package.'
bash scripts/dem-communication-test.sh 2>&1 | tee "$report/log.communication-wrapper"
stage=complete
echo 'Dependency build and communication checks completed; full CFD coupling remains pending.'
