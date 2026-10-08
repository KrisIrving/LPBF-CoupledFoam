#!/usr/bin/env bash
# Source this file; OpenFOAM setup is not compatible with nounset.
COUPLED_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if source "${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2506/etc/bashrc}"; then
    :
else
    echo 'Failed to load OpenFOAM'; return 1
fi
test "${WM_PROJECT:-}" = OpenFOAM || { echo 'OpenFOAM environment unavailable'; return 1; }
export FOAM_USER_APPBIN="$COUPLED_ROOT/.build/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$COUPLED_ROOT/.build/$WM_OPTIONS/lib"
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:${LD_LIBRARY_PATH:-}"
# WSL graphics probing may block CPU-only OpenMPI startup on an X server.
if grep -qi microsoft /proc/sys/kernel/osrelease; then
    export HWLOC_COMPONENTS="${HWLOC_COMPONENTS:--gl}"
fi
mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN"
