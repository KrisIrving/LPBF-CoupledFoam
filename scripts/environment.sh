#!/usr/bin/env bash
# Source this file after selecting an OpenFOAM environment, or set OPENFOAM_BASHRC.
COUPLED_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -n "${OPENFOAM_BASHRC:-}" ]; then
    if source "$OPENFOAM_BASHRC"; then :; else
        echo 'Failed to load OPENFOAM_BASHRC'; return 1
    fi
elif [ "${WM_PROJECT:-}" != OpenFOAM ] || [ -z "${WM_PROJECT_DIR:-}" ]; then
    if [ -f /usr/lib/openfoam/openfoam2506/etc/bashrc ]; then
        if source /usr/lib/openfoam/openfoam2506/etc/bashrc; then :; else
            echo 'Failed to load OpenFOAM v2506'; return 1
        fi
    else
        echo 'Load OpenFOAM first (for example: of2512), or set OPENFOAM_BASHRC.'
        return 1
    fi
fi
test "${WM_PROJECT:-}" = OpenFOAM || { echo 'OpenFOAM environment unavailable'; return 1; }
test -n "${WM_OPTIONS:-}" && test -n "${WM_PROJECT_VERSION:-}" || {
    echo 'Incomplete OpenFOAM environment'; return 1;
}
export FOAM_USER_APPBIN="$COUPLED_ROOT/.build/$WM_PROJECT_VERSION/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$COUPLED_ROOT/.build/$WM_PROJECT_VERSION/$WM_OPTIONS/lib"
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:${LD_LIBRARY_PATH:-}"
# WSL graphics probing may block CPU-only OpenMPI startup on an X server.
if grep -qi microsoft /proc/sys/kernel/osrelease; then
    export HWLOC_COMPONENTS="${HWLOC_COMPONENTS:--gl}"
fi
mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN"
