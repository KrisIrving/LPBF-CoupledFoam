#!/usr/bin/env bash
# Two otherwise identical cases; do not loosen the T tolerance.
set -eo pipefail
export M1_T_SOLVERS='upstream PBiCGStab'
export M1_AVERAGING_MODES=false
export M1_SUBCYCLES=1
export M1_REFERENCE_FLAG=--allow-upstream-reference-failure
bash "$(dirname "$0")/m1-compressible-test.sh"
