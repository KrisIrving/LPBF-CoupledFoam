#!/usr/bin/env bash
set -eo pipefail
export M1_T_SOLVERS=PBiCGStab
export M1_PHASE_MODES='evaporation condensation'
export M1_AVERAGING_MODES='false true'
export M1_SUBCYCLES='1 2 4'
unset M1_REFERENCE_FLAG
bash "$(dirname "$0")/m1-compressible-test.sh"
