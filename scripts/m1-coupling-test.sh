#!/usr/bin/env bash
set -eo pipefail
export M1_T_SOLVERS=PBiCGStab
export M1_PHASE_MODES='evaporation condensation'
export M1_AVERAGING_MODES=true
export M1_SUBCYCLES=1
export M1_OUTER_CORRECTORS='1 3'
unset M1_REFERENCE_FLAG
bash "$(dirname "$0")/m1-compressible-test.sh"
