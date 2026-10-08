#!/usr/bin/env bash
# Same 10 ns duration; compare temporal/subcycle sensitivity after coupling gate.
set -eo pipefail
export M1_T_SOLVERS=PBiCGStab
export M1_PHASE_MODES='evaporation condensation'
export M1_AVERAGING_MODES=true
export M1_SUBCYCLES='1 2 4'
export M1_OUTER_CORRECTORS=3
export M1_TIME_STEPS='1e-9 5e-10'
unset M1_REFERENCE_FLAG
bash "$(dirname "$0")/m1-compressible-test.sh"
