#!/usr/bin/env bash
set -eo pipefail
export M1_T_SOLVERS=PBiCGStab
export M1_PHASE_MODES='equilibrium evaporation condensation'
export M1_AVERAGING_MODES=true
export M1_SUBCYCLES=1
export M1_OUTER_CORRECTORS=3
export M1_TIME_STEPS=1e-9
export M1_PROTECT_HISTORIES=true
export M1_ENERGY_AUDIT=true
unset M1_REFERENCE_FLAG
bash "$(dirname "$0")/m1-compressible-test.sh"
