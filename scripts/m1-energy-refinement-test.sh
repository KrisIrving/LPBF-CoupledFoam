#!/usr/bin/env bash
# Thirteen small cases: one equilibrium reference and twelve exchange runs.
set -eo pipefail
export M1_T_SOLVERS=PBiCGStab
export M1_PHASE_MODES='equilibrium evaporation condensation'
export M1_THERMO_FIXTURE=consistent
export M1_COMMON_LATENT=true
export M1_REFERENCE_ENERGY=true
export M1_AVERAGING_MODES=true
export M1_SUBCYCLES='1 2 4'
export M1_OUTER_CORRECTORS=5
export M1_TIME_STEPS='1e-10 5e-11'
export M1_END_TIME=1e-9
export M1_PROTECT_HISTORIES=true
export M1_ENERGY_AUDIT=true
export M1_ENERGY_GATE=true
export M1_SINGLE_EQUILIBRIUM=true
unset M1_REFERENCE_FLAG
bash "$(dirname "$0")/m1-compressible-test.sh"
