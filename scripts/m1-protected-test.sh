#!/usr/bin/env bash
# Regression of the candidate, plus native thermo-energy observations.
set -eo pipefail
export M1_PROTECT_HISTORIES=true
bash "$(dirname "$0")/m1-refinement-test.sh"
