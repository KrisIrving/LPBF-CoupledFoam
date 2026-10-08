#!/usr/bin/env bash
set -eo pipefail
export M1_REFERENCE_ENERGY='false true'
bash "$(dirname "$0")/m1-consistent-test.sh"
