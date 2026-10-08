#!/usr/bin/env bash
set -eo pipefail
export M1_COMMON_LATENT='false true'
bash "$(dirname "$0")/m1-energy-test.sh"
