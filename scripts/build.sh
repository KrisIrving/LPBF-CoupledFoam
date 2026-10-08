#!/usr/bin/env bash
set -eo pipefail
source "$(dirname "$0")/environment.sh"
cd "$COUPLED_ROOT"
./Allwmake -j "${BUILD_JOBS:-4}" 2>&1 | tee .build/build.log
for executable in laserbeamFoam laserMeltFoam compressibleLaserbeamFoam setSolidFraction; do
    test -x "$FOAM_USER_APPBIN/$executable" || { echo "Missing: $executable"; exit 1; }
done
