#!/usr/bin/env bash
set -eo pipefail
source "$(dirname "$0")/environment.sh"
cd "$COUPLED_ROOT"
signature="$WM_PROJECT_VERSION|$WM_PROJECT_DIR|$WM_OPTIONS"
if [ -f .build/environment.signature ] && [ "$(cat .build/environment.signature)" != "$signature" ]; then
    echo 'Build environment changed. Use a fresh checkout for each OpenFOAM version.'
    exit 1
fi
printf '%s\n' "$signature" > .build/environment.signature
./Allwmake -j "${BUILD_JOBS:-4}" 2>&1 | tee .build/build.log
for executable in laserbeamFoam laserMeltFoam compressibleLaserbeamFoam setSolidFraction; do
    test -x "$FOAM_USER_APPBIN/$executable" || { echo "Missing: $executable"; exit 1; }
done
