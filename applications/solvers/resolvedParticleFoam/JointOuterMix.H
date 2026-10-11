// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LPBF_JOINT_OUTER_MIX_H
#define LPBF_JOINT_OUTER_MIX_H
#include <algorithm>
#include <cmath>
namespace lpbfM2
{
// Minimise ||R_previous + alpha (R_current-R_previous)||_2. The
// caller supplies global dot products and must recheck the actual candidate.
inline double jointOuterFactor(double previousDotDifference,double differenceSquared)
{
    if(!std::isfinite(previousDotDifference)||!std::isfinite(differenceSquared)
       ||differenceSquared<=0) return 1;
    const double alpha=-previousDotDifference/differenceSquared;
    return std::isfinite(alpha) ? std::clamp(alpha,1.0,4.0):1;
}
}
#endif
