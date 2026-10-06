"""Shared-floor Laplace scoring for mathematical checks and later comparisons.

The caller supplies a positive lower bound on the common noise scale, in the
units of weighted residuals. No floor estimator or default calibration is
implied. This module does not change ASV's production scoring rule.
"""

import math


def profile_loss(error_sum, n, noise_floor):
    """Return the profiled per-observation loss above the exact-fit baseline.

    Minimize log(b) + error_sum/(n*b) over b >= noise_floor, then subtract
    log(noise_floor), a constant shared by every candidate for this history.
    For t = error_sum/(n*noise_floor), the result is t if t <= 1, and
    1 + log(t) otherwise. Log-space evaluation avoids overflowing t or
    the product n*noise_floor.

    error_sum can be independent weighted absolute error, or the minimized
    conditional AR(1) error. The latter does not give an independent-error
    Potts path a global optimization guarantee.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError('n must be a positive integer')
    if not math.isfinite(error_sum) or error_sum < 0:
        raise ValueError('error_sum must be finite and nonnegative')
    if not math.isfinite(noise_floor) or noise_floor <= 0:
        raise ValueError('noise_floor must be finite and positive')
    if error_sum == 0:
        return 0.0

    log_ratio = math.log(error_sum) - math.log(n) - math.log(noise_floor)
    return math.exp(log_ratio) if log_ratio <= 0 else 1 + log_ratio


def selection_score(error_sum, n, segments, noise_floor, beta):
    """Add a declared per-segment penalty to the shared-floor profile loss."""
    loss = profile_loss(error_sum, n, noise_floor)
    if not isinstance(segments, int) or isinstance(segments, bool) or not 1 <= segments <= n:
        raise ValueError('segments must be an integer between 1 and n')
    if not math.isfinite(beta) or beta < 0:
        raise ValueError('beta must be finite and nonnegative')
    return beta * segments + loss
