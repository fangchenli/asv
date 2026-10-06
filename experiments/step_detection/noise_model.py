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


def correlation_cap(max_half_life):
    """Convert a maximum shock half-life in observations to a correlation cap.

    Zero denotes independent noise. Positive finite H gives 2**(-1/H).
    Reject bounds rounded to one, since they would lose strict contraction.
    """
    if not math.isfinite(max_half_life) or max_half_life < 0:
        raise ValueError('max_half_life must be finite and nonnegative')
    if max_half_life == 0:
        return 0.0
    cap = math.exp(-math.log(2) / max_half_life)
    if cap == 1:
        raise ValueError('max_half_life is too large to resolve a cap below one')
    return cap


def fit_correlated_score(residuals, weights, segments, noise_floor, beta, rho_max):
    """Score a supplied candidate with an explicit bound on noise persistence.

    residuals and weights are finite, same-length, nonempty sequences, and
    weights are positive. The levels and boundaries of the candidate are fixed.
    The production correlation helper supplies the exact conditional fit;
    production penalty selection does not use this experimental score.
    """
    from asv.step_detect import _fit_ar1

    n = len(residuals)
    if n == 0 or len(weights) != n:
        raise ValueError('residuals and weights must have the same nonzero length')
    if any(not math.isfinite(e) for e in residuals):
        raise ValueError('residuals must be finite')
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError('weights must be finite and positive')
    rho = _fit_ar1(residuals, weights, rho_max=rho_max)
    error_sum = weights[0] * abs(residuals[0]) + math.fsum(
        weight * abs(current - rho * previous)
        for previous, current, weight in zip(residuals, residuals[1:], weights[1:])
    )
    return {
        'rho': rho,
        'rho_max': rho_max,
        'error_sum': error_sum,
        'score': selection_score(error_sum, n, segments, noise_floor, beta),
    }
