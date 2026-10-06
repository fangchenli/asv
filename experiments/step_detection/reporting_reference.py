"""Conservative median evidence for independent, at-most-one-change histories.

This mathematical reference does not change ASV's detector or reporter.
See reporting_uncertainty.rst for the model, coverage proof, and limitations.
"""

import bisect
import math
from fractions import Fraction

MAX_POINTS = 200


def median_rank(n, tail):
    """Largest 1-based k with P(Binomial(n, 1/2) < k) <= tail; 0 if none.

    Integer sums and rational comparisons avoid rounding a tail probability
    across its requested bound. Rank 0 represents an unbounded endpoint.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError('n must be a positive integer')
    tail = Fraction(tail)
    if not 0 < tail < Fraction(1, 2):
        raise ValueError('tail must be between 0 and 1/2')
    count = 0
    rank = 0
    for k in range(1, n + 1):
        count += math.comb(n, k - 1)
        if count * tail.denominator > tail.numerator * (1 << n):
            break
        rank = k
    return rank


def _inputs(values, alpha, threshold):
    y = list(values)
    if not 1 <= len(y) <= MAX_POINTS:
        raise ValueError(f'Use 1 to {MAX_POINTS} retained observations')
    if any(not math.isfinite(value) for value in y):
        raise ValueError('Observations must be finite; remove missing readings explicitly')
    alpha = Fraction(alpha)
    if not 0 < alpha < 1:
        raise ValueError('alpha must be between 0 and 1')
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError('threshold must be finite and nonnegative')
    return y, alpha


def known_boundary_lower(values, split, *, alpha=0.05, threshold=0.05):
    """Lower bound on mu_after - (1 + threshold)*mu_before at a FIXED split.

    Two one-sided order-statistic bounds each spend alpha/2. The split must
    be chosen independently of these data. Positive population medians are
    assumed; using zero as a lower bound adds that parameter-space constraint.
    """
    y, alpha = _inputs(values, alpha, threshold)
    if not isinstance(split, int) or isinstance(split, bool) or not 0 < split < len(y):
        raise ValueError('split must leave observations on both sides')
    before, after = sorted(y[:split]), sorted(y[split:])
    k_before = median_rank(len(before), alpha / 2)
    k_after = median_rank(len(after), alpha / 2)
    upper_before = before[-k_before] if k_before else math.inf
    lower_after = max(0, after[k_after - 1]) if k_after else 0
    return lower_after - (1 + threshold) * upper_before


def single_change_evidence(values, *, alpha=0.05, threshold=0.05):
    """Invert simultaneous sign bounds over ALL splits and constant models.

    Target: final median minus (1 + threshold)*initial median in an independent
    history with positive medians and at most one change. No distribution,
    variance, or minimum segment length is fitted. This does not certify the
    boundary returned by a separate detector, or handle multiple true changes.
    """
    y, alpha = _inputs(values, alpha, threshold)
    n = len(y)
    intervals = n * (n + 1) // 2
    ranks = {m: median_rank(m, alpha / (2 * intervals)) for m in range(1, n + 1)}
    bounds = {}

    # C[l,r] is a simultaneous median interval if [l,r) has a common median.
    for left in range(n):
        ordered = []
        for right in range(left + 1, n + 1):
            bisect.insort(ordered, y[right - 1])
            k = ranks[right - left]
            bounds[left, right] = (max(0, ordered[k - 1]), ordered[-k]) if k else (0, math.inf)

    # J[l,r] intersects C over every subinterval contained in [l,r).
    # A hypothesized constant level must satisfy ALL those constraints.
    for length in range(2, n + 1):
        for left in range(n - length + 1):
            right = left + length
            children = [bounds[left, right], bounds[left + 1, right], bounds[left, right - 1]]
            bounds[left, right] = (
                max(item[0] for item in children),
                min(item[1] for item in children),
            )

    def feasible(bound):
        return bound[0] <= bound[1] and bound[1] > 0

    candidates = []
    constant = feasible(bounds[0, n])
    if constant:
        upper = bounds[0, n][1]
        candidates.append(0 if threshold == 0 else -threshold * upper)
    splits = []
    for split in range(1, n):
        before, after = bounds[0, split], bounds[split, n]
        if feasible(before) and feasible(after):
            splits.append(split)
            candidates.append(after[0] - (1 + threshold) * before[1])

    # An empty confidence set flags incompatibility; it must not imply an alert
    # through the mathematical convention inf(empty) = +infinity.
    lower = min(candidates) if candidates else None
    alert = lower is not None and lower > 0
    return {
        'status': 'incompatible_model'
        if lower is None
        else ('evidence' if alert else 'insufficient'),
        'has_alert': alert,
        'lower_excess': lower,
        'constant_feasible': constant,
        'feasible_splits': splits,
        'observations': n,
        'intervals': intervals,
    }
