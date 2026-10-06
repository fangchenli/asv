"""Jointly calibrate interval median constraints using independent fair signs.

See reporting_protocol.rst for the calibration and coverage arguments.
This experiment assumes independent observations and at most one true change.
"""

import bisect
import math
import random
from fractions import Fraction
from functools import cache

from .reporting_reference import MAX_POINTS, _inputs


@cache
def score_table(n):
    """Integer ranks of exact two-sided binomial tails; larger is more extreme."""
    if not isinstance(n, int) or not 1 <= n <= MAX_POINTS:
        raise ValueError(f'n must be between 1 and {MAX_POINTS}')
    tails = {}
    for m in range(1, n + 1):
        total = 0
        for c in range(m // 2 + 1):
            total += math.comb(m, c)
            tails[m, c] = min(Fraction(1), Fraction(2 * total, 1 << m))
    probabilities = sorted(set(tails.values()), reverse=True)
    rank = {p: i for i, p in enumerate(probabilities)}
    table = [()]
    for m in range(1, n + 1):
        table.append(tuple(rank[tails[m, min(c, m - c)]] for c in range(m + 1)))
    return tuple(table)


def maximum_score(signs):
    """Largest interval score in a binary sequence, including all lengths."""
    signs = list(signs)
    if any(x not in (0, 1) for x in signs):
        raise ValueError('signs must be zero or one')
    table = score_table(len(signs))
    prefix = [0]
    for sign in signs:
        prefix.append(prefix[-1] + sign)
    largest = 0
    for m in range(1, len(signs) + 1):
        scores = table[m]
        for right in range(m, len(signs) + 1):
            largest = max(largest, scores[prefix[right] - prefix[right - m]])
    return largest


def tolerance_rank(draws, alpha, failure):
    """Choose k with P(Binomial(draws, 1-alpha) >= k) <= failure, exactly.

    For the kth sorted calibration maximum, this bounds the probability
    (over calibration) that its conditional future exceedance risk is > alpha.
    """
    if not isinstance(draws, int) or draws < 1:
        raise ValueError('draws must be a positive integer')
    alpha, failure = Fraction(alpha), Fraction(failure)
    if not 0 < alpha < 1 or not 0 < failure < 1:
        raise ValueError('alpha and failure must be between zero and one')
    p = 1 - alpha
    a, b = p.numerator, p.denominator
    denominator = b**draws
    term = a**draws
    total = 0
    best, best_tail = draws + 1, Fraction(0)
    for k in range(draws, 0, -1):
        total += term
        if total * failure.denominator > denominator * failure.numerator:
            break
        best, best_tail = k, Fraction(total, denominator)
        term = term * k * (b - a) // (a * (draws - k + 1))
    return best, best_tail


def calibrate(n, *, draws, seed, alpha=Fraction(1, 20), failure=Fraction(1, 200)):
    """Calibrate using sign sequences only, never benchmark histories."""
    table = score_table(n)
    k, tail = tolerance_rank(draws, alpha, failure)
    rng = random.Random(seed)
    maxima = [maximum_score([rng.getrandbits(1) for _ in range(n)]) for _ in range(draws)]
    # If no finite order statistic can meet the tolerance requirement, accepting
    # every sign sequence is the safe fallback.
    critical = sorted(maxima)[k - 1] if k <= draws else max(max(row) for row in table[1:])
    return {
        'n': n,
        'draws': draws,
        'seed': seed,
        'alpha': str(Fraction(alpha)),
        'calibration_failure': str(Fraction(failure)),
        'order_rank': k,
        'order_rank_failure_bound': float(tail),
        'critical_score': critical,
        'maxima': maxima,
    }


def interval_ranks(n, critical):
    """Invert score <= critical into symmetric order-statistic endpoints."""
    table = score_table(n)
    if not isinstance(critical, int) or critical < 0:
        raise ValueError('critical must be a nonnegative integer score')
    return {m: sum(score > critical for score in table[m][: m // 2 + 1]) for m in range(1, n + 1)}


def evidence(values, calibration, *, threshold=0.05):
    """History-level evidence with all boundaries retained; no location claim.

    Calibration must be generated independently for exactly this retained length.
    Both sign orientations use the same two-sided interval score.
    """
    y, _ = _inputs(values, Fraction(calibration['alpha']), threshold)
    n = len(y)
    if calibration['n'] != n:
        raise ValueError('Calibration length must equal retained observation count')
    ranks = interval_ranks(n, calibration['critical_score'])
    bounds = {}
    for left in range(n):
        ordered = []
        for right in range(left + 1, n + 1):
            bisect.insort(ordered, y[right - 1])
            k = ranks[right - left]
            bounds[left, right] = (max(0, ordered[k - 1]), ordered[-k]) if k else (0, math.inf)
    for size in range(2, n + 1):
        for left in range(n - size + 1):
            right = left + size
            parts = [bounds[left, right], bounds[left + 1, right], bounds[left, right - 1]]
            bounds[left, right] = max(p[0] for p in parts), min(p[1] for p in parts)

    def feasible(bound):
        return bound[0] <= bound[1] and bound[1] > 0

    constant = feasible(bounds[0, n])
    lower = []
    if constant:
        lower.append(-threshold * bounds[0, n][1] if threshold else 0)
    splits = []
    for split in range(1, n):
        before, after = bounds[0, split], bounds[split, n]
        if feasible(before) and feasible(after):
            splits.append(split)
            lower.append(after[0] - (1 + threshold) * before[1])
    bound = min(lower) if lower else None
    alert = bound is not None and bound > 0
    return {
        'status': 'incompatible_model'
        if bound is None
        else ('evidence' if alert else 'insufficient'),
        'has_alert': alert,
        'lower_excess': bound,
        'constant_feasible': constant,
        'feasible_splits': splits,
        'observations': n,
    }
