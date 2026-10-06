"""Exact small-history reference for the independent shared-floor score.

Compute minimum weighted absolute error at every feasible segment count,
then select the count with the declared floor and complexity penalty.
No gamma search, approximation, or data-dependent calibration is used.
"""

import math

from .harness import load_detector
from .noise_model import selection_score

MAX_POINTS = 200


def fit_independent(y, w, noise_floor, beta, *, min_size=1, max_size=None, backend='python'):
    """Return the global independent-score optimum and the segment-count frontier.

    Inputs are prepared, finite observations and positive finite weights.
    The floor uses these weights' scale; no filtering or normalization occurs.
    Segment endpoints are exclusive and refer to this input sequence.
    Each feasible count gets a fit, its independent error, and its score.

    The dynamic program takes O(n**3) time and O(n**2) space; interval
    preprocessing also depends on the selected ASV backend. Limit inputs
    to MAX_POINTS so this correctness reference is not used for long histories.
    Exactness is combinatorial, with floating-point interval costs and scores.
    """
    n = len(y)
    if n == 0 or len(w) != n:
        raise ValueError('y and w must have the same nonzero length')
    if n > MAX_POINTS:
        raise ValueError(f'The exact reference is limited to {MAX_POINTS} observations')
    if any(not math.isfinite(value) for value in y):
        raise ValueError('observations must be finite')
    if any(not math.isfinite(weight) or weight <= 0 for weight in w):
        raise ValueError('weights must be finite and positive')
    if max_size is None:
        max_size = n
    for name, value in [('min_size', min_size), ('max_size', max_size)]:
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f'{name} must be a positive integer')
    if min_size > max_size:
        raise ValueError('min_size must not exceed max_size')
    if min_size > n:
        raise ValueError('No feasible partition for the segment-size constraints')
    # Validate the declared model before building interval tables.
    selection_score(0, n, 1, noise_floor, beta)

    module = load_detector(backend)
    # The native interval backend requires lists, even for tuple inputs.
    mu_dist = module.get_mu_dist(list(y), list(w))
    intervals = {}
    for left in range(n):
        for right in range(left + min_size, min(n, left + max_size) + 1):
            level = mu_dist.mu(left, right - 1)
            cost = mu_dist.dist(left, right - 1)
            if not math.isfinite(cost) or cost < 0 or not math.isfinite(level):
                raise ValueError('Interval levels and nonnegative costs must be representable')
            intervals[left, right] = (level, cost)

    # E[k,t] is the minimum error for t observations in exactly k segments.
    previous = [0.0] + [math.inf] * n
    parents = [[None] * (n + 1)]
    frontier = []
    for count in range(1, n // min_size + 1):
        current = [math.inf] * (n + 1)
        parent = [None] * (n + 1)
        for right in range(count * min_size, min(n, count * max_size) + 1):
            first = max((count - 1) * min_size, right - max_size)
            last = min((count - 1) * max_size, right - min_size)
            for left in range(first, last + 1):
                if previous[left] == math.inf:
                    continue
                cost = previous[left] + intervals[left, right][1]
                if not math.isfinite(cost):
                    raise ValueError('Accumulated interval costs must be representable')
                # Strict comparison retains the earliest start when costs tie.
                if cost < current[right]:
                    current[right] = cost
                    parent[right] = left
        parents.append(parent)

        if current[n] != math.inf:
            right = n
            rights, values, costs = [], [], []
            for k in range(count, 0, -1):
                left = parents[k][right]
                level, cost = intervals[left, right]
                rights.append(right)
                values.append(level)
                costs.append(cost)
                right = left
            rights.reverse()
            values.reverse()
            costs.reverse()
            score = selection_score(current[n], n, count, noise_floor, beta)
            if not math.isfinite(score):
                raise ValueError('Selection scores must be representable')
            frontier.append(
                {
                    'segments': count,
                    'right': rights,
                    'values': values,
                    'costs': costs,
                    'error_sum': current[n],
                    'score': score,
                }
            )
        previous = current

    if not frontier:
        raise ValueError('No feasible partition for the segment-size constraints')
    return {
        'selected': min(frontier, key=lambda fit: (fit['score'], fit['segments'])),
        'frontier': frontier,
        'noise_floor': noise_floor,
        'beta': beta,
        'size': n,
        'min_size': min_size,
        'max_size': max_size,
        'backend': backend,
    }
