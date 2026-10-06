"""Test every below-threshold one-change explanation under Gaussian noise.

Independent observations and one common positive noise variance are assumed.
See direct_protocol.rst for the finite-sample argument and model limitations.
"""

import math


def critical_values(n, *, alpha=0.05, size_share=0.8, threshold=0.05):
    """Analytical critical values; no benchmark data enter calibration."""
    from scipy.stats import f, t

    if not isinstance(n, int) or isinstance(n, bool) or not 4 <= n <= 200:
        raise ValueError('n must be an integer from 4 to 200')
    if not 0 < alpha < 1 or not 0 < size_share < 1:
        raise ValueError('alpha and size_share must be between zero and one')
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError('threshold must be finite and nonnegative')
    size_alpha = alpha * size_share
    fit_alpha = alpha * (1 - size_share)
    return {
        'n': n,
        'alpha': alpha,
        'size_share': size_share,
        'threshold': threshold,
        'size_alpha': size_alpha,
        'fit_alpha': fit_alpha,
        'size_t_critical': float(t.isf(size_alpha, n - 2)),
        'fit_f_critical': float(f.isf(fit_alpha / (n - 2), 1, n - 3)),
    }


def interval_moments(values):
    """Means and squared residual sums, using stable incremental updates."""
    means, costs = {}, {}
    for left in range(len(values)):
        mean = total = 0.0
        for right in range(left + 1, len(values) + 1):
            value = values[right - 1]
            delta = value - mean
            mean += delta / (right - left)
            total += delta * (value - mean)
            means[left, right] = mean
            costs[left, right] = max(0.0, total)
    return means, costs


def _ratio(numerator, denominator):
    if denominator:
        return numerator / denominator
    return math.copysign(math.inf, numerator) if numerator else 0.0


def evidence(values, calibration):
    """Reject only if every fixed-boundary null is rejected.

    At each split, either a one-sided t contrast rejects the size null, or
    a Bonferroni-corrected extra-boundary F scan rejects its two-level shape.
    Array entry i corresponds to split i+1, in retained-observation coordinates.
    The surviving null splits are not a confidence interval for the boundary.
    """
    y = list(values)
    n = len(y)
    if not 4 <= n <= 200 or calibration['n'] != n:
        raise ValueError('Use 4 to 200 observations with matching calibration length')
    if any(not math.isfinite(x) for x in y):
        raise ValueError('Observations must be finite')
    factor = 1 + calibration['threshold']
    if not math.isfinite(factor) or factor < 1:
        raise ValueError('threshold must be finite and nonnegative')
    for key in ('size_t_critical', 'fit_f_critical'):
        if not math.isfinite(calibration[key]) or calibration[key] <= 0:
            raise ValueError('Critical values must be finite and positive')
    # Ratios are invariant to the timing unit. Scaling avoids square overflow.
    scale = max(map(abs, y)) or 1.0
    means, costs = interval_moments([x / scale for x in y])
    t_scores, f_scores, extra_splits, surviving = [], [], [], []
    size_rejected, fit_rejected = [], []
    for split in range(1, n):
        two = costs[0, split] + costs[split, n]
        excess = means[split, n] - factor * means[0, split]
        variance = two / (n - 2) * (1 / (n - split) + factor**2 / split)
        size_t = _ratio(excess, math.sqrt(variance))
        best_cost, best_split = math.inf, None
        for extra in range(1, n):
            if extra == split:
                continue
            left, right = sorted((extra, split))
            three = costs[0, left] + costs[left, right] + costs[right, n]
            if three < best_cost:
                best_cost, best_split = three, extra
        fit_f = _ratio(max(0.0, two - best_cost), best_cost / (n - 3))
        size_reject = size_t > calibration['size_t_critical']
        fit_reject = fit_f > calibration['fit_f_critical']
        if size_reject:
            size_rejected.append(split)
        if fit_reject:
            fit_rejected.append(split)
        if not (size_reject or fit_reject):
            surviving.append(split)
        t_scores.append(size_t)
        f_scores.append(fit_f)
        extra_splits.append(best_split)
    return {
        'has_alert': not surviving,
        'null_splits': surviving,
        'size_rejected_splits': size_rejected,
        'fit_rejected_splits': fit_rejected,
        'size_t': t_scores,
        'lack_of_fit_f': f_scores,
        'best_extra_split': extra_splits,
    }
