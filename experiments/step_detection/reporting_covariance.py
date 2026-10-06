"""Direct Gaussian threshold test with a supplied covariance shape.

Cov(error) = sigma**2 * covariance; only sigma is estimated from this history.
See known_covariance.rst for both tests, their joint bound, and numerical limits.
"""

import math

import numpy as np
from scipy.linalg import solve_triangular


class _InsufficientPrecision(Exception):
    pass


def _whiten(values, covariance):
    """Normalize units, validate the shape matrix, and transform data/design."""
    n = len(values)
    matrix = np.asarray(covariance, dtype=float)
    if matrix.shape != (n, n) or not np.isfinite(matrix).all():
        raise ValueError('Covariance must be a finite n-by-n matrix')
    scale = float(np.max(np.abs(matrix)))
    if scale == 0:
        raise ValueError('Covariance must be positive definite')
    matrix = matrix / scale
    if np.max(np.abs(matrix - matrix.T)) > 64 * np.finfo(float).eps:
        raise ValueError('Covariance must be symmetric')
    matrix = (matrix + matrix.T) / 2
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues[0] <= 0:
        raise ValueError('Covariance must be positive definite')
    condition = float(eigenvalues[-1] / eigenvalues[0])
    if condition > 1e12:
        raise ValueError('Covariance condition number must not exceed 1e12')
    lower = np.linalg.cholesky(matrix)
    # Column zero is the constant; column t is the indicator for positions >= t.
    steps = (np.arange(n)[:, None] >= np.arange(n)[None, :]).astype(float)
    y = values / (float(np.max(np.abs(values))) or 1.0)
    return (
        solve_triangular(lower, y, lower=True),
        solve_triangular(lower, steps, lower=True),
        condition,
    )


def _fit(design, values, resolution):
    basis, triangular = np.linalg.qr(design, mode='reduced')
    coefficients = solve_triangular(triangular, basis.T @ values)
    residual = values - design @ coefficients
    norm = float(np.linalg.norm(residual))
    if not math.isfinite(norm) or norm <= resolution:
        raise _InsufficientPrecision('Residual variation is below numerical resolution')
    return coefficients, triangular, norm * norm


def evidence(values, covariance, calibration):
    """Test all one-change nulls using one fixed, known covariance shape.

    Array entry i corresponds to split i+1. Calibration comes from
    reporting_direct.critical_values; the same t/F cutoffs apply after GLS.
    No covariance is estimated or changed between candidate models here.

    status='ok' carries all statistics. status='insufficient_precision'
    abstains (has_alert=False) and leaves diagnostic arrays as None; that
    status is not statistical evidence that any particular split is plausible.
    Malformed or ill-conditioned inputs raise ValueError.
    """
    y = np.asarray(list(values), dtype=float)
    n = len(y)
    if y.ndim != 1 or not 4 <= n <= 200 or calibration['n'] != n:
        raise ValueError('Use 4 to 200 observations with matching calibration length')
    if not np.isfinite(y).all():
        raise ValueError('Observations must be finite')
    factor = 1 + calibration['threshold']
    if not math.isfinite(factor) or factor < 1:
        raise ValueError('threshold must be finite and nonnegative')
    for key in ('size_t_critical', 'fit_f_critical'):
        if not math.isfinite(calibration[key]) or calibration[key] <= 0:
            raise ValueError('Critical values must be finite and positive')
    z, steps, condition = _whiten(y, covariance)
    resolution = 64 * np.finfo(float).eps * n * float(np.linalg.norm(z))
    contrast = np.array([-factor, 1.0])
    t_scores, f_scores, extra_splits, surviving = [], [], [], []
    size_rejected, fit_rejected, three_costs = [], [], {}
    try:
        for split in range(1, n):
            design = np.column_stack((steps[:, 0] - steps[:, split], steps[:, split]))
            coefficients, triangular, two = _fit(design, z, resolution)
            # If B=U*A is a QR decomposition, G=(A'*A)^-1 and
            # h'*G*h = ||solve(A', h)||^2. Avoid forming a normal-matrix inverse.
            direction = solve_triangular(triangular, contrast, trans='T')
            variance_factor = float(direction @ direction)
            size_t = float(contrast @ coefficients) / math.sqrt(two / (n - 2) * variance_factor)
            best_cost, best_split = math.inf, None
            for extra in range(1, n):
                if extra == split:
                    continue
                left, right = sorted((split, extra))
                key = (left, right)
                if key not in three_costs:
                    extended = np.column_stack(
                        (
                            steps[:, 0] - steps[:, left],
                            steps[:, left] - steps[:, right],
                            steps[:, right],
                        )
                    )
                    three_costs[key] = _fit(extended, z, resolution)[2]
                if three_costs[key] < best_cost:
                    best_cost, best_split = three_costs[key], extra
            fit_f = max(0.0, two - best_cost) / (best_cost / (n - 3))
            if not (math.isfinite(size_t) and math.isfinite(fit_f)):
                raise _InsufficientPrecision('Test statistics exceed numerical range')
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
    except _InsufficientPrecision as exc:
        return {
            'status': 'insufficient_precision',
            'reason': str(exc),
            'covariance_condition': condition,
            'has_alert': False,
            'null_splits': None,
            'size_rejected_splits': None,
            'fit_rejected_splits': None,
            'size_t': None,
            'lack_of_fit_f': None,
            'best_extra_split': None,
        }
    return {
        'status': 'ok',
        'covariance_condition': condition,
        'has_alert': not surviving,
        'null_splits': surviving,
        'size_rejected_splits': size_rejected,
        'fit_rejected_splits': fit_rejected,
        'size_t': t_scores,
        'lack_of_fit_f': f_scores,
        'best_extra_split': extra_splits,
    }
