"""Joint location/correlation confidence set and continuous AR(1) null checks.

Stationary Gaussian AR(1), one common marginal variance, at most one mean step.
See unknown_correlation.rst. This is a mathematical reference, not a default.
"""

import math
from fractions import Fraction

import numpy as np
from scipy.special import logsumexp

from . import rational_polynomial as p
from .reporting_direct import critical_values


def calibration(n, *, alpha=0.05, confidence_alpha=0.01, threshold=0.05):
    if not isinstance(n, int) or isinstance(n, bool) or not 8 <= n <= 200:
        raise ValueError('Use 8 to 200 observations')
    if not 0 < confidence_alpha < alpha < 1:
        raise ValueError('Require 0 < confidence_alpha < alpha < 1')
    result = critical_values(n, alpha=alpha - confidence_alpha, threshold=threshold)
    if result['size_t_critical'] <= 0 or result['fit_f_critical'] <= 0:
        raise ValueError('Reporting cutoffs must be positive')
    result['total_alpha'] = alpha
    result['confidence_alpha'] = confidence_alpha
    return result


def predictor(training):
    """Choose a proper predictor using only a normalized training prefix."""
    y = np.asarray(training, dtype=float)
    mean = float(np.mean(y))
    lag, current = y[:-1] - mean, y[1:] - mean
    denominator = float(lag @ lag)
    rho = float(np.clip(float(lag @ current) / denominator, -0.95, 0.95)) if denominator else 0.0
    variance = max(float(np.mean((current - rho * lag) ** 2)), 1e-12)
    jump_variance = max(mean * mean, float(np.var(y)), 1e-12)
    return {'mean': mean, 'rho': rho, 'variance': variance, 'jump_variance': jump_variance}


def predictive_log_density(values, training_count, fitted):
    """Equal mixture of no jump and one Gaussian-prior jump in validation."""
    y = np.asarray(values, dtype=float)
    rho, variance = fitted['rho'], fitted['variance']
    z = y[training_count:] - rho * y[training_count - 1 : -1] - (1 - rho) * fitted['mean']
    m = len(z)
    common = m * math.log(2 * math.pi * variance)
    components = [-0.5 * (common + float(z @ z) / variance)]
    for j in range(m):
        vector = np.zeros(m)
        vector[j], vector[j + 1 :] = 1, 1 - rho
        norm = float(np.linalg.norm(vector))
        direction = vector / norm
        projection = float(direction @ z)
        residual = z - projection * direction
        inflation = fitted['jump_variance'] * norm**2 / variance
        quadratic = (float(residual @ residual) + projection**2 / (1 + inflation)) / variance
        components.append(-0.5 * (common + math.log1p(inflation) + quadratic))
    return float(logsumexp(components) - math.log(m + 1))


def confidence_threshold(values, confidence_alpha):
    training_count = len(values) // 2
    fitted = predictor(values[:training_count])
    with np.errstate(over='ignore', invalid='ignore'):
        try:
            log_q = predictive_log_density(values, training_count, fitted)
        except OverflowError:
            log_q = -math.inf
    if not math.isfinite(log_q):
        # Dropping the confidence exclusion is conservative; it never drops rho.
        log_q = -math.inf
    m = len(values) - training_count
    log_threshold = (
        math.log(m) - math.log(2 * math.pi) - 1 + 2 / m * (-math.log(confidence_alpha) - log_q)
    )
    # An infinite bound retains every pair, so it cannot cause a false rejection.
    bound = math.inf if log_threshold > 709 else math.nextafter(math.exp(log_threshold), math.inf)
    return {
        'training_count': training_count,
        'validation_count': m,
        'predictor': fitted,
        'log_predictive_density': log_q,
        'residual_bound': bound,
    }


def confidence_residual(values, split, training_count):
    """Exact quadratic after profiling two intercepts and a free transition."""
    values = tuple(Fraction(x) for x in values)
    result = p.polynomial([0])
    groups = (range(training_count, split), range(max(training_count, split + 1), len(values)))
    for indices in groups:
        indices = list(indices)
        if not indices:
            continue
        current = [values[i] for i in indices]
        lag = [values[i - 1] for i in indices]
        mean_current = sum(current) / len(indices)
        mean_lag = sum(lag) / len(indices)
        current = [x - mean_current for x in current]
        lag = [x - mean_lag for x in lag]
        result = p.add(
            result,
            p.polynomial(
                [
                    sum(x * x for x in current),
                    -2 * sum(x * z for x, z in zip(current, lag, strict=True)),
                    sum(z * z for z in lag),
                ]
            ),
        )
    return result


def precision_product(a, b):
    """a' P(rho) b for the stationary AR(1) tridiagonal precision numerator."""
    return p.polynomial(
        [
            sum(x * y for x, y in zip(a, b, strict=True)),
            -sum(a[i] * b[i + 1] + a[i + 1] * b[i] for i in range(len(a) - 1)),
            sum(a[i] * b[i] for i in range(1, len(a) - 1)),
        ]
    )


class _DegenerateResidual(Exception):
    pass


class Models:
    """Cache exact GLS residual and contrast polynomials for plateau designs."""

    def __init__(self, values):
        self.values = tuple(Fraction(x) for x in values)
        self.n = len(values)
        self.total = precision_product(self.values, self.values)
        self.cache = {}

    def fit(self, boundaries):
        boundaries = tuple(sorted(boundaries))
        if boundaries not in self.cache:
            edges = (0, *boundaries, self.n)
            columns = [
                tuple(int(a <= i < b) for i in range(self.n))
                for a, b in zip(edges[:-1], edges[1:], strict=True)
            ]
            gram = [[precision_product(a, b) for b in columns] for a in columns]
            linear = [precision_product(a, self.values) for a in columns]
            determinant, adjugate = p.determinant(gram), p.adjugate(gram)
            explained = p.polynomial([0])
            for i in range(len(columns)):
                for j in range(len(columns)):
                    explained = p.add(
                        explained, p.multiply(p.multiply(linear[i], adjugate[i][j]), linear[j])
                    )
            numerator = p.subtract(p.multiply(self.total, determinant), explained)
            if numerator == (0,):
                raise _DegenerateResidual('Exactly zero residual variation in a candidate model')
            self.cache[boundaries] = (determinant, numerator, adjugate, linear)
        return self.cache[boundaries]

    def size(self, split, config):
        _, numerator, adjugate, linear = self.fit((split,))
        contrast = (-Fraction(1 + config['threshold']), Fraction(1))
        excess, variance = p.polynomial([0]), p.polynomial([0])
        for i in range(2):
            for j in range(2):
                excess = p.add(excess, p.scale(p.multiply(adjugate[i][j], linear[j]), contrast[i]))
                variance = p.add(variance, p.scale(adjugate[i][j], contrast[i] * contrast[j]))
        comparison = p.subtract(
            p.scale(p.multiply(excess, excess), self.n - 2),
            p.scale(p.multiply(numerator, variance), Fraction(config['size_t_critical']) ** 2),
        )
        return tuple(p.remove_stationary_endpoint_factors(x) for x in (excess, comparison))

    def shape(self, split, extra, config):
        d2, n2, _, _ = self.fit((split,))
        d3, n3, _, _ = self.fit((split, extra))
        comparison = p.subtract(
            p.scale(p.multiply(n2, d3), self.n - 3),
            p.scale(p.multiply(n3, d2), self.n - 3 + Fraction(config['fit_f_critical'])),
        )
        return p.remove_stationary_endpoint_factors(comparison)


def evidence(values, config=None, *, max_cells=4096, max_depth=16):
    """Certify elimination of every (split,rho) for -1<rho<1, or abstain.

    A surviving rational point gives status='surviving_explanation'. Exhausted
    certification budgets give 'unresolved', never an alert. Successful results
    include exact interval endpoints and the route certifying each interval.
    """
    y = list(values)
    n = len(y)
    if not 8 <= n <= 200 or any(not math.isfinite(x) for x in y):
        raise ValueError('Use 8 to 200 finite observations')
    if config is None:
        config = calibration(n)
    expected = calibration(
        n,
        alpha=config['total_alpha'],
        confidence_alpha=config['confidence_alpha'],
        threshold=config['threshold'],
    )
    if config != expected:
        raise ValueError('Use a matching unmodified calibration')
    if not isinstance(max_cells, int) or isinstance(max_cells, bool) or max_cells < 1:
        raise ValueError('max_cells must be a positive integer')
    if not isinstance(max_depth, int) or isinstance(max_depth, bool) or max_depth < 0:
        raise ValueError('max_depth must be a nonnegative integer')
    # Scaling is chosen from training alone, preserving conditional-density logic.
    scale = max(map(abs, y[: n // 2])) or 1.0
    y = [x / scale for x in y]
    if any(not math.isfinite(x) for x in y):
        raise ValueError('Observations exceed numerical range after training normalization')
    confidence = confidence_threshold(y, config['confidence_alpha'])
    confidence['normalization_scale'] = scale
    models = Models(y)
    bound = confidence['residual_bound']
    certificates = []
    residuals = {
        split: confidence_residual(models.values, split, confidence['training_count'])
        for split in range(1, n)
    }
    confidence_polynomials = {
        str(split): [str(x) for x in poly] for split, poly in residuals.items()
    }
    visited = 0

    def result(status, witness=None):
        return {
            'has_alert': status == 'certified_alert',
            'status': status,
            'witness': witness,
            'cells_visited': visited,
            'certificate': certificates,
            'confidence': confidence,
            'confidence_polynomials': confidence_polynomials,
            'calibration': config,
        }

    try:
        for split in range(1, n):
            residual = residuals[split]
            excluded = (
                p.polynomial([-1])
                if math.isinf(bound)
                else p.subtract(residual, p.polynomial([bound]))
            )
            size = None
            shapes = {}
            pending = [(Fraction(-1), Fraction(1), 0)]
            while pending:
                left, right, depth = pending.pop()
                if visited >= max_cells:
                    return result(
                        'unresolved', {'split': split, 'left': str(left), 'right': str(right)}
                    )
                visited += 1
                midpoint = (left + right) / 2
                candidates = []
                if p.evaluate(excluded, midpoint) > 0:
                    candidates.append(('confidence', None, (excluded,)))
                if size is None:
                    size = models.size(split, config)
                if all(p.evaluate(x, midpoint) > 0 for x in size):
                    candidates.append(('size', None, size))
                for extra in range(1, n):
                    if extra == split:
                        continue
                    if extra not in shapes:
                        shapes[extra] = models.shape(split, extra, config)
                    if p.evaluate(shapes[extra], midpoint) > 0:
                        candidates.append(('shape', extra, (shapes[extra],)))
                if not candidates:
                    return result(
                        'surviving_explanation',
                        {'split': split, 'rho': str(midpoint), 'rho_float': float(midpoint)},
                    )
                for route, extra, polynomials in candidates:
                    if all(p.positive_on(x, left, right) for x in polynomials):
                        certificates.append(
                            {
                                'split': split,
                                'left': str(left),
                                'right': str(right),
                                'route': route,
                                'extra': extra,
                            }
                        )
                        break
                else:
                    if depth >= max_depth:
                        return result(
                            'unresolved', {'split': split, 'left': str(left), 'right': str(right)}
                        )
                    pending.extend(((midpoint, right, depth + 1), (left, midpoint, depth + 1)))
                    continue
    except _DegenerateResidual:
        return result('insufficient_variation')
    return result('certified_alert')
