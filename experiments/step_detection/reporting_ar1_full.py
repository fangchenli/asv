"""Full-history AR(1) confidence reference with a proper predictive mixture.

See full_history_confidence.rst. Units must be specified independently of the
observations. This module does not alter production or earlier frozen studies.
"""

import math
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np
from scipy.special import gammaln, logsumexp

from . import rational_polynomial as p
from . import reporting_ar1 as a

PRIOR = {
    'mean': 0.0,
    'mean_precision': 1e-4,
    'shape': 0.5,
    'scale_exponents': list(range(-12, 13)),
}


def component_log_density(values, split, beta):
    """Integrate independent plateau means and one shared noise variance."""
    y = np.asarray(values, dtype=float)
    n = len(y)
    kappa, alpha = PRIOR['mean_precision'], PRIOR['shape']
    blocks = [y] if split is None else [y[:split], y[split:]]
    cost, determinant = 0.0, 0.0
    for block in blocks:
        center = float(np.mean(block))
        residual = block - center
        cost += (
            float(residual @ residual)
            + kappa * len(block) / (kappa + len(block)) * (center - PRIOR['mean']) ** 2
        )
        determinant += math.log(kappa / (kappa + len(block)))
    return float(
        -n / 2 * math.log(2 * math.pi)
        + determinant / 2
        + gammaln(alpha + n / 2)
        - gammaln(alpha)
        + alpha * math.log(beta)
        - (alpha + n / 2) * math.log(beta + cost / 2)
    )


def predictive_log_density(values):
    """Proper density in externally normalized units, without parameter fitting."""
    n = len(values)
    components = []
    exponents = PRIOR['scale_exponents']
    with np.errstate(over='ignore', invalid='ignore'):
        try:
            for split in [None, *range(1, n)]:
                noise_mixture = logsumexp(
                    [component_log_density(values, split, 0.5 * 2.0 ** (2 * j)) for j in exponents]
                ) - math.log(len(exponents))
                weight = 0.5 if split is None else 0.5 / (n - 1)
                components.append(noise_mixture + math.log(weight))
            result = float(logsumexp(components))
        except (OverflowError, ValueError):
            return -math.inf
    return result if math.isfinite(result) else -math.inf


def confidence_coefficient(log_q, n, delta):
    """Lower exponential of the computed log coefficient, as a rational.

    Density and transcendental arithmetic preceding this conversion remain
    floating point. Subsequent interval checks are exact for this coefficient.
    """
    log_value = 2 * (math.log(delta) + log_q) + n * (math.log(2 * math.pi) + 1 - math.log(n))
    if not math.isfinite(log_value) or abs(log_value) > 10000:
        return Fraction(0)  # Disable confidence exclusion on numerical failure.
    with localcontext() as context:
        context.prec = 80
        return Fraction(Decimal.from_float(log_value).exp().next_minus())


def _divide_endpoint(poly, endpoint):
    quotient = [Fraction(0)] * (len(poly) - 1)
    quotient[-1] = poly[-1]
    for i in range(len(quotient) - 2, -1, -1):
        quotient[i] = poly[i + 1] + endpoint * quotient[i + 1]
    return p.scale(quotient, -1 if endpoint == 1 else 1)


def residual_ratio(models, split):
    """Cancel only COMMON endpoint factors, preserving the value Q=N/d."""
    denominator, numerator, _, _ = models.fit((split,))
    for endpoint in (Fraction(1), Fraction(-1)):
        while (
            len(numerator) > 1
            and len(denominator) > 1
            and p.evaluate(numerator, endpoint) == p.evaluate(denominator, endpoint) == 0
        ):
            numerator = _divide_endpoint(numerator, endpoint)
            denominator = _divide_endpoint(denominator, endpoint)
    return numerator, denominator


def excluded_at(numerator, denominator, coefficient, n, rho):
    num, den = p.evaluate(numerator, rho), p.evaluate(denominator, rho)
    return coefficient * num**n > (1 - rho**2) * den**n


def excluded_on(numerator, denominator, coefficient, n, left, right):
    """Sufficient exact confidence certificate without expanding degree-O(n) powers.

    The denominator is positive inside (-1,1) by the GLS Gram determinant
    identity. Bernstein bounds enclose the remaining low-degree polynomials.
    """
    lower = min(p.bernstein(numerator, left, right))
    upper = max(p.bernstein(denominator, left, right))
    stationary = Fraction(1) if left <= 0 <= right else max(1 - left**2, 1 - right**2)
    return lower > 0 and upper > 0 and coefficient * lower**n > stationary * upper**n


def evidence(values, *, unit, config=None, max_cells=4096, max_depth=16):
    """Combine full-history confidence exclusion and existing t/F certificates."""
    return _evidence(
        values,
        unit=unit,
        config=config,
        max_cells=max_cells,
        max_depth=max_depth,
        predictor=predictive_log_density,
        prior=PRIOR,
    )


def _evidence(values, *, unit, config, max_cells, max_depth, predictor, prior, indexed=False):
    """Shared engine; indexed predictors must be proper for each fixed split."""
    raw = list(values)
    n = len(raw)
    if not 8 <= n <= 200 or any(not math.isfinite(x) for x in raw):
        raise ValueError('Use 8 to 200 finite observations')
    if not math.isfinite(unit) or unit <= 0:
        raise ValueError('Specify a positive finite external unit')
    y = [x / unit for x in raw]
    if any(not math.isfinite(x) for x in y):
        raise ValueError('External units overflow the observations')
    if config is None:
        config = a.calibration(n)
    expected = a.calibration(
        n,
        alpha=config['total_alpha'],
        confidence_alpha=config['confidence_alpha'],
        threshold=config['threshold'],
    )
    if config != expected:
        raise ValueError('Use an unmodified matching calibration')
    for count, minimum in [(max_cells, 1), (max_depth, 0)]:
        if not isinstance(count, int) or isinstance(count, bool) or count < minimum:
            raise ValueError('Invalid certification budget')
    log_q = None if indexed else predictor(y)
    global_coefficient = (
        None if indexed else confidence_coefficient(log_q, n, config['confidence_alpha'])
    )
    location_densities = {}
    models = a.Models(y)
    certificates, residuals = [], {}
    visited = 0

    def result(status, witness=None):
        confidence = {'unit': unit, 'prior': prior, 'residual_ratios': residuals}
        if indexed:
            confidence['by_split'] = location_densities
        else:
            confidence.update(
                {'log_predictive_density': log_q, 'coefficient': str(global_coefficient)}
            )
        return {
            'has_alert': status == 'certified_alert',
            'status': status,
            'witness': witness,
            'cells_visited': visited,
            'certificate': certificates,
            'confidence': confidence,
            'calibration': config,
        }

    try:
        for split in range(1, n):
            if indexed:
                location_log_q = predictor(y, split)
                coefficient = confidence_coefficient(location_log_q, n, config['confidence_alpha'])
                location_densities[str(split)] = {
                    'log_predictive_density': location_log_q,
                    'coefficient': str(coefficient),
                }
            else:
                coefficient = global_coefficient
            numerator, denominator = residual_ratio(models, split)
            residuals[str(split)] = {
                'numerator': [str(x) for x in numerator],
                'denominator': [str(x) for x in denominator],
            }
            size = models.size(split, config)
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
                if excluded_on(numerator, denominator, coefficient, n, left, right):
                    certificates.append(
                        {
                            'split': split,
                            'left': str(left),
                            'right': str(right),
                            'route': 'confidence',
                            'extra': None,
                        }
                    )
                    continue
                candidates = []
                if all(p.evaluate(poly, midpoint) > 0 for poly in size):
                    candidates.append(('size', None, size))
                for extra in range(1, n):
                    if extra == split:
                        continue
                    if extra not in shapes:
                        shapes[extra] = models.shape(split, extra, config)
                    if p.evaluate(shapes[extra], midpoint) > 0:
                        candidates.append(('shape', extra, (shapes[extra],)))
                if not candidates and not excluded_at(
                    numerator, denominator, coefficient, n, midpoint
                ):
                    return result(
                        'surviving_explanation',
                        {
                            'split': split,
                            'rho': str(midpoint),
                            'rho_float': float(midpoint),
                        },
                    )
                for route, extra, polynomials in candidates:
                    if all(p.positive_on(poly, left, right) for poly in polynomials):
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
    except a._DegenerateResidual:
        return result('insufficient_variation')
    return result('certified_alert')
