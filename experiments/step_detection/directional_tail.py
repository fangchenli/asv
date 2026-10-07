"""Exact upper tail certificates near rho=1; no complete reporting search."""

import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from . import rational_polynomial as p
from . import residual_direction as direction

HERE = Path(__file__).parent
TILT = F(1, 16)
BITS = 32


def ceil_dyadic(value, bits=BITS):
    """Outward rational rounding to limit subsequent determinant sizes."""
    value = F(value)
    if not isinstance(bits, int) or bits < 0:
        raise ValueError('Require nonnegative integer precision')
    scale = 1 << bits
    return F(-(-value.numerator * scale // value.denominator), scale)


def gram_eigenvalue_upper(n, split):
    """Upper bound on four distinct eigenvalues of the difference Gram matrix."""
    if n < 6 or not 1 <= split < n:
        raise ValueError('Require at least four residual coordinates')
    bounds = sorted(
        F(22 * j, 7 * length) ** 2
        for length in (split, n - split)
        for j in range(1, min(length, 5))
    )
    return bounds[3]


def determinant_lower(n, split, diagonal_shift, gram_scale):
    """det(gram_scale*G + diagonal_shift*I) / det(G), using two recurrences."""
    diagonal_shift, gram_scale = F(diagonal_shift), F(gram_scale)
    if n < 3 or not 1 <= split < n or diagonal_shift <= 0 or gram_scale <= 0:
        raise ValueError('Require an interior split and positive matrix coefficients')
    diagonal = 2 * gram_scale + diagonal_shift
    result = F(1)
    for length in (split, n - split):
        previous, current = F(1), diagonal
        for _ in range(2, length):
            previous, current = current, diagonal * current - gram_scale**2 * previous
        if length > 1:
            result *= current / length
    return result


def tilted_bound_squared(weights, tilt):
    """Exact bound for explicit rational weights, used in analytical tests."""
    weights, tilt = [F(x) for x in weights], F(tilt)
    if not weights or tilt <= 0 or any(1 + 2 * tilt * a <= 0 for a in weights):
        raise ValueError('Require a finite exponential tilt')
    determinant = math.prod(1 + 2 * tilt * a for a in weights)
    positive = sorted((a / (1 + 2 * tilt * a) for a in weights if a > 0), reverse=True)
    correction = max(F(1), 4 * tilt * positive[3]) if len(positive) >= 4 else F(1)
    return min(F(1), 1 / (determinant * correction**2))


def certify_interval(model, left, right, tilt=TILT, delta=direction.DELTA):
    """Uniform upper bound on the directional lower-tail probability.

    The scope is 0<=left<=right<=1, where 1 denotes a limiting endpoint.
    An unhelpful bound returns one, never a surviving-explanation claim.
    """
    left, right, tilt, delta = map(F, (left, right, tilt, delta))
    if not 0 <= left <= right <= 1 or not 0 < tilt < F(1, 2) or not 0 < delta < 1:
        raise ValueError('Require a nonnegative correlation interval and valid tilt/budget')
    result = {
        'left': str(left),
        'right': str(right),
        'tilt': str(tilt),
        'delta': str(delta),
        'status': 'not_certified',
        'p_upper': '1',
        'p_upper_squared': '1',
    }
    n, split = model['n'], model['split']
    kappa = 2 - (n - 3) * (1 - left)
    result['covariance_lower'] = str(kappa)
    if kappa <= 0:
        result['reason'] = 'Covariance lower bound is nonpositive'
        return result

    def coefficients(poly):
        return (p.evaluate(poly, left),) if left == right else p.bernstein(poly, left, right)

    lower_n, upper_h = min(coefficients(model['num'])), max(coefficients(model['den']))
    if lower_n <= 0 or upper_h <= 0:
        result['reason'] = 'Residual ratio enclosure is inconclusive'
        return result
    radius_upper = ceil_dyadic((1 + right) * model['ss'] * upper_h / lower_n)
    det = determinant_lower(n, split, 2 * tilt * kappa / radius_upper, 1 - 2 * tilt)
    correction = F(1)
    result.update({'radius_upper': str(radius_upper), 'determinant_lower': str(det)})
    if n >= 6:
        gram_upper = gram_eigenvalue_upper(n, split)
        weight_lower = kappa / (radius_upper * gram_upper) - 1
        result.update(
            {
                'four_gram_eigenvalues_upper': str(gram_upper),
                'four_weights_lower': str(weight_lower),
            }
        )
        if weight_lower > 0:
            tilted_lower = weight_lower / (1 + 2 * tilt * weight_lower)
            correction = max(F(1), 4 * tilt * tilted_lower)
            result['four_tilted_weights_lower'] = str(tilted_lower)
    squared = min(F(1), 1 / (det * correction**2))
    upper = min(F(1), direction.sqrt_lower(squared, bits=48) + F(1, 2**48))
    result.update(
        {
            'density_correction': str(correction),
            'chernoff_upper_squared': str(min(F(1), 1 / det)),
            'p_upper_squared': str(squared),
            'p_upper': str(upper),
            'status': 'certified_excluded' if squared < delta**2 else 'not_certified',
        }
    )
    return result


def approximate_tail(values, split, rho):
    """Imhof inversion as a floating diagnostic, never an interval certificate."""
    if not -1 < rho < 1:
        raise ValueError('Require stationary correlation')
    y = np.array(values, dtype=float)
    basis = direction.residual_basis(len(y), split)
    y[:split] -= y[:split].mean()
    y[split:] -= y[split:].mean()
    z = basis.T @ y
    if np.linalg.norm(z) == 0:
        raise ValueError('Require nonzero residuals')
    s = z / np.linalg.norm(z)
    distance = np.abs(np.arange(len(y))[:, None] - np.arange(len(y))[None, :])
    covariance = basis.T @ (rho**distance) @ basis
    radius = 1 / float(s @ np.linalg.solve(covariance, s))
    weights = np.linalg.eigvalsh(covariance) / radius - 1
    if rho == 0 or len(z) == 1:
        return {'probability': 1.0, 'quadrature_error_estimate': 0.0}

    def integrand(u):
        if u == 0:
            return float(weights.sum())
        terms = 2 * u * weights
        return float(
            np.sin(0.5 * np.arctan(terms).sum()) * np.exp(-0.25 * np.log1p(terms**2).sum()) / u
        )

    integral, error = quad(integrand, 0, np.inf, epsabs=1e-10, limit=500)
    return {
        'probability': float(0.5 - integral / math.pi),
        'quadrature_error_estimate': float(error / math.pi),
    }


def diagnose():
    independent, independent_source = direction.info.load_example()
    positive, positive_source = direction.load_positive_example()
    cases = []
    for case, source, witness in (
        (independent, independent_source, F(4095, 4096)),
        (positive, positive_source, F(16383, 16384)),
    ):
        model = direction.state(case['values'], case['position'])
        points = sorted({F(0), F(str(case['correlation'])), F(9, 10), F(99, 100), witness})
        cases.append(
            {
                'id': case['id'],
                'archive': source,
                'split': case['position'],
                'witness': str(witness),
                'pointwise_diagnostics': [
                    {
                        'rho': str(rho),
                        **approximate_tail(case['values'], case['position'], float(rho)),
                    }
                    for rho in points
                ],
                'interval_certificates': [
                    certify_interval(model, left, F(1))
                    for left in (F(99, 100), F(999, 1000), witness)
                ],
            }
        )
    return {
        'purpose': 'Archived mathematical diagnosis, not a new detection study',
        'tilt': str(TILT),
        'radius_upper_bits': BITS,
        'probability_display_bits': 48,
        'cases': cases,
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in (
                'directional_tail.py',
                'residual_direction.py',
                'rational_polynomial.py',
                'reporting_ar1.py',
                'reporting_ar1_full.py',
                'ar1_information.py',
            )
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = diagnose()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
