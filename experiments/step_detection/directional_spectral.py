"""Directional tail certificates using several residual covariance directions.

All tilts and direction counts bound the same p-value. This standalone
mathematical reference leaves the existing reporting searches unchanged.
"""

import math
from fractions import Fraction as F

from . import directional_determinant as determinant
from . import residual_direction as direction

TILTS = (F(1, 16), F(1, 8))
COUNTS = (4, 8, 12, 16)


def density_constant(count):
    """If count=4m tilted weights exceed b>0, density <= constant/b."""
    if not isinstance(count, int) or isinstance(count, bool) or count < 4 or count % 4:
        raise ValueError('Require a positive multiple of four directions')
    m = count // 4
    return F(math.comb(2 * m - 2, m - 1), 4**m)


def eigenvalue_lower(n, count, left, right, *, stationary_normalized=False):
    """Uniform lower bound on the count-th largest eigenvalue of U' R U.

    U removes two plateau levels. Interlacing costs two eigenvalue positions;
    an AR precision comparison and sin(x)<=x remove trigonometric arithmetic.
    In stationary_normalized mode, bound U' R U/(1-rho**2) instead.
    """
    left, right = F(left), F(right)
    density_constant(count)
    if (
        not isinstance(n, int)
        or isinstance(n, bool)
        or n < count + 2
        or not -1 < left <= right < 1
    ):
        raise ValueError('Require enough residual directions and a stationary interval')
    lower = F(0) if left <= 0 <= right else min(abs(left), abs(right))
    upper = max(abs(left), abs(right))
    frequency_squared = F(22 * (count + 2), 7 * (n + 1)) ** 2
    # This quadratic is convex, so its maximum is at one of the endpoints.
    denominator = max((1 - value) ** 2 + value * frequency_squared for value in (lower, upper))
    numerator = F(1) if stationary_normalized else 1 - upper**2
    return numerator / denominator


def coefficient_bound_squared(weights, tilt, count):
    """Exact scalar reference for analytical checks with supplied eigenweights."""
    weights, tilt = [F(x) for x in weights], F(tilt)
    constant = density_constant(count)
    if not weights or tilt <= 0 or any(1 + 2 * tilt * a <= 0 for a in weights):
        raise ValueError('Require a finite exponential tilt')
    det = math.prod(1 + 2 * tilt * a for a in weights)
    positive = sorted((a / (1 + 2 * tilt * a) for a in weights if a > 0), reverse=True)
    correction = (
        max(F(1), tilt * positive[count - 1] / constant) if len(positive) >= count else F(1)
    )
    return min(F(1), 1 / (det * correction**2))


def at_tilt(model, left, right, tilt, *, delta=direction.DELTA, bits=determinant.BITS):
    """Refine one full-determinant interval certificate at a fixed tilt."""
    tilt, delta = F(tilt), F(delta)
    base = determinant.certify_interval(model, left, right, tilt=tilt, delta=delta, bits=bits)
    result = {
        'tilt': str(tilt),
        'base': base,
        'directions': [],
        'density_correction': '1',
        'selected_count': None,
        'p_upper_squared': '1',
        'p_upper': '1',
        'status': 'not_certified',
    }
    if 'determinant_enclosure' not in base or F(base['determinant_enclosure'][0]) <= 0:
        return result
    radius = F(base['radius_upper'])
    correction = F(1)
    for count in COUNTS:
        if count > model['n'] - 2:
            continue
        eigen_lower = eigenvalue_lower(model['n'], count, left, right)
        weight = eigen_lower / radius - 1
        tilted_weight = weight / (1 + 2 * tilt * weight) if weight > 0 else F(0)
        candidate = max(F(1), tilt * tilted_weight / density_constant(count))
        result['directions'].append(
            {
                'count': count,
                'eigenvalue_lower': str(eigen_lower),
                'weight_lower': str(weight),
                'positive_tilted_weight_lower': str(tilted_weight),
                'density_constant': str(density_constant(count)),
                'correction': str(candidate),
            }
        )
        if candidate > correction:
            correction, result['selected_count'] = candidate, count
    # Refine the untruncated Chernoff bound; clipping it at one first is invalid.
    squared = min(F(1), 1 / (F(base['determinant_enclosure'][0]) * correction**2))
    upper = min(F(1), direction.sqrt_lower(squared, bits=48) + F(1, 2**48))
    result.update(
        density_correction=str(correction),
        p_upper_squared=str(squared),
        p_upper=str(upper),
        status='certified_excluded' if squared < delta**2 else 'not_certified',
    )
    return result


def certify_interval(model, left, right, *, delta=direction.DELTA, bits=determinant.BITS):
    """Minimum of two valid interval bounds; no extra probability allocation."""
    attempts = [at_tilt(model, left, right, tilt, delta=delta, bits=bits) for tilt in TILTS]
    best = min(attempts, key=lambda result: F(result['p_upper_squared']))
    return {
        'left': str(F(left)),
        'right': str(F(right)),
        'delta': str(F(delta)),
        'bits': bits,
        'tilts': [str(tilt) for tilt in TILTS],
        'counts': list(COUNTS),
        'selected_tilt': best['tilt'],
        'selected_count': best['selected_count'],
        'p_upper_squared': best['p_upper_squared'],
        'p_upper': best['p_upper'],
        'status': best['status'],
        'attempts': attempts,
    }
