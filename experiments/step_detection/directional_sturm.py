"""Split-aware eigenvalue certificates from exact tridiagonal Sturm counts.

This mathematical reference keeps the existing tilts and density inequality.
It tightens only the eigenvalue lower bounds used by the spectral certificate.
"""

from fractions import Fraction as F
from functools import lru_cache

from . import directional_determinant as determinant
from . import directional_spectral as spectral
from . import residual_direction as direction

EIGEN_BITS = 32


def precision_count_below(n, correlation, threshold):
    """Count eigenvalues strictly below threshold in the AR precision numerator.

    Leading-principal-minor signs form a Sturm sequence. Interior zeros are
    skipped; the adjacent nonzero signs are opposite for nonzero off-diagonals.
    A final zero is skipped too, counting strict inequality at an eigenvalue.
    """
    correlation, threshold = F(correlation), F(threshold)
    if not isinstance(n, int) or isinstance(n, bool) or n < 2 or not 0 <= correlation < 1:
        raise ValueError('Require n>=2 and a stationary nonnegative correlation')
    if correlation == 0:
        return n if threshold > 1 else 0
    previous, current = F(1), 1 - threshold
    last_sign, count = 1, 0
    for i in range(1, n + 1):
        if current:
            sign = 1 if current > 0 else -1
            count += sign != last_sign
            last_sign = sign
        if i < n:
            diagonal = 1 if i == n - 1 else 1 + correlation**2
            previous, current = (
                current,
                (diagonal - threshold) * current - correlation**2 * previous,
            )
    return count


def precision_eigenvalue_upper(n, index, correlation):
    """Rational upper bound on the index-th smallest AR precision eigenvalue."""
    correlation = F(correlation)
    if (
        not isinstance(n, int)
        or isinstance(n, bool)
        or n < 2
        or not isinstance(index, int)
        or isinstance(index, bool)
        or not 1 <= index <= n
        or not 0 <= correlation < 1
    ):
        raise ValueError('Invalid precision eigenvalue request')
    return _precision_eigenvalue_upper(n, index, correlation)


@lru_cache(maxsize=4096)
def _precision_eigenvalue_upper(n, index, correlation):
    """Cache only canonical inputs, after public argument validation."""
    if correlation == 0:
        return F(1)
    lower, upper = F(0), (1 + correlation) ** 2
    # Gershgorin supplies the initial upper bound; all eigenvalues are positive.
    for _ in range(EIGEN_BITS):
        midpoint = (lower + upper) / 2
        if precision_count_below(n, correlation, midpoint) >= index:
            upper = midpoint
        else:
            lower = midpoint
    return upper


def covariance_eigenvalue_lower(n, index, left, right):
    """Uniform lower bound on the index-th largest AR covariance eigenvalue."""
    left, right = F(left), F(right)
    if not -1 < left <= right < 1:
        raise ValueError('Require a stationary interval')
    low = F(0) if left <= 0 <= right else min(abs(left), abs(right))
    high = max(abs(left), abs(right))
    midpoint = (low + high) / 2
    # A(v)-A(midpoint) has at most two off-diagonal changes per row.
    # Its operator norm is bounded by the maximum absolute row sum.
    perturbation = high**2 - midpoint**2 + 2 * (high - midpoint)
    precision_upper = precision_eigenvalue_upper(n, index, midpoint)
    return (1 - high**2) / (precision_upper + perturbation)


def residual_eigenvalue_bounds(n, split, count, left, right):
    """Compare full-history and single-plateau subspaces with old bounds."""
    if not isinstance(split, int) or isinstance(split, bool) or not 1 <= split < n:
        raise ValueError('Require an interior split')
    original = spectral.eigenvalue_lower(n, count, left, right)
    full = covariance_eigenvalue_lower(n, count + 2, left, right)
    length = max(split, n - split)
    block = covariance_eigenvalue_lower(length, count + 1, left, right) if count < length else F(0)
    return {'toeplitz': original, 'full_sturm': full, 'plateau_sturm': block}


def _floor_nth_root(value, degree):
    """Exact floor of a nonnegative integer nth root."""
    low, high = 0, 1 << ((value.bit_length() + degree - 1) // degree)
    while high**degree <= value:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if middle**degree <= value:
            low = middle
        else:
            high = middle
    return low


def grouped_density_bounds(tilted_weight_lowers, *, bits=64):
    """Fourier density bounds retaining successive groups of four weights.

    Entry j is a lower bound on the j-th selected positive tilted coefficient,
    for ranks 4, 8, 12, ... . Hölder combines the first m groups into a bound
    ``c_m / geometric_mean(weights)``. Roots are rounded down exactly, so the
    returned density bounds are rational upper bounds.
    """
    weights = [F(value) for value in tilted_weight_lowers]
    if not isinstance(bits, int) or isinstance(bits, bool) or bits < 1:
        raise ValueError('Require positive integer root precision')
    if any(value <= 0 for value in weights):
        raise ValueError('Require positive tilted-weight lower bounds')
    scale = 1 << bits
    product = F(1)
    results = []
    for groups, weight in enumerate(weights, start=1):
        product *= weight
        degree = groups
        scaled_floor = product.numerator * scale**degree // product.denominator
        root_lower = F(_floor_nth_root(scaled_floor, degree), scale)
        if root_lower <= 0:
            break
        density_upper = spectral.density_constant(4 * groups) / root_lower
        results.append(
            {
                'count': 4 * groups,
                'geometric_mean_lower': str(root_lower),
                'density_upper': str(density_upper),
            }
        )
    return results


def at_tilt(model, left, right, tilt, *, delta=direction.DELTA, bits=determinant.BITS):
    """Same determinant and density correction, with stronger eigenvalue inputs."""
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
    correction = F(1)
    tilted_weight_lowers = []
    for count in spectral.COUNTS:
        if count > model['n'] - 2:
            continue
        bounds = residual_eigenvalue_bounds(model['n'], model['split'], count, left, right)
        selected = max(bounds, key=bounds.get)
        eigen_lower = bounds[selected]
        weight = eigen_lower / F(base['radius_upper']) - 1
        tilted = weight / (1 + 2 * tilt * weight) if weight > 0 else F(0)
        candidate = max(F(1), tilt * tilted / spectral.density_constant(count))
        if tilted > 0:
            tilted_weight_lowers.append(tilted)
        result['directions'].append(
            {
                'count': count,
                'eigenvalue_bounds': {key: str(value) for key, value in bounds.items()},
                'selected_eigenvalue_bound': selected,
                'weight_lower': str(weight),
                'positive_tilted_weight_lower': str(tilted),
                'correction': str(candidate),
            }
        )
        if candidate > correction:
            correction, result['selected_count'] = candidate, count
    grouped = grouped_density_bounds(tilted_weight_lowers)
    result['grouped_density_bounds'] = grouped
    for item in grouped:
        candidate = max(F(1), tilt / F(item['density_upper']))
        if candidate > correction:
            correction, result['selected_count'] = candidate, item['count']
    squared = min(F(1), 1 / (F(base['determinant_enclosure'][0]) * correction**2))
    result.update(
        density_correction=str(correction),
        p_upper_squared=str(squared),
        p_upper=str(min(F(1), direction.sqrt_lower(squared, bits=48) + F(1, 2**48))),
        status='certified_excluded' if squared < delta**2 else 'not_certified',
    )
    return result


def certify_interval(model, left, right, *, delta=direction.DELTA, bits=determinant.BITS):
    """Strengthen the two-tilt certificate without changing its probability rule."""
    attempts = [
        at_tilt(model, left, right, tilt, delta=delta, bits=bits) for tilt in spectral.TILTS
    ]
    best = min(attempts, key=lambda result: F(result['p_upper_squared']))
    return {
        'left': str(F(left)),
        'right': str(F(right)),
        'delta': str(F(delta)),
        'bits': bits,
        'eigen_bits': EIGEN_BITS,
        'tilts': [str(tilt) for tilt in spectral.TILTS],
        'counts': list(spectral.COUNTS),
        'selected_tilt': best['tilt'],
        'selected_count': best['selected_count'],
        'p_upper_squared': best['p_upper_squared'],
        'p_upper': best['p_upper'],
        'status': best['status'],
        'attempts': attempts,
    }
