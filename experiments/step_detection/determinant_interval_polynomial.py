"""Outward-rounded polynomial recurrence for projected AR(1) determinants."""

from fractions import Fraction as F
from math import comb, lcm


def _floor_ratio(numerator, denominator):
    return numerator // denominator


def _ceil_ratio(numerator, denominator):
    return -((-numerator) // denominator)


def _constant(value, scale):
    value = F(value)
    numerator = value.numerator * scale
    return (
        _floor_ratio(numerator, value.denominator),
        _ceil_ratio(numerator, value.denominator),
    )


def _trim(polynomial):
    result = list(polynomial)
    while len(result) > 1 and result[-1] == (0, 0):
        result.pop()
    return tuple(result)


def _add(first, second):
    size = max(len(first), len(second))
    return _trim(
        (
            (first[i][0] if i < len(first) else 0) + (second[i][0] if i < len(second) else 0),
            (first[i][1] if i < len(first) else 0) + (second[i][1] if i < len(second) else 0),
        )
        for i in range(size)
    )


def _scale(polynomial, factor):
    factor = F(factor)
    if factor >= 0:
        return _trim(
            (
                _floor_ratio(value[0] * factor.numerator, factor.denominator),
                _ceil_ratio(value[1] * factor.numerator, factor.denominator),
            )
            for value in polynomial
        )
    return _trim(
        (
            _floor_ratio(value[1] * factor.numerator, factor.denominator),
            _ceil_ratio(value[0] * factor.numerator, factor.denominator),
        )
        for value in polynomial
    )


def _subtract(first, second):
    return _add(first, _scale(second, -1))


def _multiply(first, second, scale):
    result = []
    for degree in range(len(first) + len(second) - 1):
        lower, upper = 0, 0
        for i in range(max(0, degree - len(second) + 1), min(len(first) - 1, degree) + 1):
            left, right = first[i], second[degree - i]
            products = (
                left[0] * right[0],
                left[0] * right[1],
                left[1] * right[0],
                left[1] * right[1],
            )
            lower += min(products)
            upper += max(products)
        result.append((_floor_ratio(lower, scale), _ceil_ratio(upper, scale)))
    return _trim(result)


def _continuants(n, split, a, b, center, scale):
    zero, one = (_constant(0, scale),), (_constant(1, scale),)
    rho = (_constant(center, scale), _constant(1, scale))
    rho2 = _multiply(rho, rho, scale)
    q = _subtract(one, rho2)
    end = _add((_constant(a, scale),), _scale(q, b))
    inner = _add(_scale(_add(one, rho2), a), _scale(q, b))
    off2 = _scale(rho2, F(a) ** 2)

    leading = [one, end]
    for size in range(2, n + 1):
        diagonal = end if size == n else inner
        leading.append(
            _subtract(
                _multiply(diagonal, leading[-1], scale),
                _multiply(off2, leading[-2], scale),
            )
        )

    trailing = [zero for _ in range(n + 1)]
    trailing[n], trailing[n - 1] = one, end
    for start in reversed(range(n - 1)):
        diagonal = end if start == 0 else inner
        trailing[start] = _subtract(
            _multiply(diagonal, trailing[start + 1], scale),
            _multiply(off2, trailing[start + 2], scale),
        )

    a_rho = _scale(rho, a)
    sums = []
    for lo, hi in ((0, split), (split, n)):
        partial, prefix_sum, diagonal_sum = zero, zero, zero
        for j in range(lo, hi):
            partial = _add(_multiply(a_rho, partial, scale), leading[j])
            prefix_sum = _add(prefix_sum, _multiply(partial, trailing[j + 1], scale))
            diagonal_sum = _add(diagonal_sum, _multiply(leading[j], trailing[j + 1], scale))
        sums.append(_subtract(_scale(prefix_sum, 2), diagonal_sum))

    partial, cross_sum = zero, zero
    for j in range(n):
        partial = _multiply(a_rho, partial, scale)
        if j < split:
            partial = _add(partial, leading[j])
        else:
            cross_sum = _add(cross_sum, _multiply(partial, trailing[j + 1], scale))

    full_det = leading[n]
    shift = _scale(q, b)
    h00 = _subtract(_scale(full_det, split), _multiply(shift, sums[0], scale))
    h11 = _subtract(_scale(full_det, n - split), _multiply(shift, sums[1], scale))
    h01 = _scale(_multiply(shift, cross_sum, scale), -1)
    numerator = _subtract(_multiply(h00, h11, scale), _multiply(h01, h01, scale))
    denominator = _scale(_multiply(q, full_det, scale), F(a * a * split * (n - split)))
    return numerator, denominator


def _bernstein(polynomial, left, right, scale):
    degree = len(polynomial) - 1
    width = right - left
    left_denominator, width_denominator = left.denominator, width.denominator
    common_denominator = left_denominator**degree * width_denominator**degree
    left_powers = [left.numerator**i for i in range(degree + 1)]
    width_powers = [width.numerator**i for i in range(degree + 1)]
    left_denominator_powers = [left_denominator**i for i in range(degree + 1)]
    width_denominator_powers = [width_denominator**i for i in range(degree + 1)]
    transformed = []
    for i in range(degree + 1):
        lower, upper = 0, 0
        for j in range(i, degree + 1):
            weight = (
                comb(j, i)
                * left_powers[j - i]
                * width_powers[i]
                * left_denominator_powers[degree - (j - i)]
                * width_denominator_powers[degree - i]
            )
            lo, hi = polynomial[j]
            lower += weight * (lo if weight >= 0 else hi)
            upper += weight * (hi if weight >= 0 else lo)
        transformed.append((lower, upper))

    bernstein_denominator = lcm(*(comb(degree, i) for i in range(degree + 1)))
    result = []
    for k in range(degree + 1):
        lower, upper = 0, 0
        for i in range(k + 1):
            weight = comb(k, i) * (bernstein_denominator // comb(degree, i))
            lower += weight * transformed[i][0]
            upper += weight * transformed[i][1]
        denominator = common_denominator * bernstein_denominator * scale
        result.append(
            (_floor_ratio(lower * scale, denominator), _ceil_ratio(upper * scale, denominator))
        )
    return tuple(result)


def determinant_interval(n, split, left, right, a, b, *, bits=192):
    """Return an outward-rounded Bernstein lower bound, or ``None``."""
    left, right, a, b = map(F, (left, right, a, b))
    if (
        not isinstance(n, int)
        or isinstance(n, bool)
        or not isinstance(split, int)
        or isinstance(split, bool)
        or n < 3
        or not 1 <= split < n
        or a <= 0
        or b <= 0
    ):
        raise ValueError('Require a valid residual split and positive coefficients')
    if not 0 <= left < right < 1:
        raise ValueError('Require 0 <= left < right < 1')
    if not isinstance(bits, int) or isinstance(bits, bool) or bits < 1:
        raise ValueError('Require positive interval precision')
    scale = 1 << bits
    center, radius = (left + right) / 2, (right - left) / 2
    numerator, denominator = _continuants(n, split, a, b, center, scale)
    numerator = _bernstein(numerator, -radius, radius, scale)
    denominator = _bernstein(denominator, -radius, radius, scale)
    degree = max(len(numerator), len(denominator)) - 1
    numerator = _elevate(numerator, degree, scale)
    denominator = _elevate(denominator, degree, scale)
    if any(lower <= 0 for lower, _ in denominator) or any(lower <= 0 for lower, _ in numerator):
        return None
    return min(F(lower, upper) for (lower, _), (_, upper) in zip(numerator, denominator))


def _elevate(coefficients, target_degree, scale):
    result = tuple(coefficients)
    while len(result) - 1 < target_degree:
        degree = len(result) - 1
        elevated = [result[0], result[-1]]
        for index in range(1, degree + 1):
            lower = index * result[index - 1][0] + (degree + 1 - index) * result[index][0]
            upper = index * result[index - 1][1] + (degree + 1 - index) * result[index][1]
            elevated.insert(-1, (_floor_ratio(lower, degree + 1), _ceil_ratio(upper, degree + 1)))
        result = tuple(elevated)
    return result
