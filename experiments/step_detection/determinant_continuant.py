"""Centered polynomial enclosure for the projected AR(1) determinant.

This is an experimental alternative to interval LDL arithmetic. It uses the
tridiagonal continuant formula and exact rational polynomials in correlation.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, gcd, lcm

from . import rational_polynomial as poly


def _scaled_polynomial(coefficients, denominator=1):
    """Canonical integer coefficients with one shared positive denominator."""
    coefficients = list(coefficients)
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    divisor = denominator
    for coefficient in coefficients:
        divisor = gcd(divisor, abs(coefficient))
    divisor = max(divisor, 1)
    return tuple(value // divisor for value in coefficients), denominator // divisor


def _sp_constant(value):
    value = F(value)
    return _scaled_polynomial((value.numerator,), value.denominator)


def _sp_add(first, second):
    first_coefficients, first_denominator = first
    second_coefficients, second_denominator = second
    common = gcd(first_denominator, second_denominator)
    first_scale = second_denominator // common
    second_scale = first_denominator // common
    length = max(len(first_coefficients), len(second_coefficients))
    coefficients = [0] * length
    for index in range(length):
        left = first_coefficients[index] if index < len(first_coefficients) else 0
        right = second_coefficients[index] if index < len(second_coefficients) else 0
        coefficients[index] = left * first_scale + right * second_scale
    return _scaled_polynomial(coefficients, first_denominator * first_scale)


def _sp_scale(polynomial, factor):
    factor = F(factor)
    return _scaled_polynomial(
        (value * factor.numerator for value in polynomial[0]),
        polynomial[1] * factor.denominator,
    )


def _sp_multiply(first, second):
    left, right = first[0], second[0]
    product = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            product[i + j] += x * y
    return _scaled_polynomial(product, first[1] * second[1])


def _sp_subtract(first, second):
    return _sp_add(first, _sp_scale(second, -1))


def _sp_continuants(n, a, b, center):
    center = F(center)
    rho = _scaled_polynomial((center.numerator, center.denominator), center.denominator)
    rho2 = _sp_multiply(rho, rho)
    q = _sp_subtract(_sp_constant(1), rho2)
    end = _sp_add(_sp_constant(a), _sp_scale(q, b))
    inner = _sp_add(
        _sp_scale(_sp_add(_sp_constant(1), rho2), a),
        _sp_scale(q, b),
    )
    off2 = _sp_scale(rho2, F(a) ** 2)

    leading = [_sp_constant(1), end]
    for size in range(2, n + 1):
        diagonal = end if size == n else inner
        leading.append(
            _sp_subtract(
                _sp_multiply(diagonal, leading[-1]),
                _sp_multiply(off2, leading[-2]),
            )
        )

    trailing = [_sp_constant(0) for _ in range(n + 1)]
    trailing[n] = _sp_constant(1)
    trailing[n - 1] = end
    for start in reversed(range(n - 1)):
        diagonal = end if start == 0 else inner
        trailing[start] = _sp_subtract(
            _sp_multiply(diagonal, trailing[start + 1]),
            _sp_multiply(off2, trailing[start + 2]),
        )
    return q, leading, trailing, rho


def _sp_within_block_sum(lo, hi, a_rho, leading, trailing):
    partial = _sp_constant(0)
    prefix_sum = _sp_constant(0)
    diagonal_sum = _sp_constant(0)
    for j in range(lo, hi):
        partial = _sp_add(_sp_multiply(a_rho, partial), leading[j])
        prefix_sum = _sp_add(prefix_sum, _sp_multiply(partial, trailing[j + 1]))
        diagonal_sum = _sp_add(diagonal_sum, _sp_multiply(leading[j], trailing[j + 1]))
    return _sp_subtract(_sp_scale(prefix_sum, 2), diagonal_sum)


def _sp_between_block_sum(lo_a, hi_a, lo_b, hi_b, a_rho, leading, trailing):
    partial = _sp_constant(0)
    total = _sp_constant(0)
    for j in range(lo_a, hi_b):
        partial = _sp_multiply(a_rho, partial)
        if lo_a <= j < hi_a:
            partial = _sp_add(partial, leading[j])
        if lo_b <= j < hi_b:
            total = _sp_add(total, _sp_multiply(partial, trailing[j + 1]))
    return total


def _sp_to_fraction(polynomial):
    return poly.polynomial(F(value, polynomial[1]) for value in polynomial[0])


def _tm_constant(value, degree):
    return (F(value),) + (F(0),) * degree, F(0)


def _tm_add(left, right):
    return (
        tuple(a + b for a, b in zip(left[0], right[0])),
        left[1] + right[1],
    )


def _tm_scale(value, factor):
    factor = F(factor)
    return tuple(coefficient * factor for coefficient in value[0]), value[1] * abs(factor)


def _tm_multiply(left, right, radius):
    degree = len(left[0]) - 1
    coefficients = [F(0)] * (degree + 1)
    tail = F(0)
    for i, first in enumerate(left[0]):
        for j, second in enumerate(right[0]):
            product = first * second
            if i + j <= degree:
                coefficients[i + j] += product
            else:
                tail += abs(product) * radius ** (i + j)
    left_bound = sum(abs(value) * radius**i for i, value in enumerate(left[0]))
    right_bound = sum(abs(value) * radius**i for i, value in enumerate(right[0]))
    error = tail + left_bound * right[1] + right_bound * left[1] + left[1] * right[1]
    return tuple(coefficients), error


def _tm_range(value, radius):
    variation = sum(abs(coefficient) * radius**i for i, coefficient in enumerate(value[0]) if i)
    return value[0][0] - variation - value[1], value[0][0] + variation + value[1]


def _tm_determinant_polynomials(n, split, a, b, center, radius, degree):
    zero, one = _tm_constant(0, degree), _tm_constant(1, degree)
    rho = ((F(center), F(1)) + (F(0),) * (degree - 1), F(0))
    rho2 = _tm_multiply(rho, rho, radius)
    q = _tm_add(one, _tm_scale(rho2, -1))
    end = _tm_add(_tm_constant(a, degree), _tm_scale(q, b))
    inner = _tm_add(
        _tm_scale(_tm_add(one, rho2), a),
        _tm_scale(q, b),
    )
    off2 = _tm_scale(rho2, F(a) ** 2)

    leading = [one, end]
    for size in range(2, n + 1):
        diagonal = end if size == n else inner
        leading.append(
            _tm_add(
                _tm_multiply(diagonal, leading[-1], radius),
                _tm_scale(_tm_multiply(off2, leading[-2], radius), -1),
            )
        )

    trailing = [zero for _ in range(n + 1)]
    trailing[n] = one
    trailing[n - 1] = end
    for start in reversed(range(n - 1)):
        diagonal = end if start == 0 else inner
        trailing[start] = _tm_add(
            _tm_multiply(diagonal, trailing[start + 1], radius),
            _tm_scale(_tm_multiply(off2, trailing[start + 2], radius), -1),
        )

    a_rho = _tm_scale(rho, a)
    blocks = ((0, split), (split, n))
    block_sums = []
    for lo, hi in blocks:
        partial, prefix_sum, diagonal_sum = zero, zero, zero
        for j in range(lo, hi):
            partial = _tm_add(_tm_multiply(a_rho, partial, radius), leading[j])
            prefix_sum = _tm_add(prefix_sum, _tm_multiply(partial, trailing[j + 1], radius))
            diagonal_sum = _tm_add(diagonal_sum, _tm_multiply(leading[j], trailing[j + 1], radius))
        block_sums.append(_tm_add(_tm_scale(prefix_sum, 2), _tm_scale(diagonal_sum, -1)))

    partial, cross_sum = zero, zero
    for j in range(n):
        partial = _tm_multiply(a_rho, partial, radius)
        if j < split:
            partial = _tm_add(partial, leading[j])
        elif j >= split:
            cross_sum = _tm_add(cross_sum, _tm_multiply(partial, trailing[j + 1], radius))

    full_det = leading[n]
    shift = _tm_scale(q, b)
    h00 = _tm_add(
        _tm_scale(full_det, split), _tm_scale(_tm_multiply(shift, block_sums[0], radius), -1)
    )
    h11 = _tm_add(
        _tm_scale(full_det, n - split), _tm_scale(_tm_multiply(shift, block_sums[1], radius), -1)
    )
    h01 = _tm_scale(_tm_multiply(shift, cross_sum, radius), -1)
    numerator = _tm_add(
        _tm_multiply(h00, h11, radius),
        _tm_scale(_tm_multiply(h01, h01, radius), -1),
    )
    denominator = _tm_scale(
        _tm_multiply(q, full_det, radius),
        F(a * a * split * (n - split)),
    )
    return numerator, denominator


def determinant_interval_taylor(n, split, left, right, a, b, degree=4):
    """Return a rigorous positive determinant-ratio lower bound using a
    degree-limited Taylor model, or ``None`` when the enclosure is inconclusive.
    """
    left, right, a, b = F(left), F(right), F(a), F(b)
    if n < 3 or not 1 <= split < n or a <= 0 or b <= 0:
        raise ValueError('Require a valid residual split and positive coefficients')
    if left < 0 or right >= 1 or left > right:
        raise ValueError('Require 0 <= left <= right < 1')
    if not isinstance(degree, int) or degree < 1:
        raise ValueError('Taylor degree must be a positive integer')
    center, radius = (left + right) / 2, (right - left) / 2
    numerator, denominator = _tm_determinant_polynomials(n, split, a, b, center, radius, degree)
    n_lower, _ = _tm_range(numerator, radius)
    d_lower, d_upper = _tm_range(denominator, radius)
    if n_lower <= 0 or d_lower <= 0:
        return None
    return n_lower / d_upper


@lru_cache(maxsize=32)
def _determinant_polynomials_scaled(n, split, a, b, center):
    if n < 3 or not 1 <= split < n or a <= 0 or b <= 0:
        raise ValueError('Require a valid residual split and positive coefficients')
    a, b = F(a), F(b)
    q, leading, trailing, rho = _sp_continuants(n, a, b, F(center))
    full_det = leading[n]
    blocks = ((0, split), (split, n))
    a_rho = _sp_scale(rho, a)
    sums = [_sp_within_block_sum(*block, a_rho, leading, trailing) for block in blocks]
    cross_sum = _sp_between_block_sum(*blocks[0], *blocks[1], a_rho, leading, trailing)

    shift = _sp_scale(q, b)
    h00 = _sp_subtract(_sp_scale(full_det, split), _sp_multiply(shift, sums[0]))
    h11 = _sp_subtract(_sp_scale(full_det, n - split), _sp_multiply(shift, sums[1]))
    h01 = _sp_scale(_sp_multiply(shift, cross_sum), -1)
    numerator = _sp_subtract(_sp_multiply(h00, h11), _sp_multiply(h01, h01))
    denominator = _sp_scale(_sp_multiply(q, full_det), F(a * a * split * (n - split)))
    return numerator, denominator


@lru_cache(maxsize=32)
def determinant_polynomial(n, split, a, b, center=F(0)):
    """Return numerator/denominator polynomials for det(aI+b U'R U)."""
    scaled = _determinant_polynomials_scaled(n, split, F(a), F(b), F(center))
    return tuple(_sp_to_fraction(polynomial) for polynomial in scaled)


def _centered_bounds_scaled(polynomial, radius):
    coefficients, denominator = polynomial
    variation = sum(
        F(abs(value), denominator) * radius**degree
        for degree, value in enumerate(coefficients)
        if degree
    )
    constant = F(coefficients[0], denominator)
    return constant - variation, constant + variation


def _bernstein_integer_coefficients(integers, coefficient_denominator, left, right):
    degree = len(integers) - 1
    width = right - left
    left_denominator, width_denominator = left.denominator, width.denominator
    common_denominator = (
        coefficient_denominator * left_denominator**degree * width_denominator**degree
    )
    left_powers = [left.numerator**i for i in range(degree + 1)]
    width_powers = [width.numerator**i for i in range(degree + 1)]
    left_denominator_powers = [left_denominator**i for i in range(degree + 1)]
    width_denominator_powers = [width_denominator**i for i in range(degree + 1)]
    power_numerators = []
    for i in range(degree + 1):
        total = 0
        for j in range(i, degree + 1):
            total += (
                integers[j]
                * comb(j, i)
                * left_powers[j - i]
                * width_powers[i]
                * left_denominator_powers[degree - (j - i)]
                * width_denominator_powers[degree - i]
            )
        power_numerators.append(total)

    bernstein_denominator = lcm(*(comb(degree, i) for i in range(degree + 1)))
    bernstein_numerators = [
        sum(
            power_numerators[i] * comb(k, i) * (bernstein_denominator // comb(degree, i))
            for i in range(k + 1)
        )
        for k in range(degree + 1)
    ]
    denominator = common_denominator * bernstein_denominator
    return tuple(F(value, denominator) for value in bernstein_numerators)


def _bernstein_coefficients(coefficients, left, right):
    """Convert power coefficients to Bernstein coefficients on [left,right]."""
    coefficient_denominator = lcm(*(value.denominator for value in coefficients))
    integers = [int(value * coefficient_denominator) for value in coefficients]
    return _bernstein_integer_coefficients(integers, coefficient_denominator, left, right)


def _split_bernstein(coefficients):
    """Bisect Bernstein coefficients by de Casteljau's exact scheme."""
    denominator = lcm(*(value.denominator for value in coefficients))
    work = [int(value * denominator) for value in coefficients]
    degree = len(work) - 1
    left, right = [work[0] * (1 << degree)], [work[-1] * (1 << degree)]
    for level in range(1, degree + 1):
        work = [work[i] + work[i + 1] for i in range(len(work) - 1)]
        scale = 1 << (degree - level)
        left.append(work[0] * scale)
        right.append(work[-1] * scale)
    denominator *= 1 << degree
    return (
        tuple(F(value, denominator) for value in left),
        tuple(F(value, denominator) for value in reversed(right)),
    )


def _elevate_bernstein(coefficients, degree):
    """Elevate a Bernstein sequence without changing its polynomial."""
    result = tuple(coefficients)
    while len(result) - 1 < degree:
        old_degree = len(result) - 1
        new_degree = old_degree + 1
        result = tuple(
            result[0]
            if index == 0
            else result[-1]
            if index == new_degree
            else F(index, new_degree) * result[index - 1]
            + F(new_degree - index, new_degree) * result[index]
            for index in range(new_degree + 1)
        )
    return result


def _bernstein_ratio_lower(numerator, denominator):
    """Lower-bound a polynomial ratio using paired Bernstein coefficients."""
    degree = max(len(numerator), len(denominator)) - 1
    numerator = _elevate_bernstein(numerator, degree)
    denominator = _elevate_bernstein(denominator, degree)
    if any(value <= 0 for value in denominator):
        return None
    return min(n / d for n, d in zip(numerator, denominator))


def _bernstein_bounds(centered_coefficients, radius):
    """Bound a centered polynomial by Bernstein coefficients on [-radius,radius]."""
    if radius == 0:
        return centered_coefficients[0], centered_coefficients[0]
    bernstein = _bernstein_coefficients(centered_coefficients, -radius, radius)
    return min(bernstein), max(bernstein)


def determinant_interval(n, split, left, right, a, b):
    """Return a rigorous centered-polynomial lower bound, or ``None``."""
    left, right = F(left), F(right)
    if left < 0 or right >= 1 or left > right:
        raise ValueError('Require 0 <= left <= right < 1')
    center = (left + right) / 2
    radius = (right - left) / 2
    numerator, denominator = _determinant_polynomials_scaled(n, split, F(a), F(b), center)
    numerator_bernstein = _bernstein_integer_coefficients(
        numerator[0], numerator[1], -radius, radius
    )
    denominator_bernstein = _bernstein_integer_coefficients(
        denominator[0], denominator[1], -radius, radius
    )
    n_lower = max(_centered_bounds_scaled(numerator, radius)[0], min(numerator_bernstein))
    d_upper = min(_centered_bounds_scaled(denominator, radius)[1], max(denominator_bernstein))
    if n_lower <= 0 or d_upper <= 0:
        return None
    independent_ratio = n_lower / d_upper
    paired_ratio = _bernstein_ratio_lower(numerator_bernstein, denominator_bernstein)
    return max(independent_ratio, paired_ratio or F(0))
