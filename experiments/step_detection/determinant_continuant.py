"""Centered polynomial enclosure for the projected AR(1) determinant.

This is an experimental alternative to interval LDL arithmetic. It uses the
tridiagonal continuant formula and exact rational polynomials in correlation.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, lcm

from . import rational_polynomial as poly


def _multiply(first, second):
    """Multiply rational polynomials with integer convolution and shared scales."""
    first_scale = lcm(*(value.denominator for value in first))
    second_scale = lcm(*(value.denominator for value in second))
    left = [int(value * first_scale) for value in first]
    right = [int(value * second_scale) for value in second]
    product = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            product[i + j] += x * y
    denominator = first_scale * second_scale
    return poly.polynomial(F(value, denominator) for value in product)


def _continuants(n, a, b, center):
    rho = poly.polynomial((center, 1))
    rho2 = _multiply(rho, rho)
    q = poly.subtract(poly.polynomial((1,)), rho2)
    end = poly.add(poly.polynomial((a,)), poly.scale(q, b))
    inner = poly.add(poly.scale(poly.add(poly.polynomial((1,)), rho2), a), poly.scale(q, b))
    off2 = poly.scale(rho2, a * a)

    leading = [poly.polynomial((1,)), end]
    for size in range(2, n + 1):
        diagonal = end if size == n else inner
        leading.append(
            poly.subtract(_multiply(diagonal, leading[-1]), _multiply(off2, leading[-2]))
        )

    trailing = [poly.polynomial((0,)) for _ in range(n + 1)]
    trailing[n] = poly.polynomial((1,))
    trailing[n - 1] = end
    for start in reversed(range(n - 1)):
        diagonal = end if start == 0 else inner
        trailing[start] = poly.subtract(
            _multiply(diagonal, trailing[start + 1]),
            _multiply(off2, trailing[start + 2]),
        )
    return q, leading, trailing, rho


def _within_block_sum(lo, hi, a_rho, leading, trailing):
    """Numerator of the sum of inverse entries within one contiguous block."""
    partial = poly.polynomial((0,))
    prefix_sum = poly.polynomial((0,))
    diagonal_sum = poly.polynomial((0,))
    for j in range(lo, hi):
        partial = poly.add(_multiply(a_rho, partial), leading[j])
        prefix_sum = poly.add(prefix_sum, _multiply(partial, trailing[j + 1]))
        diagonal_sum = poly.add(diagonal_sum, _multiply(leading[j], trailing[j + 1]))
    return poly.subtract(poly.scale(prefix_sum, 2), diagonal_sum)


def _between_block_sum(lo_a, hi_a, lo_b, hi_b, a_rho, leading, trailing):
    """Numerator of the inverse-entry sum between ordered disjoint blocks."""
    partial = poly.polynomial((0,))
    total = poly.polynomial((0,))
    for j in range(lo_a, hi_b):
        partial = _multiply(a_rho, partial)
        if lo_a <= j < hi_a:
            partial = poly.add(partial, leading[j])
        if lo_b <= j < hi_b:
            total = poly.add(total, _multiply(partial, trailing[j + 1]))
    return total


@lru_cache(maxsize=32)
def determinant_polynomial(n, split, a, b, center=F(0)):
    """Return numerator/denominator polynomials for det(aI+b U'R U)."""
    if n < 3 or not 1 <= split < n or a <= 0 or b <= 0:
        raise ValueError('Require a valid residual split and positive coefficients')
    q, leading, trailing, rho = _continuants(n, F(a), F(b), F(center))
    full_det = leading[n]
    blocks = ((0, split), (split, n))
    a_rho = poly.scale(rho, F(a))
    sums = [_within_block_sum(*block, a_rho, leading, trailing) for block in blocks]
    cross_sum = _between_block_sum(*blocks[0], *blocks[1], a_rho, leading, trailing)

    shift = poly.scale(q, F(b))
    h00 = poly.subtract(poly.scale(full_det, split), poly.multiply(shift, sums[0]))
    h11 = poly.subtract(poly.scale(full_det, n - split), poly.multiply(shift, sums[1]))
    h01 = poly.scale(poly.multiply(shift, cross_sum), -1)
    numerator = poly.subtract(_multiply(h00, h11), _multiply(h01, h01))
    denominator = poly.scale(_multiply(q, full_det), F(a * a * split * (n - split)))
    return numerator, denominator


def _centered_bounds(coefficients, radius):
    """Bound a polynomial on ``[-radius, radius]`` by its Taylor coefficients."""
    variation = sum(abs(value) * radius**k for k, value in enumerate(coefficients) if k)
    return coefficients[0] - variation, coefficients[0] + variation


def _bernstein_coefficients(coefficients, left, right):
    """Convert power coefficients to Bernstein coefficients on [left,right]."""
    degree = len(coefficients) - 1
    width = right - left
    power_coefficients = [
        sum(
            coefficients[j] * comb(j, i) * left ** (j - i) * width**i for j in range(i, degree + 1)
        )
        for i in range(degree + 1)
    ]
    return tuple(
        sum(power_coefficients[i] * F(comb(k, i), comb(degree, i)) for i in range(k + 1))
        for k in range(degree + 1)
    )


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
    numerator, denominator = determinant_polynomial(n, split, F(a), F(b), center)
    n_lower = max(_centered_bounds(numerator, radius)[0], _bernstein_bounds(numerator, radius)[0])
    d_upper = min(
        _centered_bounds(denominator, radius)[1], _bernstein_bounds(denominator, radius)[1]
    )
    if n_lower <= 0 or d_upper <= 0:
        return None
    return n_lower / d_upper
