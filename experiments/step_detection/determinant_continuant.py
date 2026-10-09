"""Centered polynomial enclosure for the projected AR(1) determinant.

This is an experimental alternative to interval LDL arithmetic. It uses the
tridiagonal continuant formula and exact rational polynomials in correlation.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb

from . import rational_polynomial as poly


def _continuants(n, a, b):
    q = poly.polynomial((1, 0, -1))
    rho2 = poly.polynomial((0, 0, 1))
    end = poly.add(poly.polynomial((a,)), poly.scale(q, b))
    inner = poly.add(poly.scale(poly.add(poly.polynomial((1,)), rho2), a), poly.scale(q, b))
    off2 = poly.scale(rho2, a * a)

    leading = [poly.polynomial((1,)), end]
    for size in range(2, n + 1):
        diagonal = end if size == n else inner
        leading.append(
            poly.subtract(poly.multiply(diagonal, leading[-1]), poly.multiply(off2, leading[-2]))
        )

    trailing = [poly.polynomial((0,)) for _ in range(n + 1)]
    trailing[n] = poly.polynomial((1,))
    trailing[n - 1] = end
    for start in reversed(range(n - 1)):
        diagonal = end if start == 0 else inner
        trailing[start] = poly.subtract(
            poly.multiply(diagonal, trailing[start + 1]),
            poly.multiply(off2, trailing[start + 2]),
        )
    return q, leading, trailing


def _inverse_numerator(i, j, a, leading, trailing):
    """Numerator of T^-1[i,j], assuming nonnegative correlation and i <= j."""
    distance = j - i
    factor = poly.polynomial((0,) * distance + (a**distance,))
    return poly.multiply(poly.multiply(factor, leading[i]), trailing[j + 1])


@lru_cache(maxsize=32)
def determinant_polynomial(n, split, a, b):
    """Return numerator/denominator polynomials for det(aI+b U'R U)."""
    if n < 3 or not 1 <= split < n or a <= 0 or b <= 0:
        raise ValueError('Require a valid residual split and positive coefficients')
    q, leading, trailing = _continuants(n, F(a), F(b))
    full_det = leading[n]
    blocks = ((0, split), (split, n))
    sums = [[poly.polynomial((0,)) for _ in range(2)] for _ in range(2)]
    for x, (lo_x, hi_x) in enumerate(blocks):
        for y, (lo_y, hi_y) in enumerate(blocks):
            if x > y:
                sums[x][y] = sums[y][x]
                continue
            total = poly.polynomial((0,))
            if x == y:
                for i in range(lo_x, hi_x):
                    total = poly.add(total, _inverse_numerator(i, i, F(a), leading, trailing))
                    for j in range(i + 1, hi_x):
                        total = poly.add(
                            total, poly.scale(_inverse_numerator(i, j, F(a), leading, trailing), 2)
                        )
            else:
                for i in range(lo_x, hi_x):
                    for j in range(lo_y, hi_y):
                        total = poly.add(total, _inverse_numerator(i, j, F(a), leading, trailing))
            sums[x][y] = total

    shift = poly.scale(q, F(b))
    h00 = poly.subtract(poly.scale(full_det, split), poly.multiply(shift, sums[0][0]))
    h11 = poly.subtract(poly.scale(full_det, n - split), poly.multiply(shift, sums[1][1]))
    h01 = poly.scale(poly.multiply(shift, sums[0][1]), -1)
    numerator = poly.subtract(poly.multiply(h00, h11), poly.multiply(h01, h01))
    denominator = poly.scale(poly.multiply(q, full_det), F(a * a * split * (n - split)))
    return numerator, denominator


def _centered_bounds(coefficients, center, radius):
    """Bound a polynomial on center +/- radius by its exact Taylor coefficients."""
    shifted = [
        sum(coefficients[j] * comb(j, k) * center ** (j - k) for j in range(k, len(coefficients)))
        for k in range(len(coefficients))
    ]
    variation = sum(abs(value) * radius**k for k, value in enumerate(shifted) if k)
    return shifted[0] - variation, shifted[0] + variation


def determinant_interval(n, split, left, right, a, b):
    """Return a rigorous centered-polynomial lower bound, or ``None``."""
    left, right = F(left), F(right)
    if left < 0 or right >= 1 or left > right:
        raise ValueError('Require 0 <= left <= right < 1')
    numerator, denominator = determinant_polynomial(n, split, F(a), F(b))
    center = (left + right) / 2
    radius = (right - left) / 2
    n_lower, _ = _centered_bounds(numerator, center, radius)
    _, d_upper = _centered_bounds(denominator, center, radius)
    if n_lower <= 0 or d_upper <= 0:
        return None
    return n_lower / d_upper
