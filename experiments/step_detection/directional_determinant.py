"""Interior-correlation directional certificates from a full determinant.

This is a standalone mathematical reference. It does not change the frozen
reporting search. All interval arithmetic rounds outward to dyadic rationals.
"""

from fractions import Fraction as F

from . import directional_tail as tail
from . import rational_polynomial as p
from . import residual_direction as direction

BITS = 192


class Arithmetic:
    """Rational interval operations with bounded denominator sizes."""

    def __init__(self, bits=BITS):
        if not isinstance(bits, int) or isinstance(bits, bool) or bits < 1:
            raise ValueError('Require positive integer precision')
        self.scale = 1 << bits

    def interval(self, lower, upper=None):
        lower = F(lower)
        upper = lower if upper is None else F(upper)
        if lower > upper:
            raise ValueError('Reversed interval')
        return (
            F(lower.numerator * self.scale // lower.denominator, self.scale),
            F(-(-upper.numerator * self.scale // upper.denominator), self.scale),
        )

    def add(self, x, y):
        return self.interval(x[0] + y[0], x[1] + y[1])

    def sub(self, x, y):
        return self.interval(x[0] - y[1], x[1] - y[0])

    def mul(self, x, y):
        products = [a * b for a in x for b in y]
        return self.interval(min(products), max(products))

    def div(self, x, y):
        if y[0] <= 0 <= y[1]:
            raise ArithmeticError('Denominator enclosure contains zero')
        return self.mul(x, self.interval(1 / y[1], 1 / y[0]))

    def square(self, x):
        lower = 0 if x[0] <= 0 <= x[1] else min(x[0] ** 2, x[1] ** 2)
        return self.interval(lower, max(x[0] ** 2, x[1] ** 2))


def determinant_interval(n, split, left, right, a, b, *, bits=BITS):
    """Enclose det(a I + b U' R_rho U) for every rho in [left, right].

    U is any orthonormal complement of the two plateau indicators. A failed
    pivot enclosure raises ArithmeticError; subdivision can resolve it.
    """
    left, right, a, b = map(F, (left, right, a, b))
    if (
        not isinstance(n, int)
        or isinstance(n, bool)
        or not isinstance(split, int)
        or isinstance(split, bool)
        or n < 3
        or not 1 <= split < n
        or not -1 < left <= right < 1
        or a <= 0
        or b <= 0
    ):
        raise ValueError('Require an interior split/interval and positive coefficients')
    ar = Arithmetic(bits)
    zero, one = ar.interval(0), ar.interval(1)
    ai, bi, rho = ar.interval(a), ar.interval(b), ar.interval(left, right)
    rho2 = ar.square(rho)
    stationary = ar.sub(one, rho2)
    if stationary[0] <= 0:
        raise ArithmeticError('Stationary factor enclosure is nonpositive')
    shift = ar.mul(bi, stationary)
    endpoint = ar.add(ai, shift)
    interior = ar.add(ar.mul(ai, ar.add(one, rho2)), shift)
    off = ar.mul(ar.sub(zero, ai), rho)

    # T = a A_rho + b (1-rho**2) I; A_rho is the AR(1) precision numerator.
    pivots, multipliers = [endpoint], []
    determinant = endpoint
    for i in range(1, n):
        if pivots[-1][0] <= 0:
            raise ArithmeticError('LDL pivot enclosure is nonpositive')
        multiplier = ar.div(off, pivots[-1])
        pivot = ar.sub(endpoint if i == n - 1 else interior, ar.mul(multiplier, off))
        multipliers.append(multiplier)
        pivots.append(pivot)
        determinant = ar.mul(determinant, pivot)
    if pivots[-1][0] <= 0:
        raise ArithmeticError('LDL pivot enclosure is nonpositive')

    def solve(start, stop):
        forward = []
        for i in range(n):
            rhs = one if start <= i < stop else zero
            forward.append(rhs if i == 0 else ar.sub(rhs, ar.mul(multipliers[i - 1], forward[-1])))
        result = [zero] * n
        for i in reversed(range(n)):
            value = ar.div(forward[i], pivots[i])
            result[i] = (
                value if i == n - 1 else ar.sub(value, ar.mul(multipliers[i], result[i + 1]))
            )
        return result

    blocks = ((0, split), (split, n))
    gram = [[zero, zero], [zero, zero]]
    for j, (start, stop) in enumerate(blocks):
        solved = solve(start, stop)
        for i, (lo, hi) in enumerate(blocks):
            value = zero
            for item in solved[lo:hi]:
                value = ar.add(value, item)
            diagonal = ar.interval(hi - lo) if i == j else zero
            # M^-1 = (I - b (1-rho**2) T^-1) / a.
            gram[i][j] = ar.div(ar.sub(diagonal, ar.mul(shift, value)), ai)
    gram_det = ar.sub(ar.mul(gram[0][0], gram[1][1]), ar.mul(gram[0][1], gram[1][0]))
    return ar.div(
        ar.mul(determinant, gram_det),
        ar.mul(stationary, ar.interval(split * (n - split))),
    )


def certify_interval(model, left, right, *, tilt=tail.TILT, delta=direction.DELTA, bits=BITS):
    """Bound the same directional p-value using the complete covariance.

    Only stationary interior intervals are supported. An inconclusive bound
    returns one and never establishes that an explanation truly survives.
    """
    left, right, tilt, delta = map(F, (left, right, tilt, delta))
    if not -1 < left <= right < 1 or not 0 < tilt < F(1, 2) or not 0 < delta < 1:
        raise ValueError('Require an interior interval and valid tilt/budget')
    Arithmetic(bits)
    result = {
        'left': str(left),
        'right': str(right),
        'tilt': str(tilt),
        'delta': str(delta),
        'bits': bits,
        'status': 'not_certified',
        'p_upper_squared': '1',
        'p_upper': '1',
    }

    def enclosure(poly):
        coefficients = (
            (p.evaluate(poly, left),) if left == right else p.bernstein(poly, left, right)
        )
        return min(coefficients), max(coefficients)

    numerator, denominator = enclosure(model['num']), enclosure(model['den'])
    if numerator[0] <= 0 or denominator[0] <= 0 or model['ss'] <= 0:
        result['reason'] = 'Residual ratio enclosure is inconclusive'
        return result
    stationary_upper = 1 if left <= 0 <= right else 1 - min(left**2, right**2)
    radius_upper = tail.ceil_dyadic(
        stationary_upper * model['ss'] * denominator[1] / numerator[0], bits=bits
    )
    result['radius_upper'] = str(radius_upper)
    try:
        lower, upper = determinant_interval(
            model['n'],
            model['split'],
            left,
            right,
            1 - 2 * tilt,
            2 * tilt / radius_upper,
            bits=bits,
        )
    except ArithmeticError as error:
        result['reason'] = str(error)
        return result
    result['determinant_enclosure'] = [str(lower), str(upper)]
    if lower <= 0:
        result['reason'] = 'Determinant enclosure is nonpositive'
        return result
    squared = min(F(1), 1 / lower)
    upper_p = min(F(1), direction.sqrt_lower(squared, bits=48) + F(1, 2**48))
    result.update(
        p_upper_squared=str(squared),
        p_upper=str(upper_p),
        status='certified_excluded' if squared < delta**2 else 'not_certified',
    )
    return result
