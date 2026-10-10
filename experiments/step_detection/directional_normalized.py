"""Opt-in polynomial certificate after cancelling the stationary AR(1) scale.

The radius and covariance are both divided by q = 1-rho**2. This preserves
their ratio while avoiding independent interval bounds on a shared factor.
"""

from fractions import Fraction as F

from . import determinant_interval_polynomial as polynomial
from . import directional_determinant as determinant
from . import directional_sturm as sturm
from . import directional_tail as tail
from . import rational_polynomial as rational
from . import residual_direction as direction

BITS = 128
RADIUS_BITS = determinant.BITS


def radius_upper(model, left, right):
    """Bound r/q directly, retaining cancellation of q in C/r."""

    def enclosure(poly):
        return (
            (rational.evaluate(poly, left),)
            if left == right
            else rational.bernstein(poly, left, right)
        )

    numerator = min(enclosure(model['num']))
    denominator = enclosure(model['den'])
    if numerator <= 0 or min(denominator) <= 0 or model['ss'] <= 0:
        raise ArithmeticError('Normalized residual ratio enclosure is inconclusive')
    return tail.ceil_dyadic(model['ss'] * max(denominator) / numerator, bits=RADIUS_BITS)


def certify_interval(model, left, right, *, delta=direction.DELTA):
    """Try the two existing tilts on a positive, interior interval.

    A failed enclosure returns a non-certificate. No subdivision, tilt search,
    floating-point acceptance, or change in the confidence allocation occurs.
    """
    left, right, delta = F(left), F(right), F(delta)
    if not 0 <= left < right < 1 or not 0 < delta < 1:
        raise ValueError('Require a positive interior interval and a probability cutoff')
    try:
        radius = radius_upper(model, left, right)
    except ArithmeticError:
        radius = None

    def certifier(_model, lower, upper, *, tilt, delta, bits):
        assert (lower, upper, bits) == (left, right, BITS)
        base = {
            'left': str(left),
            'right': str(right),
            'tilt': str(tilt),
            'delta': str(delta),
            'bits': bits,
            'radius_bits': RADIUS_BITS,
            'stationary_normalized': True,
        }
        if radius is None:
            return {**base, 'status': 'radius_inconclusive'}
        base['radius_upper'] = str(radius)
        lower_bound = polynomial.determinant_interval(
            model['n'],
            model['split'],
            left,
            right,
            1 - 2 * tilt,
            2 * tilt / radius,
            bits=bits,
            stationary_normalized=True,
        )
        if lower_bound is None:
            return {**base, 'status': 'determinant_inconclusive'}
        return {
            **base,
            'status': 'rounded_polynomial_lower_bound',
            'determinant_enclosure': [str(lower_bound), str(lower_bound)],
        }

    proof = sturm.certify_interval(
        model,
        left,
        right,
        delta=delta,
        bits=BITS,
        determinant_certifier=certifier,
        rank_groups_on_intervals=True,
    )
    return {**proof, 'stationary_normalized': True, 'radius_bits': RADIUS_BITS}
