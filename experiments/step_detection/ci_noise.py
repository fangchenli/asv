"""Exact algebra and fixed-sample sign bounds for the CI noise investigation.

This is an analytical reference, not a CI runner or a production decision rule.
See ci_noise.rst for the measurement model and assumptions behind the bounds.
"""

from fractions import Fraction
from math import comb
from operator import index


def balanced_log_contrast(log_times, order='ABBA'):
    """Return mean(B log times) - mean(A log times) for four fixed slots.

    Inputs are already log times. Fractions allow exact checks of the linear
    identities without claiming that logarithms of actual timings are rational.
    A is the baseline; B is the candidate. Both orientations get equal weight.
    """
    if order not in ('ABBA', 'BAAB'):
        raise ValueError('order must be ABBA or BAAB')
    values = tuple(Fraction(value) for value in log_times)
    if len(values) != 4:
        raise ValueError('exactly four log times are required')
    return (
        sum(
            (value if version == 'B' else -value for value, version in zip(values, order)),
            Fraction(0),
        )
        / 2
    )


def sign_tail_bound(jobs, exceedances, contamination=Fraction(0)):
    """Bound P(K >= exceedances) under a fixed-sample one-sided null.

    Independent jobs have clean exceedance probability at most 1/2 and each
    is contaminated with probability at most epsilon. Contamination can force
    an exceedance, giving p <= (1 + epsilon)/2. Ties count as non-exceedances;
    they remain in jobs. This is not a bound for arbitrarily dependent jobs
    or for an unknown contamination probability estimated from the same data.
    """
    jobs, exceedances = index(jobs), index(exceedances)
    epsilon = Fraction(contamination)
    if jobs < 1 or not 0 <= exceedances <= jobs + 1:
        raise ValueError('require jobs >= 1 and 0 <= exceedances <= jobs + 1')
    if not 0 <= epsilon <= 1:
        raise ValueError('contamination must be between zero and one')
    p = (1 + epsilon) / 2
    return sum(
        (comb(jobs, k) * p**k * (1 - p) ** (jobs - k) for k in range(exceedances, jobs + 1)),
        Fraction(0),
    )


def required_exceedances(jobs, alpha=Fraction(1, 20), contamination=Fraction(0)):
    """Smallest attainable count with tail <= alpha, or None if unattainable."""
    alpha = Fraction(alpha)
    if not 0 < alpha < 1:
        raise ValueError('alpha must lie strictly between zero and one')
    # Validate even if no count will meet the requested budget.
    sign_tail_bound(jobs, 0, contamination)
    for k in range(1, jobs + 1):
        if sign_tail_bound(jobs, k, contamination) <= alpha:
            return k
    return None


def main():
    print('Fixed one-sided 5% budget; counts are independent jobs, not samples.')
    print('jobs  clean required  10% contamination required  all-positive bound (10%)')
    for jobs in (3, 5, 6, 10, 20):
        clean = required_exceedances(jobs)
        robust = required_exceedances(jobs, contamination=Fraction(1, 10))
        bound = sign_tail_bound(jobs, jobs, Fraction(1, 10))
        print(f'{jobs:4}  {clean!s:>14}  {robust!s:>26}  {float(bound):.8f}')


if __name__ == '__main__':
    main()
