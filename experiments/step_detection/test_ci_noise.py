"""Deterministic identities and exhaustive probability checks; no simulation."""

from fractions import Fraction
from itertools import product

import pytest

from experiments.step_detection.ci_noise import (
    balanced_log_contrast,
    required_exceedances,
    sign_tail_bound,
)


@pytest.mark.parametrize('order', ['ABBA', 'BAAB'])
def test_balanced_order_cancels_shared_level_and_linear_log_drift(order):
    delta = Fraction(3, 50)
    values = [100 + Fraction(7, 3) * t + delta * (v == 'B') for t, v in enumerate(order)]
    assert balanced_log_contrast(values, order) == delta


def test_shared_multiplicative_slowdown_cancels_in_ratio():
    base, candidate, slowdown = Fraction(10), Fraction(21, 2), Fraction(3, 2)
    assert candidate * slowdown / (base * slowdown) == Fraction(21, 20)
    assert candidate * slowdown / base == Fraction(63, 40)
    # A disturbance affecting only B does not cancel.
    assert candidate * slowdown / base != candidate / base


def test_random_orientation_symmetrizes_arbitrary_slot_noise():
    delta = Fraction(2, 25)
    for nuisance in product((-3, 0, 7), repeat=4):
        contrasts = []
        for order in ('ABBA', 'BAAB'):
            values = [x + delta * (v == 'B') for x, v in zip(nuisance, order)]
            contrasts.append(balanced_log_contrast(values, order))
        assert sum(contrasts) == 2 * delta
        # At the null boundary, at most one of two equiprobable orders exceeds.
        assert sum(value > delta for value in contrasts) <= 1


def test_curvature_and_unequal_spacing_do_not_cancel_with_fixed_order():
    assert balanced_log_contrast([t * t for t in range(4)]) == -2
    assert balanced_log_contrast([0, 1, 2, 10]) == Fraction(-7, 2)


def test_repeating_correlated_signs_does_not_create_independent_jobs():
    # All repeats share one fair environment sign. Five repeats still yield
    # an all-positive probability of 1/2, not the independent-job bound 1/32.
    cases = [(-1,) * 5, (1,) * 5]
    actual = Fraction(sum(all(x > 0 for x in case) for case in cases), len(cases))
    assert actual == Fraction(1, 2) > sign_tail_bound(5, 5)


@pytest.mark.parametrize('epsilon', [Fraction(0), Fraction(1, 10), Fraction(1)])
def test_bound_against_exhaustive_contamination_outcomes(epsilon):
    # Per job: clean negative, clean positive, contaminated positive.
    weights = ((1 - epsilon) / 2, (1 - epsilon) / 2, epsilon)
    jobs = 4
    mass_by_count = [Fraction(0)] * (jobs + 1)
    for outcomes in product(range(3), repeat=jobs):
        mass = Fraction(1)
        for outcome in outcomes:
            mass *= weights[outcome]
        mass_by_count[sum(outcome != 0 for outcome in outcomes)] += mass
    for k in range(jobs + 2):
        assert sign_tail_bound(jobs, k, epsilon) == sum(mass_by_count[k:])


def test_two_of_three_is_not_five_percent_evidence():
    assert sign_tail_bound(3, 2) == Fraction(1, 2)
    assert sign_tail_bound(3, 3) == Fraction(1, 8)
    assert required_exceedances(3) is None


def test_contamination_changes_minimum_job_count():
    assert required_exceedances(5) == 5
    assert sign_tail_bound(5, 5, Fraction(1, 10)) == Fraction(11, 20) ** 5
    assert sign_tail_bound(5, 5, Fraction(1, 10)) > Fraction(1, 20)
    assert required_exceedances(5, contamination=Fraction(1, 10)) is None
    assert required_exceedances(6, contamination=Fraction(1, 10)) == 6
    assert required_exceedances(20, contamination=1) is None


@pytest.mark.parametrize('jobs', [1, 5, 6, 10, 20])
@pytest.mark.parametrize('epsilon', [Fraction(0), Fraction(1, 10), Fraction(1, 2)])
def test_critical_count_is_minimal(jobs, epsilon):
    alpha = Fraction(1, 20)
    count = required_exceedances(jobs, alpha, epsilon)
    if count is None:
        assert sign_tail_bound(jobs, jobs, epsilon) > alpha
    else:
        assert sign_tail_bound(jobs, count, epsilon) <= alpha
        assert sign_tail_bound(jobs, count - 1, epsilon) > alpha


@pytest.mark.parametrize('jobs,k,epsilon', [(0, 0, 0), (3, -1, 0), (3, 5, 0), (3, 2, -1)])
def test_invalid_bound_inputs(jobs, k, epsilon):
    with pytest.raises(ValueError):
        sign_tail_bound(jobs, k, epsilon)


def test_invalid_design_and_budget():
    with pytest.raises(ValueError):
        balanced_log_contrast([0, 1, 2])
    with pytest.raises(ValueError):
        balanced_log_contrast([0, 1, 2, 3], 'AABB')
    with pytest.raises(ValueError):
        required_exceedances(5, alpha=0)
