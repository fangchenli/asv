"""Analytical and exhaustive checks; no generated performance experiment."""

import itertools
import math
from fractions import Fraction

import pytest

from experiments.step_detection.reporting_reference import (
    MAX_POINTS,
    known_boundary_lower,
    median_rank,
    single_change_evidence,
)


@pytest.mark.parametrize('tail', [Fraction(1, 40), Fraction(1, 8), Fraction(1, 4)])
def test_rank_coverage_by_exhaustive_sign_patterns(tail):
    # Each pattern is equally probable for independent continuous observations
    # at their common median. Check the actual order-statistic exclusion event.
    for n in range(1, 11):
        k = median_rank(n, tail)
        upper_failures = lower_failures = 0
        for signs in itertools.product((-1, 1), repeat=n):
            ordered = sorted(signs)
            if k:
                lower_failures += ordered[k - 1] > 0
                upper_failures += ordered[-k] < 0
        assert Fraction(lower_failures, 2**n) <= tail
        assert Fraction(upper_failures, 2**n) <= tail
        # A tighter rank would exceed the requested error probability.
        assert Fraction(sum(math.comb(n, j) for j in range(k + 1)), 2**n) > tail


def test_exact_tail_boundary_and_short_segment_limit():
    assert median_rank(3, Fraction(1, 8)) == 1
    assert median_rank(3, Fraction(1, 8) - Fraction(1, 10**20)) == 0
    assert median_rank(5, Fraction(1, 40)) == 0
    assert median_rank(6, Fraction(1, 40)) == 1
    for n, first_useful in [(40, 16), (100, 18)]:
        tail = Fraction(1, 20) / (n * (n + 1))
        assert median_rank(first_useful - 1, tail) == 0
        assert median_rank(first_useful, tail) == 1


def test_known_boundary_requires_evidence_above_threshold():
    assert known_boundary_lower([10] * 6 + [12] * 6, 6) == pytest.approx(1.5)
    assert known_boundary_lower([10] * 6 + [10.4] * 6, 6) < 0
    assert known_boundary_lower([10] * 6 + [10.5] * 6, 6) == 0
    # A huge point estimate cannot replace the missing lower confidence bound.
    assert known_boundary_lower([10] * 28 + [100] * 3, 28) < 0


def test_unknown_boundary_keeps_all_plausible_splits():
    result = single_change_evidence([10] * 40 + [12] * 40)
    assert result['status'] == 'evidence'
    assert result['has_alert']
    assert result['lower_excess'] == pytest.approx(1.5)
    assert not result['constant_feasible']
    assert 40 in result['feasible_splits']
    assert len(result['feasible_splits']) > 1


def test_short_unknown_boundary_example_is_conservative():
    values = [10] * 20 + [12] * 20
    assert known_boundary_lower(values, 20) > 0
    result = single_change_evidence(values)
    assert not result['has_alert']
    assert 20 in result['feasible_splits']


@pytest.mark.parametrize('change', [0, 0.04, 0.05])
def test_null_median_model_remains_feasible(change):
    result = single_change_evidence([10] * 40 + [10 * (1 + change)] * 40)
    assert not result['has_alert']
    assert result['lower_excess'] <= 0
    assert 40 in result['feasible_splits']
    assert result['constant_feasible'] == (change == 0)


def test_no_change_model_and_unbounded_case():
    result = single_change_evidence([10])
    assert result['constant_feasible']
    assert result['lower_excess'] == -math.inf
    assert not result['has_alert']
    assert single_change_evidence([10], threshold=0)['lower_excess'] == 0


def test_incompatible_model_never_becomes_vacuous_alert():
    result = single_change_evidence([10] * 40 + [12] * 40 + [10] * 40)
    assert result['status'] == 'incompatible_model'
    assert not result['has_alert']
    assert result['lower_excess'] is None
    assert not result['constant_feasible']
    assert result['feasible_splits'] == []


def test_simultaneous_constraints_against_direct_level_enumeration():
    # Independent oracle: enumerate every possible positive level cell and
    # directly check every subinterval, with no recursive intersections.
    values = [10] * 16 + [11] * 8 + [12] * 16
    alpha = Fraction(1, 4)
    n = len(values)
    tail = alpha / (n * (n + 1))
    ranks = {m: median_rank(m, tail) for m in range(1, n + 1)}
    levels = [1, 10, 10.5, 11, 11.5, 12, 13]

    def accepted(left, right, level):
        for start in range(left, right):
            for end in range(start + 1, right + 1):
                k = ranks[end - start]
                if k:
                    ordered = sorted(values[start:end])
                    if not ordered[k - 1] <= level <= ordered[-k]:
                        return False
        return True

    expected = [
        split
        for split in range(1, n)
        if any(accepted(0, split, level) for level in levels)
        and any(accepted(split, n, level) for level in levels)
    ]
    result = single_change_evidence(values, alpha=alpha)
    assert result['feasible_splits'] == expected
    assert result['constant_feasible'] == any(accepted(0, n, level) for level in levels)


@pytest.mark.parametrize('scale', [0.001, 1000])
def test_units_do_not_change_evidence(scale):
    values = [10] * 40 + [12] * 40
    before = single_change_evidence(values)
    after = single_change_evidence([scale * value for value in values])
    assert before['feasible_splits'] == after['feasible_splits']
    assert before['has_alert'] == after['has_alert']
    assert after['lower_excess'] == pytest.approx(scale * before['lower_excess'])


@pytest.mark.parametrize(
    'values,kwargs',
    [
        ([], {}),
        ([1] * (MAX_POINTS + 1), {}),
        ([math.nan], {}),
        ([math.inf], {}),
        ([1], {'alpha': 0}),
        ([1], {'alpha': 1}),
        ([1], {'threshold': -1}),
        ([1], {'threshold': math.inf}),
    ],
)
def test_invalid_inputs(values, kwargs):
    with pytest.raises(ValueError):
        single_change_evidence(values, **kwargs)


@pytest.mark.parametrize('split', [0, 2, 0.5, True])
def test_invalid_fixed_split(split):
    with pytest.raises(ValueError):
        known_boundary_lower([10, 12], split)
