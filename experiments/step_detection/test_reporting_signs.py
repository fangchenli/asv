"""Exact checks for the sign statistic, Monte Carlo tolerance, and inversion."""

import itertools
import math
from fractions import Fraction

import pytest

from experiments.step_detection.reporting_reference import single_change_evidence
from experiments.step_detection.reporting_signs import (
    calibrate,
    evidence,
    interval_ranks,
    maximum_score,
    score_table,
    tolerance_rank,
)


@pytest.mark.parametrize('draws', [1, 10, 50])
@pytest.mark.parametrize('alpha', [Fraction(1, 20), Fraction(1, 4)])
def test_tolerance_rank_against_direct_binomial_sum(draws, alpha):
    failure = Fraction(1, 10)
    k, bound = tolerance_rank(draws, alpha, failure)
    p = 1 - alpha

    def tail(start):
        return sum(
            (math.comb(draws, j) * p**j * (1 - p) ** (draws - j) for j in range(start, draws + 1)),
            Fraction(0),
        )

    assert tail(k) == bound <= failure
    assert tail(k - 1) > failure


def test_exhaustive_sign_statistic_and_interval_inversion():
    n = 8
    table = score_table(n)
    sequences = list(itertools.product((0, 1), repeat=n))
    maxima = [maximum_score(s) for s in sequences]
    critical = sorted(maxima)[3 * len(maxima) // 4 - 1]
    assert sum(s > critical for s in maxima) / len(maxima) <= 0.25
    ranks = interval_ranks(n, critical)
    for signs, score in zip(sequences, maxima):
        direct = max(
            table[right - left][sum(signs[left:right])]
            for left in range(n)
            for right in range(left + 1, n + 1)
        )
        assert score == direct
        excluded = False
        for left in range(n):
            for right in range(left + 1, n + 1):
                values = sorted(2 * x - 1 for x in signs[left:right])
                k = ranks[len(values)]
                if k:
                    excluded |= not values[k - 1] <= 0 <= values[-k]
        assert excluded == (score > critical)


def test_matches_union_bound_reference_when_given_its_cutoff():
    n = 40
    table = score_table(n)
    alpha = Fraction(1, 20)
    reject = []
    for m in range(1, n + 1):
        for c in range(m // 2 + 1):
            p = min(Fraction(1), Fraction(2 * sum(math.comb(m, j) for j in range(c + 1)), 2**m))
            if p <= alpha / (n * (n + 1) // 2):
                reject.append(table[m][c])
    calibration = {'n': n, 'alpha': str(alpha), 'critical_score': min(reject) - 1}
    for values in ([10] * 40, [10] * 20 + [12] * 20, list(range(1, 41))):
        expected = single_change_evidence(values, alpha=alpha)
        expected.pop('intervals')
        assert evidence(values, calibration) == expected


def test_calibration_is_reproducible_and_rejects_wrong_length():
    result = calibrate(8, draws=63, seed=123, failure=Fraction(1, 10))
    assert result == calibrate(8, draws=63, seed=123, failure=Fraction(1, 10))
    assert result['order_rank_failure_bound'] <= 0.1
    with pytest.raises(ValueError, match='length'):
        evidence([10] * 7, result)


def test_insufficient_calibration_draws_accept_all_signs():
    result = calibrate(8, draws=1, seed=123)
    assert result['order_rank'] == 2
    assert all(k == 0 for k in interval_ranks(8, result['critical_score']).values())
    assert not evidence([10] * 4 + [100] * 4, result)['has_alert']


@pytest.mark.parametrize('n', [0, 201])
def test_bad_length(n):
    with pytest.raises(ValueError):
        score_table(n)


def test_invalid_sign_and_cutoff():
    with pytest.raises(ValueError):
        maximum_score([0, 2])
    with pytest.raises(ValueError):
        interval_ranks(8, -1)
