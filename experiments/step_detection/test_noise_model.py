"""Analytical checks for the proposed score, without a benchmark campaign."""

import importlib.util
import itertools
import math
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    'noise_model', Path(__file__).with_name('noise_model.py')
)
model = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(model)


@pytest.mark.parametrize(
    'error_sum,n,floor,expected',
    [(0, 10, 0.5, 0), (2, 10, 0.5, 0.4), (5, 10, 0.5, 1), (10, 10, 0.5, 1 + math.log(2))],
)
def test_profile_loss_analytical(error_sum, n, floor, expected):
    assert model.profile_loss(error_sum, n, floor) == pytest.approx(expected)


@pytest.mark.parametrize('error_sum', [0, 2, 5, 20])
def test_profile_loss_matches_constrained_likelihood(error_sum):
    n, floor = 10, 0.5
    b_hat = max(error_sum / n, floor)

    def likelihood(b):
        return math.log(b) + error_sum / (n * b) - math.log(floor)

    result = model.profile_loss(error_sum, n, floor)
    assert result == pytest.approx(likelihood(b_hat))
    for b in [floor, 0.75, 1, 2, 4, 16]:
        assert result <= likelihood(b) + 1e-14


@pytest.mark.parametrize('error_sum', [0, 1, 10])
@pytest.mark.parametrize('scale', [1e-200, 1, 1e200])
def test_units_and_weight_scaling(error_sum, scale):
    expected = model.selection_score(error_sum, 10, 2, 0.5, 0.2)
    actual = model.selection_score(error_sum * scale, 10, 2, 0.5 * scale, 0.2)
    assert actual == pytest.approx(expected, abs=1e-12)


def test_extreme_ratios_stay_finite():
    assert model.profile_loss(1e308, 10, 1e-300) == pytest.approx(1 + 607 * math.log(10))
    assert model.profile_loss(1e308, 4, 1e308) == pytest.approx(0.25)
    assert model.profile_loss(0, 100, 1e-300) == 0


@pytest.mark.parametrize('n', [0, -1, 1.5, True])
def test_invalid_observation_count(n):
    with pytest.raises(ValueError, match='positive integer'):
        model.profile_loss(1, n, 0.5)


@pytest.mark.parametrize('error_sum', [-1, math.nan, math.inf])
def test_invalid_error(error_sum):
    with pytest.raises(ValueError, match='nonnegative'):
        model.profile_loss(error_sum, 10, 0.5)


@pytest.mark.parametrize('floor', [0, -1, math.nan, math.inf])
def test_invalid_floor(floor):
    with pytest.raises(ValueError, match='positive'):
        model.profile_loss(1, 10, floor)


@pytest.mark.parametrize('segments', [0, 11, 1.5, True])
def test_invalid_segment_count(segments):
    with pytest.raises(ValueError, match='segments'):
        model.selection_score(1, 10, segments, 0.5, 0.2)


@pytest.mark.parametrize('beta', [-1, math.nan, math.inf])
def test_invalid_penalty(beta):
    with pytest.raises(ValueError, match='beta'):
        model.selection_score(1, 10, 1, 0.5, beta)


def test_exact_fits_prefer_fewer_segments():
    assert model.selection_score(0, 10, 1, 0.5, 0.2) < model.selection_score(0, 10, 3, 0.5, 0.2)


@pytest.mark.parametrize('floor,beta', list(itertools.product([0.1, 1, 10], [0.01, 0.2, 1])))
def test_supporting_penalty_for_profile_score(floor, beta):
    # An arbitrary finite set of achievable (error, count) pairs suffices to
    # check the tangent argument; neither gamma selection nor data are fitted.
    candidates = [(20, 1), (9, 2), (3, 3), (0, 4)]
    n = 10
    best_error, best_count = min(
        candidates, key=lambda p: model.selection_score(p[0], n, p[1], floor, beta)
    )
    gamma = beta * max(best_error, n * floor)
    best_linear = best_error + gamma * (best_count - 1)
    for error, count in candidates:
        assert best_linear <= error + gamma * (count - 1) + 1e-12


@pytest.mark.parametrize('n,split_wins', [(32, True), (4096, False)])
def test_correlation_can_hide_a_persistent_step(n, split_wins):
    # Two equal noiseless blocks at 0 and 2: a one-level median fit at 1
    # has optimal AR(1) error 3 for n >= 4; the two-level fit has error 0.
    beta = 4 * math.log(n) / n
    merged = model.selection_score(3, n, 1, 0.1, beta)
    split = model.selection_score(0, n, 2, 0.1, beta)
    assert (split < merged) == split_wins
    # Independent error for the same one-level fit is n, so it still splits.
    independent = model.selection_score(n, n, 1, 0.1, beta)
    assert split < independent


@pytest.mark.parametrize(
    'half_life,expected', [(0, 0), (1, 0.5), (2, math.sqrt(0.5)), (4, 2**-0.25)]
)
def test_half_life_cap(half_life, expected):
    cap = model.correlation_cap(half_life)
    assert cap == pytest.approx(expected)
    if half_life:
        assert cap**half_life == pytest.approx(0.5)


@pytest.mark.parametrize('half_life', [-1, math.inf, math.nan, 1e300])
def test_invalid_or_unresolvable_half_life(half_life):
    with pytest.raises(ValueError, match='max_half_life'):
        model.correlation_cap(half_life)


@pytest.mark.parametrize('cap', [0, 0.5, 0.9, 1])
def test_bounded_step_score(cap):
    residuals = [-1] * 4 + [1] * 4
    result = model.fit_correlated_score(residuals, [1] * 8, 1, 0.1, 0.2, cap)
    assert result['rho'] == cap
    assert result['error_sum'] == pytest.approx(8 - 5 * cap)
    assert result['score'] == pytest.approx(model.selection_score(8 - 5 * cap, 8, 1, 0.1, 0.2))


@pytest.mark.parametrize('cap', [0, 0.5, 0.9])
@pytest.mark.parametrize('weights', [[1, 1, 1], [1, 4, 2]])
def test_contraction_bound(cap, weights):
    for residuals in itertools.product((-2, -1, 0, 1, 2), repeat=3):
        result = model.fit_correlated_score(residuals, weights, 1, 0.1, 0.2, cap)
        independent = sum(w * abs(e) for w, e in zip(weights, residuals))
        lower = (1 - cap) * min(weights) / max(weights) * independent
        assert lower - 1e-12 <= result['error_sum'] <= independent + 1e-12
        assert abs(result['rho']) <= cap


def test_zero_residuals_choose_zero_correlation():
    result = model.fit_correlated_score([0, 0, 0], [1, 2, 3], 1, 0.1, 0.2, 0.5)
    assert result == {'rho': 0, 'rho_max': 0.5, 'error_sum': 0, 'score': 0.2}


def test_half_life_bound_restores_lasting_step_preference():
    n = 4096
    residuals = [-1] * (n // 2) + [1] * (n // 2)
    beta = 4 * math.log(n) / n
    split = model.selection_score(0, n, 2, 0.1, beta)
    unbounded = model.fit_correlated_score(residuals, [1] * n, 1, 0.1, beta, 1)
    cap = model.correlation_cap(4)
    bounded = model.fit_correlated_score(residuals, [1] * n, 1, 0.1, beta, cap)
    assert unbounded['score'] < split < bounded['score']


@pytest.mark.parametrize(
    'residuals,weights',
    [([], []), ([1], [1, 2]), ([math.nan], [1]), ([math.inf], [1]), ([1], [0]), ([1], [math.inf])],
)
def test_invalid_residual_score_inputs(residuals, weights):
    with pytest.raises(ValueError):
        model.fit_correlated_score(residuals, weights, 1, 0.1, 0.2, 0.5)
