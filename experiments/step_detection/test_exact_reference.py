"""Verify the independent reference against every partition of tiny inputs."""

import itertools
import math
from fractions import Fraction

import pytest

from experiments.step_detection.exact_reference import MAX_POINTS, fit_independent
from experiments.step_detection.harness import load_detector


@pytest.fixture(params=['python', 'native'])
def backend(request):
    if request.param == 'native':
        try:
            load_detector('native')
        except RuntimeError:
            pytest.skip('Build the C++ extension to check the native backend')
    return request.param


def exhaustive_frontier(y, w, min_size=1, max_size=None):
    """Enumerate partitions and observed segment levels using rational costs."""
    n = len(y)
    max_size = n if max_size is None else max_size
    frontier = {}
    for mask in range(1 << (n - 1)):
        rights = [i for i in range(1, n) if mask & (1 << (i - 1))] + [n]
        left, total = 0, Fraction(0)
        for right in rights:
            if not min_size <= right - left <= max_size:
                break
            total += min(
                sum(
                    Fraction(w[i]) * abs(Fraction(y[i]) - Fraction(level))
                    for i in range(left, right)
                )
                for level in set(y[left:right])
            )
            left = right
        else:
            count = len(rights)
            frontier[count] = min(frontier.get(count, total), total)
    return frontier


@pytest.mark.parametrize('min_size,max_size', [(1, None), (1, 2), (2, None), (2, 3)])
@pytest.mark.parametrize('floor,beta', [(0.1, 0.2), (1, 0.5)])
def test_against_every_partition(backend, min_size, max_size, floor, beta):
    # 81 sequences, two weight patterns, four constraints, two scoring models:
    # 1296 input/model cases per backend, with every feasible count checked.
    for y in itertools.product((-1, 0, 2), repeat=4):
        for w in ([1, 1, 1, 1], [1, 3, 2, 4]):
            expected = exhaustive_frontier(y, w, min_size, max_size)
            result = fit_independent(
                y, w, floor, beta, min_size=min_size, max_size=max_size, backend=backend
            )
            assert {fit['segments'] for fit in result['frontier']} == expected.keys()
            scores = []
            for fit in result['frontier']:
                count = fit['segments']
                assert fit['error_sum'] == float(expected[count])
                assert len(fit['right']) == len(fit['values']) == len(fit['costs']) == count
                left = 0
                actual_error = Fraction(0)
                for right, level, cost in zip(fit['right'], fit['values'], fit['costs']):
                    assert min_size <= right - left <= (max_size or len(y))
                    actual = sum(
                        Fraction(w[i]) * abs(Fraction(y[i]) - Fraction(level))
                        for i in range(left, right)
                    )
                    assert float(actual) == cost
                    actual_error += actual
                    left = right
                assert left == len(y)
                assert actual_error == expected[count]
                ratio = float(expected[count]) / (len(y) * floor)
                score = beta * count + (ratio if ratio <= 1 else 1 + math.log(ratio))
                assert fit['score'] == pytest.approx(score)
                scores.append(score)
            assert result['selected']['score'] == pytest.approx(min(scores))
            best = result['selected']
            gamma = beta * max(best['error_sum'], len(y) * floor)
            one_level_error = float(exhaustive_frontier(y, w, len(y), len(y))[1])
            assert beta * len(y) * floor <= gamma + 1e-12
            assert gamma <= beta * max(one_level_error, len(y) * floor) + 1e-12
            linear = best['error_sum'] + gamma * (best['segments'] - 1)
            for fit in result['frontier']:
                assert linear <= fit['error_sum'] + gamma * (fit['segments'] - 1) + 1e-12


def test_known_step_and_frontier(backend):
    result = fit_independent([10, 10, 10, 12, 12, 12], [1] * 6, 0.1, 0.2, backend=backend)
    assert result['selected']['right'] == [3, 6]
    assert result['selected']['values'] == [10, 12]
    assert result['selected']['error_sum'] == 0
    assert result['selected']['score'] == 0.4
    assert [fit['segments'] for fit in result['frontier']] == list(range(1, 7))


def test_ties_prefer_fewer_segments_and_earlier_starts(backend):
    result = fit_independent([0, 0, 1, 1], [1] * 4, 1, 0.5, backend=backend)
    assert result['frontier'][0]['score'] == result['frontier'][1]['score'] == 1
    assert result['selected']['segments'] == 1
    constant = fit_independent([3] * 4, [1] * 4, 1, 0, backend=backend)
    assert constant['selected']['segments'] == 1
    assert constant['frontier'][1]['right'] == [1, 4]
    assert constant['frontier'][2]['right'] == [1, 2, 4]


@pytest.mark.parametrize(
    'scale,offset,weight_scale', [(1, 1000, 1), (1e-100, 0, 1), (1e100, 0, 1), (1, 0, 10)]
)
def test_transformation_invariance(backend, scale, offset, weight_scale):
    y = [10, 10.2, 10.1, 12, 12.2, 12.1]
    w = [1, 2, 1, 1, 2, 1]
    before = fit_independent(y, w, 0.1, 0.2, backend=backend)
    after = fit_independent(
        [scale * value + offset for value in y],
        [weight_scale * weight for weight in w],
        0.1 * scale * weight_scale,
        0.2,
        backend=backend,
    )
    assert after['selected']['right'] == before['selected']['right']
    assert after['selected']['score'] == pytest.approx(before['selected']['score'], abs=1e-10)
    for first, second in zip(before['frontier'], after['frontier']):
        assert second['error_sum'] == pytest.approx(first['error_sum'] * scale * weight_scale)


def test_single_observation(backend):
    result = fit_independent([4], [2], 0.1, 0.2, backend=backend)
    assert result['selected'] == {
        'segments': 1,
        'right': [1],
        'values': [4],
        'costs': [0],
        'error_sum': 0,
        'score': 0.2,
    }


def test_no_feasible_partition():
    with pytest.raises(ValueError, match='No feasible partition'):
        fit_independent([1] * 5, [1] * 5, 1, 0.2, min_size=3, max_size=3)


def test_reference_size_limit():
    with pytest.raises(ValueError, match='limited to'):
        fit_independent([1] * (MAX_POINTS + 1), [1] * (MAX_POINTS + 1), 1, 0.2)


@pytest.mark.parametrize(
    'y,w',
    [
        ([], []),
        ([1], [1, 2]),
        ([math.nan], [1]),
        ([math.inf], [1]),
        ([1], [0]),
        ([1], [-1]),
        ([1], [math.inf]),
    ],
)
def test_invalid_observations_and_weights(y, w):
    with pytest.raises(ValueError):
        fit_independent(y, w, 0.1, 0.2)


@pytest.mark.parametrize(
    'kwargs',
    [{'min_size': 0}, {'max_size': 1.5}, {'min_size': True}, {'min_size': 3, 'max_size': 2}],
)
def test_invalid_constraints(kwargs):
    with pytest.raises(ValueError):
        fit_independent([1] * 4, [1] * 4, 0.1, 0.2, **kwargs)


@pytest.mark.parametrize('floor,beta', [(0, 1), (-1, 1), (math.inf, 1), (1, -1), (1, math.nan)])
def test_invalid_model(floor, beta):
    with pytest.raises(ValueError):
        fit_independent([1, 2], [1, 1], floor, beta)
