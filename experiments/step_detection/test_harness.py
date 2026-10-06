"""Checks for experiment validity, independent of the full ASV test fixtures."""

import importlib.util
import math
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    'step_experiment', Path(__file__).with_name('harness.py')
)
h = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(h)


@pytest.fixture(params=['python', 'native'])
def backend(request):
    if request.param == 'native':
        try:
            h.load_detector('native')
        except RuntimeError:
            pytest.skip('Build the C++ extension to check the native backend')
    return request.param


def test_interval_oracle():
    assert h.direct_interval([10, 10, 100], [1, 1, 1]) == (10, 90)
    assert h.direct_interval([10, 11, 20], [1, 4, 1]) == (11, 10)
    assert h.direct_interval([1, 3], [1, 1])[1] == 2


def test_oracle_known_partition_and_constraints():
    result = h.exhaustive_fit([1, 1, 3, 3], [1] * 4, 1)
    assert result['right'] == [2, 4]
    assert result['objective'] == 1
    assert h.exhaustive_fit([1, 1, 3, 3], [1] * 4, 1, min_size=3)['objective'] == 4
    with pytest.raises(ValueError, match='No feasible'):
        h.exhaustive_fit([1, 1, 3], [1] * 3, 1, min_size=2, max_size=2)


def test_exact_against_exhaustive_oracle(backend):
    report = h.verify_oracle(backend)
    assert report['checked'] == 200
    assert report['failures'] == []


def test_matching_maximizes_count_before_distance():
    metrics = h.match_boundaries([3, 5], [1, 4], tolerance=2)
    assert metrics['matched'] == 2
    assert metrics['pairs'] == [[3, 1], [5, 4]]
    assert metrics['location_error_sum'] == 3
    assert h.match_boundaries([3], [1, 4], tolerance=2)['pairs'] == [[3, 4]]
    assert h.match_boundaries([], [1, 4])['false'] == 2
    assert h.match_boundaries([3], [])['missed'] == 1
    assert h.match_boundaries(None, [3]) is None


def test_prepare_uses_production_filtering():
    module = h.load_detector()
    case = {
        'values': [1.0, None, 3.0, 4.0, math.nan, 6.0],
        'weights': [1.0, 2.0, 0.0, None, 1.0, math.nan],
    }
    assert h.prepare(module, case) == ([1.0, 4.0, 6.0], [1.0, 1.0, 1.0], [0, 3, 5])
    assert module.solve_potts_autogamma.__name__ == 'solve_potts_autogamma'


@pytest.mark.parametrize('scenario', h.SCENARIOS)
def test_current_experiment_matches_production(backend, scenario):
    case = h.make_case(scenario, 40, 3)
    module = h.load_detector(backend)
    expected_steps = module.detect_steps(case['values'], case['weights'])
    y, w, _ = h.prepare(module, case)
    right, values, costs, gamma = module.solve_potts_autogamma(y, w)
    result = h.run_case(case, backend, 'current')
    assert result['steps'] == expected_steps
    assert result['right'] == right
    assert result['values'] == values
    assert result['costs'] == costs
    assert result['gamma'] == gamma


@pytest.mark.parametrize('budget', [4, 10])
def test_grid_uses_production_score_and_full_bracket(backend, budget):
    case = h.make_case('step', 50, 2)
    grid = h.run_case(case, backend, 'grid', budget=budget)
    assert len(grid['trace']) == budget
    assert len({row['log_gamma_ratio'] for row in grid['trace']}) == budget
    lo, hi = grid['search_bounds']
    assert min(row['log_gamma_ratio'] for row in grid['trace']) == lo
    assert max(row['log_gamma_ratio'] for row in grid['trace']) == hi
    # Force unmodified production autogamma to evaluate each recorded gamma alone.
    for sample in grid['trace']:
        module = h.load_detector(backend)
        y, w, _ = h.prepare(module, case)
        golden = module.golden_search
        scores = []
        active = False

        def force(f, a, b, _golden=golden, _sample=sample, _scores=scores, **kwargs):
            nonlocal active
            if active:
                return _golden(f, a, b, **kwargs)
            active = True
            x = _sample['log_gamma_ratio']
            try:
                _scores.append(f(x))
            finally:
                active = False
            return x

        module.golden_search = force
        right, values, costs, gamma = module.solve_potts_autogamma(y, w)
        assert right == sample['right']
        assert values == sample['values']
        assert costs == sample['costs']
        assert gamma == sample['gamma']
        assert scores == [sample['score']]


def test_hybrid_retains_current_candidate(backend):
    case = h.make_case('weak', 60, 7)
    current = h.run_case(case, backend, 'current')
    budget = len(current['trace'])
    hybrid = h.run_case(case, backend, 'hybrid', budget=budget)
    assert hybrid['selection_score'] <= current['selection_score']
    assert len(hybrid['trace']) == 2 * budget
    assert [r['gamma'] for r in hybrid['trace'][:budget]] == [r['gamma'] for r in current['trace']]
    assert hybrid['counts']['rho_cache_hits'] > 0


def test_patches_restored_on_error():
    module = h.load_detector()
    original = module.golden_search
    with pytest.raises(ValueError, match='at least 4'), h.SearchExperiment(module, 'grid', 1):
        module.solve_potts_autogamma([1, 2, 3], [1, 1, 1])
    assert module.golden_search is original


@pytest.mark.parametrize(
    'scenario',
    ['flat', 'step', 'weak', 'recent', 'recovered', 'partial', 'dip', 'multiple', 'sparse'],
)
def test_reporting_labels_match_noiseless_policy(scenario):
    case = h.make_case(scenario, 60, 0)
    case['values'] = case['true_levels']
    module = h.load_detector()
    y = case['values']
    starts = [0] + [i for i in range(1, len(y)) if y[i] != y[i - 1]]
    rights = starts[1:] + [len(y)]
    steps = [
        (case['revisions'][l], case['revisions'][r - 1] + 1, y[l], y[l], 0)
        for l, r in zip(starts, rights)
    ]
    _, _, alerts = module.detect_regressions(steps, threshold=0.05)
    observed = [case['revisions'].index(a[1]) for a in alerts or []]
    assert observed == case['true_alerts']


def test_known_approximation_gap_and_offset():
    import random

    module = h.load_detector()
    rng = random.Random(7)
    y = [round(rng.gauss(1 if i < 35 else 1.3, 0.3), 2) for i in range(70)]
    exact = module.solve_potts(y, [1] * 70, 1)
    approximate = module.solve_potts_approx(y, [1] * 70, 1)
    assert h.check_solution(y, [1] * 70, 1, exact) == pytest.approx(16.05)
    assert h.check_solution(y, [1] * 70, 1, approximate) == pytest.approx(16.52)
    outputs = []
    for offset in (0, 1000):
        y = [offset + 1] * 30 + [offset + 1.02] * 30
        outputs.append(module.solve_potts_autogamma(y, [1] * 60)[0])
    assert outputs == [[30, 60], [60]]


def test_missing_truth_maps_to_observable_boundary():
    assert h.observed_truth([5], [0, 1, 3, 8, 9]) == [8]
    assert h.observed_truth([0, 20], [3, 8]) == []
    assert h.observed_truth(None, [1, 2]) is None


def test_aggregate_excludes_drift_and_provides_interval():
    rows = [
        h.run_case(h.make_case(name, 30, 0), 'python', 'current') for name in ('flat', 'drift')
    ]
    result = h.summarize(rows)['groups'][0]
    assert result['runs'] == 2
    assert result['boundary_metrics']['histories'] == 1
    low, high = result['boundary_metrics']['false_history_rate_ci95']
    assert 0 <= low <= high <= 1
