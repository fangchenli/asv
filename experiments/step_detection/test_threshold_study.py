"""Check the experimental comparison, especially scoring and reporting controls."""

import copy
import math
import statistics

import pytest

from experiments.step_detection import harness as h
from experiments.step_detection import threshold_study as t
from experiments.step_detection.exact_reference import fit_independent
from experiments.step_detection.noise_model import correlation_cap, selection_score


def test_design_is_paired_and_has_disjoint_seeds():
    dev, heldout = list(t.cases('development')), list(t.cases('heldout'))
    assert len(dev) == 360 and len(heldout) == 600
    assert not {c['seed'] for c in dev} & {c['seed'] for c in heldout}
    flat = t.make_case(40, 0, 0.02, 0.7, 'middle', 0)
    step = t.make_case(40, 0.06, 0.02, 0.7, 'middle', 0)
    assert [b - a for a, b in zip(flat['values'], step['values'])] == pytest.approx(
        [0] * 20 + [0.6] * 20
    )
    assert len(t.configurations()) == 27


@pytest.mark.parametrize('rho', [0, 0.7])
@pytest.mark.parametrize('location', ['middle', 'recent'])
def test_candidate_capture_preserves_production(rho, location):
    module = h.load_detector()
    case = t.make_case(40, 0.06, 0.05, rho, location, 2)
    expected = module.detect_steps(case['values'])
    original = module.solve_potts_approx
    production, pool = t.candidate_pool(module, case['values'], 'python')
    assert module.solve_potts_approx is original
    assert production in pool
    result = t.assess(module, case, production)
    assert result['boundaries'] == [s[0] for s in expected[1:]]
    alerts = module.detect_regressions(expected, threshold=0.05)[2] or []
    assert result['alerts'] == [a[1] for a in alerts]


def test_independent_pool_contains_global_shared_floor_optimum():
    module = h.load_detector()
    case = t.make_case(40, 0.06, 0.02, 0.7, 'middle', 1)
    y, n = case['values'], case['n']
    _, pool = t.candidate_pool(module, y, 'python')
    floor = 0.5 * statistics.median(abs(a - b) for a, b in zip(y, y[1:]))
    beta = 4 * math.log(n) / n
    score = min(
        selection_score(t.innovations(module, y, fit, 0), n, len(fit['right']), floor, beta)
        for fit in pool
    )
    exact = fit_independent(y, [1] * n, floor, beta)['selected']
    assert score == pytest.approx(exact['score'], abs=1e-12)


@pytest.mark.parametrize('half_life', [0, 1, 4])
def test_innovation_score_matches_exhaustive_convex_knots(half_life):
    module = h.load_detector()
    y = [-2, 1, 0, 3, -1]
    fit = {'right': [5], 'values': [0], 'costs': [7]}
    cap = correlation_cap(half_life)
    knots = [-cap, 0, cap] + [max(-cap, min(cap, b / a)) for a, b in zip(y, y[1:]) if a]
    expected = min(abs(y[0]) + sum(abs(b - rho * a) for a, b in zip(y, y[1:])) for rho in knots)
    assert t.innovations(module, y, fit, cap) == pytest.approx(expected)


def test_exact_threshold_is_excluded_from_error_rates():
    module = h.load_detector()
    rows = []
    for change in (0.04, 0.05, 0.06):
        case = t.make_case(40, change, 0, 0, 'middle', 0)
        fit = {'right': [20, 40], 'values': [10, 10 * (1 + change)], 'costs': [0, 0]}
        rows.append({'case': case, 'methods': {'test': t.assess(module, case, fit)}})
    result = t.summarize(rows)['test']
    assert result['positive_histories'] == result['negative_histories'] == 1
    assert result['at_threshold_histories'] == 1
    assert result['false_rate'] == result['miss_rate'] == 0


def test_reporting_retains_noise_gate():
    module = h.load_detector()
    case = t.make_case(40, 0.06, 0.05, 0, 'middle', 0)
    quiet = {'right': [20, 40], 'values': [10, 10.6], 'costs': [2, 2]}
    noisy = dict(quiet, costs=[16, 16])
    assert t.assess(module, case, quiet)['has_alert']
    assert not t.assess(module, case, noisy)['has_alert']


def test_frozen_settings_reject_changed_code_or_protocol():
    frozen = {
        'source_hashes': {'code': 'abc'},
        'design': copy.deepcopy(t.DESIGN),
        'backend': 'python',
    }
    t.check_frozen(frozen, {'code': 'abc'}, 'python')
    with pytest.raises(ValueError, match='Frozen'):
        t.check_frozen(frozen, {'code': 'changed'}, 'python')
    frozen['design']['sizes'] = [12]
    with pytest.raises(ValueError, match='Frozen'):
        t.check_frozen(frozen, {'code': 'abc'}, 'python')


def test_tuning_uses_both_errors_and_declared_ties():
    configs = t.configurations()
    summary = {cfg['name']: {'false_rate': 0.1, 'miss_rate': 0.2} for cfg in configs}
    chosen = t.select_settings(summary, configs)
    assert [(c['penalty'], c['floor_factor'], c['half_life']) for c in chosen] == [
        (8, 1, 0),
        (8, 1, 1),
    ]
    alternative = next(c for c in configs if c['family'] == 'independent')
    summary[alternative['name']] = {'false_rate': 0, 'miss_rate': 0.5}
    assert t.select_settings(summary, configs) == chosen
