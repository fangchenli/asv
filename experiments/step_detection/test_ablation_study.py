"""Verify scoring controls and calibration without using evaluation outcomes."""

import copy
import math

import pytest

from experiments.step_detection import ablation_study as a
from experiments.step_detection import harness as h
from experiments.step_detection import threshold_study as t
from experiments.step_detection.inspect_comparison import score_candidate


@pytest.mark.parametrize('rho', [0, 0.7])
@pytest.mark.parametrize('noise', [0.005, 0.05])
def test_current_score_replays_production_and_penalty_change(rho, noise):
    module = h.load_detector()
    case = t.make_case(40, 0.06, noise, rho, 'recent', 200)
    y = case['values']
    _, pool = t.candidate_pool(module, y, 'python')
    cfg = next(c for c in a.configurations() if c['name'] == 'current-r1-c4')
    for fit in pool:
        error = t.innovations(module, y, fit, 1)
        actual = a.candidate_score(fit, error, len(y), 1, cfg)
        assert actual == score_candidate(module, y, [1] * len(y), fit)
        reduced = a.candidate_score(fit, error, len(y), 1, dict(cfg, penalty=2))
        assert actual - reduced == pytest.approx(2 * math.log(len(y)) / len(y) * len(fit['right']))


@pytest.mark.parametrize('rho', [0, 0.7])
def test_previous_bounded_setting_is_preserved(rho):
    case = t.make_case(40, 0.06, 0.02, rho, 'middle', 200)
    cfg = next(c for c in a.configurations() if c['name'] == 'shared-r0.5-c2')
    old = next(c for c in t.configurations() if c['name'] == 'bounded-c2-f0.5-h1')
    assert (
        a.evaluate(case, [cfg], 'python')['methods'][cfg['name']]
        == t.evaluate(case, [old], 'python')['methods'][old['name']]
    )


def test_calibration_obeys_budget_and_can_be_infeasible():
    configs = [c for c in a.configurations() if c['family'] == 'shared-r0.5']
    summary = {c['name']: {'false_rate': 0.02, 'miss_rate': 0.5} for c in configs}
    summary['shared-r0.5-c2'] = {'false_rate': 0.04, 'miss_rate': 0.1}
    selected = a.calibrate(summary, configs)
    assert selected[0]['config'] is None
    assert selected[1]['config'] is None
    assert selected[2]['config']['penalty'] == 2
    summary['shared-r0.5-c8'] = {'false_rate': 0, 'miss_rate': 0.2}
    selected = a.calibrate(summary, configs)
    assert [r['config']['penalty'] for r in selected] == [8, 8, 2]


def test_frozen_design_checks_and_fixed_factorial_membership():
    frozen = {
        'source_hashes': {'code': 'abc'},
        'design': copy.deepcopy(a.DESIGN),
        'backend': 'python',
        'penalties': a.PENALTIES,
        'budgets': a.BUDGETS,
    }
    a.check_frozen(frozen, {'code': 'abc'}, 'python')
    frozen['budgets'] = [0.5]
    with pytest.raises(ValueError, match='Frozen'):
        a.check_frozen(frozen, {'code': 'abc'}, 'python')
    assert len(a.evaluation_configs([])) == 8
    assert not set(a.DESIGN['development_seeds']) & set(a.DESIGN['heldout_seeds'])
    assert not (set(a.DESIGN['development_seeds']) | set(a.DESIGN['heldout_seeds'])) & set(
        t.DESIGN['development_seeds'] + t.DESIGN['heldout_seeds']
    )


def test_native_and_python_agree():
    try:
        h.load_detector('native')
    except RuntimeError:
        pytest.skip('Build the C++ extension to compare backends')
    case = t.make_case(40, 0.04, 0.02, 0.7, 'recent', 200)
    configs = a.evaluation_configs([])
    assert (
        a.evaluate(case, configs, 'native')['methods']
        == a.evaluate(case, configs, 'python')['methods']
    )
