"""Check the frozen evaluator without generating evaluation-seed histories."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.step_detection import ar1_study as s
from experiments.step_detection import covariance_study as c
from experiments.step_detection import reporting_ar1 as a
from experiments.step_detection import reporting_covariance as gls


def inherited():
    return json.loads((Path(__file__).parent / 'data/covariance_v1_frozen.json').read_text())


def test_design_and_matrix_registry_without_evaluation_data():
    base = math.prod(len(s.DESIGN[k]) for k in ('sizes', 'changes', 'noise', 'locations', 'seeds'))
    assert base == 120 and base * len(s.CONDITIONS) == 360
    assert set(s.DESIGN['seeds']).isdisjoint(c.DESIGN['seeds'])
    assert len(s.shapes()) == 6
    for matrix in s.shapes().values():
        assert np.linalg.eigvalsh(matrix)[0] > 0


@pytest.mark.parametrize('location', ['early', 'middle', 'recent'])
def test_paired_stationary_innovations(location):
    cases = {
        condition: s.make_case(condition, 40, 0.05, 0.02, location, 12345)
        for condition in s.CONDITIONS
    }
    control = cases['independent_constant']
    center = np.array([10 * (1 + 0.05 * (i >= control['position'])) for i in range(40)])
    independent = (np.array(control['values']) - center) / 0.2
    for condition, rho in s.CONDITIONS.items():
        case = cases[condition]
        residual = (np.array(case['values']) - center) / 0.2
        assert residual[0] == pytest.approx(independent[0])
        assert residual[1:] - rho * residual[:-1] == pytest.approx(
            math.sqrt(1 - rho * rho) * independent[1:]
        )
        assert case['pair_id'] == control['pair_id']
        assert not case['positive']
    assert control['position'] == {'early': 10, 'middle': 20, 'recent': 36}[location]


@pytest.mark.parametrize(
    'coefficients,bound,expected',
    [
        ([0, 0, 1], 0.25, [-0.5, 0.5]),
        ([2], 1, None),
        ([1], 1, [-1, 1]),
        ([0, 1], 0, [-1, 0]),
        ([0, -1], 0, [0, 1]),
        ([0, 0, 1], 0, [0, 0]),
        ([1, -2, 1], 0, None),
        ([1, 2, 1], 0, None),
        ([0, 0, 1], math.inf, [-1, 1]),
    ],
)
def test_quadratic_region(coefficients, bound, expected):
    actual = s.quadratic_interval(coefficients, bound)
    if expected is None:
        assert actual is None
    else:
        assert actual == pytest.approx(expected)


@pytest.mark.parametrize(
    'coefficients,bound,expected',
    [
        ([1, -2, 1], 0, False),
        ([0, 0, 1], 1, True),
        ([1], 1, True),
        ([0, 0, 1], 0.99, False),
        ([0, 0, 1], math.inf, True),
    ],
)
def test_near_one_membership_is_exact(coefficients, bound, expected):
    assert s.permits_near_one(coefficients, bound) is expected


def test_rethreshold_matches_recomputed_oracle():
    case = s.make_case('positive', 40, 0.06, 0.02, 'early', 12345)
    matrix = c.covariance(case)
    config = a.calibration(40)
    old = gls.evidence(case['values'], matrix, c.s.d.critical_values(40))
    assert s.rethreshold(old, config) == gls.evidence(case['values'], matrix, config)


def test_complete_pipeline_reused_and_true_pair_checked():
    parent = inherited()
    case = s.make_case('negative', 40, 0.04, 0.02, 'early', 12345)
    expected = c.evaluate(case, parent, 'native')
    actual = s.evaluate(
        case,
        {'inherited': parent, 'calibrations': {'40': a.calibration(40)}, 'work': s.WORK},
        'native',
    )
    unknown = actual.pop('ar1')
    reporting = actual.pop('oracle_reporting')
    diag = actual.pop('ar1_diagnostics')
    assert actual.pop('ar1_seconds') >= 0
    for prefix, evidence in [('ar1', unknown), ('oracle_reporting', reporting)]:
        assert actual['methods'].pop(prefix + '_alone') == evidence['has_alert']
        assert actual['methods'].pop(prefix + '_gate') == (
            evidence['has_alert'] and expected['methods']['shared_existing']
        )
    assert actual == expected
    assert 0 <= diag['union_width'] <= 2
    assert 0 <= diag['generating_split_width'] <= diag['union_width']
    if unknown['has_alert']:
        assert diag['generating_pair_route'] not in ('neither', 'abstain')


def test_frozen_budget_change_rejected():
    frozen = s.freeze(inherited(), 'native')
    s.check_frozen(frozen, 'native')
    frozen['work'] = dict(frozen['work'], max_cells=100)
    with pytest.raises(ValueError, match='Frozen'):
        s.check_frozen(frozen, 'native')
