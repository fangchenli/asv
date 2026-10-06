"""Validate the new evaluator using only test seed labels."""

import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.step_detection import ar1_study as previous
from experiments.step_detection import indexed_study as s
from experiments.step_detection import reporting_ar1 as a
from experiments.step_detection import reporting_ar1_full as full
from experiments.step_detection import reporting_ar1_jump as jump


def inherited():
    return json.loads((Path(__file__).parent / 'data/ar1_v1_frozen.json').read_text())


def test_design_without_evaluation_observations():
    bases = math.prod(
        len(s.DESIGN[key]) for key in ('sizes', 'changes', 'noise', 'locations', 'seeds')
    )
    assert bases == 120 and bases * len(s.CONDITIONS) == 360
    assert set(s.DESIGN['seeds']).isdisjoint(previous.DESIGN['seeds'])
    assert s.DESIGN['stream_prefix'] != 'ar1-study'
    assert sum(x > 0.05 for x in s.DESIGN['changes']) / len(s.DESIGN['changes']) * 360 == 144
    assert s.UNIT == 1 and s.WORK == previous.WORK


@pytest.mark.parametrize('location', ['early', 'middle', 'recent'])
def test_new_stream_pairing_and_covariance(location):
    control = s.make_case('independent_constant', 40, 0.05, 0.02, location, 12345)
    mean = 10 * (1 + 0.05 * (np.arange(40) >= control['position']))
    independent = (np.array(control['values']) - mean) / 0.2
    for condition, rho in s.CONDITIONS.items():
        case = s.make_case(condition, 40, 0.05, 0.02, location, 12345)
        residual = (np.array(case['values']) - mean) / 0.2
        assert residual[0] == pytest.approx(independent[0])
        assert residual[1:] - rho * residual[:-1] == pytest.approx(
            math.sqrt(1 - rho**2) * independent[1:]
        )
        assert case['pair_id'] == control['pair_id']
        assert case['covariance_id'] in previous.shapes()
        assert not case['positive']
    assert (
        control['values']
        != previous.make_case('independent_constant', 40, 0.05, 0.02, location, 12345)['values']
    )


def test_pipeline_reuse_and_true_pair_after_early_stop():
    parent = inherited()
    case = s.make_case('negative', 40, 0, 0.02, 'middle', 12345)
    expected = previous.evaluate(case, parent, 'native')
    frozen = {
        'inherited': parent,
        'calibrations': {'40': a.calibration(40)},
        'unit': s.UNIT,
        'work': s.WORK,
    }
    actual = s.evaluate(case, frozen, 'native')
    for name in s.NEW_METHODS:
        evidence = actual.pop(name)
        diag = actual.pop(name + '_diagnostics')
        assert actual.pop(name + '_seconds') >= 0
        assert actual['methods'].pop(name + '_alone') == evidence['has_alert']
        assert actual['methods'].pop(name + '_gate') == (
            evidence['has_alert'] and expected['methods']['shared_existing']
        )
        if name != 'full_history':
            assert str(case['position']) not in evidence['confidence']['by_split']
            log_q = jump._indexed_log_density(
                case['values'],
                case['position'],
                full.predictive_log_density(case['values']),
                original_levels=name == 'indexed_original',
            )
            assert diag['generating_split_log_density'] == pytest.approx(log_q)
        assert diag['true_pair_retained']
        assert not diag['confidence_disabled_at_true_split']
    actual.pop('ar1_seconds')
    expected.pop('ar1_seconds')
    assert actual == expected


@pytest.mark.parametrize('field,value', [('unit', 2), ('work', {'max_cells': 1, 'max_depth': 16})])
def test_freeze_rejects_changes(field, value):
    frozen = s.freeze(inherited(), 'native')
    s.check_frozen(frozen, 'native')
    changed = copy.deepcopy(frozen)
    changed[field] = value
    with pytest.raises(ValueError, match='Frozen indexed-study'):
        s.check_frozen(changed, 'native')


def test_comparison_counts():
    rows = []
    for positive, before, after in [
        (True, False, True),
        (True, True, False),
        (False, False, True),
        (False, True, False),
    ]:
        methods = {
            name + '_' + suffix: False
            for name in ('ar1', 'oracle_reporting', *s.NEW_METHODS)
            for suffix in ('alone', 'gate')
        }
        for suffix in ('alone', 'gate'):
            methods['indexed_original_' + suffix] = before
            methods['indexed_jump_' + suffix] = after
        rows.append({'case': {'positive': positive, 'condition': 'negative'}, 'methods': methods})
    counts = s.comparisons(rows)['all']['indexed_jump_vs_indexed_original_gate']
    assert counts == {'positive_gained': 1, 'positive_lost': 1, 'null_added': 1, 'null_removed': 1}


def test_diagnostics_keep_unresolved_and_false_alerts_visible():
    rows = []
    for i, (positive, status) in enumerate(
        [(True, 'certified_alert'), (True, 'unresolved'), (False, 'certified_alert')], 1
    ):
        row = {
            'case': {'condition': 'negative', 'n': 40, 'location': 'recent', 'positive': positive},
            'methods': {},
        }
        for method in ('ar1', *s.NEW_METHODS):
            row[method] = {'status': status, 'cells_visited': i * 10}
            row[method + '_seconds'] = i
            row[method + '_diagnostics'] = {
                'true_pair_retained': positive,
                'confidence_disabled_at_true_split': status == 'unresolved',
                'generating_pair_route': 'both' if positive else 'confidence',
            }
            row['methods'][method + '_gate'] = status == 'certified_alert'
        rows.append(row)
    result = s.diagnostics(rows)
    for method, summary in result['all'].items():
        assert summary['statuses']['certified_alert'] == 2
        assert summary['statuses']['unresolved'] == 1
        assert summary['true_pair_excluded'] == 1
        assert summary['false_alert_routes']['confidence'] == 1
        assert summary['median_seconds'] == 2 and summary['max_cells'] == 30
        assert summary['confidence_disabled_at_true_split'] == (None if method == 'ar1' else 1)
