"""Study validation using test streams, never the fresh evaluation seeds."""

import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.step_detection import directional_artifacts as artifacts
from experiments.step_detection import directional_study as s
from experiments.step_detection import indexed_study as previous
from experiments.step_detection import reporting_ar1 as ar1
from experiments.step_detection import reporting_directional as directional


def inherited():
    return json.loads((Path(__file__).parent / 'data/indexed_v1_frozen.json').read_text())


def test_design_without_evaluation_observations():
    count = math.prod(
        len(s.DESIGN[key]) for key in ('sizes', 'changes', 'noise', 'locations', 'seeds')
    )
    assert count == 120 and count * len(s.CONDITIONS) == 360
    assert set(s.DESIGN['seeds']).isdisjoint(previous.DESIGN['seeds'])
    assert set(s.DESIGN['seeds']).isdisjoint(previous.previous.DESIGN['seeds'])
    assert s.DESIGN['stream_prefix'] != previous.DESIGN['stream_prefix']
    assert sum(change > 0.05 for change in s.DESIGN['changes']) * 72 == 144
    assert s.WORK == previous.WORK and s.UNIT == previous.UNIT


@pytest.mark.parametrize('location', ['early', 'middle', 'recent'])
def test_paired_innovations_and_fresh_stream(location):
    control = s.make_case('independent_constant', 40, 0.05, 0.02, location, 12345)
    mean = 10 * (1 + 0.05 * (np.arange(40) >= control['position']))
    innovations = (np.array(control['values']) - mean) / 0.2
    for condition, rho in s.CONDITIONS.items():
        case = s.make_case(condition, 40, 0.05, 0.02, location, 12345)
        residual = (np.array(case['values']) - mean) / 0.2
        assert residual[0] == pytest.approx(innovations[0])
        assert residual[1:] - rho * residual[:-1] == pytest.approx(
            math.sqrt(1 - rho**2) * innovations[1:]
        )
        assert case['pair_id'] == control['pair_id'] and not case['positive']
    assert (
        control['values']
        != previous.make_case('independent_constant', 40, 0.05, 0.02, location, 12345)['values']
    )


def test_full_pipeline_and_verification_on_test_stream():
    parent = inherited()
    case = s.make_case('negative', 40, 0, 0.02, 'middle', 12345)
    frozen = {'inherited': parent, 'calibrations': {'40': ar1.calibration(40)}, 'work': s.WORK}
    actual = s.evaluate(case, frozen, 'native')
    counts = artifacts.verify_case((case, actual))
    assert sum(counts.values()) == 6
    expected = previous.evaluate(case, parent, 'native')
    for name in s.NEW_METHODS:
        result = actual.pop(name)
        assert result == directional.evidence(case['values'], use_tail=name == 'directional_tail')
        diagnostic = actual.pop(name + '_diagnostics')
        assert diagnostic == s.true_pair_diagnostic(
            name, case, actual['oracle_reporting'], result['calibration']
        )
        assert diagnostic['true_pair_retained']
        assert str(case['position']) not in result['confidence']['by_split']
        assert actual.pop(name + '_seconds') >= 0
        assert actual['methods'].pop(name + '_alone') == result['has_alert']
        assert actual['methods'].pop(name + '_gate') == (
            result['has_alert'] and actual['methods']['shared_existing']
        )
    for name in ('ar1', *previous.NEW_METHODS):
        actual.pop(name + '_seconds')
        expected.pop(name + '_seconds')
    assert actual == expected


@pytest.mark.parametrize(
    'field,value',
    [
        ('work', {'max_cells': 1, 'max_depth': 16}),
        ('confidence', {'tilt': '1/8', 'radius_upper_bits': 32}),
        ('design', {}),
        ('source_hashes', {}),
    ],
)
def test_freeze_detects_changes(field, value):
    frozen = s.freeze(inherited(), 'native')
    s.check_frozen(frozen, 'native')
    changed = copy.deepcopy(frozen)
    changed[field] = value
    with pytest.raises(ValueError, match='Frozen directional-study'):
        s.check_frozen(changed, 'native')


def test_paired_counts_and_incomplete_search_diagnostics():
    rows = []
    for i, (positive, before, after) in enumerate(
        [
            (True, False, True),
            (True, True, False),
            (False, False, True),
            (False, True, False),
        ],
        1,
    ):
        row = {
            'case': {'positive': positive, 'condition': 'negative', 'n': 40, 'location': 'recent'},
            'methods': {},
        }
        for method in (*s.ALL_UNKNOWN, 'oracle_reporting'):
            decision = before if method == 'directional_uniform' else after
            for suffix in ('alone', 'gate'):
                row['methods'][method + '_' + suffix] = decision
            if method == 'oracle_reporting':
                continue
            row[method] = {
                'status': 'certified_alert' if decision else 'unresolved',
                'cells_visited': i * 10,
            }
            row[method + '_seconds'] = i
            row[method + '_diagnostics'] = {
                'true_pair_retained': positive,
                'confidence_disabled_at_true_split': False,
                'generating_pair_route': 'confidence',
            }
        rows.append(row)
    counts = s.comparisons(rows)['all']['directional_tail_vs_directional_uniform_gate']
    assert counts == {'positive_gained': 1, 'positive_lost': 1, 'null_added': 1, 'null_removed': 1}
    for result in s.diagnostics(rows)['all'].values():
        assert result['statuses']['unresolved'] == 2
        assert result['statuses']['certified_alert'] == 2
        assert result['true_pair_excluded'] == 2
        assert result['false_alert_routes']['confidence'] == 1
        assert result['max_cells'] == 40
