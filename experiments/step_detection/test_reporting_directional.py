"""Complete-search checks for directional confidence and reporting certificates."""

import copy
import math
from fractions import Fraction as F

import pytest

from experiments.step_detection import directional_diagnosis as diagnosis
from experiments.step_detection import directional_tail as tail
from experiments.step_detection import reporting_ar1 as ar1
from experiments.step_detection import reporting_directional as reporting
from experiments.step_detection import residual_direction as direction


def history(change=0.2):
    return [10 * (1 + change * (i >= 12)) + 0.1 * math.sin(1.7 * i) for i in range(24)]


@pytest.mark.parametrize('split', [1, 7, 12, 23])
def test_reporting_fit_recovers_independently_centered_confidence_state(split):
    y = history()
    assert reporting.confidence_state(ar1.Models(y), split) == direction.state(y, split)


@pytest.fixture(scope='module')
def alert():
    y = history()
    return y, reporting.evidence(y)


def test_complete_search_covers_every_location_and_correlation(alert):
    y, result = alert
    assert result['status'] == 'certified_alert'
    verified = diagnosis.verify(y, result)
    assert verified['complete_coverage']
    assert set(verified['certificate_routes']) == {'size', 'shape', 'direction_uniform'}


@pytest.mark.parametrize('tampering', ['missing_interval', 'wrong_route', 'wrong_state'])
def test_verifier_rejects_incomplete_or_invalid_proofs(alert, tampering):
    y, original = alert
    result = copy.deepcopy(original)
    if tampering == 'missing_interval':
        result['certificate'].pop()
    elif tampering == 'wrong_route':
        result['certificate'][0]['route'] = 'direction_uniform'
    else:
        result['confidence']['by_split']['1']['ss'] = '0'
    with pytest.raises(AssertionError):
        diagnosis.verify(y, result)


@pytest.fixture(scope='module')
def correlated():
    case, _ = direction.load_positive_example()
    y = case['values']
    return y, reporting.evidence(y, use_tail=False), reporting.evidence(y)


def test_tail_bound_removes_saved_witness_and_certifies_complete_history(correlated):
    y, control, combined = correlated
    assert control['status'] == 'surviving_explanation'
    assert combined['status'] == 'certified_alert'
    assert diagnosis.verify(y, control)['verified']
    verified = diagnosis.verify(y, combined)
    assert verified['complete_coverage']
    assert verified['certificate_routes']['direction_tail'] == 2
    assert 'direction_tail' not in diagnosis.verify(y, control)['certificate_routes']
    witness = control['witness']
    state = direction.state(y, witness['split'])
    delta = F(control['confidence']['delta'])
    rho = F(witness['rho'])
    assert not reporting.point_bounds(state, rho, delta, use_tail=False)['excluded']
    assert reporting.point_bounds(state, rho, delta)['excluded']


def test_tail_interval_proof_agrees_with_point_membership_and_rejects_tampering(correlated):
    y, _, result = correlated
    delta = F(result['confidence']['delta'])
    for cell in result['certificate']:
        if cell['route'] != 'direction_tail':
            continue
        model = direction.state(y, cell['split'])
        left, right = F(cell['left']), F(cell['right'])
        for rho in [left, (left + right) / 2, (left + 3 * right) / 4]:
            point = reporting.point_bounds(model, rho, delta)
            assert point['excluded']
            assert F(point['tail']['p_upper_squared']) <= F(cell['proof']['p_upper_squared'])
    damaged = copy.deepcopy(result)
    cell = next(c for c in damaged['certificate'] if c['route'] == 'direction_tail')
    cell['proof']['determinant_lower'] = '1'
    with pytest.raises(AssertionError):
        diagnosis.verify(y, damaged)


@pytest.mark.parametrize('use_tail', [False, True])
def test_below_threshold_history_keeps_a_verified_explanation(use_tail):
    y = history(0.04)
    result = reporting.evidence(y, use_tail=use_tail)
    assert result['status'] == 'surviving_explanation'
    assert not result['has_alert']
    assert diagnosis.verify(y, result)['verified']


@pytest.mark.parametrize(
    'budget,reason', [({'max_cells': 1}, 'max_cells'), ({'max_depth': 0}, 'max_depth')]
)
def test_incomplete_search_abstains(budget, reason):
    y = history()
    result = reporting.evidence(y, **budget)
    assert result['status'] == 'unresolved' and not result['has_alert']
    assert result['witness']['reason'] == reason
    assert diagnosis.verify(y, result)['verified']


def test_constant_history_abstains():
    result = reporting.evidence([10] * 24)
    assert result['status'] == 'insufficient_variation' and not result['has_alert']
    assert diagnosis.verify([10] * 24, result)['verified']


def test_unit_conversion_preserves_reporting_decision_and_interval_partition(alert):
    y, original = alert
    converted = reporting.evidence([1024 * value for value in y])
    assert converted['status'] == original['status']
    assert converted['cells_visited'] == original['cells_visited']
    assert converted['certificate'] == original['certificate']
    assert diagnosis.verify([1024 * value for value in y], converted)['verified']


def test_confidence_removes_both_levels_and_noise_scale():
    y = [F(10 + (i >= 5) + (-1) ** i) for i in range(16)]
    transformed = [3 * value + (10 if i < 5 else 25) for i, value in enumerate(y)]
    models = [reporting.confidence_state(ar1.Models(values), 5) for values in [y, transformed]]
    bounds = [reporting.point_bounds(model, F(999, 1000), F(1, 100)) for model in models]
    assert bounds[0] == bounds[1]
    assert bounds[0]['tail']['tilt'] == str(tail.TILT)


@pytest.mark.parametrize('values', [[1] * 7, [1] * 201, [float('nan')] * 8, [float('inf')] * 8])
def test_invalid_histories(values):
    with pytest.raises(ValueError, match='finite observations'):
        reporting.evidence(values)


@pytest.mark.parametrize(
    'kwargs', [{'max_cells': 0}, {'max_depth': -1}, {'max_cells': True}, {'use_tail': 1}]
)
def test_invalid_work_controls(kwargs):
    with pytest.raises(ValueError):
        reporting.evidence(history(), **kwargs)


def test_modified_calibration_is_rejected():
    config = ar1.calibration(24)
    config['size_t_critical'] *= 0.9
    with pytest.raises(ValueError, match='matching calibration'):
        reporting.evidence(history(), config)
