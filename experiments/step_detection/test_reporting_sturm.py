"""Coverage, conservative stopping, and certificate integrity for the new route."""

import copy
import math
from fractions import Fraction as F

import pytest

from experiments.step_detection import directional_sturm as sturm
from experiments.step_detection import reporting_directional as inherited
from experiments.step_detection import reporting_sturm as reporting
from experiments.step_detection import residual_direction as direction
from experiments.step_detection import sturm_replay as replay
from experiments.step_detection import sturm_verification as verification


def history(change=0.2):
    return [10 * (1 + change * (i >= 12)) + 0.1 * math.sin(1.7 * i) for i in range(24)]


@pytest.mark.parametrize('change', [0.2, 0.04])
def test_disabling_new_bound_exactly_recovers_inherited_search(change):
    values = history(change)
    previous = inherited.evidence(values)
    result = reporting.evidence(values, use_sturm=False)
    assert verification.verify(values, result)['verified']
    comparable = copy.deepcopy(result)
    comparable['confidence']['rule'] = previous['confidence']['rule']
    del comparable['confidence']['sturm_enabled']
    del comparable['confidence']['sturm_bits']
    del comparable['confidence']['sturm_eigen_bits']
    del comparable['confidence']['sturm_tilts']
    del comparable['confidence']['sturm_counts']
    if comparable['status'] == 'surviving_explanation':
        assert comparable['witness']['confidence_bounds'].pop('sturm') is None
    assert comparable == previous


@pytest.fixture(scope='module')
def alert():
    values = history()
    return values, reporting.evidence(values)


def test_complete_search_and_unit_conversion(alert):
    values, result = alert
    assert result['has_alert']
    assert verification.verify(values, result)['complete_coverage']
    converted = reporting.evidence([1024 * value for value in values])
    assert converted['status'] == result['status']
    assert converted['certificate'] == result['certificate']
    assert verification.verify([1024 * value for value in values], converted)['verified']


@pytest.mark.parametrize('tampering', ['gap', 'overlap', 'state', 'rule', 'budget'])
def test_full_verifier_rejects_damaged_certificates(alert, tampering):
    values, original = alert
    result = copy.deepcopy(original)
    if tampering == 'gap':
        result['certificate'].pop()
    elif tampering == 'overlap':
        result['certificate'].append(result['certificate'][0])
    elif tampering == 'state':
        result['confidence']['by_split']['1']['ss'] = '0'
    elif tampering == 'rule':
        result['confidence']['rule'] = 'other'
    else:
        result['confidence']['delta'] = '1/2'
    with pytest.raises(AssertionError):
        verification.verify(values, result)


@pytest.fixture(scope='module')
def new_interval():
    item = next(item for item in replay.examples() if '-early-seed1601' in item['case']['id'])
    values = item['case']['values']
    result = reporting.evidence(values, max_cells=1)
    model = direction.state(values, 1)
    left, right = F(52223, 65536), F(52225, 65536)
    proof = sturm.certify_interval(model, left, right, delta=F(result['confidence']['delta']))
    assert proof['status'] == 'certified_excluded'
    result['certificate'] = [
        {
            'split': 1,
            'left': str(left),
            'right': str(right),
            'extra': None,
            'route': 'direction_sturm',
            'proof': proof,
        }
    ]
    return values, result


def test_new_route_verifies_partial_coverage_and_changes_old_point_decision(new_interval):
    values, result = new_interval
    checked = verification.verify(values, result)
    assert checked['certificate_routes'] == {'direction_sturm': 1}
    assert not checked['complete_coverage']
    state = direction.state(values, 1)
    delta = F(result['confidence']['delta'])
    assert not inherited.point_bounds(state, F(51, 64), delta)['excluded']
    assert reporting.point_bounds(state, F(51, 64), delta)['excluded']


@pytest.mark.parametrize(
    'tampering', ['proof', 'endpoint', 'disabled', 'precision', 'eigen_precision', 'false_alert']
)
def test_new_route_rejects_tampering(new_interval, tampering):
    values, original = new_interval
    result = copy.deepcopy(original)
    if tampering == 'proof':
        result['certificate'][0]['proof']['p_upper_squared'] = '0'
    elif tampering == 'endpoint':
        result['certificate'][0]['right'] = '1'
    elif tampering == 'disabled':
        result['confidence']['sturm_enabled'] = False
    elif tampering == 'precision':
        result['confidence']['sturm_bits'] = 10
    elif tampering == 'eigen_precision':
        result['confidence']['sturm_eigen_bits'] = 16
    else:
        result['status'], result['has_alert'] = 'certified_alert', True
    with pytest.raises(AssertionError):
        verification.verify(values, result)


def test_below_threshold_survivor_includes_the_new_bound():
    values = history(0.04)
    result = reporting.evidence(values)
    assert result['status'] == 'surviving_explanation'
    assert result['witness']['confidence_bounds']['sturm'] is not None
    assert verification.verify(values, result)['verified']
    result['witness']['confidence_bounds']['sturm']['p_upper_squared'] = '0'
    with pytest.raises(AssertionError):
        verification.verify(values, result)


@pytest.mark.parametrize(
    'budget,reason', [({'max_cells': 1}, 'max_cells'), ({'max_depth': 0}, 'max_depth')]
)
def test_search_limits_abstain_and_verify(budget, reason):
    values = history()
    result = reporting.evidence(values, **budget)
    assert result['status'] == 'unresolved' and not result['has_alert']
    assert result['witness']['reason'] == reason
    assert verification.verify(values, result)['verified']


def test_constant_history_abstains():
    values = [10] * 24
    result = reporting.evidence(values)
    assert result['status'] == 'insufficient_variation' and not result['has_alert']
    assert verification.verify(values, result)['verified']


@pytest.mark.parametrize(
    'kwargs', [{'use_sturm': 1}, {'use_tail': 0}, {'max_cells': True}, {'max_depth': -1}]
)
def test_invalid_controls(kwargs):
    with pytest.raises(ValueError):
        reporting.evidence(history(), **kwargs)


def test_modified_cutoff_rejected():
    from experiments.step_detection import reporting_ar1 as ar1

    config = ar1.calibration(24)
    config['size_t_critical'] *= 0.9
    with pytest.raises(ValueError):
        reporting.evidence(history(), config)
