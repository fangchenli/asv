"""Check the opt-in route, finite work budget, and reconstructed certificates."""

import copy
import json
from fractions import Fraction as F

import pytest

from experiments.step_detection import directional_normalized as normalized
from experiments.step_detection import interval_polynomial_diagnosis as diagnosis
from experiments.step_detection import reporting_sturm as reporting
from experiments.step_detection import residual_direction as direction
from experiments.step_detection import sturm_verification as verification
from experiments.step_detection.test_reporting_sturm import history


@pytest.fixture(scope='module')
def saved_interval():
    saved = json.loads(
        (
            diagnosis.HERE / 'data/sturm_fresh_v2_interval_poly_n40_early_normalized.json'
        ).read_text()
    )
    case, _ = diagnosis.saved_case(saved['case_id'])
    return case['values'], saved


@pytest.fixture(scope='module')
def certificate(saved_interval):
    values, saved = saved_interval
    model = direction.state(values, saved['split'])
    proof = normalized.certify_interval(model, F(saved['left']), F(saved['right']))
    assert proof['status'] == 'certified_excluded'
    assert proof['p_upper_squared'] == saved['p_upper_squared']
    for attempt, expected in zip(proof['attempts'], saved['cells'][0]['attempts'], strict=True):
        assert attempt['base']['determinant_enclosure'][0] == expected['determinant_lower']
        assert attempt['density_correction'] == expected['density_correction']
    return proof


def test_reusable_kernel_reproduces_saved_diagnostic(certificate):
    assert certificate['stationary_normalized'] is True
    assert certificate['bits'] == 128 and certificate['radius_bits'] == 192


@pytest.mark.parametrize('failure', ['radius', 'determinant'])
def test_inconclusive_enclosure_abstains(monkeypatch, failure):
    model = direction.state(history(), 1)
    if failure == 'radius':
        model['num'] = [F(0)]
    else:
        monkeypatch.setattr(normalized.polynomial, 'determinant_interval', lambda *a, **k: None)
    proof = normalized.certify_interval(model, F(1, 2), F(3, 4))
    assert proof['status'] == 'not_certified' and proof['p_upper_squared'] == '1'
    assert all(a['base']['status'] == f'{failure}_inconclusive' for a in proof['attempts'])


@pytest.mark.parametrize('left,right', [(-1, 0), (0, 1), (F(1, 2), F(1, 2))])
def test_unsupported_interval_rejected(left, right):
    with pytest.raises(ValueError):
        normalized.certify_interval({}, left, right)


def test_default_never_uses_fallback(monkeypatch, saved_interval):
    def unexpected(*args, **kwargs):
        pytest.fail('Default search invoked the normalized route')

    monkeypatch.setattr(normalized, 'certify_interval', unexpected)
    values, _ = saved_interval
    result = reporting.evidence(values)
    assert result['status'] == 'unresolved' and result['witness']['reason'] == 'max_depth'
    assert 'normalized' not in result['confidence']
    assert reporting.evidence(values, use_normalized=False, max_normalized_cells=0) == result


@pytest.mark.parametrize('change', [0.2, 0.04])
def test_opt_in_preserves_completed_search_when_fallback_is_unneeded(change):
    values = history(change)
    original = reporting.evidence(values)
    result = reporting.evidence(values, use_normalized=True)
    assert verification.verify(values, result)['verified']
    assert result['confidence'].pop('normalized')['calls'] == 0
    assert result == original


@pytest.mark.parametrize('budget', [0, 1, 2])
def test_fallback_runs_only_at_depth_limit_and_obeys_budget(monkeypatch, saved_interval, budget):
    values, saved = saved_interval
    seen = []

    def fake_certificate(model, left, right, *, delta):
        # This test isolates control flow; the real proof is tested separately.
        assert 0 <= left < right < 1
        assert right - left == F(1, 65536)
        assert delta == F(reporting.ar1.calibration(len(values))['confidence_alpha'])
        seen.append((left, right))
        return {'status': 'certified_excluded'}

    monkeypatch.setattr(normalized, 'certify_interval', fake_certificate)
    result = reporting.evidence(values, use_normalized=True, max_normalized_cells=budget)
    assert len(seen) == budget == result['confidence']['normalized']['calls']
    if seen:
        assert seen[0] == (F(saved['left']), F(saved['right']))
    assert result['status'] == 'unresolved'
    assert result['witness']['reason'] == 'max_normalized_cells'
    assert sum(c['route'] == 'direction_normalized' for c in result['certificate']) == budget


def test_failed_fallback_stops_conservatively(monkeypatch, saved_interval):
    values, _ = saved_interval
    monkeypatch.setattr(
        normalized, 'certify_interval', lambda *a, **k: {'status': 'not_certified'}
    )
    result = reporting.evidence(values, use_normalized=True)
    assert result['status'] == 'unresolved' and result['witness']['reason'] == 'max_depth'
    assert result['confidence']['normalized']['calls'] == 1
    assert verification.verify(values, result)['verified']


def test_cell_limit_does_not_invoke_fallback(monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail('Cell budget bypassed')

    monkeypatch.setattr(normalized, 'certify_interval', unexpected)
    result = reporting.evidence(history(), use_normalized=True, max_cells=1)
    assert result['witness']['reason'] == 'max_cells'
    assert result['confidence']['normalized']['calls'] == 0
    assert verification.verify(history(), result)['verified']


@pytest.fixture(scope='module')
def partial_result(saved_interval, certificate):
    values, saved = saved_interval
    result = reporting.evidence(values, use_normalized=True, max_cells=1)
    model = direction.state(values, saved['split'])
    result['confidence']['by_split'][str(saved['split'])] = {
        'ss': str(model['ss']),
        **{key: [str(x) for x in model[key]] for key in ('num', 'den', 'b')},
    }
    # Use the exact binary-float allocation employed by the search.
    proof = normalized.certify_interval(
        model,
        F(saved['left']),
        F(saved['right']),
        delta=F(result['confidence']['delta']),
    )
    result['confidence']['normalized']['calls'] = 1
    result['certificate'] = [
        {
            'split': saved['split'],
            'left': saved['left'],
            'right': saved['right'],
            'extra': None,
            'route': 'direction_normalized',
            'proof': proof,
        }
    ]
    return values, result


def test_real_normalized_route_verifies_partial_coverage(partial_result):
    values, result = partial_result
    checked = verification.verify(values, result)
    assert checked['verified'] and not checked['complete_coverage']
    assert checked['certificate_routes'] == {'direction_normalized': 1}


@pytest.mark.parametrize(
    'tampering',
    ['proof', 'units', 'endpoint', 'disabled', 'precision', 'calls', 'budget', 'alert'],
)
def test_verifier_rejects_normalized_tampering(partial_result, tampering):
    values, original = partial_result
    result = copy.deepcopy(original)
    cell, config = result['certificate'][0], result['confidence']['normalized']
    if tampering == 'proof':
        cell['proof']['p_upper_squared'] = '0'
    elif tampering == 'units':
        cell['proof']['attempts'][0]['base']['stationary_normalized'] = False
    elif tampering == 'endpoint':
        cell['right'] = '1'
    elif tampering == 'disabled':
        del result['confidence']['normalized']
    elif tampering == 'precision':
        config['bits'] = 192
    elif tampering == 'calls':
        config['calls'] = 0
    elif tampering == 'budget':
        config['max_cells'] = 0
    else:
        result['status'], result['has_alert'] = 'certified_alert', True
    with pytest.raises(AssertionError):
        verification.verify(values, result)


@pytest.mark.parametrize(
    'kwargs',
    [
        {'use_normalized': 1},
        {'max_normalized_cells': True},
        {'max_normalized_cells': -1},
    ],
)
def test_invalid_controls(kwargs):
    with pytest.raises(ValueError):
        reporting.evidence(history(), **kwargs)
