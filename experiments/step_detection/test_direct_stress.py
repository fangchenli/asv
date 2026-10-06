"""Check covariance mathematics, paired noise, and the unchanged reporting test."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.step_detection import direct_stress as s
from experiments.step_detection import direct_study as d


@pytest.mark.parametrize('rho', [0, 0.7, -0.3])
@pytest.mark.parametrize('profile', ['constant', 'up', 'down'])
def test_covariance_diagnostic_against_innovation_matrix(rho, profile):
    n, split = 9, 6
    # Independently construct the linear mapping from independent N(0,1)
    # innovations (including the stationary initial draw) to observations.
    transform = np.zeros((n, n))
    transform[:, 0] = [rho**i for i in range(n)]
    for i in range(n):
        for j in range(1, i + 1):
            transform[i, j] = math.sqrt(1 - rho * rho) * rho ** (i - j)
    transform = np.diag(s.amplitudes(n, split, profile)) @ transform
    covariance = transform @ transform.T
    design = np.column_stack((np.arange(n) < split, np.arange(n) >= split)).astype(float)
    projection = design @ np.linalg.inv(design.T @ design) @ design.T
    contrast = np.array([-1.05 / split] * split + [1 / (n - split)] * (n - split))
    actual = contrast @ covariance @ contrast
    residual = np.trace((np.eye(n) - projection) @ covariance)
    result = s.variance_diagnostic(n, split, rho, profile)
    assert result['contrast_variance'] == pytest.approx(actual)
    assert result['expected_residual_sum'] == pytest.approx(residual)
    assert result['expected_reported_variance'] == pytest.approx(
        residual / (n - 2) * (contrast @ contrast)
    )


def test_independent_variance_change_closed_form():
    n, before, after = 40, 36, 4
    for profile, va, vb in [('constant', 1, 1), ('up', 1, 9), ('down', 9, 1)]:
        actual = 1.05**2 * va / before + vb / after
        expected = (
            ((before - 1) * va + (after - 1) * vb) / (n - 2) * (1.05**2 / before + 1 / after)
        )
        assert s.variance_diagnostic(n, before, 0, profile)['variance_ratio'] == pytest.approx(
            actual / expected
        )
    assert s.variance_diagnostic(40, 36, 0, 'constant')['variance_ratio'] == pytest.approx(1)
    assert s.variance_diagnostic(40, 36, 0, 'up')['variance_ratio'] == pytest.approx(
        4.9810397, rel=1e-5
    )


@pytest.mark.parametrize('location', ['middle', 'recent'])
def test_paired_amplitudes_and_ar_innovations(location):
    cases = {c: s.make_case(c, 40, 0.04, 0.02, location, 12345) for c in s.CONDITIONS}
    control = cases['independent_constant']
    pos = control['position']
    center = np.array([10 * (1 + 0.04 * (i >= pos)) for i in range(40)])
    independent = (np.array(control['values']) - center) / 0.2
    correlated = (np.array(cases['correlated_constant']['values']) - center) / 0.2
    assert correlated[0] == pytest.approx(independent[0])
    assert correlated[1:] - 0.7 * correlated[:-1] == pytest.approx(
        math.sqrt(1 - 0.7**2) * independent[1:]
    )
    for name, (rho, profile) in s.CONDITIONS.items():
        base = correlated if rho else independent
        error = (np.array(cases[name]['values']) - center) / 0.2
        assert error == pytest.approx(base * s.amplitudes(40, pos, profile))
        assert cases[name]['pair_id'] == control['pair_id']
        assert not cases[name]['positive']
    assert (
        s.make_case('independent_constant', 40, 0.05, 0.02, location, 12345)['positive'] is False
    )


def test_design_counts_without_generating_evaluation_histories():
    base = math.prod(len(s.DESIGN[k]) for k in ('sizes', 'changes', 'noise', 'locations', 'seeds'))
    assert base == 400
    assert base * len(s.CONDITIONS) == 2400
    assert base // len(s.DESIGN['changes']) == 80
    assert set(s.DESIGN['seeds']).isdisjoint(d.DESIGN['seeds'])


@pytest.mark.parametrize(
    'size,shape,route',
    [
        (False, False, 'neither'),
        (True, False, 'size_only'),
        (False, True, 'shape_only'),
        (True, True, 'both'),
    ],
)
def test_rejection_routes(size, shape, route):
    direct = {
        'size_rejected_splits': [3] if size else [],
        'fit_rejected_splits': [3] if shape else [],
    }
    assert s.rejection_route(direct, 3) == route


@pytest.mark.parametrize('backend', ['python', 'native'])
def test_entire_frozen_pipeline_is_reused(backend):
    parent = json.loads((Path(__file__).parent / 'data/direct_v1_frozen.json').read_text())
    case = s.make_case('correlated_up', 40, 0.06, 0.02, 'recent', 12345)
    expected = d.evaluate(case, parent, backend)
    actual = s.evaluate(case, {'inherited': parent}, backend)
    route = actual.pop('generating_split_route')
    assert actual == expected
    if actual['direct']['has_alert']:
        assert route != 'neither'
