"""Deterministic GLS and Gaussian-projection checks; no simulation study."""

import math

import numpy as np
import pytest

from experiments.step_detection import reporting_covariance as gls
from experiments.step_detection import reporting_direct as direct


def covariance(n, rho=0.7, amplitude=3):
    positions = np.arange(n)
    scales = np.where(positions >= n // 2, amplitude, 1.0)
    return np.outer(scales, scales) * rho ** np.abs(positions[:, None] - positions[None, :])


def fixture(n):
    return np.array([10 + 0.13 * math.sin(2.1 * i) + 0.7 * (i >= n // 2) for i in range(n)])


def independent_fit(values, matrix, design):
    # Raw precision-matrix GLS with an intercept/jump parameterization,
    # independently of the implementation's whitened plateau QR calculation.
    precision = np.linalg.inv(matrix)
    gram_inverse = np.linalg.inv(design.T @ precision @ design)
    coefficients = gram_inverse @ design.T @ precision @ values
    residual = values - design @ coefficients
    return coefficients, gram_inverse, float(residual @ precision @ residual)


@pytest.mark.parametrize(
    'n,rho,amplitude', [(4, 0.7, 3), (9, 0, 3), (13, 0.7, 3), (17, -0.4, 0.3)]
)
def test_every_statistic_against_independent_gls(n, rho, amplitude):
    y, matrix = fixture(n), covariance(n, rho, amplitude)
    calibration = direct.critical_values(n, threshold=0.07)
    result = gls.evidence(y, matrix, calibration)
    assert result['status'] == 'ok'
    for split in range(1, n):
        design = np.column_stack((np.ones(n), np.arange(n) >= split))
        beta, gram, two = independent_fit(y, matrix, design)
        contrast = np.array([-0.07, 1.0])
        expected_t = contrast @ beta / math.sqrt(two / (n - 2) * (contrast @ gram @ contrast))
        costs = [
            (independent_fit(y, matrix, np.column_stack((design, np.arange(n) >= j)))[2], j)
            for j in range(1, n)
            if j != split
        ]
        best, extra = min(costs)
        expected_f = (two - best) / (best / (n - 3))
        assert result['size_t'][split - 1] == pytest.approx(expected_t, rel=1e-8, abs=1e-8)
        assert result['lack_of_fit_f'][split - 1] == pytest.approx(expected_f, rel=1e-8, abs=1e-8)
        assert result['best_extra_split'][split - 1] == extra
    size = [i + 1 for i, x in enumerate(result['size_t']) if x > calibration['size_t_critical']]
    shape = [
        i + 1 for i, x in enumerate(result['lack_of_fit_f']) if x > calibration['fit_f_critical']
    ]
    survivors = sorted(set(range(1, n)) - set(size) - set(shape))
    assert result['size_rejected_splits'] == size
    assert result['fit_rejected_splits'] == shape
    assert result['null_splits'] == survivors
    assert result['has_alert'] == (not survivors)


@pytest.mark.parametrize('n,threshold', [(4, 0), (17, 0.05), (40, 0.1), (100, 0.05)])
def test_identity_covariance_recovers_previous_test(n, threshold):
    y = fixture(n)
    calibration = direct.critical_values(n, threshold=threshold)
    old = direct.evidence(y, calibration)
    new = gls.evidence(y, np.eye(n), calibration)
    assert new['status'] == 'ok'
    for key in old:
        if key in ('size_t', 'lack_of_fit_f'):
            assert new[key] == pytest.approx(old[key], rel=1e-8, abs=1e-8)
        else:
            assert new[key] == old[key]


@pytest.mark.parametrize('rho', [0, 0.7, -0.4])
def test_projection_identities_underpin_t_and_f_laws(rho):
    n, split, extra = 9, 4, 7
    matrix = covariance(n, rho)
    # Use a symmetric spectral square root, different from Cholesky whitening.
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    whiten = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T
    np.testing.assert_allclose(whiten @ matrix @ whiten.T, np.eye(n), atol=1e-13)
    two = np.column_stack((np.ones(n), np.arange(n) >= split))
    three = np.column_stack((two, np.arange(n) >= extra))
    b2, b3 = whiten @ two, whiten @ three
    g2 = np.linalg.inv(b2.T @ b2)
    p2 = b2 @ g2 @ b2.T
    p3 = b3 @ np.linalg.inv(b3.T @ b3) @ b3.T
    added, remaining = p3 - p2, np.eye(n) - p3
    for projection, rank in [(p2, 2), (added, 1), (remaining, n - 3)]:
        np.testing.assert_allclose(projection @ projection, projection, atol=1e-13)
        np.testing.assert_allclose(projection.T, projection, atol=1e-13)
        assert np.trace(projection) == pytest.approx(rank)
    np.testing.assert_allclose(added @ remaining, np.zeros((n, n)), atol=1e-13)
    np.testing.assert_allclose(added @ b2, np.zeros((n, 2)), atol=1e-13)
    np.testing.assert_allclose(remaining @ b2, np.zeros((n, 2)), atol=1e-13)
    h = np.array([-0.05, 1.0])
    direction = b2 @ g2 @ h / math.sqrt(h @ g2 @ h)
    assert direction @ direction == pytest.approx(1)
    np.testing.assert_allclose(direction @ (np.eye(n) - p2), np.zeros(n), atol=1e-13)


def test_general_covariance_with_varying_diagonal_and_dense_correlations():
    n, split = 9, 4
    positions = np.arange(n)
    loading = np.column_stack((np.sin(positions), np.cos(positions / 2)))
    matrix = np.diag(1 + positions / 3) + loading @ loading.T
    y = fixture(n)
    design = np.column_stack((np.ones(n), positions >= split))
    beta, gram, two = independent_fit(y, matrix, design)
    h = np.array([-0.05, 1.0])
    result = gls.evidence(y, matrix, direct.critical_values(n))
    assert result['size_t'][split - 1] == pytest.approx(
        h @ beta / math.sqrt(two / (n - 2) * (h @ gram @ h))
    )
    costs = [
        independent_fit(y, matrix, np.column_stack((design, positions >= extra)))[2]
        for extra in range(1, n)
        if extra != split
    ]
    assert result['lack_of_fit_f'][split - 1] == pytest.approx(
        (two - min(costs)) / (min(costs) / (n - 3))
    )


def test_known_unequal_variances_give_correct_contrast_and_scale():
    n, split = 40, 36
    matrix = np.diag([1.0] * split + [9.0] * (n - split))
    design = np.column_stack((np.arange(n) < split, np.arange(n) >= split)).astype(float)
    y = fixture(n)
    beta, gram, cost = independent_fit(y, matrix, design)
    np.testing.assert_allclose(beta, [np.mean(y[:split]), np.mean(y[split:])])
    h = np.array([-1.05, 1.0])
    assert h @ gram @ h == pytest.approx(2.280625)
    expected_cost = sum((y[:split] - beta[0]) ** 2) + sum((y[split:] - beta[1]) ** 2) / 9
    assert cost == pytest.approx(expected_cost)
    precision = np.linalg.inv(matrix)
    residual_operator = np.eye(n) - design @ gram @ design.T @ precision
    # E(Q)/sigma^2 = trace(A' R^-1 A R) = n-2 at a true split.
    assert np.trace(residual_operator.T @ precision @ residual_operator @ matrix) == pytest.approx(
        n - 2
    )
    result = gls.evidence(y, matrix, direct.critical_values(n))
    assert result['size_t'][split - 1] == pytest.approx(h @ beta / math.sqrt(cost / 38 * 2.280625))


@pytest.mark.parametrize(
    'scale,which', [(1e-100, 'y'), (1000, 'y'), (1e100, 'y'), (1e-100, 'R'), (1e100, 'R')]
)
def test_unit_and_covariance_scale_invariance(scale, which):
    y, matrix = fixture(13), covariance(13)
    calibration = direct.critical_values(len(y))
    before = gls.evidence(y, matrix, calibration)
    after = gls.evidence(
        y * (scale if which == 'y' else 1), matrix * (scale if which == 'R' else 1), calibration
    )
    assert before['status'] == after['status'] == 'ok'
    for key in ('size_t', 'lack_of_fit_f'):
        assert before[key] == pytest.approx(after[key], rel=1e-8, abs=1e-8)
    for key in (
        'has_alert',
        'null_splits',
        'size_rejected_splits',
        'fit_rejected_splits',
        'best_extra_split',
    ):
        assert before[key] == after[key]


@pytest.mark.parametrize('change,expected', [(0.04, False), (0.08, True)])
def test_complete_decision_on_resolved_history(change, expected):
    n = 16
    y = [10 + 0.001 * math.sin(i * 2.1) + 10 * change * (i >= 8) for i in range(n)]
    result = gls.evidence(y, covariance(n), direct.critical_values(n))
    assert result['status'] == 'ok'
    assert result['has_alert'] is expected
    if not expected:
        assert 8 in result['null_splits']
        assert result['fit_rejected_splits']


@pytest.mark.parametrize('values', [[0] * 8, [10] * 8, [10] * 4 + [10.8] * 4])
def test_unresolved_residual_variation_abstains(values):
    result = gls.evidence(values, covariance(8), direct.critical_values(8))
    assert result['status'] == 'insufficient_precision'
    assert result['has_alert'] is False
    assert result['null_splits'] is None
    assert result['size_t'] is None
    assert result['lack_of_fit_f'] is None


@pytest.mark.parametrize(
    'matrix',
    [
        np.eye(3),
        np.zeros((4, 4)),
        np.ones((4, 4)),
        np.diag([1, 1, 1, -1]),
        np.diag([1, 1, 1, 1e-14]),
        np.diag([1, 1, 1, math.nan]),
        np.diag([1, 1, 1, math.inf]),
        np.eye(4) + np.triu(np.ones((4, 4)), 1),
    ],
)
def test_bad_covariance(matrix):
    with pytest.raises(ValueError):
        gls.evidence(fixture(4), matrix, direct.critical_values(4))


@pytest.mark.parametrize(
    'values', [[], [10] * 3, [10] * 201, [[10, 11]] * 4, [10, 11, 12, math.nan]]
)
def test_bad_observations(values):
    with pytest.raises(ValueError):
        gls.evidence(values, np.eye(4), direct.critical_values(4))


@pytest.mark.parametrize(
    'key,value',
    [
        ('n', 5),
        ('threshold', -1),
        ('threshold', math.inf),
        ('size_t_critical', 0),
        ('fit_f_critical', math.nan),
    ],
)
def test_bad_calibration(key, value):
    calibration = direct.critical_values(4)
    calibration[key] = value
    with pytest.raises(ValueError):
        gls.evidence(fixture(4), np.eye(4), calibration)
