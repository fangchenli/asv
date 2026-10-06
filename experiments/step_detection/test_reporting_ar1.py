"""Check conditional likelihoods, exact algebra, and continuous certificates."""

import math
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from experiments.step_detection import rational_polynomial as p
from experiments.step_detection import reporting_ar1 as a
from experiments.step_detection import reporting_covariance as gls


def values(n=12, change=0.08):
    return [10 + 0.1 * (-1) ** i + 10 * change * (i >= n // 2) for i in range(n)]


@pytest.mark.parametrize('left,right', [(F(-1), F(1)), (F(-3, 4), F(1, 8)), (F(1, 3), F(3, 4))])
def test_bernstein_reconstructs_polynomial_exactly(left, right):
    poly = p.polynomial([3, -2, F(1, 7), 5, -3])
    coefficients = p.bernstein(poly, left, right)
    degree = len(poly) - 1
    for u in [F(0), F(1, 7), F(1, 2), F(4, 5), F(1)]:
        reconstructed = sum(
            coefficients[k] * math.comb(degree, k) * u**k * (1 - u) ** (degree - k)
            for k in range(degree + 1)
        )
        assert reconstructed == p.evaluate(poly, left + (right - left) * u)


def test_grid_can_miss_narrow_surviving_region_but_certificate_cannot():
    center, radius = F(1, 32), F(1, 4096)
    poly = p.polynomial([center**2 - radius**2, -2 * center, 1])
    assert all(p.evaluate(poly, F(i, 4)) > 0 for i in range(-4, 5))
    assert p.evaluate(poly, center) < 0
    assert not p.positive_on(poly, F(-1), F(1))


def test_subdivision_can_certify_when_whole_interval_is_inconclusive():
    poly = p.polynomial([F(1, 10), 0, 1])
    assert not p.positive_on(poly, F(-1), F(1))
    assert p.positive_on(poly, F(-1), F(0))
    assert p.positive_on(poly, F(0), F(1))
    # An internal zero is not a strict-positive certificate.
    assert not p.positive_on(p.polynomial([0, 0, 1]), F(0), F(1))
    assert p.positive_on(p.polynomial([1, 1]), F(-1), F(1))


def test_removing_endpoint_factors_preserves_sign_and_exact_value():
    base = p.polynomial([2, 1])
    poly = p.multiply(
        p.multiply(p.polynomial([1, -1]), p.polynomial([1, -1])),
        p.multiply(p.polynomial([1, 1]), base),
    )
    assert p.remove_stationary_endpoint_factors(poly) == base
    for rho in [F(-7, 8), F(0), F(7, 8)]:
        assert p.evaluate(poly, rho) == (1 - rho) ** 2 * (1 + rho) * p.evaluate(base, rho)


def test_predictive_mixture_matches_direct_gaussian_densities():
    y = np.array(values()) / 10
    training = 6
    fitted = a.predictor(y[:training])
    rho, variance, prior = (fitted[k] for k in ('rho', 'variance', 'jump_variance'))
    z = y[training:] - rho * y[training - 1 : -1] - (1 - rho) * fitted['mean']
    m = len(z)
    log_densities = [multivariate_normal.logpdf(z, mean=np.zeros(m), cov=variance * np.eye(m))]
    for j in range(m):
        vector = np.zeros(m)
        vector[j], vector[j + 1 :] = 1, 1 - rho
        matrix = variance * np.eye(m) + prior * np.outer(vector, vector)
        log_densities.append(multivariate_normal.logpdf(z, mean=np.zeros(m), cov=matrix))
    assert a.predictive_log_density(y, training, fitted) == pytest.approx(
        float(logsumexp(log_densities) - math.log(m + 1)), rel=1e-8, abs=1e-8
    )
    changed = np.array(y)
    changed[training:] += 3
    assert a.confidence_threshold(changed, 0.01)['predictor'] == fitted
    # Check the proper density in observation coordinates, including dependence
    # on the final training value and all earlier validation observations.
    transform = np.eye(m) - rho * np.eye(m, k=-1)
    assert np.linalg.det(transform) == pytest.approx(1)
    offset = np.full(m, (1 - rho) * fitted['mean'])
    offset[0] += rho * y[training - 1]
    observation_mean = np.linalg.solve(transform, offset)
    inverse = np.linalg.inv(transform)
    observation_logs = []
    for j in range(m + 1):
        vector = np.zeros(m)
        if j < m:
            vector[j], vector[j + 1 :] = 1, 1 - rho
        innovation_covariance = variance * np.eye(m) + prior * np.outer(vector, vector)
        observation_logs.append(
            multivariate_normal.logpdf(
                y[training:],
                mean=observation_mean,
                cov=inverse @ innovation_covariance @ inverse.T,
            )
        )
    assert a.predictive_log_density(y, training, fitted) == pytest.approx(
        float(logsumexp(observation_logs) - math.log(m + 1)),
        rel=1e-8,
        abs=1e-8,
    )


@pytest.mark.parametrize('split', [2, 6, 8, 11])
@pytest.mark.parametrize('rho', [F(-2, 5), F(0), F(7, 10)])
def test_conditional_quadratic_matches_relaxed_likelihood_and_bounds_true_fit(split, rho):
    y = np.array(values()) + 0.017 * np.sin(np.arange(12) * 2.1)
    training, n = 6, len(y)
    poly = a.confidence_residual(y, split, training)
    z = y[training:] - float(rho) * y[training - 1 : -1]
    indices = np.arange(training, n)
    enlarged = np.column_stack((indices < split, indices == split, indices > split)).astype(float)
    residual = z - enlarged @ np.linalg.lstsq(enlarged, z, rcond=None)[0]
    assert float(p.evaluate(poly, rho)) == pytest.approx(residual @ residual, abs=1e-11)
    original = np.column_stack((np.arange(n) < split, np.arange(n) >= split)).astype(float)
    conditional = original[training:] - float(rho) * original[training - 1 : -1]
    actual_residual = z - conditional @ np.linalg.lstsq(conditional, z, rcond=None)[0]
    assert residual @ residual <= actual_residual @ actual_residual + 1e-11


def test_residual_threshold_equals_likelihood_ratio_threshold():
    y = np.array(values()) / 10
    info = a.confidence_threshold(y, 0.01)
    m, log_q, bound = (
        info[k] for k in ('validation_count', 'log_predictive_density', 'residual_bound')
    )
    log_u = -m / 2 * (math.log(2 * math.pi) + 1 + math.log(bound / m))
    assert log_q - log_u == pytest.approx(math.log(100), abs=1e-10)
    for factor in [0.9, 1.1]:
        log_u = -m / 2 * (math.log(2 * math.pi) + 1 + math.log(factor * bound / m))
        assert (log_q - log_u > math.log(100)) == (factor > 1)


def test_unrepresentable_density_retains_all_correlations(monkeypatch):
    monkeypatch.setattr(a, 'predictive_log_density', lambda *args: -math.inf)
    result = a.confidence_threshold(values(), 0.01)
    assert result['residual_bound'] == math.inf


def test_error_budget_and_modified_cutoffs():
    config = a.calibration(12)
    assert config['size_alpha'] == pytest.approx(0.032)
    assert config['fit_alpha'] == pytest.approx(0.008)
    assert config['size_alpha'] + config['fit_alpha'] + config[
        'confidence_alpha'
    ] == pytest.approx(0.05)
    config['fit_f_critical'] *= 0.9
    with pytest.raises(ValueError, match='unmodified'):
        a.evidence(values(), config)


@pytest.mark.parametrize('rho', [F(-4, 5), F(0), F(7, 10)])
def test_gls_residuals_and_all_rejection_polynomials(rho):
    n = 9
    y = [10 + 0.17 * math.sin(2.1 * i) + 0.6 * (i >= 4) for i in range(n)]
    models = a.Models(y)
    config = a.calibration(n)
    indices = np.arange(n)
    matrix = float(rho) ** np.abs(indices[:, None] - indices[None, :])
    reference = gls.evidence(y, matrix, config)
    precision = np.linalg.inv(matrix) * (1 - float(rho) ** 2)
    for split in range(1, n):
        d, numerator, adjugate, linear = models.fit((split,))
        design = np.column_stack((indices < split, indices >= split)).astype(float)
        gram = design.T @ precision @ design
        beta = np.linalg.solve(gram, design.T @ precision @ y)
        residual = y - design @ beta
        assert float(p.evaluate(d, rho)) == pytest.approx(np.linalg.det(gram))
        assert float(p.evaluate(numerator, rho) / p.evaluate(d, rho)) == pytest.approx(
            residual @ precision @ residual, abs=1e-10
        )
        size = models.size(split, config)
        assert all(p.evaluate(poly, rho) > 0 for poly in size) == (
            split in reference['size_rejected_splits']
        )
        any_shape = False
        for extra in range(1, n):
            if extra == split:
                continue
            poly = models.shape(split, extra, config)
            extended = np.column_stack((design, indices >= extra))
            coefficients = np.linalg.solve(
                extended.T @ precision @ extended, extended.T @ precision @ y
            )
            remaining = y - extended @ coefficients
            q3 = remaining @ precision @ remaining
            statistic = (residual @ precision @ residual - q3) / (q3 / (n - 3))
            assert (p.evaluate(poly, rho) > 0) == (statistic > config['fit_f_critical'])
            any_shape |= p.evaluate(poly, rho) > 0
        assert any_shape == (split in reference['fit_rejected_splits'])


def test_alert_certificate_covers_every_split_and_interval():
    y = values()
    result = a.evidence(y, max_cells=256, max_depth=8)
    assert result['status'] == 'certified_alert'
    assert result['has_alert']
    assert len(result['confidence_polynomials']) == len(y) - 1
    scale = result['confidence']['normalization_scale']
    normalized = [v / scale for v in y]
    models = a.Models(normalized)
    for split in range(1, len(y)):
        cells = sorted(
            (F(row['left']), F(row['right']), row)
            for row in result['certificate']
            if row['split'] == split
        )
        assert cells[0][0] == -1 and cells[-1][1] == 1
        assert all(cells[i][1] == cells[i + 1][0] for i in range(len(cells) - 1))
        for left, right, row in cells:
            if row['route'] == 'confidence':
                residual = a.confidence_residual(models.values, split, len(y) // 2)
                polynomials = (
                    p.subtract(residual, p.polynomial([result['confidence']['residual_bound']])),
                )
            elif row['route'] == 'size':
                polynomials = models.size(split, result['calibration'])
            else:
                polynomials = (models.shape(split, row['extra'], result['calibration']),)
            assert all(p.positive_on(poly, left, right) for poly in polynomials)


def test_surviving_pair_is_valid_in_independent_gls_check():
    y = values(16, 0.04)
    result = a.evidence(y)
    assert result['status'] == 'surviving_explanation'
    assert not result['has_alert']
    split, rho = result['witness']['split'], result['witness']['rho_float']
    indices = np.arange(len(y))
    reference = gls.evidence(
        y, rho ** np.abs(indices[:, None] - indices[None, :]), result['calibration']
    )
    assert split in reference['null_splits']
    poly = tuple(F(x) for x in result['confidence_polynomials'][str(split)])
    assert p.evaluate(poly, F(result['witness']['rho'])) <= F(
        result['confidence']['residual_bound']
    )


@pytest.mark.parametrize('kwargs', [{'max_cells': 1}, {'max_depth': 0}])
def test_work_limits_abstain_instead_of_using_unchecked_gaps(kwargs):
    result = a.evidence(values(), **kwargs)
    assert result['status'] == 'unresolved'
    assert not result['has_alert']


@pytest.mark.parametrize('y', [[10] * 12, [10] * 6 + [10.8] * 6])
def test_exactly_zero_candidate_variation_abstains(y):
    result = a.evidence(y)
    assert result['status'] == 'insufficient_variation'
    assert not result['has_alert']


@pytest.mark.parametrize('scale', [1e-100, 1000, 1e100])
def test_timing_units_preserve_decisions(scale):
    original = a.evidence(values())
    changed = a.evidence([x * scale for x in values()])
    assert changed['status'] == original['status'] == 'certified_alert'
    assert changed['confidence']['residual_bound'] == pytest.approx(
        original['confidence']['residual_bound']
    )


@pytest.mark.parametrize(
    'kwargs',
    [
        {'alpha': 0.01},
        {'confidence_alpha': 0},
        {'threshold': -1},
        {'alpha': 0.99, 'confidence_alpha': 0.01},
    ],
)
def test_invalid_calibration(kwargs):
    with pytest.raises(ValueError):
        a.calibration(12, **kwargs)


@pytest.mark.parametrize('kwargs', [{'max_cells': 0}, {'max_cells': True}, {'max_depth': -1}])
def test_invalid_certification_budget(kwargs):
    with pytest.raises(ValueError):
        a.evidence(values(), **kwargs)


@pytest.mark.parametrize('y', [[10] * 7, [10] * 201, [10] * 11 + [math.nan]])
def test_invalid_observations(y):
    with pytest.raises(ValueError):
        a.evidence(y)
