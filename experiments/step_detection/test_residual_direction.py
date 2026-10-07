"""Check the nuisance-free density, endpoint limits, and proposed bound."""

import math
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import gamma

from experiments.step_detection import rational_polynomial as p
from experiments.step_detection import residual_direction as direction


def values(n=12):
    return [10 + (i >= n // 2) + math.sin(1.7 * i) + 0.2 * math.cos(0.3 * i) for i in range(n)]


@pytest.mark.parametrize('split', [1, 4, 7])
@pytest.mark.parametrize('rho', [F(-9, 10), F(0), F(7, 10), F(9999, 10000)])
def test_exact_scalar_identity_matches_independent_dense_projection(split, rho):
    y = values(8)
    model = direction.state(y, split)
    actual = direction.log_fraction(direction.squared_density(model, rho)) / 2
    assert actual == pytest.approx(direction.matrix_log_density(y, split, float(rho)), abs=2e-9)
    basis = direction.residual_basis(8, split)
    design = np.column_stack((np.arange(8) < split, np.arange(8) >= split))
    assert basis.T @ basis == pytest.approx(np.eye(6), abs=1e-14)
    assert basis.T @ design == pytest.approx(np.zeros((6, 2)), abs=1e-14)
    distance = np.abs(np.arange(8)[:, None] - np.arange(8)[None, :])
    covariance = basis.T @ (float(rho) ** distance) @ basis
    determinant = (
        (1 - float(rho) ** 2) ** 5
        * (1 - float(rho))
        * float(p.evaluate(model['b'], rho))
        / (split * (8 - split))
    )
    assert np.linalg.det(covariance) == pytest.approx(determinant, rel=2e-9, abs=1e-30)


@pytest.mark.parametrize('rho', [-0.8, 0.0, 0.85])
def test_circle_density_integrates_to_one(rho):
    basis = direction.residual_basis(4, 2)

    def density(angle):
        y = basis @ np.array([math.cos(angle), math.sin(angle)])
        return math.exp(direction.matrix_log_density(y, 2, rho)) / (2 * math.pi)

    integral, error = quad(density, 0, 2 * math.pi, epsabs=1e-11)
    assert integral == pytest.approx(1, abs=1e-10)
    assert error < 1e-8


def test_polar_integration_of_gaussian_recovers_density_and_removes_scale():
    y, split, rho = np.array(values(5)), 2, 0.6
    basis = direction.residual_basis(len(y), split)
    z = basis.T @ y
    s = z / np.linalg.norm(z)
    distance = np.abs(np.arange(len(y))[:, None] - np.arange(len(y))[None, :])
    covariance = basis.T @ (rho**distance) @ basis
    d = len(s)
    quadratic = float(s @ np.linalg.solve(covariance, s))
    area = 2 * math.pi ** (d / 2) / gamma(d / 2)
    expected = math.exp(direction.matrix_log_density(y, split, rho))
    for sigma in [0.1, 1, 10]:
        constant = (2 * math.pi * sigma**2) ** (-d / 2) / math.sqrt(np.linalg.det(covariance))
        integrated, _ = quad(
            lambda radius, constant=constant, sigma=sigma: area
            * constant
            * radius ** (d - 1)
            * math.exp(-(radius**2) * quadratic / (2 * sigma**2)),
            0,
            np.inf,
            epsabs=1e-10,
        )
        assert integrated == pytest.approx(expected, rel=1e-9)


def test_exact_invariance_to_levels_scale_and_reversal():
    y = [F(x) for x in values(8)]
    original = direction.state(y, 3)
    shifted = direction.state([7 * x + (F(19, 7) if i < 3 else -100) for i, x in enumerate(y)], 3)
    reversed_model = direction.state(y[::-1], 5)
    for rho in [F(-9, 10), F(0), F(9, 10), F(1)]:
        density = direction.squared_density(original, rho)
        assert density == direction.squared_density(shifted, rho)
        assert density == direction.squared_density(reversed_model, rho)
    assert direction.squared_density(original, 0) == 1
    assert direction.predictor_lower(original) == direction.predictor_lower(shifted)


@pytest.mark.parametrize('split', [1, 5, 11])
def test_positive_endpoint_limit_is_within_plateau_differences(split):
    y = [F(x) for x in values()]
    model = direction.state(y, split)
    cost = sum((y[i] - y[i - 1]) ** 2 for i in range(1, len(y)) if i != split)
    assert p.evaluate(model['num'], F(1)) / p.evaluate(model['den'], F(1)) == cost
    expected = split * (len(y) - split) * (model['ss'] / cost) ** (len(y) - 2)
    assert direction.squared_density(model, F(1)) == expected
    assert float(direction.squared_density(model, F(1) - F(1, 2**24))) == pytest.approx(
        float(expected), rel=1e-4
    )


def test_negative_endpoint_generic_density_vanishes_but_alternating_density_diverges():
    generic = direction.state(values(6), 3)
    alternating = direction.state([(-1) ** i + 10 + (i >= 3) for i in range(6)], 3)
    assert p.evaluate(generic['num'], F(-1)) > 0
    assert p.evaluate(alternating['num'], F(-1)) == 0
    eps = F(1, 2**20)
    generic_ratio = direction.squared_density(generic, -1 + eps / 2) / direction.squared_density(
        generic, -1 + eps
    )
    pole_ratio = direction.squared_density(alternating, -1 + eps / 2) / direction.squared_density(
        alternating, -1 + eps
    )
    assert float(generic_ratio) == pytest.approx(0.5, rel=1e-4)
    assert float(pole_ratio) == pytest.approx(2 ** (6 - 3), rel=1e-4)
    assert not direction.interval_excluded(alternating, -1, -1 + eps, F(1))


@pytest.mark.parametrize('split', [1, 2])
def test_three_readings_have_no_directional_correlation_information(split):
    model = direction.state([1, 3, 2], split)
    for rho in [F(-999, 1000), F(0), F(3, 4), F(1)]:
        assert direction.squared_density(model, rho) == 1
    assert direction.predictor_lower(model) == 1


@pytest.mark.parametrize('value', [F(0), F(2), F(9, 4), F(1, 10**100), F(10**100)])
def test_integer_square_root_proves_two_sided_enclosure(value):
    lower = direction.sqrt_lower(value)
    assert lower**2 <= value < (lower + F(1, 2**128)) ** 2


def test_interval_bound_is_nonvacuous_and_implies_pointwise_exclusion():
    model = direction.state([(-1) ** i for i in range(16)], 8)
    left, right = F(99, 100), F(1)
    assert direction.interval_excluded(model, left, right, F(1))
    for j in range(21):
        rho = left + (right - left) * F(j, 20)
        assert direction.DELTA**2 > direction.squared_density(model, rho)
    assert not direction.interval_excluded(model, left, right, F(0))


def test_invalid_and_degenerate_inputs_do_not_create_evidence():
    with pytest.raises(ValueError, match='nonzero'):
        direction.state([10, 10, 11, 11], 2)
    with pytest.raises(ValueError, match='interior'):
        direction.state([1, 2, 3, 4], 0)
    with pytest.raises(ValueError, match='stationary'):
        direction.matrix_log_density(values(), 6, 1)
    with pytest.raises(ValueError, match='nonnegative'):
        direction.sqrt_lower(-1)
    with pytest.raises(ValueError, match='stationary'):
        direction.predictor_lower(direction.state(values(), 6), [F(1)])


def test_directional_statistic_is_a_weighted_chi_square_ratio():
    basis = direction.residual_basis(8, 3)
    distance = np.abs(np.arange(8)[:, None] - np.arange(8)[None, :])
    covariance = basis.T @ (0.7**distance) @ basis
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    xi = np.array([1, -2, 0.5, 3, -0.25, 1.5])
    for sigma in [0.01, 1, 100]:
        z = sigma * eigenvectors @ (np.sqrt(eigenvalues) * xi)
        s = z / np.linalg.norm(z)
        observed = 1 / float(s @ np.linalg.solve(covariance, s))
        expected = float(eigenvalues @ (xi * xi) / (xi @ xi))
        assert observed == pytest.approx(expected, rel=1e-12)


def test_archived_density_comparisons_and_decisions_are_unambiguous():
    result = direction.diagnose()
    for case in result['cases']:
        assert F(case['q_upper']) - F(case['q_lower']) == F(1, 2**128)
        for row in case['rows']:
            assert row['mixture_excluded_exact'] != row['mixture_retained_exact']
            if not row['positive_endpoint_limit']:
                assert row['log_direction_density'] == pytest.approx(
                    row['dense_log_density'], abs=2e-9
                )
    independent, positive = result['cases']
    assert independent['witness_to_one_interval']['mixture_excluded_exact']
    assert not positive['witness_to_one_interval']['mixture_excluded_exact']
    row = next(row for row in positive['rows'] if row['rho'] == positive['witness'])
    assert row['mixture_retained_exact']
