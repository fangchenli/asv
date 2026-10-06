"""Independent algebra and likelihood checks for the information diagnosis."""

import math
from fractions import Fraction

import numpy as np
import pytest
from scipy.stats import multivariate_normal

from experiments.step_detection import ar1_information as info
from experiments.step_detection import rational_polynomial as p
from experiments.step_detection import reporting_ar1 as a
from experiments.step_detection import reporting_covariance as gls


def fixture_values():
    indices = np.arange(14)
    return 10 + 0.6 * (indices >= 9) + 0.13 * np.sin(1.7 * indices)


@pytest.mark.parametrize('rho', [-0.7, 0, 0.6, 0.99])
@pytest.mark.parametrize('split', [3, 7, 9, 13])
def test_conditional_reparameterization_matches_original_design(rho, split):
    y = fixture_values()
    k = len(y) // 2
    design = np.column_stack((np.arange(len(y)) < split, np.arange(len(y)) >= split))
    conditional = design[k:].astype(float) - rho * design[k - 1 : -1]
    z = y[k:] - rho * y[k - 1 : -1]
    residual = z - conditional @ np.linalg.lstsq(conditional, z, rcond=None)[0]
    assert info.conditional_sse(y, split, k, rho) == pytest.approx(residual @ residual, abs=1e-12)


@pytest.mark.parametrize('rho', [-0.7, 0, 0.6, 0.99, 1])
def test_relaxation_gap_formula_and_profiled_limit(rho):
    y = fixture_values()
    k, split = 7, 9
    relaxed = float(p.evaluate(a.confidence_residual(y, split, k), Fraction(rho)))
    exact = info.conditional_sse(y, split, k, rho)
    assert exact - relaxed == pytest.approx(info.relaxation_gap(y, split, k, rho), abs=1e-12)
    if rho == 1:
        differences = np.delete(np.diff(y)[k - 1 :], split - k)
        assert exact == pytest.approx(sum((differences - differences.mean()) ** 2))
        assert info.conditional_sse(y, split, k, 1 - 1e-9) == pytest.approx(exact, rel=1e-7)


def test_profiling_then_taking_limit_differs_from_fitting_design_at_one():
    y = fixture_values() + 0.1 * np.arange(14)
    k, split = 7, 9
    design = np.column_stack((np.arange(len(y)) < split, np.arange(len(y)) >= split))
    conditional = design[k:].astype(float) - design[k - 1 : -1]
    z = np.diff(y)[k - 1 :]
    residual = z - conditional @ np.linalg.lstsq(conditional, z, rcond=None)[0]
    profiled_limit = info.conditional_sse(y, split, k, 1)
    assert residual @ residual > profiled_limit + 0.01


@pytest.mark.parametrize('rho', [-0.7, 0, 0.6, 0.99])
def test_full_likelihood_matches_stationary_gaussian_density_and_reversal(rho):
    y = fixture_values()
    split, n = 9, len(y)
    matrix = rho ** np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    precision = np.linalg.inv(matrix)
    design = np.column_stack((np.arange(n) < split, np.arange(n) >= split)).astype(float)
    levels = np.linalg.solve(design.T @ precision @ design, design.T @ precision @ y)
    residual = y - design @ levels
    sigma_squared = float(residual @ precision @ residual) / n
    expected = multivariate_normal.logpdf(y, mean=design @ levels, cov=sigma_squared * matrix)
    result = info.full_profile(y, split, rho)
    assert result['log_u'] == pytest.approx(expected, abs=1e-9)
    assert result['levels'] == pytest.approx(levels)
    reverse = info.full_profile(y[::-1], n - split, rho)
    assert reverse['log_u'] == pytest.approx(result['log_u'], abs=1e-9)
    assert reverse['levels'][::-1] == pytest.approx(result['levels'])


def test_full_likelihood_retains_stationary_endpoint_factor():
    y = fixture_values()
    split, n = 9, len(y)
    differences = np.delete(np.diff(y), split - 1)
    limit = float(differences @ differences)
    for rho in [1 - 1e-6, 1 - 1e-8, 1 - 1e-10]:
        result = info.full_profile(y, split, rho)
        determinant_term = 0.5 * (math.log1p(-rho) + math.log1p(rho))
        assert result['sse'] == pytest.approx(limit, rel=1e-3)
        assert result['log_u'] - determinant_term == pytest.approx(
            info.conditional_log_maximum(limit, n), abs=1e-3
        )


def test_archived_witness_survives_both_directions_exact_fit_and_average():
    diagnosis = info.diagnose()
    for direction in ['forward', 'reverse']:
        result = diagnosis[direction]
        for row in result['rows'][-2:]:
            assert row['exact_log_e'] < math.log(200)
            assert row['relaxed_log_e'] < math.log(200)
        witness = result['rows'][-2]
        assert witness['nonnegative_sse'] == pytest.approx(witness['exact_sse'])
        assert witness['nonnegative_log_e'] < math.log(200)
        assert result['exact_near_one_likelihood_gain_from_rho_zero'] - result[
            'prediction_regret_against_rho_zero_profile'
        ] == pytest.approx(result['rows'][-1]['exact_log_e'])
        assert result['generating_density_near_one_log_e_diagnostic_only'] > math.log(200)
    for row in diagnosis['averaged_evidence'][-2:]:
        assert row['exact_log_e_average'] < math.log(100)
        assert row['relaxed_log_e_average'] < math.log(100)
    shared = diagnosis['shared_nuisance_witness']
    assert shared['forward_log_e'] < math.log(200)
    assert shared['reverse_log_e'] < math.log(200)
    assert shared['average_log_e'] < math.log(100)
    case, _ = info.load_example()
    n, rho = len(case['values']), 4095 / 4096
    matrix = rho ** np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    reporting = gls.evidence(case['values'], matrix, a.calibration(n))
    assert reporting['status'] == 'ok'
    assert case['position'] not in reporting['size_rejected_splits']
    assert case['position'] not in reporting['fit_rejected_splits']
