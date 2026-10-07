"""Analytical probability and matrix checks for the spectral refinement."""

import gzip
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import f as f_distribution

from experiments.step_detection import directional_determinant as determinant
from experiments.step_detection import directional_spectral as spectral
from experiments.step_detection import residual_direction as direction


@pytest.mark.parametrize(
    'count,constant', [(4, F(1, 4)), (8, F(1, 8)), (12, F(3, 32)), (16, F(5, 64))]
)
def test_density_constants_from_integral_recurrence(count, constant):
    assert spectral.density_constant(count) == constant
    if count > 4:
        m = count // 4
        assert constant == spectral.density_constant(count - 4) * F(2 * m - 3, 2 * m - 2)


@pytest.mark.parametrize('count', [4, 8, 12, 16])
@pytest.mark.parametrize('negative_count', [1, 4, 20])
def test_bound_dominates_closed_form_gaussian_square_probability(count, negative_count):
    for positive in (F(1, 2), F(4), F(20)):
        negative = F(3, 4)
        probability = f_distribution.cdf(
            float(negative * negative_count / (positive * count)), count, negative_count
        )
        bound = spectral.coefficient_bound_squared(
            [positive] * count + [-negative] * negative_count, F(1, 8), count
        )
        assert probability**2 <= float(bound)
        assert 0 < bound <= 1


def test_do_not_refine_a_chernoff_bound_after_clipping_it_at_one():
    # Its exponential moment is huge despite a useful tilted-density factor.
    # Clipping that moment first would manufacture a probability bound below 1.
    weights = [F(100)] * 8 + [F(-9, 10)] * 1000
    assert spectral.coefficient_bound_squared(weights, F(1, 4), 8) == 1


@pytest.mark.parametrize('n,split,count', [(6, 1, 4), (14, 7, 12), (24, 23, 16), (100, 1, 8)])
@pytest.mark.parametrize(
    'left,right', [(F(-7, 8), F(-13, 16)), (F(-1, 8), F(1, 8)), (F(13, 16), F(7, 8))]
)
def test_uniform_eigenvalue_bound_against_projected_covariance(n, split, count, left, right):
    lower = spectral.eigenvalue_lower(n, count, left, right)
    assert lower > 0
    assert lower == spectral.eigenvalue_lower(n, count, -right, -left)
    basis = direction.residual_basis(n, split)
    distance = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    for rho in (left, (3 * left + right) / 4, (left + right) / 2, right):
        covariance = basis.T @ (float(rho) ** distance) @ basis
        assert float(lower) <= np.linalg.eigvalsh(covariance)[-count] + 1e-13


def saved_case():
    path = Path(__file__).parent / 'data' / 'determinant_reporting_v1_diagnosis.json.gz'
    return json.loads(gzip.decompress(path.read_bytes()))['rows'][0]['case']


def test_new_saved_witness_and_its_neighborhood_are_certified():
    case = saved_case()
    model = direction.state(case['values'], 1)
    rho, width = F(13, 16), F(1, 4096)
    proof = spectral.certify_interval(model, rho - width, rho + width)
    assert proof['status'] == 'certified_excluded'
    assert F(proof['p_upper']) < F(736, 100000)
    assert proof['selected_tilt'] == '1/8' and proof['selected_count'] == 8
    assert determinant.certify_interval(model, rho, rho)['status'] == 'not_certified'
    for attempt in proof['attempts']:
        assert F(attempt['p_upper_squared']) <= F(attempt['base']['p_upper_squared'])
    for candidate in (rho - width, rho, rho + width):
        point = spectral.certify_interval(model, candidate, candidate)
        assert F(point['p_upper_squared']) <= F(proof['p_upper_squared'])


def test_scale_and_plateau_level_invariance():
    values = [F(10 + (i >= 5) + (-1) ** i) for i in range(16)]
    transformed = [7 * value + (3 if i < 5 else -5) for i, value in enumerate(values)]
    proofs = [
        spectral.certify_interval(direction.state(y, 5), F(3, 4), F(3, 4))
        for y in (values, transformed)
    ]
    assert proofs[0] == proofs[1]


def test_uninformative_cases_and_wide_intervals_remain_non_exclusions():
    model = direction.state([1, 4, 2, 5, 3, 8], 2)
    assert spectral.certify_interval(model, 0, 0)['p_upper_squared'] == '1'
    assert spectral.certify_interval(model, F(-99, 100), F(99, 100))['status'] == 'not_certified'
    model = direction.state([1, 4, 2], 1)
    assert spectral.certify_interval(model, F(7, 8), F(7, 8))['p_upper_squared'] == '1'


@pytest.mark.parametrize('count', [0, 3, 5, True])
def test_invalid_direction_counts(count):
    with pytest.raises(ValueError):
        spectral.density_constant(count)


def test_invalid_spectral_domain():
    with pytest.raises(ValueError):
        spectral.eigenvalue_lower(8, 8, F(1, 2), F(3, 4))
    with pytest.raises(ValueError):
        spectral.eigenvalue_lower(100, 8, 0, 1)
