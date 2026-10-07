"""Independent checks of the tilted tail bound and its interval certificate."""

import math
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.stats import f as f_distribution

from experiments.step_detection import directional_tail as tail
from experiments.step_detection import residual_direction as direction


def differences(n, split):
    edges = [i for i in range(n - 1) if i != split - 1]
    matrix = np.zeros((n, n - 2))
    for j, i in enumerate(edges):
        matrix[i, j], matrix[i + 1, j] = 1, -1
    return matrix, edges


@pytest.mark.parametrize('n,split', [(6, 1), (9, 4), (12, 9), (100, 25)])
def test_difference_covariance_and_four_eigenvalue_enclosure(n, split):
    matrix, edges = differences(n, split)
    gram = matrix.T @ matrix
    eigenvalues = np.linalg.eigvalsh(gram)
    assert eigenvalues[3] <= float(tail.gram_eigenvalue_upper(n, split))
    assert np.linalg.det(gram) == pytest.approx(split * (n - split), rel=1e-11)
    for rho in [0, 0.7, 0.999]:
        distance = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        actual = matrix.T @ (rho**distance) @ matrix / (1 - rho)
        expected = np.array(
            [[2 if i == j else -(1 - rho) * rho ** (abs(i - j) - 1) for j in edges] for i in edges]
        )
        assert actual == pytest.approx(expected, abs=2e-12)
        lower = 2 - (n - 3) * (1 - rho)
        assert np.linalg.eigvalsh(actual)[0] >= lower - 1e-10


@pytest.mark.parametrize('n,split', [(3, 1), (4, 2), (12, 1), (12, 7), (100, 25)])
def test_tridiagonal_recurrence_matches_dense_determinant(n, split):
    matrix, _ = differences(n, split)
    gram = matrix.T @ matrix
    a, c = F(7, 8), F(3, 17)
    expected = np.linalg.det(float(a) * gram + float(c) * np.eye(n - 2)) / np.linalg.det(gram)
    assert float(tail.determinant_lower(n, split, c, a)) == pytest.approx(expected, rel=1e-11)


@pytest.mark.parametrize('positive', [F(4), F(20), F(100)])
@pytest.mark.parametrize('negative_count', [1, 4, 10])
def test_tilted_bound_dominates_closed_form_f_tail(positive, negative_count):
    negative = F(1, 2)
    exact_probability = f_distribution.cdf(
        float(negative * negative_count / (4 * positive)), 4, negative_count
    )
    bound = tail.tilted_bound_squared([positive] * 4 + [-negative] * negative_count, F(1, 16))
    assert exact_probability**2 <= float(bound)
    assert 0 < bound <= 1


def test_density_refinement_improves_chernoff_without_changing_the_event():
    weights, tilt = [F(40)] * 4 + [F(-1, 2)] * 20, F(1, 16)
    chernoff_squared = 1 / math.prod(1 + 2 * tilt * a for a in weights)
    assert tail.tilted_bound_squared(weights, tilt) < chernoff_squared < 1
    assert tail.tilted_bound_squared([F(0)] * 4, tilt) == 1


@pytest.mark.parametrize('value', [F(-1, 3), F(0), F(1, 3), F(1234567, 13)])
def test_outward_rounding_is_rigorous(value):
    rounded = tail.ceil_dyadic(value)
    assert value <= rounded < value + F(1, 2**tail.BITS)


def test_saved_correlated_interval_has_an_exact_certificate():
    case, _ = direction.load_positive_example()
    model = direction.state(case['values'], case['position'])
    left = F(16383, 16384)
    certificate = tail.certify_interval(model, left, F(1))
    assert certificate['status'] == 'certified_excluded'
    assert F(certificate['p_upper']) < F(82, 10000)
    assert F(certificate['chernoff_upper_squared']) > direction.DELTA**2
    determinant = F(certificate['determinant_lower'])
    correction = F(certificate['density_correction'])
    assert determinant * correction**2 > 1 / direction.DELTA**2

    y = np.array(case['values'])
    n, split = len(y), case['position']
    basis = direction.residual_basis(n, split)
    y[:split] -= y[:split].mean()
    y[split:] -= y[split:].mean()
    s = basis.T @ y
    s /= np.linalg.norm(s)
    distance = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    radius_upper = float(F(certificate['radius_upper']))
    lower_weight = float(F(certificate['four_weights_lower']))
    for rho in [float(left), float((left + 1) / 2), 1 - 1e-7]:
        covariance = basis.T @ (rho**distance) @ basis / (1 - rho)
        radius = 1 / float(s @ np.linalg.solve(covariance, s))
        assert radius <= radius_upper
        weights = np.linalg.eigvalsh(covariance) / radius - 1
        assert weights[-4] >= lower_weight
        actual_det = math.prod(1 + 2 * float(tail.TILT) * a for a in weights)
        assert actual_det >= float(determinant)
        probability = tail.approximate_tail(case['values'], split, rho)['probability']
        assert probability < float(F(certificate['p_upper']))


def test_independent_witness_and_positive_endpoint_limit_also_certify():
    case, _ = direction.info.load_example()
    model = direction.state(case['values'], case['position'])
    assert tail.certify_interval(model, F(4095, 4096), 1)['status'] == 'certified_excluded'
    assert tail.certify_interval(model, 1, 1)['status'] == 'certified_excluded'


def test_unhelpful_intervals_and_constant_direction_return_no_certificate():
    case, _ = direction.load_positive_example()
    model = direction.state(case['values'], case['position'])
    certificate = tail.certify_interval(model, F(7, 10), F(9, 10))
    assert certificate['status'] == 'not_certified' and certificate['p_upper'] == '1'
    short = direction.state([1, 3, 2], 1)
    assert tail.certify_interval(short, F(999, 1000), 1)['status'] == 'not_certified'
    assert tail.approximate_tail([1, 3, 2], 1, 0.7)['probability'] == 1
    assert tail.approximate_tail(case['values'], case['position'], 0)['probability'] == 1


def test_interval_certificates_preserve_levels_scale_and_reversal():
    y = [F(10 + (i >= 5) + (-1) ** i) for i in range(16)]
    models = [
        direction.state(y, 5),
        direction.state(y[::-1], 11),
        direction.state([3 * x + (10 if i < 5 else 25) for i, x in enumerate(y)], 5),
    ]
    results = [tail.certify_interval(model, F(999, 1000), 1) for model in models]
    assert results[0] == results[1] == results[2]


def test_invalid_domains_are_rejected():
    model = direction.state([1, 3, 2], 1)
    with pytest.raises(ValueError, match='correlation interval'):
        tail.certify_interval(model, -1, 1)
    with pytest.raises(ValueError, match='correlation interval'):
        tail.certify_interval(model, 0, 1, tilt=F(1, 2))
    with pytest.raises(ValueError, match='exponential tilt'):
        tail.tilted_bound_squared([-1, 2], F(1, 2))
