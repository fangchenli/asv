"""Independent algebra and enclosure checks for the full determinant bound."""

import gzip
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest

from experiments.step_detection import determinant_continuant as continuant
from experiments.step_detection import directional_determinant as bound
from experiments.step_detection import directional_tail as old
from experiments.step_detection import residual_direction as direction


def rational_determinant(matrix):
    """Small dense elimination oracle, independent of the AR(1) recurrence."""
    matrix = [list(row) for row in matrix]
    determinant = F(1)
    for i in range(len(matrix)):
        pivot = matrix[i][i]
        assert pivot > 0  # All oracle matrices here are positive definite.
        determinant *= pivot
        for j in range(i + 1, len(matrix)):
            factor = matrix[j][i] / pivot
            for k in range(i + 1, len(matrix)):
                matrix[j][k] -= factor * matrix[i][k]
    return determinant


def dense_rational_determinant(n, split, rho, a, b):
    # Difference columns span the residual space; their Gram determinant is t(n-t).
    edges = [i for i in range(n - 1) if i != split - 1]
    matrix = []
    for i in edges:
        row = []
        for j in edges:
            gram = 2 * (i == j) - (abs(i - j) == 1)
            covariance = 2 * rho ** abs(i - j) - rho ** abs(i - j - 1) - rho ** abs(i - j + 1)
            row.append(a * gram + b * covariance)
        matrix.append(row)
    return rational_determinant(matrix) / (split * (n - split))


@pytest.mark.parametrize('bits', [8, 64, 192])
def test_outward_arithmetic_contains_exact_endpoint_operations(bits):
    ar = bound.Arithmetic(bits)
    for lo, hi in [(F(-5, 7), F(3, 11)), (F(2, 3), F(9, 7)), (F(-9, 7), F(-2, 3))]:
        x = ar.interval(lo, hi)
        assert x[0] <= lo <= hi <= x[1]
        y = ar.interval(F(3, 13), F(8, 9))
        for operation, exact in [
            (ar.add, lambda a, b: a + b),
            (ar.sub, lambda a, b: a - b),
            (ar.mul, lambda a, b: a * b),
            (ar.div, lambda a, b: a / b),
        ]:
            lower, upper = operation(x, y)
            for first in (lo, (lo + hi) / 2, hi):
                for second in (F(3, 13), F(8, 9)):
                    assert lower <= exact(first, second) <= upper
        lower, upper = ar.square(x)
        assert lower <= min(lo**2, hi**2) <= max(lo**2, hi**2) <= upper
        if lo <= 0 <= hi:
            assert lower == 0
    with pytest.raises(ArithmeticError):
        ar.div(ar.interval(1), ar.interval(-1, 1))


@pytest.mark.parametrize('n,split', [(3, 1), (6, 1), (6, 5), (9, 4)])
@pytest.mark.parametrize('rho', [F(-7, 8), F(0), F(7, 8)])
def test_interval_contains_independent_exact_projected_determinants(n, split, rho):
    a, b, width = F(7, 8), F(3, 7), F(1, 4096)
    point = bound.determinant_interval(n, split, rho, rho, a, b)
    exact = dense_rational_determinant(n, split, rho, a, b)
    assert point[0] <= exact <= point[1]
    assert point[1] - point[0] < F(1, 10**40)
    interval = bound.determinant_interval(n, split, rho - width, rho + width, a, b)
    for offset in (-width, -width / 2, F(0), width / 2, width):
        exact = dense_rational_determinant(n, split, rho + offset, a, b)
        assert interval[0] <= exact <= interval[1]


@pytest.mark.parametrize('n,split', [(6, 1), (6, 3), (9, 4)])
def test_centered_continuant_bound_contains_exact_interval_values(n, split):
    a, b = F(7, 8), F(3, 7)
    left, right = F(7, 8) - F(1, 4096), F(7, 8) + F(1, 4096)
    lower = continuant.determinant_interval(n, split, left, right, a, b)
    assert lower is not None and lower > 0
    for rho in (left, (left + right) / 2, right):
        assert lower <= dense_rational_determinant(n, split, rho, a, b)


def test_de_casteljau_bisection_preserves_exact_bernstein_bounds():
    left, right = continuant._split_bernstein((F(0), F(0), F(1, 3)))
    assert left == (F(0), F(0), F(1, 12))
    assert right == (F(1, 12), F(1, 6), F(1, 3))


def saved_losses():
    root = Path(__file__).parent / 'data'
    diagnosis = json.loads((root / 'directional_v1_loss_diagnosis.json').read_text())
    for item in diagnosis['lost_histories']:
        archive = item['archives']['inputs']
        raw = (root / archive['file']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == archive['sha256']
        rows = [json.loads(line) for line in gzip.decompress(raw).splitlines()]
        values = next(row['values'] for row in rows if row['id'] == item['case']['id'])
        yield item, values


def test_archived_witness_is_excluded_over_a_nonzero_interval():
    for item, values in saved_losses():
        split = item['witness']['split']
        rho = F(item['witness']['rho'])
        model = direction.state(values, split)
        half_width = F(1, 4096)
        proof = bound.certify_interval(model, rho - half_width, rho + half_width)
        if item['case']['location'] == 'early':
            assert proof['status'] == 'certified_excluded'
            assert F(proof['p_upper']) < F(71, 10000)
            assert old.certify_interval(model, rho, rho)['p_upper_squared'] == '1'
        else:
            assert proof['status'] == 'not_certified'
        # Independent dense calculation checks radius and determinant at interior points.
        basis = direction.residual_basis(len(values), split)
        residual = basis.T @ np.array(values)
        residual /= np.linalg.norm(residual)
        distance = np.abs(np.arange(len(values))[:, None] - np.arange(len(values))[None, :])
        lower = float(F(proof['determinant_enclosure'][0]))
        radius_upper = float(F(proof['radius_upper']))
        for candidate in (rho - half_width, rho, rho + half_width):
            covariance = basis.T @ (float(candidate) ** distance) @ basis
            radius = 1 / (residual @ np.linalg.solve(covariance, residual))
            assert radius <= radius_upper
            eigenvalues = np.linalg.eigvalsh(covariance)
            determinant = np.prod(1 + 2 * float(old.TILT) * (eigenvalues / radius - 1))
            assert determinant >= lower


def test_independence_and_one_dimensional_directions_cannot_be_excluded():
    model = direction.state([1, 4, 2, 5, 3, 8], 2)
    assert bound.certify_interval(model, 0, 0)['p_upper_squared'] == '1'
    model = direction.state([1, 4, 2], 1)
    for rho in (F(-7, 8), F(0), F(7, 8)):
        assert bound.certify_interval(model, rho, rho)['p_upper_squared'] == '1'


def test_level_scale_and_reversal_invariance():
    values = [F(10 + (i >= 4) + (-1) ** i) for i in range(12)]
    models = [
        direction.state(values, 4),
        direction.state([7 * x + (3 if i < 4 else -5) for i, x in enumerate(values)], 4),
        direction.state(values[::-1], 8),
    ]
    proofs = [bound.certify_interval(model, F(-7, 8), F(-7, 8)) for model in models]
    assert len({proof['radius_upper'] for proof in proofs}) == 1
    for proof in proofs[1:]:
        assert float(F(proof['p_upper_squared'])) == pytest.approx(
            float(F(proofs[0]['p_upper_squared'])), abs=1e-14
        )


def test_wide_interval_and_low_precision_fail_conservatively():
    model = direction.state([1, 4, 2, 5, 3, 8], 2)
    proof = bound.certify_interval(model, F(-99, 100), F(99, 100))
    assert proof['status'] == 'not_certified' and proof['p_upper_squared'] == '1'
    proof = bound.certify_interval(model, F(999, 1000), F(999, 1000), bits=1)
    assert proof['status'] == 'not_certified' and proof['p_upper_squared'] == '1'


@pytest.mark.parametrize(
    'left,right,tilt', [(-1, 0, F(1, 16)), (0, 1, F(1, 16)), (1, 0, F(1, 16)), (0, 0, F(1, 2))]
)
def test_invalid_domains_raise(left, right, tilt):
    model = direction.state([1, 4, 2, 5], 2)
    with pytest.raises(ValueError):
        bound.certify_interval(model, left, right, tilt=tilt)
