"""Exact Sturm edge cases and independent checks of split-aware bounds."""

import gzip
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest

from experiments.step_detection import directional_spectral as previous
from experiments.step_detection import directional_sturm as sturm
from experiments.step_detection import residual_direction as direction


def precision(n, rho):
    matrix = np.diag([1] + [1 + rho**2] * (n - 2) + [1])
    return matrix + np.diag([-rho] * (n - 1), 1) + np.diag([-rho] * (n - 1), -1)


def fraction_precision_count_reference(n, correlation, threshold):
    previous, current = F(1), F(1) - threshold
    last_sign, count = 1, 0
    for i in range(1, n + 1):
        if current:
            sign = 1 if current > 0 else -1
            count += sign != last_sign
            last_sign = sign
        if i < n:
            diagonal = 1 if i == n - 1 else 1 + correlation**2
            previous, current = (
                current,
                (diagonal - threshold) * current - correlation**2 * previous,
            )
    return count


def test_exact_roots_and_interior_zeros_in_sturm_sequence():
    # For n=2 the eigenvalues are exactly 1-rho and 1+rho.
    for threshold, count in [(F(1, 4), 0), (F(1, 2), 0), (F(1), 1), (F(3, 2), 1), (F(2), 2)]:
        assert sturm.precision_count_below(2, F(1, 2), threshold) == count
    # n=3 has an eigenvalue exactly 1. The first and final minors vanish there.
    assert sturm.precision_count_below(3, F(1, 2), 1) == 1
    # A zero first minor at threshold 1 must not invalidate the next sign change.
    assert sturm.precision_count_below(4, F(1, 2), 1) == 2
    assert sturm.precision_count_below(9, 0, 1) == 0
    assert sturm.precision_count_below(9, 0, F(1001, 1000)) == 9


@pytest.mark.parametrize(
    'n,rho,threshold',
    [
        (2, F(1, 2), F(1, 2)),
        (4, F(1, 2), F(1)),
        (7, F(7, 9), F(5, 4)),
        (20, F(99, 100), F(3, 2)),
        (100, F(101, 128), F(1, 10)),
    ],
)
def test_integer_sturm_signs_match_fraction_recurrence(n, rho, threshold):
    assert sturm.precision_count_below(n, rho, threshold) == fraction_precision_count_reference(
        n, rho, threshold
    )


@pytest.mark.parametrize('n,rho', [(3, F(1, 4)), (8, F(3, 4)), (25, F(99, 100))])
def test_counts_and_eigenvalue_upper_bounds_match_dense_precision(n, rho):
    eigenvalues = np.linalg.eigvalsh(precision(n, float(rho)))
    for threshold in (F(1, 1000), F(2, 7), F(7, 6), F(8, 3), F(4)):
        assert sturm.precision_count_below(n, rho, threshold) == int(
            np.count_nonzero(eigenvalues < float(threshold))
        )
    for index in (1, (n + 1) // 2, n):
        upper = sturm.precision_eigenvalue_upper(n, index, rho)
        actual = eigenvalues[index - 1]
        assert float(upper) >= actual - 1e-14
        assert float(upper) - actual <= float((1 + rho) ** 2 / 2**sturm.EIGEN_BITS) + 1e-14
        assert sturm.precision_count_below(n, rho, upper) >= index


@pytest.mark.parametrize(
    'left,right', [(F(-7, 8), F(-13, 16)), (F(-1, 8), F(1, 8)), (F(13, 16), F(7, 8))]
)
def test_precision_perturbation_and_covariance_interval_enclosure(left, right):
    n, index = 20, 6
    low = F(0) if left <= 0 <= right else min(abs(left), abs(right))
    high = max(abs(left), abs(right))
    mid = (low + high) / 2
    perturbation = high**2 - mid**2 + 2 * (high - mid)
    bound = sturm.covariance_eigenvalue_lower(n, index, left, right)
    distance = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    for rho in (left, (left + right) / 2, right):
        difference = precision(n, float(abs(rho))) - precision(n, float(mid))
        assert np.linalg.eigvalsh(difference)[-1] <= float(perturbation) + 1e-14
        actual = np.linalg.eigvalsh(float(rho) ** distance)[-index]
        assert float(bound) <= actual + 1e-13
    assert bound == sturm.covariance_eigenvalue_lower(n, index, -right, -left)


@pytest.mark.parametrize('n,split,count', [(8, 1, 4), (16, 8, 12), (24, 23, 16), (100, 1, 12)])
@pytest.mark.parametrize('rho', [F(-51, 64), F(0), F(51, 64)])
@pytest.mark.parametrize('stationary_normalized', [False, True])
def test_each_subspace_bound_is_below_projected_eigenvalue(
    n, split, count, rho, stationary_normalized
):
    width = F(1, 65536)
    bounds = sturm.residual_eigenvalue_bounds(
        n, split, count, rho - width, rho + width, stationary_normalized=stationary_normalized
    )
    assert bounds == sturm.residual_eigenvalue_bounds(
        n, n - split, count, rho - width, rho + width, stationary_normalized=stationary_normalized
    )
    basis = direction.residual_basis(n, split)
    distance = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    for candidate in (rho - width, rho, rho + width):
        covariance = basis.T @ (float(candidate) ** distance) @ basis
        if stationary_normalized:
            covariance /= float(1 - candidate**2)
        actual = np.linalg.eigvalsh(covariance)[-count]
        assert all(float(value) <= actual + 1e-12 for value in bounds.values())


@pytest.mark.parametrize('stationary_normalized', [False, True])
def test_rank_groups_respect_plateau_dimension_and_dense_fourier_weights(stationary_normalized):
    n, split, left, right, radius, tilt = 30, 15, F(4, 5), F(801, 1000), F(1, 1000), F(1, 8)
    groups = sturm.rank_group_weights(
        n, split, left, right, radius, tilt, stationary_normalized=stationary_normalized
    )
    assert len(groups) == (max(split, n - split) - 1) // 4
    basis = direction.residual_basis(n, split)
    distance = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    for rho in (left, (left + right) / 2, right):
        covariance = basis.T @ (float(rho) ** distance) @ basis
        if stationary_normalized:
            covariance /= float(1 - rho**2)
        weights = np.linalg.eigvalsh(covariance)[::-1] / float(radius) - 1
        tilted = weights / (1 + 2 * float(tilt) * weights)
        for index, group in enumerate(groups):
            assert float(group) ** 4 <= np.prod(tilted[4 * index : 4 * index + 4]) + 1e-12


def saved_model():
    path = Path(__file__).parent / 'data/spectral_reporting_v1_diagnosis.json.gz'
    row = json.loads(gzip.decompress(path.read_bytes()))['rows'][0]
    return direction.state(row['case']['values'], 1)


def test_saved_witness_and_neighboring_interval_cross_cutoff():
    model, rho, width = saved_model(), F(51, 64), F(1, 65536)
    point = sturm.certify_interval(model, rho, rho)
    interval = sturm.certify_interval(model, rho - width, rho + width)
    assert previous.certify_interval(model, rho, rho)['status'] == 'not_certified'
    assert point['status'] == interval['status'] == 'certified_excluded'
    assert F(point['p_upper']) < F(957, 100000)
    assert F(interval['p_upper']) < F(964, 100000)
    assert interval['selected_tilt'] == '1/8' and interval['selected_count'] == 16
    best = next(item for item in interval['attempts'] if item['tilt'] == '1/8')
    selected = next(item for item in best['directions'] if item['count'] == 12)
    assert selected['selected_eigenvalue_bound'] == 'plateau_sturm'
    for candidate in (rho - width, rho, rho + width):
        assert F(sturm.certify_interval(model, candidate, candidate)['p_upper_squared']) <= F(
            interval['p_upper_squared']
        )


def test_grouped_fourier_density_bounds_are_conservative():
    weights = [F(3), F(2), F(1)]
    bounds = sturm.grouped_density_bounds(weights, bits=16)
    product = F(1)
    for groups, item in enumerate(bounds, start=1):
        product *= weights[groups - 1]
        root_lower = F(item['geometric_mean_lower'])
        density_upper = F(item['density_upper'])
        constant = sturm.spectral.density_constant(4 * groups)
        assert root_lower**groups <= product
        # This comparison avoids approximating the irrational exact root.
        assert density_upper**groups * product >= constant**groups
        assert density_upper <= constant / weights[groups - 1]
    with pytest.raises(ValueError):
        sturm.grouped_density_bounds([F(1), F(0)])
    with pytest.raises(ValueError):
        sturm.grouped_density_bounds([F(1)], bits=0)


def test_heterogeneous_fourier_integral_is_exact_and_handles_equal_weights():
    bounds = sturm.heterogeneous_density_bounds([F(3), F(2), F(1)])
    assert F(bounds[0]['density_upper']) == F(1, 12)
    assert F(bounds[1]['density_upper']) == F(1, 20)
    equal = sturm.heterogeneous_density_bounds([F(2), F(2), F(1)])
    assert len(equal) == 3
    assert all(F(item['group_weight_lowers'][i]) <= F(2) for i, item in enumerate(equal) if i < 2)
    for groups, item in enumerate(equal, start=1):
        weights = [F(value) for value in item['group_weight_lowers']]
        holder_upper = F(sturm.grouped_density_bounds(weights)[-1]['density_upper'])
        assert F(item['density_upper']) <= holder_upper


def test_individual_rank_groups_bound_the_corresponding_fourier_factors():
    n, split, rho = 30, 9, F(3, 4)
    radius, tilt = F(1, 2), F(1, 8)
    groups = sturm.rank_group_weights(n, split, rho, rho, radius, tilt, bits=64)
    assert groups
    length = max(split, n - split)
    for group_index, mean_lower in enumerate(groups):
        start = 4 * group_index + 1
        product = F(1)
        for rank in range(start, start + 4):
            eigen_lower = sturm.covariance_eigenvalue_lower(length, rank + 1, rho, rho)
            weight = eigen_lower / radius - 1
            assert weight > 0
            product *= weight / (1 + 2 * tilt * weight)
        assert mean_lower**4 <= product


def test_rank_group_density_sharpens_the_early_survivor_without_claiming_alert():
    path = Path(__file__).parent / 'data/sturm_reporting_v3_depth17_diagnosis.json.gz'
    archive = json.loads(gzip.decompress(path.read_bytes()))
    row = next(
        row
        for row in archive['rows']
        if row['case']['id'].endswith('early-seed1601')
        and row['result']['status'] == 'surviving_explanation'
    )
    witness = row['result']['witness']
    rho = F(witness['rho'])
    model = direction.state(row['case']['values'], witness['split'])
    old = F(witness['confidence_bounds']['sturm']['p_upper'])
    improved = sturm.certify_interval(model, rho, rho)
    assert F(1, 100) < F(improved['p_upper']) < old
    best = next(attempt for attempt in improved['attempts'] if attempt['tilt'] == '1/8')
    assert best['selected_count'] == 28
    assert len(best['rank_group_weights']) == 7


def test_grouped_density_certificate_covers_early_replay_witness():
    model, rho, width = saved_model(), F(101, 128), F(1, 65536)
    point = sturm.certify_interval(model, rho, rho)
    interval = sturm.certify_interval(model, rho - width, rho + width)
    assert point['status'] == interval['status'] == 'certified_excluded'
    assert F(point['p_upper']) < F(1, 100)
    assert F(interval['p_upper']) < F(1, 100)
    assert F(point['p_upper']) < F(9, 1000)
    assert interval['selected_tilt'] == '1/8' and interval['selected_count'] == 16
    best = next(item for item in interval['attempts'] if item['tilt'] == '1/8')
    assert len(best['grouped_density_bounds']) == 4
    assert len(best['heterogeneous_density_bounds']) == 4
    ranked = sturm.certify_interval(model, rho - width, rho + width, rank_groups_on_intervals=True)
    assert ranked['status'] == 'certified_excluded'
    ranked_best = next(item for item in ranked['attempts'] if item['tilt'] == '1/8')
    assert ranked_best['rank_group_weights']


def test_refinement_preserves_old_bound_and_invariances():
    values = [F(10 + (i >= 5) + (-1) ** i) for i in range(16)]
    changed = [7 * value + (5 if i < 5 else -3) for i, value in enumerate(values)]
    proofs = [
        sturm.certify_interval(direction.state(y, 5), F(3, 4), F(3, 4)) for y in (values, changed)
    ]
    assert proofs[0] == proofs[1]
    old = previous.certify_interval(direction.state(values, 5), F(3, 4), F(3, 4))
    assert F(proofs[0]['p_upper_squared']) <= F(old['p_upper_squared'])


def test_uninformative_and_wide_intervals_fail_conservatively():
    model = direction.state([1, 4, 2, 5, 3, 8], 2)
    assert sturm.certify_interval(model, 0, 0)['p_upper_squared'] == '1'
    assert sturm.certify_interval(model, F(-99, 100), F(99, 100))['status'] == 'not_certified'
    assert sturm.precision_eigenvalue_upper(10, 5, 0) == 1


@pytest.mark.parametrize(
    'args', [(1, 1, F(1, 2)), (4, 0, F(1, 2)), (4, 5, F(1, 2)), (4, 1, 1), (True, 1, 0)]
)
def test_invalid_eigenvalue_requests(args):
    with pytest.raises(ValueError):
        sturm.precision_eigenvalue_upper(*args)


def test_cache_does_not_bypass_argument_validation():
    sturm.precision_eigenvalue_upper(4, 1, F(1, 2))
    for args in ((4, True, F(1, 2)), (4.0, 1, F(1, 2))):
        with pytest.raises(ValueError):
            sturm.precision_eigenvalue_upper(*args)
