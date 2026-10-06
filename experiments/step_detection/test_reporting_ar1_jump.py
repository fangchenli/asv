"""Check proper indexed densities, prior equivalence, and full certificates."""

import math
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import multivariate_t

from experiments.step_detection import ar1_information as info
from experiments.step_detection import jump_prior_diagnosis as diagnosis
from experiments.step_detection import rational_polynomial as p
from experiments.step_detection import reporting_ar1 as ar1
from experiments.step_detection import reporting_ar1_full as full
from experiments.step_detection import reporting_ar1_jump as jump


def values():
    return np.array([10 + 0.8 * (i >= 5) + 0.1 * math.sin(i * 1.7) for i in range(12)])


@pytest.mark.parametrize('split', [1, 5, 11])
@pytest.mark.parametrize('beta', [0.002, 0.5, 32])
def test_reparameterization_alone_recovers_original_density(split, beta):
    y = values()
    old_jump_precision = full.PRIOR['mean_precision'] / 2
    assert jump.component_log_density(y, split, beta, old_jump_precision) == pytest.approx(
        full.component_log_density(y, split, beta), abs=1e-10
    )


@pytest.mark.parametrize('split', [1, 5, 11])
@pytest.mark.parametrize('jump_precision', [1, 1 / 64, 1 / 4096])
def test_new_component_matches_independent_student_density(split, jump_precision):
    y = values()
    design = np.column_stack((np.ones(len(y)), (np.arange(len(y)) >= split) - 0.5))
    prior_variance = np.diag([1 / jump.PRIOR['midpoint_precision'], 1 / jump_precision])
    alpha, beta = full.PRIOR['shape'], 0.002
    shape = beta / alpha * (np.eye(len(y)) + design @ prior_variance @ design.T)
    expected = multivariate_t.logpdf(y, loc=np.zeros(len(y)), shape=shape, df=2 * alpha)
    assert jump.component_log_density(y, split, beta, jump_precision) == pytest.approx(
        expected, abs=1e-8
    )


@pytest.mark.parametrize('split', [1, 5, 11])
def test_indexed_density_is_a_fixed_mixture_and_reversal_maps_location(split):
    y = values()
    old = full.predictive_log_density(y)
    alternatives = [
        jump.component_log_density(y, split, 0.5 * 2.0 ** (2 * j), 1 / scale**2)
        for j in full.PRIOR['scale_exponents']
        for scale in jump.JUMP_SCALES
    ]
    proper_alternative = float(logsumexp(alternatives) - math.log(len(alternatives)))
    expected = float(logsumexp([old, proper_alternative]) - math.log(2))
    result = jump.predictive_log_density(y, split)
    assert result == pytest.approx(expected)
    assert result >= old - math.log(2)
    assert result == pytest.approx(jump.predictive_log_density(y[::-1], len(y) - split))


def test_global_control_and_indexed_densities_keep_original_weight():
    for y in [values(), values() * 10, np.zeros(12)]:
        old = full.predictive_log_density(y)
        assert jump.global_log_density(y) >= old - math.log(2) - 1e-12
        for split in [1, 6, 11]:
            assert jump.predictive_log_density(y, split) >= old - math.log(2) - 1e-12


def test_cost_decomposition_is_an_exact_accounting_identity():
    case, _ = info.load_example()
    result = diagnosis.cost_decomposition(case['values'], case['position'])
    assert result['sum'] == pytest.approx(result['total_gap'], abs=1e-10)
    assert all(value > 0 for value in result['parts'].values())
    assert result['parts']['location_mixture'] == pytest.approx(
        math.log(2 * (len(case['values']) - 1))
    )


def test_indexed_engine_records_each_candidates_density_and_preserves_units():
    y = values()
    result = jump.evidence(y, unit=1, max_cells=1)
    assert result['status'] == 'unresolved'
    assert not result['has_alert']
    for split, item in result['confidence']['by_split'].items():
        assert item['log_predictive_density'] == pytest.approx(
            jump.predictive_log_density(y, int(split))
        )
        assert F(item['coefficient']) == full.confidence_coefficient(
            item['log_predictive_density'], len(y), result['calibration']['confidence_alpha']
        )
    scaled = jump.evidence(y * 1024, unit=1024, max_cells=1)
    assert scaled['confidence']['by_split'] == result['confidence']['by_split']
    assert scaled['certificate'] == result['certificate']


def test_saved_slowdown_alert_has_complete_reconstructible_certificate():
    case, _ = info.load_example()
    y = case['values']
    result = jump.evidence(y, unit=1)
    assert result['status'] == 'certified_alert'
    assert result['witness'] is None
    assert set(result['confidence']['by_split']) == {str(t) for t in range(1, len(y))}
    models = ar1.Models(y)
    coverage = {split: [] for split in range(1, len(y))}
    for cell in result['certificate']:
        split, left, right = cell['split'], F(cell['left']), F(cell['right'])
        coverage[split].append((left, right))
        if cell['route'] == 'confidence':
            coefficient = F(result['confidence']['by_split'][str(split)]['coefficient'])
            numerator, denominator = full.residual_ratio(models, split)
            assert full.excluded_on(numerator, denominator, coefficient, len(y), left, right)
        else:
            polynomials = (
                models.size(split, result['calibration'])
                if cell['route'] == 'size'
                else (models.shape(split, cell['extra'], result['calibration']),)
            )
            assert all(p.positive_on(poly, left, right) for poly in polynomials)
    for intervals in coverage.values():
        intervals.sort()
        assert intervals[0][0] == -1 and intervals[-1][1] == 1
        assert all(a[1] == b[0] for a, b in zip(intervals[:-1], intervals[1:], strict=True))


def test_degenerate_data_and_invalid_units_remain_conservative():
    assert jump.evidence([10.0] * 8, unit=1)['status'] == 'insufficient_variation'
    with pytest.raises(ValueError, match='external unit'):
        jump.evidence(values(), unit=0)
