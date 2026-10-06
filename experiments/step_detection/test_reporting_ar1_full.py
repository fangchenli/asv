"""Proper predictive densities, rational enclosures, and reporting certificates."""

import math
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import multivariate_t

from experiments.step_detection import ar1_information as info
from experiments.step_detection import rational_polynomial as p
from experiments.step_detection import reporting_ar1 as old
from experiments.step_detection import reporting_ar1_full as full
from experiments.step_detection import reporting_covariance as gls


def strong_values():
    return [10 + 2 * (i >= 50) + 0.1 * math.sin(i * 1.7) for i in range(100)]


@pytest.mark.parametrize('split', [None, 1, 4, 7])
@pytest.mark.parametrize('beta', [0.002, 0.5, 32.0])
def test_integrated_density_matches_multivariate_student_distribution(split, beta):
    y = np.array(strong_values()[:8])
    design = (
        np.ones((8, 1))
        if split is None
        else np.column_stack((np.arange(8) < split, np.arange(8) >= split)).astype(float)
    )
    alpha, kappa = full.PRIOR['shape'], full.PRIOR['mean_precision']
    shape = beta / alpha * (np.eye(8) + design @ design.T / kappa)
    expected = multivariate_t.logpdf(y, loc=np.zeros(8), shape=shape, df=2 * alpha)
    assert full.component_log_density(y, split, beta) == pytest.approx(expected, abs=1e-8)


def test_mixture_weights_and_reversal():
    y = strong_values()[:8]
    terms = []
    for split in [None, *range(1, len(y))]:
        weight = 0.5 if split is None else 0.5 / (len(y) - 1)
        for j in full.PRIOR['scale_exponents']:
            terms.append(
                full.component_log_density(y, split, 0.5 * 2.0 ** (2 * j)) + math.log(weight / 25)
            )
    assert full.predictive_log_density(y) == pytest.approx(logsumexp(terms))
    assert full.predictive_log_density(y[::-1]) == pytest.approx(full.predictive_log_density(y))


@pytest.mark.parametrize('split', [1, 5, 11])
def test_common_endpoint_cancellation_preserves_profile_likelihood(split):
    y = strong_values()[:12]
    models = old.Models(y)
    num, den = full.residual_ratio(models, split)
    d, n, _, _ = models.fit((split,))
    for rho in [F(-7, 8), F(0), F(7, 8), F(4095, 4096)]:
        ratio = p.evaluate(num, rho) / p.evaluate(den, rho)
        assert ratio == p.evaluate(n, rho) / p.evaluate(d, rho)
        assert float(ratio) == pytest.approx(info.full_profile(y, split, float(rho))['sse'])


def test_confidence_interval_bound_uses_stationary_factor_at_endpoint():
    num, den = (F(1),), (F(1),)
    assert full.excluded_on(num, den, F(1, 5), 8, F(15, 16), F(1))
    assert not full.excluded_on(num, den, F(1, 5), 8, F(-1), F(1))
    for i in range(15, 17):
        assert full.excluded_at(num, den, F(1, 5), 8, F(i, 16))


def test_midpoint_grid_cannot_certify_a_narrow_surviving_region():
    num = p.polynomial([F(1, 64) + F(1, 1000), F(-1, 4), 1])
    den = p.polynomial([1])
    assert all(full.excluded_at(num, den, F(10000), 2, F(i, 4)) for i in range(-4, 5))
    assert not full.excluded_at(num, den, F(10000), 2, F(1, 8))
    assert not full.excluded_on(num, den, F(10000), 2, F(-1), F(1))


def test_numerical_failure_disables_confidence_exclusion():
    assert full.confidence_coefficient(-math.inf, 12, 0.01) == 0
    assert full.confidence_coefficient(1e6, 12, 0.01) == 0
    assert not full.excluded_on((F(1),), (F(1),), F(0), 12, F(0), F(1))


def test_saved_witness_survives_but_closer_endpoint_is_excluded():
    case, _ = info.load_example()
    y, split = case['values'], case['position']
    log_q = full.predictive_log_density(y)
    coefficient = full.confidence_coefficient(log_q, len(y), 0.01)
    num, den = full.residual_ratio(old.Models(y), split)
    for rho in [F(4095, 4096), F(9999, 10000)]:
        log_e = log_q - info.full_profile(y, split, float(rho))['log_u']
        assert full.excluded_at(num, den, coefficient, len(y), rho) == (log_e > math.log(100))
    assert not full.excluded_at(num, den, coefficient, len(y), F(4095, 4096))
    assert full.excluded_on(num, den, coefficient, len(y), F(99999, 100000), F(1))
    rho = 4095 / 4096
    matrix = rho ** np.abs(np.arange(len(y))[:, None] - np.arange(len(y))[None, :])
    report = gls.evidence(y, matrix, old.calibration(len(y)))
    assert split not in report['size_rejected_splits']
    assert split not in report['fit_rejected_splits']


def test_full_alert_certificate_covers_continuum_and_can_be_reconstructed():
    y = strong_values()
    result = full.evidence(y, unit=1)
    assert result['status'] == 'certified_alert'
    coefficient = F(result['confidence']['coefficient'])
    models = old.Models(y)
    by_split = {i: [] for i in range(1, len(y))}
    for cell in result['certificate']:
        split, left, right = cell['split'], F(cell['left']), F(cell['right'])
        by_split[split].append((left, right))
        if cell['route'] == 'confidence':
            num, den = full.residual_ratio(models, split)
            assert full.excluded_on(num, den, coefficient, len(y), left, right)
        else:
            polynomials = (
                models.size(split, result['calibration'])
                if cell['route'] == 'size'
                else (models.shape(split, cell['extra'], result['calibration']),)
            )
            assert all(p.positive_on(poly, left, right) for poly in polynomials)
    for cells in by_split.values():
        cells.sort()
        assert cells[0][0] == -1 and cells[-1][1] == 1
        assert all(a[1] == b[0] for a, b in zip(cells[:-1], cells[1:], strict=True))


def test_units_and_work_limits():
    y = strong_values()
    reference = full.evidence(y, unit=1, max_cells=1)
    scaled = full.evidence([x * 1024 for x in y], unit=1024, max_cells=1)
    assert reference['status'] == scaled['status'] == 'unresolved'
    assert not reference['has_alert']
    assert reference['certificate'] == scaled['certificate']
    assert reference['confidence']['coefficient'] == scaled['confidence']['coefficient']
    with pytest.raises(ValueError, match='external unit'):
        full.evidence(y, unit=0)
    with pytest.raises(ValueError, match='budget'):
        full.evidence(y, unit=1, max_depth=-1)
    config = old.calibration(len(y))
    config['size_t_critical'] *= 0.9
    with pytest.raises(ValueError, match='calibration'):
        full.evidence(y, unit=1, config=config)


def test_degenerate_history_is_not_an_alert():
    assert full.evidence([10.0] * 8, unit=1)['status'] == 'insufficient_variation'
