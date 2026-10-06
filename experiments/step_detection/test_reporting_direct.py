"""Independent least-squares checks of the direct composite-null test."""

import math
from fractions import Fraction

import numpy as np
import pytest
from scipy.stats import f, t

from experiments.step_detection.reporting_direct import critical_values, evidence, interval_moments


def test_interval_costs_against_rational_calculation():
    values = [3, 1, 4, 1, 5, 9, 2]
    means, costs = interval_moments(values)
    for left in range(len(values)):
        for right in range(left + 1, len(values) + 1):
            mean = Fraction(sum(values[left:right]), right - left)
            total = sum((Fraction(v) - mean) ** 2 for v in values[left:right])
            assert means[left, right] == pytest.approx(float(mean))
            assert costs[left, right] == pytest.approx(float(total))


@pytest.mark.parametrize('n', [8, 17, 40])
def test_every_statistic_against_independent_linear_regressions(n):
    values = np.array([10 + 0.17 * math.sin(i * 3) + 0.6 * (i >= n // 2) for i in range(n)])
    calibration = critical_values(n)
    result = evidence(values, calibration)
    for split in range(1, n):
        design = np.column_stack((np.ones(n), np.arange(n) >= split))
        coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
        residuals = values - design @ coefficients
        two = residuals @ residuals
        contrast = np.array([-0.05, 1])
        variance = two / (n - 2) * (contrast @ np.linalg.inv(design.T @ design) @ contrast)
        expected_t = contrast @ coefficients / math.sqrt(variance)
        extra_costs = []
        for extra in range(1, n):
            if extra == split:
                continue
            extended = np.column_stack((design, np.arange(n) >= extra))
            residuals = values - extended @ np.linalg.lstsq(extended, values, rcond=None)[0]
            extra_costs.append((float(residuals @ residuals), extra))
        best_cost, best_extra = min(extra_costs)
        expected_f = (two - best_cost) / (best_cost / (n - 3))
        assert result['size_t'][split - 1] == pytest.approx(expected_t, rel=1e-10, abs=1e-10)
        assert result['lack_of_fit_f'][split - 1] == pytest.approx(
            expected_f, rel=1e-10, abs=1e-10
        )
        assert result['best_extra_split'][split - 1] == best_extra
    expected_null = [
        i + 1
        for i, (ts, fs) in enumerate(zip(result['size_t'], result['lack_of_fit_f']))
        if ts <= calibration['size_t_critical'] and fs <= calibration['fit_f_critical']
    ]
    assert result['null_splits'] == expected_null
    assert result['has_alert'] == (not expected_null)


@pytest.mark.parametrize('n', [4, 40, 100])
def test_critical_values_allocate_declared_error(n):
    c = critical_values(n)
    assert t.sf(c['size_t_critical'], n - 2) == pytest.approx(0.04)
    assert (n - 2) * f.sf(c['fit_f_critical'], 1, n - 3) == pytest.approx(0.01)
    # For df=1, a squared two-sided t statistic gives the same F cutoff.
    assert t.isf(c['fit_alpha'] / (2 * (n - 2)), n - 3) ** 2 == pytest.approx(c['fit_f_critical'])


@pytest.mark.parametrize(
    'change,expected', [(0, False), (0.04, False), (0.05, False), (0.08, True)]
)
def test_noiseless_history(change, expected):
    result = evidence([10] * 36 + [10 * (1 + change)] * 4, critical_values(40))
    assert result['has_alert'] is expected
    if not expected:
        assert 36 in result['null_splits']


def test_rejecting_some_locations_is_not_enough():
    values = [10 + 0.01 * (-1) ** i for i in range(20)] + [
        10.4 + 0.01 * (-1) ** i for i in range(20)
    ]
    result = evidence(values, critical_values(40))
    assert result['fit_rejected_splits']
    assert 20 in result['null_splits']
    assert not result['has_alert']


@pytest.mark.parametrize('scale', [1e-100, 1000, 1e100])
def test_timing_units_preserve_statistics(scale):
    values = [10 + 0.01 * math.sin(i) + 0.8 * (i >= 36) for i in range(40)]
    c = critical_values(40)
    before, after = evidence(values, c), evidence([scale * x for x in values], c)
    assert before['null_splits'] == after['null_splits']
    assert before['has_alert'] == after['has_alert']
    assert after['size_t'] == pytest.approx(before['size_t'], rel=1e-9)
    assert after['lack_of_fit_f'] == pytest.approx(before['lack_of_fit_f'], rel=1e-9)


@pytest.mark.parametrize('values', [[], [10] * 3, [10] * 201, [10, 10, 10, math.nan]])
def test_bad_observations(values):
    with pytest.raises(ValueError):
        evidence(values, critical_values(4))


@pytest.mark.parametrize('kwargs', [{'alpha': 0}, {'size_share': 1}, {'threshold': -1}])
def test_bad_calibration(kwargs):
    with pytest.raises(ValueError):
        critical_values(40, **kwargs)
