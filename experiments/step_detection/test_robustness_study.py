"""Validate stress mechanisms and original-coordinate reporting after gaps."""

import copy
import json
import math

import pytest

from experiments.step_detection import ablation_study as a
from experiments.step_detection import harness as h
from experiments.step_detection import robustness_study as r
from experiments.step_detection import threshold_study as t


def test_new_seeds_and_case_count():
    assert not set(r.DESIGN['seeds']) & set(
        a.DESIGN['development_seeds'] + a.DESIGN['heldout_seeds']
    )
    dimensions = ['sizes', 'changes', 'noise', 'correlations', 'locations', 'seeds', 'conditions']
    assert math.prod(len(r.DESIGN[k]) for k in dimensions) == 7680


def test_variance_changes_share_control_noise():
    controls = r.make_case(40, 0, 0.02, 0.7, 'middle', 0, 'gaussian')['values']
    for condition in ('variance_up', 'variance_down'):
        actual = r.make_case(40, 0, 0.02, 0.7, 'middle', 0, condition)['values']
        expected = [
            10 + (x - 10) * (3 if (i >= 20) == (condition == 'variance_up') else 1)
            for i, x in enumerate(controls)
        ]
        assert actual == pytest.approx(expected)


def test_combined_masks_match_separate_stress_conditions():
    def make(condition):
        return r.make_case(100, 0.06, 0.02, 0.7, 'middle', 0, condition)['values']

    outliers, missing, combined = (
        make(c) for c in ('outliers', 'missing_random', 'outliers_missing')
    )
    assert any(v is None for v in missing)
    assert combined[0] is not None and combined[-1] is not None
    assert combined == [None if mask is None else value for mask, value in zip(missing, outliers)]
    control = make('gaussian')
    increments = [a - b for a, b in zip(outliers, control)]
    assert any(delta > 0 for delta in increments)
    assert all(delta == pytest.approx(0) or delta == pytest.approx(2) for delta in increments)


@pytest.mark.parametrize('condition', r.CONDITIONS)
def test_prepared_fit_reporting_matches_public_production(condition):
    module = h.load_detector()
    case = r.make_case(40, 0.06, 0.02, 0.7, 'recent', 0, condition)
    expected = module.detect_steps(case['values'], case['weights'])
    y, _, indices = h.prepare(module, case)
    production, _ = t.candidate_pool(module, y, 'python')
    result = r.assess(module, case, y, indices, production)
    assert result['boundaries'] == [s[0] for s in expected[1:]]
    assert result['alerts'] == [
        alert[1] for alert in (module.detect_regressions(expected, threshold=0.05)[2] or [])
    ]


def test_gap_truth_maps_to_first_usable_observation():
    module = h.load_detector()
    case = r.make_case(40, 0.1, 0, 0, 'recent', 0, 'missing_gap')
    y, _, indices = h.prepare(module, case)
    assert case['position'] == 36 and indices[-2:] == [38, 39]
    fit = {'right': [34, 36], 'values': [10, 11], 'costs': [0, 0]}
    result = r.assess(module, case, y, indices, fit)
    assert result['boundaries'] == result['alerts'] == [38]
    assert result['boundary_metrics']['pairs'] == [[38, 38]]
    assert result['alert_metrics']['matched'] == 1


def test_inherited_settings_and_freeze_guard():
    try:
        h.load_detector('native')
    except RuntimeError:
        pytest.skip('Archived calibration pins the native extension')
    calibration = json.loads((r.HERE / 'data' / 'ablation_v1_frozen.json').read_text())
    if calibration['source_hashes'] != a.source_hashes('native'):
        pytest.skip(
            'Archived calibration pins another build; generate a new calibration to replay it'
        )
    frozen = r.freeze(calibration, 'native')
    assert frozen['selections'] == [
        c for c in calibration['selections'] if c['family'] in ('current-r1', 'shared-r1')
    ]
    r.check_frozen(frozen, 'native')
    changed = copy.deepcopy(frozen)
    changed['design']['seeds'] = [0]
    with pytest.raises(ValueError, match='Frozen'):
        r.check_frozen(changed, 'native')


def test_missing_case_backend_agreement():
    try:
        h.load_detector('native')
    except RuntimeError:
        pytest.skip('Build the C++ extension to compare backends')
    case = r.make_case(40, 0.06, 0.02, 0.7, 'recent', 0, 'outliers_missing')
    configs = [c for c in a.configurations() if c['name'] == 'shared-r1-c2']
    assert (
        r.evaluate(case, configs, 'native')['methods']
        == r.evaluate(case, configs, 'python')['methods']
    )
