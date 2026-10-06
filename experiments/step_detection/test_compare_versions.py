"""Validate paired counts and that untimed diagnostics preserve both detectors."""

import copy
import math
import subprocess

import pytest

from experiments.step_detection import compare_versions as c
from experiments.step_detection import harness as h
from experiments.step_detection.inspect_comparison import score_candidate


@pytest.fixture(params=[c.BEFORE, c.AFTER])
def source(request, tmp_path):
    path = tmp_path / 'detector.py'
    path.write_bytes(
        subprocess.check_output(['git', 'show', f'{request.param}:asv/step_detect.py'], cwd=h.ROOT)
    )
    return path


@pytest.fixture(params=['python', 'native'])
def backend(request):
    if request.param == 'native':
        try:
            h.load_detector('native')
        except RuntimeError:
            pytest.skip('Build the C++ extension to check the native backend')
    return request.param


@pytest.mark.parametrize('scenario', ['weak', 'missing', 'correlated'])
def test_diagnostics_preserve_selected_fit(source, backend, scenario):
    module = c.load_source(source, backend)
    case = h.make_case(scenario, 40, 3)
    y, w, indices = h.prepare(module, case)
    expected = module.solve_potts_autogamma(y, w)
    original_search, original_fit = module.golden_search, module.solve_potts_approx
    diagnostic = c.diagnose(module, case, backend)
    assert module.golden_search is original_search
    assert module.solve_potts_approx is original_fit
    selected = diagnostic['selected']
    assert (
        selected['right'],
        selected['values'],
        selected['costs'],
        selected['gamma'],
    ) == expected
    measured = c.evaluate(module, case)
    assert measured['boundaries'] == [indices[r] for r in expected[0][:-1]]
    assert diagnostic['fixed_gamma_gap'] >= -1e-9
    assert diagnostic['same_count_error_gap'] >= -1e-9
    assert score_candidate(module, y, w, selected) == selected['score']


def test_forced_candidate_uses_the_revision_specific_correlation_score(tmp_path):
    scores = []
    for revision in (c.BEFORE, c.AFTER):
        path = tmp_path / f'{revision}.py'
        path.write_bytes(
            subprocess.check_output(['git', 'show', f'{revision}:asv/step_detect.py'], cwd=h.ROOT)
        )
        module = c.load_source(path, 'python')
        scores.append(
            score_candidate(
                module,
                [-1, 0, 1],
                [1, 1, 1],
                {'right': [3], 'values': [0], 'costs': [2]},
            )
        )
    assert scores[0] - scores[1] == pytest.approx(math.log(3 / 2))


def test_repeats_count_accuracy_once_and_pair_median_times():
    case = h.make_case('step', 40, 0)
    result = c.evaluate(h.load_detector(), case)
    records = []
    for version, times in [('before', [1, 3, 100]), ('after', [2, 4, 5])]:
        for repeat, seconds in enumerate(times):
            records.append(
                dict(
                    result,
                    case_id=case['id'],
                    corpus='test',
                    backend='python',
                    version=version,
                    repeat=repeat,
                    seconds=seconds,
                )
            )
    summary = c.summarize(records)
    group = summary['groups'][0]
    assert group['histories'] == 1
    assert group['before']['alert_metrics']['matched'] == 1
    assert group['after']['alert_metrics']['matched'] == 1
    assert group['median_paired_runtime_ratio'] == 4 / 3
    assert group['changed_fits'] == 0
    assert group['changed_alert_locations'] == 0
    broken = copy.deepcopy(records)
    broken[-1]['steps'] = []
    with pytest.raises(ValueError, match='Nondeterministic steps'):
        c.summarize(broken)


def test_numpy_snapshot_is_unlabeled_and_replayable():
    cases = c.make_cases([], include_scaling=False)
    assert len(cases) == 6
    for case in cases:
        assert case['true_boundaries'] is None
        assert case['true_alerts'] is None
        assert len(case['values']) == len(case['revisions']) == len(case['weights'])
        assert len(case['provenance']['commit_hashes']) == len(case['values'])
        assert case['revisions'] == sorted(set(case['revisions']))
        assert case['provenance']['url'].startswith('https://pv.github.io/numpy-bench/graphs/')
