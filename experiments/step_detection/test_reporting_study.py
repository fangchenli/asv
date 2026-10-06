"""Check reporting-study labels and the unchanged fitting/reporting pipeline."""

import math

import pytest

from experiments.step_detection import ablation_study as a
from experiments.step_detection import harness as h
from experiments.step_detection import reporting_study as r


def test_case_count_labels_and_distinct_streams():
    # Check design counts without generating any future evaluation histories.
    count = math.prod(len(values) for values in r.DESIGN.values())
    assert count == 2400
    assert count * sum(c > 0.05 for c in r.DESIGN['changes']) // len(r.DESIGN['changes']) == 960
    assert count // len(r.DESIGN['changes']) == 480
    cases = [
        r.make_case('gaussian', 40, change, 0.02, loc, 12345)
        for change in (0, 0.05)
        for loc in ('middle', 'recent')
    ]
    assert all(not c['positive'] for c in cases)
    assert len({tuple(c['values']) for c in cases}) == len(cases)


@pytest.mark.parametrize('backend', ['python', 'native'])
def test_same_detector_and_reporting_with_evidence_gate(backend):
    if backend == 'native':
        h.load_detector(backend)
    case = r.make_case('gaussian', 40, 0.08, 0.02, 'middle', 12345)
    cfg = {'name': 'shared-r1-c2', 'floor': 'shared', 'cap': 1, 'penalty': 2}
    # A permissive confidence set suffices to test wiring; it makes no alert.
    frozen = {
        'detector': cfg,
        'calibrations': {'40': {'n': 40, 'alpha': '1/20', 'critical_score': 10000}},
    }
    actual = r.evaluate(case, frozen, backend)
    expected = a.evaluate(case, [cfg], backend)
    assert actual['shared_report'] == expected['methods'][cfg['name']]
    assert actual['methods']['production'] == expected['methods']['production']['has_alert']
    for name in ('reference', 'joint'):
        assert actual['methods'][name + '_gate'] == (
            actual['methods']['shared_existing'] and actual[name]['has_alert']
        )


def test_infinity_serialization_is_explicit():
    assert r.json_safe({'a': [math.inf, -math.inf, None]}) == {
        'a': ['infinity', '-infinity', None]
    }
