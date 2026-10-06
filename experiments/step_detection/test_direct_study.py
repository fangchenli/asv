"""Verify the fixed pipeline and fresh design without generating evaluation data."""

import json
import math
from pathlib import Path

import pytest

from experiments.step_detection import direct_study as d
from experiments.step_detection import reporting_study as r
from experiments.step_detection.reporting_direct import critical_values


def test_fresh_design_counts():
    assert math.prod(len(v) for v in d.DESIGN.values()) == 2400
    assert set(d.DESIGN['seeds']).isdisjoint(r.DESIGN['seeds'])
    assert d.DESIGN['changes'] == [0, 0.04, 0.05, 0.06, 0.08]


@pytest.mark.parametrize('backend', ['python', 'native'])
def test_old_pipeline_is_unchanged(backend):
    inherited = json.loads((Path(__file__).parent / 'data/reporting_v1_frozen.json').read_text())
    case = r.make_case('gaussian', 40, 0.08, 0.02, 'recent', 12345)
    frozen = {'inherited': inherited, 'critical_values': {'40': critical_values(40)}}
    expected = r.evaluate(case, inherited, backend)
    actual = d.evaluate(case, frozen, backend)
    for key in ('shared_fit', 'shared_report', 'reference', 'joint'):
        assert actual[key] == expected[key]
    for name, decision in expected['methods'].items():
        assert actual['methods'][name] == decision
    assert actual['methods']['direct_gate'] == (
        actual['methods']['shared_existing'] and actual['direct']['has_alert']
    )
