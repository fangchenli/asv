"""Protect saved results when resuming an interrupted exact-arithmetic study."""

import json
import os
import subprocess
import sys

import pytest

from experiments.step_detection.directional_recovery import load_prefix


@pytest.fixture
def saved(tmp_path):
    case = {'id': 'test-stream', 'values': [1, 2]}
    row = {'case': {'id': 'test-stream'}, 'decision': False}
    (tmp_path / 'inputs.jsonl').write_text(json.dumps(case) + '\n')
    (tmp_path / 'records.jsonl').write_text(json.dumps(row) + '\n')
    return tmp_path, case, row


def test_saved_results_are_preserved(saved):
    path, case, row = saved
    assert load_prefix(path, [case, {'id': 'next', 'values': [3, 4]}]) == [row]


@pytest.mark.parametrize('damage', ['input_value', 'missing_record', 'record_identity'])
def test_mismatched_inputs_or_records_are_rejected(saved, damage):
    path, case, row = saved
    if damage == 'input_value':
        changed = {**case, 'values': [1, 3]}
        (path / 'inputs.jsonl').write_text(json.dumps(changed) + '\n')
    elif damage == 'missing_record':
        (path / 'records.jsonl').write_text('')
    else:
        row['case']['id'] = 'wrong-stream'
        (path / 'records.jsonl').write_text(json.dumps(row) + '\n')
    with pytest.raises(ValueError, match='Saved'):
        load_prefix(path, [case])


def test_large_fraction_serialization_preserves_its_exact_value():
    code = (
        'from fractions import Fraction; '
        'value = Fraction(10**5000 + 1, 10**4999 + 3); '
        'encoded = str(value); '
        'assert len(encoded) > 5000 and Fraction(encoded) == value'
    )
    subprocess.run(
        [sys.executable, '-c', code],
        env={**os.environ, 'PYTHONINTMAXSTRDIGITS': '0'},
        check=True,
    )
