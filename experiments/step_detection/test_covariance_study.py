"""Verify oracle inputs and pipeline reuse without evaluating study seeds."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from experiments.step_detection import covariance_study as c
from experiments.step_detection import direct_stress as s


def inherited():
    return json.loads((Path(__file__).parent / 'data/direct_stress_v1_frozen.json').read_text())


@pytest.mark.parametrize('condition', s.CONDITIONS)
def test_supplied_covariance_matches_innovation_transform(condition):
    case = s.make_case(condition, 9, 0.04, 0.02, 'recent', 12345)
    rho = case['correlation']
    transform = np.zeros((9, 9))
    transform[:, 0] = [rho**i for i in range(9)]
    for i in range(9):
        for j in range(1, i + 1):
            transform[i, j] = math.sqrt(1 - rho**2) * rho ** (i - j)
    transform = np.diag(s.amplitudes(9, case['position'], case['variance_profile'])) @ transform
    np.testing.assert_allclose(c.covariance(case), transform @ transform.T, atol=1e-14)
    # Neither the mean-change size nor absolute noise scale enters the oracle.
    changed = dict(case, change=0.08, noise=0.005)
    np.testing.assert_array_equal(c.covariance(case), c.covariance(changed))


@pytest.mark.parametrize('backend', ['python', 'native'])
def test_complete_old_pipeline_unchanged(backend):
    parent = inherited()
    case = s.make_case('correlated_up', 40, 0.06, 0.02, 'recent', 12345)
    expected = s.evaluate(case, parent, backend)
    actual = c.evaluate(case, {'inherited': parent}, backend)
    oracle = actual.pop('oracle')
    route = actual.pop('oracle_generating_split_route')
    assert actual['methods'].pop('oracle_alone') == oracle['has_alert']
    assert actual['methods'].pop('oracle_gate') == (
        oracle['has_alert'] and actual['methods']['shared_existing']
    )
    assert actual == expected
    if oracle['has_alert']:
        assert route not in ('neither', 'abstain')


def test_identity_control():
    case = s.make_case('independent_constant', 40, 0.06, 0.02, 'recent', 12345)
    row = c.evaluate(case, {'inherited': inherited()}, 'python')
    assert row['oracle']['status'] == 'ok'
    assert row['oracle']['size_t'] == pytest.approx(row['direct']['size_t'])
    assert row['oracle']['lack_of_fit_f'] == pytest.approx(row['direct']['lack_of_fit_f'])
    assert row['methods']['direct_gate'] == row['methods']['oracle_gate']


def test_abstention_is_saved_as_non_alert(monkeypatch):
    case = s.make_case('correlated_up', 40, 0.06, 0.02, 'recent', 12345)
    monkeypatch.setattr(
        c.c, 'evidence', lambda *args: {'status': 'insufficient_precision', 'has_alert': False}
    )
    row = c.evaluate(case, {'inherited': inherited()}, 'python')
    assert row['oracle_generating_split_route'] == 'abstain'
    assert not row['methods']['oracle_alone']
    assert not row['methods']['oracle_gate']
    assert c.diagnostics([row])['correlated_up']['abstentions'] == 1


def test_design_counts_and_matrix_registry_without_evaluation_data():
    base = math.prod(len(c.DESIGN[k]) for k in ('sizes', 'changes', 'noise', 'locations', 'seeds'))
    assert base == 400
    assert base * len(s.CONDITIONS) == 2400
    assert set(c.DESIGN['seeds']).isdisjoint(s.DESIGN['seeds'])
    matrices = c.shapes()
    assert len(matrices) == 24
    for key, matrix in matrices.items():
        n = len(matrix)
        assert n in (40, 100)
        assert np.linalg.eigvalsh(matrix)[0] > 0
        if key.startswith('independent_constant'):
            np.testing.assert_array_equal(matrix, np.eye(n))


def test_paired_counting_tracks_both_directions():
    rows = []
    for positive in (True, False):
        for old, new in [(False, True), (True, False), (True, True), (False, False)]:
            rows.append(
                {
                    'case': {'condition': 'correlated_up', 'positive': positive},
                    'methods': {
                        'direct_alone': old,
                        'oracle_alone': new,
                        'direct_gate': old,
                        'oracle_gate': new,
                    },
                }
            )
    result = c.comparisons(rows)
    expected = {'positive_gained': 1, 'positive_lost': 1, 'null_added': 1, 'null_removed': 1}
    assert result['all']['alone'] == result['correlated_up']['gate'] == expected


def test_frozen_covariance_change_is_detected():
    frozen = c.freeze(inherited(), 'native')
    c.check_frozen(frozen, 'native')
    frozen['shape_hashes']['independent_constant-n40-middle'] = 'changed'
    with pytest.raises(ValueError, match='Frozen'):
        c.check_frozen(frozen, 'native')
