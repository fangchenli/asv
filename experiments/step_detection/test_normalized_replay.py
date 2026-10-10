"""Check provenance and summary consistency of the verified post-hoc replay."""

import gzip
import hashlib
import json
from collections import Counter
from fractions import Fraction as F

from experiments.step_detection import normalized_replay as replay
from experiments.step_detection import sturm_unresolved_diagnosis as diagnosis


def test_saved_replay_matches_frozen_inputs_and_full_results():
    summary = json.loads((replay.HERE / 'data/normalized_reporting_v1_summary.json').read_text())
    raw = (replay.HERE / 'data' / summary['archive']['file']).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == summary['archive']['sha256']
    archive = json.loads(gzip.decompress(raw))
    assert summary['statuses'] == {'unresolved': 8, 'certified_alert': 1}
    assert [r['id'] for r in summary['rows'] if r['has_alert']] == [
        'positive-sturm-ar1-fresh-v1-n100-d0.08-s0.005-recent-seed1603'
    ]
    assert summary['work'] == archive['work'] == replay.WORK
    assert (
        summary['input_archive_sha256']
        == hashlib.sha256(diagnosis.INPUTS.read_bytes()).hexdigest()
    )
    assert (
        summary['record_archive_sha256']
        == hashlib.sha256(diagnosis.RECORDS.read_bytes()).hexdigest()
    )
    assert summary['source_hashes'] == archive['source_hashes']
    assert set(summary['source_hashes']) == set(replay.SOURCES)
    originals = {case['id']: (case, before) for case, before in diagnosis.unresolved_cases()}
    assert len(summary['rows']) == len(archive['rows']) == len(originals) == 9
    assert summary['case_ids'] == archive['case_ids'] == [r['case']['id'] for r in archive['rows']]
    assert len(set(summary['case_ids'])) == 9
    for item, row in zip(summary['rows'], archive['rows'], strict=True):
        case_id = row['case']['id']
        assert (row['case'], row['previous']) == originals[case_id]
        result, checked = row['result'], row['validation']
        assert item['id'] == case_id and item['previous'] == row['previous']
        for key in ('status', 'has_alert', 'cells_visited', 'witness'):
            assert item[key] == result[key]
        assert item['validation'] == checked and checked['verified'] is True
        assert checked['complete_coverage'] == result['has_alert']
        assert result['cells_visited'] <= replay.WORK['max_cells']
        config = result['confidence']['normalized']
        assert config['bits'] == summary['polynomial_bits'] == 128
        assert config['radius_bits'] == summary['radius_bits'] == 192
        assert item['normalized_calls'] == config['calls'] <= config['max_cells'] == 8
        routes = Counter(c['route'] for c in result['certificate'])
        assert dict(routes) == checked['certificate_routes']
        assert 0 <= config['calls'] - routes['direction_normalized'] <= 1
        if result['has_alert']:
            for split in range(1, row['case']['n']):
                intervals = sorted(
                    (F(c['left']), F(c['right']))
                    for c in result['certificate']
                    if c['split'] == split
                )
                assert intervals and intervals[0][0] == -1 and intervals[-1][1] == 1
                assert all(left < right for left, right in intervals)
                assert all(
                    a[1] == b[0] for a, b in zip(intervals[:-1], intervals[1:], strict=True)
                )
        for cell in result['certificate']:
            if cell['route'] != 'direction_normalized':
                continue
            proof = cell['proof']
            assert proof['stationary_normalized'] and proof['status'] == 'certified_excluded'
            assert proof['left'] == cell['left'] and proof['right'] == cell['right']
            assert F(proof['p_upper_squared']) < F(result['confidence']['delta']) ** 2
            assert proof['delta'] == result['confidence']['delta']
    assert summary['statuses'] == dict(Counter(r['result']['status'] for r in archive['rows']))
