"""Validate and archive the spectral development replay and new witness."""

import argparse
import gzip
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from . import directional_spectral as spectral
from . import directional_tail as tail
from . import residual_direction as direction
from . import spectral_verification as verification

HERE = Path(__file__).parent


def archive(run, output, summary_path):
    if output.exists() or summary_path.exists():
        raise FileExistsError('Archive or summary already exists')
    result = json.loads((run / 'diagnosis.json').read_text())
    for name, digest in result['source_hashes'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Replay source changed: {name}')
    source = result['rows'][0]['source']
    original_raw = (HERE / 'data' / source['file']).read_bytes()
    assert hashlib.sha256(original_raw).hexdigest() == source['sha256']
    original = json.loads(gzip.decompress(original_raw))
    by_id = {row['case']['id']: row for row in original['rows']}
    assert result['case_ids'] == [row['case']['id'] for row in result['rows']]
    assert len(result['case_ids']) == len(set(result['case_ids'])) == len(by_id)
    rows = []
    for row in result['rows']:
        before = by_id[row['case']['id']]
        assert row['source'] == source and row['case'] == before['case']
        assert row['previous'] == before['result']
        checked = verification.verify(row['case']['values'], row['result'])
        assert checked == row['validation']
        current, previous = row['result'], row['previous']
        if previous['status'] == 'surviving_explanation':
            witness = previous['witness']
            rho = F(witness['rho'])
            model = direction.state(row['case']['values'], witness['split'])
            assert row['previous_witness_bounds']['point'] == spectral.certify_interval(
                model, rho, rho
            )
            assert row['previous_witness_bounds']['interval'] == spectral.certify_interval(
                model, rho - F(1, 4096), rho + F(1, 4096)
            )
        item = {
            'id': row['case']['id'],
            'previous_status': previous['status'],
            'new_status': current['status'],
            'previous_cells': previous['cells_visited'],
            'new_cells': current['cells_visited'],
            'verification': checked,
        }
        if current['status'] == 'surviving_explanation':
            witness, before = current['witness'], previous['witness']
            item['witness'] = {key: witness[key] for key in ('split', 'rho')}
            item['previous_witness'] = {key: before[key] for key in ('split', 'rho')}
            item['new_spectral_p_upper'] = float(
                F(witness['confidence_bounds']['spectral']['p_upper'])
            )
            if item['witness'] != item['previous_witness']:
                item['numerical_tail'] = tail.approximate_tail(
                    row['case']['values'], witness['split'], float(F(witness['rho']))
                )
        rows.append(item)
    summary = {
        'purpose': result['purpose'],
        'work': result['work'],
        'verified_cases': len(rows),
        'gained_alerts': sum(
            row['result']['has_alert'] and not row['previous']['has_alert']
            for row in result['rows']
        ),
        'lost_alerts': sum(
            row['previous']['has_alert'] and not row['result']['has_alert']
            for row in result['rows']
        ),
        'numerical_qualification': 'New-witness tail is a floating quadrature estimate, not a certificate',
        'rows': rows,
        'artifact_builder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    result['summary'] = summary
    raw = (json.dumps(result, indent=2, allow_nan=False) + '\n').encode()
    compressed = gzip.compress(raw, mtime=0)
    assert gzip.decompress(compressed) == raw
    output.write_bytes(compressed)
    summary = {
        **summary,
        'archive': {'file': output.name, 'sha256': hashlib.sha256(compressed).hexdigest()},
    }
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--summary', required=True, type=Path)
    args = parser.parse_args()
    summary = archive(args.run, args.output, args.summary)
    print(
        f"Verified {summary['verified_cases']} saved cases: "
        f"{summary['gained_alerts']} gained alerts, {summary['lost_alerts']} lost alerts",
        flush=True,
    )
