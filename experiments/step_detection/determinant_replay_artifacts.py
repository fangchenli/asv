"""Validate and archive the determinant development replay and new witness."""

import argparse
import gzip
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from . import determinant_verification as verification
from . import directional_tail as tail

HERE = Path(__file__).parent


def archive(run, output, summary_path):
    if output.exists() or summary_path.exists():
        raise FileExistsError('Archive or summary already exists')
    result = json.loads((run / 'diagnosis.json').read_text())
    for name, digest in result['source_hashes'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Replay source changed: {name}')
    rows = []
    for row in result['rows']:
        checked = verification.verify(row['case']['values'], row['result'])
        assert checked == row['validation']
        current, previous = row['result'], row['previous']
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
            item['new_determinant_p_upper'] = float(
                F(witness['confidence_bounds']['determinant']['p_upper'])
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
