"""Replay the sturm bound on the seven saved spectral development cases."""

import argparse
import gzip
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
from pathlib import Path

from . import directional_sturm as sturm
from . import directional_tail as tail
from . import reporting_sturm as reporting
from . import residual_direction as direction
from . import sturm_verification as verification

HERE = Path(__file__).parent
SOURCES = (
    'sturm_replay.py',
    'sturm_verification.py',
    'reporting_sturm.py',
    'directional_sturm.py',
    'directional_spectral.py',
    'directional_determinant.py',
    'directional_diagnosis.py',
    'reporting_directional.py',
    'directional_tail.py',
    'residual_direction.py',
    'rational_polynomial.py',
    'reporting_ar1.py',
    'reporting_ar1_full.py',
    'reporting_covariance.py',
)


def examples():
    summary = json.loads((HERE / 'data/spectral_reporting_v1_summary.json').read_text())
    source = summary['archive']
    raw = (HERE / 'data' / source['file']).read_bytes()
    if hashlib.sha256(raw).hexdigest() != source['sha256']:
        raise ValueError('Saved replay archive hash mismatch')
    saved = json.loads(gzip.decompress(raw))
    for name, digest in saved['source_hashes'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Previous replay source changed: {name}')
    return [
        {'case': row['case'], 'previous': row['result'], 'source': source} for row in saved['rows']
    ]


def evaluate(item):
    values = item['case']['values']
    result = reporting.evidence(values, item['previous']['calibration'])
    validation = verification.verify(values, result)
    row = {**item, 'result': result, 'validation': validation}
    before = item['previous']
    if before['status'] == 'surviving_explanation':
        witness = before['witness']
        split, rho = witness['split'], F(witness['rho'])
        model = direction.state(values, split)
        width = F(1, 65536)
        row['previous_witness_bounds'] = {
            'point': sturm.certify_interval(model, rho, rho),
            'interval': sturm.certify_interval(model, rho - width, rho + width),
        }
    if result['status'] == 'surviving_explanation':
        witness = result['witness']
        row['numerical_witness_tail'] = tail.approximate_tail(
            values, witness['split'], float(F(witness['rho']))
        )
    return row


def replay(output, workers):
    output.mkdir(parents=True, exist_ok=False)
    items = examples()
    metadata = {
        'purpose': 'Verified development replay; not a fresh sensitivity evaluation',
        'work': {'max_cells': 4096, 'max_depth': 17},
        'tilts': [str(tilt) for tilt in sturm.spectral.TILTS],
        'counts': list(sturm.spectral.COUNTS),
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in SOURCES
        },
        'numerical_qualification': 'Quadrature probabilities are numerical diagnostics, not certificates',
        'case_ids': [item['case']['id'] for item in items],
    }
    (output / 'manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
    rows = [None] * len(items)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(evaluate, item): i for i, item in enumerate(items)}
        for future in as_completed(futures):
            i = futures[future]
            row = future.result()
            rows[i] = row
            (output / f'{i:02}.json').write_text(json.dumps(row, allow_nan=False) + '\n')
            print(
                f"{row['case']['id']}: {row['previous']['status']} -> {row['result']['status']} "
                f"({row['result']['cells_visited']} cells; {row['validation']['certificate_routes']})",
                flush=True,
            )
    result = {**metadata, 'rows': rows}
    (output / 'diagnosis.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    replay(args.output, args.workers)
