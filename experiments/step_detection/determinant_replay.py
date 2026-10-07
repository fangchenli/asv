"""Verified development replay of saved histories, not a fresh evaluation."""

import argparse
import gzip
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from . import determinant_verification as verification
from . import reporting_determinant as reporting

HERE = Path(__file__).parent
SOURCES = (
    'determinant_replay.py',
    'determinant_verification.py',
    'reporting_determinant.py',
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


def read_archive(archive):
    raw = (HERE / 'data' / archive['file']).read_bytes()
    if hashlib.sha256(raw).hexdigest() != archive['sha256']:
        raise ValueError(f'Archive hash mismatch: {archive["file"]}')
    return [json.loads(line) for line in gzip.decompress(raw).splitlines()]


def examples():
    data = HERE / 'data'
    losses_path = data / 'directional_v1_loss_diagnosis.json'
    losses = json.loads(losses_path.read_text())
    index = data / 'directional_v1_results.json'
    if hashlib.sha256(index.read_bytes()).hexdigest() != losses['result_index_sha256']:
        raise ValueError('Saved result index hash mismatch')
    cases = []
    for item in losses['lost_histories']:
        inputs = read_archive(item['archives']['inputs'])
        records = read_archive(item['archives']['records'])
        key = item['case']['id']
        case = next(row for row in inputs if row['id'] == key)
        previous = next(row['directional_tail'] for row in records if row['case']['id'] == key)
        cases.append({'case': case, 'previous': previous, 'source': item['archives']})
    path = data / 'directional_reporting_v1_diagnosis.json.gz'
    raw = path.read_bytes()
    original = json.loads(gzip.decompress(raw))
    for name, digest in original['source_hashes'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Original source changed: {name}')
    for row in original['rows']:
        if row['method'] == 'directional_tail':
            cases.append(
                {
                    'case': row['case'],
                    'previous': row['result'],
                    'source': {'file': path.name, 'sha256': hashlib.sha256(raw).hexdigest()},
                }
            )
    assert len(cases) == 7 and len({item['case']['id'] for item in cases}) == 7
    return cases


def evaluate(item):
    result = reporting.evidence(item['case']['values'], item['previous']['calibration'])
    validation = verification.verify(item['case']['values'], result)
    return {**item, 'result': result, 'validation': validation}


def replay(output, workers):
    output.mkdir(parents=True, exist_ok=False)
    items = examples()
    metadata = {
        'purpose': 'Development replay on saved cases; not fresh sensitivity estimates',
        'work': {'max_cells': 4096, 'max_depth': 16},
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in SOURCES
        },
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
            result, previous = row['result'], row['previous']
            print(
                f"{row['case']['id']}: {previous['status']} -> {result['status']} "
                f"({result['cells_visited']} cells; {row['validation']['certificate_routes']})",
                flush=True,
            )
    report = {**metadata, 'rows': rows}
    (output / 'diagnosis.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    replay(args.output, args.workers)
