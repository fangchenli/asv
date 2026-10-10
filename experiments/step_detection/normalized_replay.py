"""Verified, bounded full-search replay of the nine frozen unresolved histories."""

import argparse
import gzip
import hashlib
import json
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from . import directional_normalized as normalized
from . import reporting_sturm as reporting
from . import sturm_replay
from . import sturm_unresolved_diagnosis as diagnosis
from . import sturm_verification as verification

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(100_000)

HERE = Path(__file__).parent
SOURCES = (
    *sturm_replay.SOURCES,
    'directional_normalized.py',
    'determinant_interval_polynomial.py',
    'sturm_unresolved_diagnosis.py',
    'normalized_replay.py',
)
WORK = {'max_cells': 4096, 'max_depth': 17, 'max_normalized_cells': 8}


def evaluate(item):
    case, previous = item
    start = time.perf_counter()
    result = reporting.evidence(case['values'], use_normalized=True, **WORK)
    search_seconds = time.perf_counter() - start
    start = time.perf_counter()
    validation = verification.verify(case['values'], result)
    return {
        'case': case,
        'previous': previous,
        'result': result,
        'validation': validation,
        'search_seconds': search_seconds,
        'verification_seconds': time.perf_counter() - start,
    }


def replay(output, workers, case_ids=()):
    items = diagnosis.unresolved_cases()
    if case_ids:
        available = {case['id'] for case, _ in items}
        if set(case_ids) - available:
            raise ValueError('Requested case is not a saved unresolved history')
        items = [(case, old) for case, old in items if case['id'] in case_ids]
    # Get short and positive-correlation cases back first; all cases are retained.
    items.sort(
        key=lambda item: (item[0]['n'], not item[0]['id'].startswith('positive'), item[0]['id'])
    )
    output.mkdir(parents=True, exist_ok=False)
    metadata = {
        'schema_version': 1,
        'purpose': 'Post-hoc full-search replay; not a fresh sensitivity evaluation',
        'work': WORK,
        'polynomial_bits': normalized.BITS,
        'radius_bits': normalized.RADIUS_BITS,
        'input_archive_sha256': hashlib.sha256(diagnosis.INPUTS.read_bytes()).hexdigest(),
        'record_archive_sha256': hashlib.sha256(diagnosis.RECORDS.read_bytes()).hexdigest(),
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in SOURCES
        },
        'case_ids': [case['id'] for case, _ in items],
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
            result = row['result']
            print(
                f"{row['case']['id']}: unresolved -> {result['status']}; "
                f"cells={result['cells_visited']}; "
                f"normalized_calls={result['confidence']['normalized']['calls']}; "
                f"verified={row['validation']['verified']}",
                flush=True,
            )
    payload = json.dumps({**metadata, 'rows': rows}, allow_nan=False).encode()
    archive = gzip.compress(payload, mtime=0)
    archive_name = 'normalized_replay.json.gz'
    (output / archive_name).write_bytes(archive)
    summary = {
        **metadata,
        'archive': {'file': archive_name, 'sha256': hashlib.sha256(archive).hexdigest()},
        'statuses': dict(Counter(row['result']['status'] for row in rows)),
        'rows': [
            {
                'id': row['case']['id'],
                'previous': row['previous'],
                'status': row['result']['status'],
                'has_alert': row['result']['has_alert'],
                'cells_visited': row['result']['cells_visited'],
                'normalized_calls': row['result']['confidence']['normalized']['calls'],
                'witness': row['result']['witness'],
                'validation': row['validation'],
                'search_seconds': row['search_seconds'],
                'verification_seconds': row['verification_seconds'],
            }
            for row in rows
        ],
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--case-id', action='append', default=[])
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    replay(args.output, args.workers, args.case_id)


if __name__ == '__main__':
    main()
