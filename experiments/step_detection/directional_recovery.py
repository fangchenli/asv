"""Resume an interrupted frozen evaluation without changing its calculations.

Run with PYTHONINTMAXSTRDIGITS=0 for the trusted, generated rational proofs.
The wrapper changes scheduling and serialization limits, not the frozen study.
"""

import argparse
import hashlib
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from . import directional_study as study


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_prefix(run, inputs):
    saved_inputs = [json.loads(line) for line in (run / 'inputs.jsonl').read_text().splitlines()]
    rows = [json.loads(line) for line in (run / 'records.jsonl').read_text().splitlines()]
    if saved_inputs != inputs[: len(saved_inputs)] or len(saved_inputs) != len(rows):
        raise ValueError('Saved input prefix differs from the frozen cases or records')
    for case, row in zip(saved_inputs, rows, strict=True):
        if row['case'] != {key: value for key, value in case.items() if key != 'values'}:
            raise ValueError('Saved record metadata differs from its input')
    return rows


def resume(run, frozen_path, workers):
    if sys.get_int_max_str_digits() != 0 or os.environ.get('PYTHONINTMAXSTRDIGITS') != '0':
        raise ValueError('Run with PYTHONINTMAXSTRDIGITS=0 for large exact certificates')
    frozen = json.loads(frozen_path.read_text())
    study.check_frozen(frozen, frozen['backend'])
    manifest = json.loads((run / 'manifest.json').read_text())
    if manifest['frozen_sha256'] != digest(frozen_path) or manifest[
        'source_hashes'
    ] != study.source_hashes(frozen['backend']):
        raise ValueError('Run manifest differs from the frozen study')
    if workers != manifest['workers']:
        raise ValueError('Retain the declared worker count')
    if any(os.environ.get(key) != value for key, value in manifest['thread_environment'].items()):
        raise ValueError('Retain the declared numerical-library thread settings')
    inputs = list(study.cases())
    rows = load_prefix(run, inputs)
    recovery = {
        'reason': 'Python 4300-digit integer-to-text limit interrupted exact bound serialization',
        'python_int_max_str_digits': 0,
        'preserved_histories': len(rows),
        'preserved_inputs_sha256': digest(run / 'inputs.jsonl'),
        'preserved_records_sha256': digest(run / 'records.jsonl'),
        'original_manifest_sha256': digest(run / 'manifest.json'),
        'frozen_sha256': digest(frozen_path),
        'recovery_source_sha256': digest(Path(__file__)),
        'frozen_sources_unchanged': True,
        'workers': workers,
    }
    recovery_path = run / 'recovery.json'
    if recovery_path.exists():
        if json.loads(recovery_path.read_text()) != recovery:
            raise ValueError('Recovery provenance or saved prefix changed')
    else:
        with recovery_path.open('x') as stream:
            json.dump(recovery, stream, indent=2)
            stream.write('\n')
    checkpoint = run / 'recovery_rows'
    checkpoint.mkdir(exist_ok=True)
    remaining = []
    for i, case in enumerate(inputs[len(rows) :], len(rows)):
        path = checkpoint / f'{i:03d}.json'
        if path.exists():
            saved = json.loads(path.read_text())
            assert saved['input'] == case
        else:
            remaining.append((i, case))
    pool = ProcessPoolExecutor(max_workers=workers)
    futures = {
        pool.submit(study.evaluate, case, frozen, frozen['backend']): (i, case)
        for i, case in remaining
    }
    try:
        for done, future in enumerate(as_completed(futures), 1):
            i, case = futures[future]
            result = study.reporting.json_safe(future.result())
            path = checkpoint / f'{i:03d}.json'
            temporary = path.with_suffix('.tmp')
            temporary.write_text(
                json.dumps({'input': case, 'row': result}, allow_nan=False) + '\n'
            )
            temporary.replace(path)
            if done % 10 == 0 or done == len(remaining):
                print(
                    f'Saved {done}/{len(remaining)} remaining histories; latest {case["id"]}',
                    flush=True,
                )
    except BaseException:
        pool.terminate_workers()
        raise
    else:
        pool.shutdown()
    for i, case in enumerate(inputs[len(rows) :], len(rows)):
        saved = json.loads((checkpoint / f'{i:03d}.json').read_text())
        assert saved['input'] == case
        rows.append(saved['row'])
    for name, items in [('inputs', inputs), ('records', rows)]:
        temporary = run / (name + '.complete.jsonl')
        with temporary.open('w') as stream:
            for item in items:
                stream.write(json.dumps(item, allow_nan=False) + '\n')
        temporary.replace(run / (name + '.jsonl'))
    study.common.dump(run / 'summary.json', study.reporting.summarize(rows))
    study.common.dump(
        run / 'breakdown.json',
        {
            key: {
                str(value): study.reporting.summarize(
                    [row for row in rows if row['case'][key] == value]
                )
                for value in sorted({row['case'][key] for row in rows})
            }
            for key in ('condition', 'n', 'change', 'noise', 'location')
        },
    )
    study.common.dump(run / 'comparisons.json', study.comparisons(rows))
    study.common.dump(run / 'diagnostics.json', study.diagnostics(rows))
    print(f'Completed all {len(rows)} histories with the frozen evaluator.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    resume(args.run, args.frozen, args.workers)
