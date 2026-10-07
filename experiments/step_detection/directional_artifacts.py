"""Verify and archive a completed directional study without changing its decisions."""

import argparse
import gzip
import hashlib
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from . import directional_diagnosis as diagnosis
from . import directional_study as study
from . import indexed_artifacts as previous


def verify_case(pair):
    case, row = pair
    counts = Counter(previous.verify_case(pair))
    for method in study.NEW_METHODS:
        result = row[method]
        assert result['confidence']['tail_enabled'] == (method == 'directional_tail')
        assert result['has_alert'] == row['methods'][method + '_alone']
        assert row['methods'][method + '_gate'] == (
            result['has_alert'] and row['methods']['shared_existing']
        )
        assert diagnosis.verify(case['values'], result)['verified']
        assert row[method + '_diagnostics'] == study.true_pair_diagnostic(
            method, case, row['oracle_reporting'], result['calibration']
        )
        counts[method + ':' + result['status']] += 1
    return dict(counts)


def verifier_hash():
    hashes = study.source_hashes('native')
    return hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()


def verify_pairs(pairs, workers, cache_path=None):
    """Reuse checks only when source and complete input/record hashes match."""
    pairs = list(pairs)
    keys = [
        hashlib.sha256(
            json.dumps(pair, sort_keys=True, separators=(',', ':')).encode()
        ).hexdigest()
        for pair in pairs
    ]
    entries = {}
    if cache_path is not None and cache_path.exists():
        cached = json.loads(cache_path.read_text())
        if cached['verifier_sha256'] == verifier_hash():
            entries = cached['entries']
    missing = [(key, pair) for key, pair in zip(keys, pairs, strict=True) if key not in entries]

    def save_cache():
        if cache_path is not None:
            temporary = cache_path.with_suffix('.tmp')
            temporary.write_text(
                json.dumps({'verifier_sha256': verifier_hash(), 'entries': entries})
            )
            temporary.replace(cache_path)

    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, ((key, _), checked) in enumerate(
            zip(missing, pool.map(verify_case, (pair for _, pair in missing)), strict=True), 1
        ):
            entries[key] = checked
            if i % 10 == 0:
                save_cache()
                print(
                    f'Verified {i}/{len(missing)} additional histories and their six references',
                    flush=True,
                )
    save_cache()
    counts = Counter()
    for key in keys:
        counts.update(entries[key])
    print(
        f'All {len(pairs)} supplied histories verified ({len(missing)} newly checked).', flush=True
    )
    return counts


def archive(run, frozen_path, destination, workers, cache_path=None):
    frozen = json.loads(frozen_path.read_text())
    study.check_frozen(frozen, frozen['backend'])
    inputs = [json.loads(line) for line in (run / 'inputs.jsonl').read_text().splitlines()]
    rows = [json.loads(line) for line in (run / 'records.jsonl').read_text().splitlines()]
    assert inputs == list(study.cases()) and len(rows) == len(inputs) == 360
    assert set(Counter(case['pair_id'] for case in inputs).values()) == {3}
    assert (
        json.loads((run / 'covariance_shapes.json').read_text())
        == study.previous.previous.shapes()
    )
    manifest = json.loads((run / 'manifest.json').read_text())
    assert manifest['frozen_sha256'] == hashlib.sha256(frozen_path.read_bytes()).hexdigest()
    assert manifest['source_hashes'] == study.source_hashes(frozen['backend'])
    assert manifest['shape_hashes'] == study.previous.previous.shape_hashes()
    expected = {
        'summary': study.reporting.summarize(rows),
        'comparisons': study.comparisons(rows),
        'diagnostics': study.diagnostics(rows),
        'breakdown': {
            key: {
                str(value): study.reporting.summarize(
                    [row for row in rows if row['case'][key] == value]
                )
                for value in sorted({row['case'][key] for row in rows})
            }
            for key in ('condition', 'n', 'change', 'noise', 'location')
        },
    }
    for name, content in expected.items():
        assert content == json.loads((run / (name + '.json')).read_text())
    counts = verify_pairs(zip(inputs, rows, strict=True), workers, cache_path)
    destination.mkdir(parents=True, exist_ok=True)
    archives = []

    def save(name, content, metadata):
        target = destination / name
        if target.exists():
            raise FileExistsError(target)
        compressed = gzip.compress(content, mtime=0)
        assert len(compressed) < 500 * 1024
        target.write_bytes(compressed)
        assert gzip.decompress(target.read_bytes()) == content
        archives.append(
            {
                'file': name,
                'sha256': hashlib.sha256(compressed).hexdigest(),
                'bytes': len(compressed),
                **metadata,
            }
        )

    for condition in study.CONDITIONS:
        for n in study.DESIGN['sizes']:
            for kind, items in [('inputs', inputs), ('records', rows)]:
                selected = [
                    item
                    for item in items
                    if (item if kind == 'inputs' else item['case'])['condition'] == condition
                    and (item if kind == 'inputs' else item['case'])['n'] == n
                ]
                pending = [selected[i : i + 10] for i in range(0, len(selected), 10)]
                shard = 0
                while pending:
                    group = pending.pop(0)
                    content = ''.join(
                        json.dumps(item, separators=(',', ':'), sort_keys=True, allow_nan=False)
                        + '\n'
                        for item in group
                    ).encode()
                    if len(gzip.compress(content, mtime=0)) >= 500 * 1024:
                        assert len(group) > 1
                        middle = len(group) // 2
                        pending[:0] = [group[:middle], group[middle:]]
                        continue
                    shard += 1
                    save(
                        f'directional_v1_{condition}_n{n}_{kind}_{shard}.jsonl.gz',
                        content,
                        {
                            'kind': kind,
                            'condition': condition,
                            'n': n,
                            'shard': shard,
                            'histories': len(group),
                        },
                    )
    save(
        'directional_v1_shapes.json.gz',
        (run / 'covariance_shapes.json').read_bytes(),
        {'kind': 'covariance_shapes', 'matrices': 6},
    )
    index = {
        'freeze_commit': manifest['git_revision'],
        'frozen_file': frozen_path.name,
        'frozen_sha256': manifest['frozen_sha256'],
        'manifest': manifest,
        **expected,
        'condition_size_location': {
            f'{condition}:n{n}:{location}': study.reporting.summarize(
                [
                    row
                    for row in rows
                    if row['case']['condition'] == condition
                    and row['case']['n'] == n
                    and row['case']['location'] == location
                ]
            )
            for condition in study.CONDITIONS
            for n in study.DESIGN['sizes']
            for location in study.DESIGN['locations']
        },
        'archives': archives,
        'validation': {
            'verifier_sha256': verifier_hash(),
            'inputs_regenerated_exactly': 360,
            'base_histories': 120,
            'covariance_matrices': 6,
            'summaries_reproduced': True,
            'archive_roundtrip': True,
            'certificates_and_witnesses': dict(counts),
        },
    }
    target = destination / 'directional_v1_results.json'
    if target.exists():
        raise FileExistsError(target)
    text = json.dumps(index, indent=2, allow_nan=False) + '\n'
    assert len(text.encode()) < 500 * 1024
    target.write_text(text)
    print(f'Archived {len(archives)} files and verified all {len(rows)} histories.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--verify-available', action='store_true')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    if args.verify_available:
        frozen = json.loads(args.frozen.read_text())
        study.check_frozen(frozen, frozen['backend'])

        def complete_lines(name):
            text = (args.run / name).read_text()
            lines = text.splitlines()
            if text and not text.endswith('\n'):
                lines.pop()
            return [json.loads(line) for line in lines]

        rows = complete_lines('records.jsonl')
        inputs = complete_lines('inputs.jsonl')[: len(rows)]
        verify_pairs(zip(inputs, rows, strict=True), args.workers, args.cache)
    else:
        archive(args.run, args.frozen, args.destination, args.workers, args.cache)
