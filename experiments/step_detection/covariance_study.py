"""Evaluate the frozen known-covariance reference on fresh paired histories."""

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import scipy

from . import direct_stress as s
from . import harness as h
from . import reporting_covariance as c
from . import reporting_study as r
from . import threshold_study as t

HERE = Path(__file__).parent
DESIGN = dict(s.DESIGN, seeds=list(range(1200, 1210)))


def covariance_id(case):
    return f'{case["condition"]}-n{case["n"]}-{case["location"]}'


def covariance(case):
    positions = np.arange(case['n'])
    scales = s.amplitudes(case['n'], case['position'], case['variance_profile'])
    distances = np.abs(positions[:, None] - positions[None, :])
    return np.outer(scales, scales) * case['correlation'] ** distances


def shapes():
    result = {}
    for condition, (rho, profile) in s.CONDITIONS.items():
        for n, location in itertools.product(DESIGN['sizes'], DESIGN['locations']):
            case = {
                'condition': condition,
                'n': n,
                'location': location,
                'position': n // 2 if location == 'middle' else n - max(4, n // 10),
                'correlation': rho,
                'variance_profile': profile,
            }
            result[covariance_id(case)] = covariance(case).tolist()
    return result


def shape_hashes():
    return {
        key: hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()
        for key, value in shapes().items()
    }


def cases():
    for condition in s.CONDITIONS:
        for args in itertools.product(
            *(DESIGN[key] for key in ('sizes', 'changes', 'noise', 'locations', 'seeds'))
        ):
            case = s.make_case(condition, *args)
            case['covariance_id'] = covariance_id(case)
            yield case


def evaluate(case, frozen, backend):
    result = s.evaluate(case, frozen['inherited'], backend)
    calibration = frozen['inherited']['inherited']['critical_values'][str(case['n'])]
    oracle = c.evidence(case['values'], covariance(case), calibration)
    result['oracle'] = oracle
    route = 'abstain' if oracle['status'] != 'ok' else s.rejection_route(oracle, case['position'])
    result['oracle_generating_split_route'] = route
    result['methods']['oracle_alone'] = oracle['has_alert']
    result['methods']['oracle_gate'] = oracle['has_alert'] and result['methods']['shared_existing']
    if oracle['has_alert'] and route in ('neither', 'abstain'):
        raise AssertionError('Global oracle alert must reject the generating split')
    if case['condition'] == 'independent_constant':
        if oracle['status'] != 'ok':
            raise AssertionError('Unexpected numerical abstention in identity control')
        for key in ('has_alert', 'null_splits', 'size_rejected_splits', 'fit_rejected_splits'):
            if oracle[key] != result['direct'][key]:
                raise AssertionError('Identity-control decisions changed')
        for key in ('size_t', 'lack_of_fit_f'):
            if not all(
                math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-8)
                for a, b in zip(oracle[key], result['direct'][key], strict=True)
            ):
                raise AssertionError('Identity-control statistics changed')
    return result


def comparisons(rows):
    result = {}
    for condition in ('all', *s.CONDITIONS):
        part = [x for x in rows if condition == 'all' or x['case']['condition'] == condition]
        result[condition] = {}
        for suffix in ('alone', 'gate'):
            counts = {'positive_gained': 0, 'positive_lost': 0, 'null_added': 0, 'null_removed': 0}
            for row in part:
                old, new = (row['methods'][prefix + suffix] for prefix in ('direct_', 'oracle_'))
                if row['case']['positive']:
                    counts['positive_gained'] += new and not old
                    counts['positive_lost'] += old and not new
                else:
                    counts['null_added'] += new and not old
                    counts['null_removed'] += old and not new
            result[condition][suffix] = counts
    return result


def diagnostics(rows):
    result = {}
    for condition in s.CONDITIONS:
        part = [row for row in rows if row['case']['condition'] == condition]
        result[condition] = {
            'abstentions': sum(row['oracle']['status'] != 'ok' for row in part),
            'false_alert_routes': {
                prefix + suffix: {
                    route: sum(
                        row[key] == route
                        for row in part
                        if not row['case']['positive'] and row['methods'][prefix + suffix]
                    )
                    for route in ('size_only', 'shape_only', 'both', 'neither', 'abstain')
                }
                for prefix, key in (
                    ('direct_', 'generating_split_route'),
                    ('oracle_', 'oracle_generating_split_route'),
                )
                for suffix in ('alone', 'gate')
            },
        }
    return result


def source_hashes(backend):
    hashes = s.source_hashes(backend)
    for name in (
        'reporting_covariance.py',
        'known_covariance.rst',
        'covariance_study.py',
        'covariance_protocol.rst',
    ):
        path = HERE / name
        hashes[str(path.relative_to(h.ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def freeze(inherited, backend):
    s.check_frozen(inherited, backend)
    return {
        'design': DESIGN,
        'backend': backend,
        'source_hashes': source_hashes(backend),
        'shape_hashes': shape_hashes(),
        'inherited': inherited,
        'python': platform.python_version(),
        'numpy': np.__version__,
        'scipy': scipy.__version__,
    }


def check_frozen(frozen, backend):
    s.check_frozen(frozen['inherited'], backend)
    if (
        frozen['design'] != DESIGN
        or frozen['backend'] != backend
        or frozen['source_hashes'] != source_hashes(backend)
        or frozen['shape_hashes'] != shape_hashes()
        or frozen['numpy'] != np.__version__
        or frozen['scipy'] != scipy.__version__
        or frozen['python'] != platform.python_version()
    ):
        raise ValueError('Frozen source, design, covariance, backend, or environment changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['freeze', 'run'])
    parser.add_argument('--inherited', type=Path)
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--backend', choices=['native', 'python'], default='native')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    if args.phase == 'freeze':
        if args.inherited is None:
            parser.error('freeze requires --inherited')
        if args.output.exists():
            raise FileExistsError(args.output)
        t.dump(args.output, freeze(json.loads(args.inherited.read_text()), args.backend))
        return
    if args.frozen is None:
        parser.error('run requires --frozen')
    frozen = json.loads(args.frozen.read_text())
    check_frozen(frozen, args.backend)
    args.output.mkdir(parents=True, exist_ok=False)
    t.dump(args.output / 'covariance_shapes.json', shapes())
    t.dump(
        args.output / 'manifest.json',
        {
            'frozen_sha256': hashlib.sha256(args.frozen.read_bytes()).hexdigest(),
            'source_hashes': source_hashes(args.backend),
            'shape_hashes': shape_hashes(),
            'python': platform.python_version(),
            'numpy': np.__version__,
            'scipy': scipy.__version__,
            'git_revision': subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], text=True
            ).strip(),
            'design': DESIGN,
            'backend': args.backend,
            'workers': args.workers,
            'thread_environment': {
                key: os.environ.get(key)
                for key in ('VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS')
            },
        },
    )
    inputs = list(cases())
    rows = []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as pool,
        (args.output / 'inputs.jsonl').open('w') as input_file,
        (args.output / 'records.jsonl').open('w') as output_file,
    ):
        evaluate_case = partial(evaluate, frozen=frozen, backend=args.backend)
        for case, row in zip(inputs, pool.map(evaluate_case, inputs), strict=True):
            input_file.write(json.dumps(case, allow_nan=False) + '\n')
            output_file.write(json.dumps(r.json_safe(row), allow_nan=False) + '\n')
            input_file.flush()
            output_file.flush()
            rows.append(row)
            if len(rows) % 100 == 0:
                print(f'Evaluated {len(rows)}: {case["condition"]}', flush=True)
    t.dump(args.output / 'summary.json', r.summarize(rows))
    t.dump(
        args.output / 'breakdown.json',
        {
            key: {
                str(value): r.summarize([row for row in rows if row['case'][key] == value])
                for value in sorted({row['case'][key] for row in rows})
            }
            for key in ('condition', 'n', 'change', 'noise', 'location')
        },
    )
    t.dump(args.output / 'comparisons.json', comparisons(rows))
    t.dump(args.output / 'paired.json', s.paired(rows))
    t.dump(args.output / 'diagnostics.json', diagnostics(rows))


if __name__ == '__main__':
    main()
