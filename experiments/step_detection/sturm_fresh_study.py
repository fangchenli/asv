"""Frozen fresh comparison of the split-aware Sturm density refinement."""

import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
import platform
import random
import statistics
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import scipy

from . import ar1_study as design_source
from . import directional_sturm as sturm
from . import reporting_ar1 as ar1
from . import reporting_sturm as reporting

HERE = Path(__file__).parent
CONDITIONS = design_source.CONDITIONS
DESIGN = {
    **design_source.DESIGN,
    'seeds': [1602, 1603],
    'stream_prefix': 'sturm-ar1-fresh-v1',
}
WORK = {'max_cells': 4096, 'max_depth': 17}
SOURCES = (
    'sturm_fresh_study.py',
    'sturm_fresh_protocol.rst',
    'reporting_sturm.py',
    'directional_sturm.py',
    'directional_spectral.py',
    'directional_determinant.py',
    'directional_tail.py',
    'residual_direction.py',
    'rational_polynomial.py',
    'reporting_ar1.py',
    'reporting_ar1_full.py',
)


def source_hashes():
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in SOURCES}


def freeze():
    return {
        'purpose': 'Independent evaluation of the rank-group Fourier refinement',
        'design': DESIGN,
        'work': WORK,
        'confidence': {
            'total_alpha': 0.05,
            'confidence_alpha': 0.01,
            'reporting_alpha': 0.04,
            'tilts': [str(x) for x in sturm.spectral.TILTS],
            'rank_groups': 'four consecutive ranks, separately bounded and geometrically combined',
        },
        'methods': {
            'control': 'same complete reporting search with split-aware Sturm route disabled',
            'candidate': 'same search with split-aware Sturm route enabled',
        },
        'source_hashes': source_hashes(),
        'calibrations': {str(n): ar1.calibration(n) for n in DESIGN['sizes']},
        'python': platform.python_version(),
        'numpy': np.__version__,
        'scipy': scipy.__version__,
    }


def check_frozen(value):
    if value != freeze():
        raise ValueError('Frozen fresh-study sources, design, settings, or environment changed')


def cases():
    for condition in CONDITIONS:
        for n, change, noise, location, seed in itertools.product(
            *(DESIGN[key] for key in ('sizes', 'changes', 'noise', 'locations', 'seeds'))
        ):
            position = {'early': n // 4, 'middle': n // 2, 'recent': n - max(4, n // 10)}[
                location
            ]
            pair_id = (
                f'{DESIGN["stream_prefix"]}-n{n}-d{change}-s{noise}-{location}-seed{seed}'
            )
            rng = random.Random(int.from_bytes(hashlib.sha256(pair_id.encode()).digest(), 'big'))
            rho = CONDITIONS[condition]
            z, values = rng.gauss(0, 1), []
            for i in range(n):
                if i:
                    z = rho * z + math.sqrt(1 - rho * rho) * rng.gauss(0, 1)
                values.append(10 * (1 + change * (i >= position)) + 10 * noise * z)
            yield {
                'id': condition + '-' + pair_id,
                'pair_id': pair_id,
                'condition': condition,
                'n': n,
                'change': change,
                'noise': noise,
                'location': location,
                'seed': seed,
                'position': position,
                'positive': change > 0.05,
                'values': values,
                'correlation': rho,
                'variance_profile': 'constant',
                'covariance_id': f'{condition}-n{n}',
            }


def evaluate(case, frozen):
    config = frozen['calibrations'][str(case['n'])]
    methods = {}
    for name, enabled in (('control', False), ('rank_group', True)):
        started = time.perf_counter()
        result = reporting.evidence(
            case['values'], config=config, use_sturm=enabled, **frozen['work']
        )
        methods[name] = {
            'status': result['status'],
            'has_alert': result['has_alert'],
            'cells_visited': result['cells_visited'],
            'seconds': time.perf_counter() - started,
            'witness': result['witness'],
        }
    return {'case': case, 'methods': methods}


def summarize(rows):
    result = {}
    keys = ('condition', 'n', 'change', 'noise', 'location', 'positive')
    groups = [('all', rows)]
    for key in keys:
        for value in sorted({row['case'][key] for row in rows}):
            groups.append((f'{key}={value}', [r for r in rows if r['case'][key] == value]))
    for label, part in groups:
        methods = {}
        for method in ('control', 'rank_group'):
            methods[method] = {
                'histories': len(part),
                'alerts': sum(row['methods'][method]['has_alert'] for row in part),
                'statuses': {
                    status: sum(row['methods'][method]['status'] == status for row in part)
                    for status in (
                        'certified_alert',
                        'surviving_explanation',
                        'unresolved',
                        'insufficient_variation',
                    )
                },
                'median_cells': statistics.median(
                    row['methods'][method]['cells_visited'] for row in part
                ),
                'median_seconds': statistics.median(
                    row['methods'][method]['seconds'] for row in part
                ),
            }
        result[label] = {
            'methods': methods,
            'paired': {
                'gained_alerts': sum(
                    r['methods']['rank_group']['has_alert']
                    and not r['methods']['control']['has_alert']
                    for r in part
                ),
                'lost_alerts': sum(
                    r['methods']['control']['has_alert']
                    and not r['methods']['rank_group']['has_alert']
                    for r in part
                ),
            },
        }
    return result


def run(frozen, output, workers):
    check_frozen(frozen)
    output.mkdir(parents=True, exist_ok=False)
    inputs = list(cases())
    if len(inputs) != 360 or sum(case['positive'] for case in inputs) != 144:
        raise AssertionError('Fresh-study design count mismatch')
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for index, row in enumerate(pool.map(partial(evaluate, frozen=frozen), inputs), start=1):
            rows.append(row)
            if index % 10 == 0:
                print(f'Evaluated {index}/360: {row["case"]["id"]}', flush=True)
    inputs_raw = ''.join(json.dumps(row['case'], allow_nan=False) + '\n' for row in rows).encode()
    records_raw = ''.join(json.dumps(row, allow_nan=False) + '\n' for row in rows).encode()
    inputs_path = output / 'inputs.jsonl.gz'
    records_path = output / 'records.jsonl.gz'
    inputs_path.write_bytes(gzip.compress(inputs_raw, mtime=0))
    records_path.write_bytes(gzip.compress(records_raw, mtime=0))
    summary = {
        'purpose': frozen['purpose'],
        'case_count': len(rows),
        'positive_count': sum(row['case']['positive'] for row in rows),
        'null_count': sum(not row['case']['positive'] for row in rows),
        'frozen_sha256': hashlib.sha256(
            json.dumps(frozen, sort_keys=True, separators=(',', ':')).encode()
        ).hexdigest(),
        'inputs': {'file': inputs_path.name, 'sha256': hashlib.sha256(inputs_path.read_bytes()).hexdigest()},
        'records': {'file': records_path.name, 'sha256': hashlib.sha256(records_path.read_bytes()).hexdigest()},
        'breakdown': summarize(rows),
        'qualification': 'Descriptive results on 120 paired base histories; correlation versions are dependent.',
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    (output / 'manifest.json').write_text(
        json.dumps(
            {
                'frozen': frozen,
                'git_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                'workers': workers,
                'thread_environment': {
                    key: os.environ.get(key)
                    for key in ('VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS')
                },
            },
            indent=2,
        )
        + '\n'
    )
    print(json.dumps({'cases': len(rows), 'summary': str(output / 'summary.json')}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('freeze', 'run'))
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    if args.phase == 'freeze':
        if args.output.exists():
            raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(freeze(), indent=2) + '\n')
    else:
        if args.frozen is None:
            parser.error('run requires --frozen')
        if args.workers < 1:
            parser.error('workers must be positive')
        run(json.loads(args.frozen.read_text()), args.output, args.workers)


if __name__ == '__main__':
    main()
