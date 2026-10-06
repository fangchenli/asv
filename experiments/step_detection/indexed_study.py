"""Frozen comparison of conditional, global, and indexed AR(1) confidence."""

import argparse
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
from fractions import Fraction
from functools import partial
from pathlib import Path

import numpy as np
import scipy

from . import ar1_study as previous
from . import harness as h
from . import reporting_ar1 as ar1
from . import reporting_ar1_full as full
from . import reporting_ar1_jump as jump
from . import reporting_study as reporting
from . import threshold_study as common

HERE = Path(__file__).parent
CONDITIONS = previous.CONDITIONS
DESIGN = {**previous.DESIGN, 'seeds': [1500, 1501], 'stream_prefix': 'indexed-ar1-study'}
WORK = dict(previous.WORK)
UNIT = 1.0
NEW_METHODS = ('full_history', 'indexed_original', 'indexed_jump')
PAIRS = (
    ('full_history', 'ar1'),
    ('indexed_original', 'full_history'),
    ('indexed_jump', 'full_history'),
    ('indexed_jump', 'indexed_original'),
    ('indexed_original', 'ar1'),
    ('indexed_jump', 'ar1'),
    *((name, 'oracle_reporting') for name in NEW_METHODS),
)


def make_case(condition, n, change, noise, location, seed):
    position = {'early': n // 4, 'middle': n // 2, 'recent': n - max(4, n // 10)}[location]
    pair_id = f'{DESIGN["stream_prefix"]}-n{n}-d{change}-s{noise}-{location}-seed{seed}'
    rng = random.Random(int.from_bytes(hashlib.sha256(pair_id.encode()).digest(), 'big'))
    rho = CONDITIONS[condition]
    z, values = rng.gauss(0, 1), []
    for i in range(n):
        if i:
            z = rho * z + math.sqrt(1 - rho * rho) * rng.gauss(0, 1)
        values.append(10 * (1 + change * (i >= position)) + 10 * noise * z)
    return {
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


def cases():
    for condition in CONDITIONS:
        for args in itertools.product(
            *(DESIGN[key] for key in ('sizes', 'changes', 'noise', 'locations', 'seeds'))
        ):
            yield make_case(condition, *args)


def true_pair_diagnostic(method, result, case, oracle, models, global_log_q, config, unit):
    split, rho = case['position'], Fraction(case['correlation'])
    if method == 'full_history':
        log_q = global_log_q
    else:
        saved = result['confidence']['by_split'].get(str(split))
        log_q = (
            saved['log_predictive_density']
            if saved is not None
            else jump._indexed_log_density(
                [x / unit for x in case['values']],
                split,
                global_log_q,
                original_levels=method == 'indexed_original',
            )
        )
    coefficient = full.confidence_coefficient(log_q, case['n'], config['confidence_alpha'])
    try:
        numerator, denominator = full.residual_ratio(models, split)
        retained = not full.excluded_at(numerator, denominator, coefficient, case['n'], rho)
    except ar1._DegenerateResidual:
        retained = True
    route = (
        'confidence'
        if not retained
        else (
            'abstain' if oracle['status'] != 'ok' else previous.c.s.rejection_route(oracle, split)
        )
    )
    if result['has_alert'] and route in ('neither', 'abstain'):
        raise AssertionError('Certified alert did not reject the generating pair')
    return {
        'true_pair_retained': retained,
        'generating_pair_route': route,
        'generating_split_log_density': log_q,
        'confidence_disabled_at_true_split': coefficient == 0,
    }


def evaluate(case, frozen, backend):
    result = previous.evaluate(case, frozen['inherited'], backend)
    config = frozen['calibrations'][str(case['n'])]
    models = ar1.Models([x / frozen['unit'] for x in case['values']])
    functions = (full.evidence, jump.original_location_evidence, jump.evidence)
    global_log_q = None
    for name, function in zip(NEW_METHODS, functions, strict=True):
        started = time.perf_counter()
        evidence = function(case['values'], unit=frozen['unit'], config=config, **frozen['work'])
        result[name + '_seconds'] = time.perf_counter() - started
        result[name] = evidence
        if name == 'full_history':
            global_log_q = evidence['confidence']['log_predictive_density']
        result[name + '_diagnostics'] = true_pair_diagnostic(
            name,
            evidence,
            case,
            result['oracle_reporting'],
            models,
            global_log_q,
            config,
            frozen['unit'],
        )
        result['methods'][name + '_alone'] = evidence['has_alert']
        result['methods'][name + '_gate'] = (
            evidence['has_alert'] and result['methods']['shared_existing']
        )
    return result


def comparisons(rows):
    result = {}
    for condition in ('all', *CONDITIONS):
        part = [row for row in rows if condition == 'all' or row['case']['condition'] == condition]
        result[condition] = {}
        for new, old in PAIRS:
            for suffix in ('alone', 'gate'):
                counts = {
                    'positive_gained': 0,
                    'positive_lost': 0,
                    'null_added': 0,
                    'null_removed': 0,
                }
                for row in part:
                    before, after = (
                        row['methods'][old + '_' + suffix],
                        row['methods'][new + '_' + suffix],
                    )
                    if row['case']['positive']:
                        counts['positive_gained'] += after and not before
                        counts['positive_lost'] += before and not after
                    else:
                        counts['null_added'] += after and not before
                        counts['null_removed'] += before and not after
                result[condition][new + '_vs_' + old + '_' + suffix] = counts
    return result


def diagnostics(rows):
    def summarize(part):
        result = {}
        for method in ('ar1', *NEW_METHODS):
            result[method] = {
                'histories': len(part),
                'statuses': {
                    status: sum(row[method]['status'] == status for row in part)
                    for status in (
                        'certified_alert',
                        'surviving_explanation',
                        'unresolved',
                        'insufficient_variation',
                    )
                },
                'true_pair_excluded': sum(
                    not row[method + '_diagnostics']['true_pair_retained'] for row in part
                ),
                'confidence_disabled_at_true_split': None
                if method == 'ar1'
                else sum(
                    row[method + '_diagnostics']['confidence_disabled_at_true_split']
                    for row in part
                ),
                'median_cells': statistics.median(row[method]['cells_visited'] for row in part),
                'max_cells': max(row[method]['cells_visited'] for row in part),
                'median_seconds': statistics.median(row[method + '_seconds'] for row in part),
                'max_seconds': max(row[method + '_seconds'] for row in part),
                'false_alert_routes': {
                    route: sum(
                        row[method + '_diagnostics']['generating_pair_route'] == route
                        for row in part
                        if not row['case']['positive'] and row['methods'][method + '_gate']
                    )
                    for route in (
                        'confidence',
                        'size_only',
                        'shape_only',
                        'both',
                        'neither',
                        'abstain',
                    )
                },
            }
        return result

    return {
        'all': summarize(rows),
        **{
            key: {
                str(value): summarize([row for row in rows if row['case'][key] == value])
                for value in sorted({row['case'][key] for row in rows})
            }
            for key in ('condition', 'n', 'location', 'positive')
        },
    }


def source_hashes(backend):
    hashes = previous.source_hashes(backend)
    for name in (
        'reporting_ar1_full.py',
        'reporting_ar1_jump.py',
        'full_history_confidence.rst',
        'baseline_jump_prior.rst',
        'indexed_study.py',
        'indexed_protocol.rst',
    ):
        path = HERE / name
        hashes[str(path.relative_to(h.ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def freeze(inherited, backend):
    previous.check_frozen(inherited, backend)
    return {
        'design': DESIGN,
        'work': WORK,
        'unit': UNIT,
        'backend': backend,
        'calibrations': {str(n): ar1.calibration(n) for n in DESIGN['sizes']},
        'priors': {
            'full_history': full.PRIOR,
            'indexed_jump': jump.PRIOR,
            'indexed_original': {'original_density_weight': 0.5, 'original_prior': full.PRIOR},
        },
        'source_hashes': source_hashes(backend),
        'shape_hashes': previous.shape_hashes(),
        'inherited': inherited,
        'python': platform.python_version(),
        'numpy': np.__version__,
        'scipy': scipy.__version__,
    }


def check_frozen(frozen, backend):
    if frozen != freeze(frozen['inherited'], backend):
        raise ValueError(
            'Frozen indexed-study source, settings, covariance, or environment changed'
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['freeze', 'run'])
    parser.add_argument('--inherited', type=Path)
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--backend', choices=['native', 'python'], default='native')
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    if args.phase == 'freeze':
        if args.inherited is None:
            parser.error('freeze requires --inherited')
        if args.output.exists():
            raise FileExistsError(args.output)
        common.dump(args.output, freeze(json.loads(args.inherited.read_text()), args.backend))
        return
    if args.frozen is None:
        parser.error('run requires --frozen')
    frozen = json.loads(args.frozen.read_text())
    check_frozen(frozen, args.backend)
    args.output.mkdir(parents=True, exist_ok=False)
    common.dump(args.output / 'covariance_shapes.json', previous.shapes())
    common.dump(
        args.output / 'manifest.json',
        {
            'frozen_sha256': hashlib.sha256(args.frozen.read_bytes()).hexdigest(),
            'source_hashes': source_hashes(args.backend),
            'shape_hashes': previous.shape_hashes(),
            'git_revision': subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], text=True
            ).strip(),
            'python': platform.python_version(),
            'numpy': np.__version__,
            'scipy': scipy.__version__,
            'design': DESIGN,
            'work': WORK,
            'unit': UNIT,
            'backend': args.backend,
            'workers': args.workers,
            'thread_environment': {
                key: os.environ.get(key)
                for key in (
                    'VECLIB_MAXIMUM_THREADS',
                    'OPENBLAS_NUM_THREADS',
                    'OMP_NUM_THREADS',
                )
            },
        },
    )
    inputs, rows = list(cases()), []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as pool,
        (args.output / 'inputs.jsonl').open('w') as input_file,
        (args.output / 'records.jsonl').open('w') as output_file,
    ):
        run_case = partial(evaluate, frozen=frozen, backend=args.backend)
        for case, row in zip(inputs, pool.map(run_case, inputs), strict=True):
            input_file.write(json.dumps(case, allow_nan=False) + '\n')
            output_file.write(json.dumps(reporting.json_safe(row), allow_nan=False) + '\n')
            input_file.flush()
            output_file.flush()
            rows.append(row)
            if len(rows) % 10 == 0:
                print(f'Evaluated {len(rows)}: {case["condition"]}, n={case["n"]}', flush=True)
    common.dump(args.output / 'summary.json', reporting.summarize(rows))
    common.dump(
        args.output / 'breakdown.json',
        {
            key: {
                str(value): reporting.summarize([row for row in rows if row['case'][key] == value])
                for value in sorted({row['case'][key] for row in rows})
            }
            for key in ('condition', 'n', 'change', 'noise', 'location')
        },
    )
    common.dump(args.output / 'comparisons.json', comparisons(rows))
    common.dump(args.output / 'diagnostics.json', diagnostics(rows))


if __name__ == '__main__':
    main()
