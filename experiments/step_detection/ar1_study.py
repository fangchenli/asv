"""Frozen evaluation of joint AR(1) correlation/location uncertainty."""

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
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import partial
from pathlib import Path

import numpy as np
import scipy

from . import covariance_study as c
from . import harness as h
from . import rational_polynomial as p
from . import reporting_ar1 as a
from . import reporting_study as r
from . import threshold_study as t

HERE = Path(__file__).parent
CONDITIONS = {'negative': -0.5, 'independent_constant': 0, 'positive': 0.7}
DESIGN = {
    'sizes': [40, 100],
    'changes': [0, 0.04, 0.05, 0.06, 0.08],
    'noise': [0.005, 0.02],
    'locations': ['early', 'middle', 'recent'],
    'seeds': [1400, 1401],
    'conditions': CONDITIONS,
}
WORK = {'max_cells': 4096, 'max_depth': 16}


def make_case(condition, n, change, noise, location, seed):
    position = {'early': n // 4, 'middle': n // 2, 'recent': n - max(4, n // 10)}[location]
    pair_id = f'ar1-study-n{n}-d{change}-s{noise}-{location}-seed{seed}'
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


def shapes():
    return {
        f'{condition}-n{n}': (
            rho ** np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        ).tolist()
        for condition, rho in CONDITIONS.items()
        for n in DESIGN['sizes']
    }


def shape_hashes():
    return {
        key: hashlib.sha256(json.dumps(matrix, separators=(',', ':')).encode()).hexdigest()
        for key, matrix in shapes().items()
    }


def rethreshold(oracle, calibration):
    result = dict(oracle)
    if oracle['status'] != 'ok':
        return result
    size = [
        i + 1 for i, value in enumerate(oracle['size_t']) if value > calibration['size_t_critical']
    ]
    shape = [
        i + 1
        for i, value in enumerate(oracle['lack_of_fit_f'])
        if value > calibration['fit_f_critical']
    ]
    surviving = sorted(set(range(1, calibration['n'])) - set(size) - set(shape))
    result.update(
        has_alert=not surviving,
        null_splits=surviving,
        size_rejected_splits=size,
        fit_rejected_splits=shape,
    )
    return result


def quadratic_interval(coefficients, bound):
    """Descriptive intersection of a convex quadratic sublevel set with [-1,1]."""
    if math.isinf(bound):
        return [-1.0, 1.0]
    coefficients = [Fraction(x) for x in coefficients]
    coefficients += [Fraction(0)] * (3 - len(coefficients))
    constant, linear, square = coefficients
    constant -= Fraction(bound)
    if square < 0:
        raise ValueError('Residual quadratic must be convex')
    if square == 0:
        if linear == 0:
            return [-1.0, 1.0] if constant <= 0 else None
        root = -constant / linear
        left, right = (
            (Fraction(-1), min(Fraction(1), root))
            if linear > 0
            else (max(Fraction(-1), root), Fraction(1))
        )
        return [float(left), float(right)] if left <= right and right > -1 and left < 1 else None
    discriminant = linear * linear - 4 * square * constant
    if discriminant < 0:
        return None
    with localcontext() as context:
        context.prec = 80

        def dec(x):
            return Decimal(x.numerator) / Decimal(x.denominator)

        radius = dec(discriminant).sqrt()
        left = max(Decimal(-1), (-dec(linear) - radius) / (2 * dec(square)))
        right = min(Decimal(1), (-dec(linear) + radius) / (2 * dec(square)))
        return [float(left), float(right)] if left <= right and right > -1 and left < 1 else None


def permits_near_one(coefficients, bound):
    if math.isinf(bound):
        return True
    coefficients = [Fraction(x) for x in coefficients]
    coefficients += [Fraction(0)] * (3 - len(coefficients))
    constant, linear, square = coefficients
    value = constant + linear + square - Fraction(bound)
    return value < 0 or (value == 0 and (linear + 2 * square > 0 or linear == square == 0))


def confidence_diagnostics(unknown, case, oracle_reporting):
    bound = unknown['confidence']['residual_bound']
    intervals = {
        key: quadratic_interval(poly, bound)
        for key, poly in unknown['confidence_polynomials'].items()
    }
    merged = []
    for left, right in sorted(value for value in intervals.values() if value is not None):
        if merged and left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    split = case['position']
    poly = tuple(Fraction(x) for x in unknown['confidence_polynomials'][str(split)])
    retained = math.isinf(bound) or p.evaluate(poly, Fraction(case['correlation'])) <= Fraction(
        bound
    )
    route = (
        'confidence'
        if not retained
        else (
            'abstain'
            if oracle_reporting['status'] != 'ok'
            else c.s.rejection_route(oracle_reporting, split)
        )
    )
    if unknown['has_alert'] and route in ('neither', 'abstain'):
        raise AssertionError('Certified alert did not reject the generating pair')
    true_interval = intervals[str(split)]
    return {
        'intervals_by_split': intervals,
        'correlation_union': merged,
        'union_width': sum(right - left for left, right in merged),
        'generating_split_width': 0
        if true_interval is None
        else true_interval[1] - true_interval[0],
        'permits_near_one': any(
            permits_near_one(poly, bound) for poly in unknown['confidence_polynomials'].values()
        ),
        'true_pair_retained': retained,
        'generating_pair_route': route,
    }


def evaluate(case, frozen, backend):
    result = c.evaluate(case, frozen['inherited'], backend)
    config = frozen['calibrations'][str(case['n'])]
    oracle_reporting = rethreshold(result['oracle'], config)
    started = time.perf_counter()
    unknown = a.evidence(case['values'], config, **frozen['work'])
    seconds = time.perf_counter() - started
    result['oracle_reporting'] = oracle_reporting
    result['ar1'] = unknown
    result['ar1_seconds'] = seconds
    result['ar1_diagnostics'] = confidence_diagnostics(unknown, case, oracle_reporting)
    for prefix, evidence in [('oracle_reporting', oracle_reporting), ('ar1', unknown)]:
        result['methods'][prefix + '_alone'] = evidence['has_alert']
        result['methods'][prefix + '_gate'] = (
            evidence['has_alert'] and result['methods']['shared_existing']
        )
    return result


def comparisons(rows):
    result = {}
    for condition in ('all', *CONDITIONS):
        part = [x for x in rows if condition == 'all' or x['case']['condition'] == condition]
        result[condition] = {}
        for old_method in ('direct', 'oracle', 'oracle_reporting'):
            for suffix in ('alone', 'gate'):
                counts = {
                    'positive_gained': 0,
                    'positive_lost': 0,
                    'null_added': 0,
                    'null_removed': 0,
                }
                for row in part:
                    old, new = (
                        row['methods'][old_method + '_' + suffix],
                        row['methods']['ar1_' + suffix],
                    )
                    if row['case']['positive']:
                        counts['positive_gained'] += new and not old
                        counts['positive_lost'] += old and not new
                    else:
                        counts['null_added'] += new and not old
                        counts['null_removed'] += old and not new
                result[condition][old_method + '_' + suffix] = counts
    return result


def diagnostics(rows):
    def part_summary(part):
        return {
            'histories': len(part),
            'statuses': {
                status: sum(row['ar1']['status'] == status for row in part)
                for status in (
                    'certified_alert',
                    'surviving_explanation',
                    'unresolved',
                    'insufficient_variation',
                )
            },
            'true_pair_excluded': sum(
                not row['ar1_diagnostics']['true_pair_retained'] for row in part
            ),
            'permits_near_one': sum(row['ar1_diagnostics']['permits_near_one'] for row in part),
            'median_union_width': statistics.median(
                row['ar1_diagnostics']['union_width'] for row in part
            ),
            'median_generating_split_width': statistics.median(
                row['ar1_diagnostics']['generating_split_width'] for row in part
            ),
            'median_cells': statistics.median(row['ar1']['cells_visited'] for row in part),
            'max_cells': max(row['ar1']['cells_visited'] for row in part),
            'median_seconds': statistics.median(row['ar1_seconds'] for row in part),
            'max_seconds': max(row['ar1_seconds'] for row in part),
            'false_alert_routes': {
                route: sum(
                    row['ar1_diagnostics']['generating_pair_route'] == route
                    for row in part
                    if not row['case']['positive'] and row['methods']['ar1_gate']
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

    return {
        'all': part_summary(rows),
        **{
            key: {
                str(value): part_summary([row for row in rows if row['case'][key] == value])
                for value in sorted({row['case'][key] for row in rows})
            }
            for key in ('condition', 'n', 'location', 'positive')
        },
    }


def source_hashes(backend):
    hashes = c.source_hashes(backend)
    for name in (
        'rational_polynomial.py',
        'reporting_ar1.py',
        'unknown_correlation.rst',
        'ar1_study.py',
        'ar1_protocol.rst',
    ):
        path = HERE / name
        hashes[str(path.relative_to(h.ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def freeze(inherited, backend):
    c.check_frozen(inherited, backend)
    return {
        'design': DESIGN,
        'work': WORK,
        'backend': backend,
        'calibrations': {str(n): a.calibration(n) for n in DESIGN['sizes']},
        'source_hashes': source_hashes(backend),
        'shape_hashes': shape_hashes(),
        'inherited': inherited,
        'python': platform.python_version(),
        'numpy': np.__version__,
        'scipy': scipy.__version__,
    }


def check_frozen(frozen, backend):
    c.check_frozen(frozen['inherited'], backend)
    current = freeze(frozen['inherited'], backend)
    if frozen != current:
        raise ValueError('Frozen AR1 configuration, source, covariance, or environment changed')


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
            'git_revision': subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], text=True
            ).strip(),
            'python': platform.python_version(),
            'numpy': np.__version__,
            'scipy': scipy.__version__,
            'design': DESIGN,
            'work': WORK,
            'backend': args.backend,
            'workers': args.workers,
            'thread_environment': {
                key: os.environ.get(key)
                for key in ('VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS')
            },
        },
    )
    inputs, rows = list(cases()), []
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
            if len(rows) % 10 == 0:
                print(f'Evaluated {len(rows)}: {case["condition"]}, n={case["n"]}', flush=True)
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
    t.dump(args.output / 'diagnostics.json', diagnostics(rows))


if __name__ == '__main__':
    main()
