"""Stress the frozen direct test with variance changes and serial correlation."""

import argparse
import hashlib
import itertools
import json
import math
import platform
import random
import subprocess
from pathlib import Path

from . import direct_study as d
from . import harness as h
from . import reporting_study as r
from . import threshold_study as t

HERE = Path(__file__).parent
CONDITIONS = {
    'independent_constant': (0, 'constant'),
    'independent_up': (0, 'up'),
    'independent_down': (0, 'down'),
    'correlated_constant': (0.7, 'constant'),
    'correlated_up': (0.7, 'up'),
    'correlated_down': (0.7, 'down'),
}
DESIGN = {
    'sizes': [40, 100],
    'changes': [0, 0.04, 0.05, 0.06, 0.08],
    'noise': [0.005, 0.02],
    'locations': ['middle', 'recent'],
    'seeds': list(range(1000, 1010)),
    'conditions': {name: list(config) for name, config in CONDITIONS.items()},
    'amplitude_multiplier': 3,
}


def amplitudes(n, position, profile):
    if profile not in ('constant', 'up', 'down'):
        raise ValueError('Unknown variance profile')
    return [
        DESIGN['amplitude_multiplier']
        if (profile == 'up' and i >= position) or (profile == 'down' and i < position)
        else 1
        for i in range(n)
    ]


def variance_diagnostic(n, position, rho, profile):
    """True contrast variance / expected variance reported at a specified split.

    Noise scale is one; it cancels from the ratio. The split is the generating
    location, not a boundary selected by ASV. This is not the variance of a t
    statistic and does not give corrected critical values.
    """
    if not 1 <= position < n or not -1 < rho < 1:
        raise ValueError('Require an internal split and stationary correlation')
    scale = amplitudes(n, position, profile)
    c = [-1.05 / position if i < position else 1 / (n - position) for i in range(n)]
    covariance = [[scale[i] * scale[j] * rho ** abs(i - j) for j in range(n)] for i in range(n)]
    actual = math.fsum(c[i] * c[j] * covariance[i][j] for i in range(n) for j in range(n))
    residual = (
        math.fsum(x * x for x in scale)
        - math.fsum(covariance[i][j] for i in range(position) for j in range(position)) / position
        - math.fsum(covariance[i][j] for i in range(position, n) for j in range(position, n))
        / (n - position)
    )
    expected = residual / (n - 2) * (1.05**2 / position + 1 / (n - position))
    return {
        'contrast_variance': actual,
        'expected_residual_sum': residual,
        'expected_reported_variance': expected,
        'variance_ratio': actual / expected,
    }


def make_case(condition, n, change, noise, location, seed):
    rho, profile = CONDITIONS[condition]
    position = n // 2 if location == 'middle' else n - max(4, n // 10)
    pair_id = f'direct-stress-n{n}-d{change}-s{noise}-{location}-seed{seed}'
    rng = random.Random(int.from_bytes(hashlib.sha256(pair_id.encode()).digest(), 'big'))
    scale = amplitudes(n, position, profile)
    z, values = rng.gauss(0, 1), []
    for i in range(n):
        if i:
            z = rho * z + math.sqrt(1 - rho * rho) * rng.gauss(0, 1)
        values.append(10 * (1 + change * (i >= position)) + 10 * noise * scale[i] * z)
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
        'variance_profile': profile,
    }


def cases():
    for condition in CONDITIONS:
        for args in itertools.product(
            *(DESIGN[key] for key in ('sizes', 'changes', 'noise', 'locations', 'seeds'))
        ):
            yield make_case(condition, *args)


def rejection_route(direct, split):
    size = split in direct['size_rejected_splits']
    shape = split in direct['fit_rejected_splits']
    return (
        'both' if size and shape else 'size_only' if size else 'shape_only' if shape else 'neither'
    )


def evaluate(case, frozen, backend):
    result = d.evaluate(case, frozen['inherited'], backend)
    result['generating_split_route'] = rejection_route(result['direct'], case['position'])
    return result


def source_hashes(backend):
    hashes = d.source_hashes(backend)
    for name in ('direct_stress.py', 'direct_stress_protocol.rst'):
        p = HERE / name
        hashes[str(p.relative_to(h.ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return hashes


def freeze(inherited, backend):
    d.check_frozen(inherited, backend)
    return {
        'design': DESIGN,
        'backend': backend,
        'source_hashes': source_hashes(backend),
        'inherited': inherited,
        'variance_diagnostics': {
            f'{condition}-n{n}-{location}': variance_diagnostic(
                n, n // 2 if location == 'middle' else n - max(4, n // 10), rho, profile
            )
            for condition, (rho, profile) in CONDITIONS.items()
            for n, location in itertools.product(DESIGN['sizes'], DESIGN['locations'])
        },
    }


def check_frozen(frozen, backend):
    d.check_frozen(frozen['inherited'], backend)
    if (
        frozen['design'] != DESIGN
        or frozen['backend'] != backend
        or frozen['source_hashes'] != source_hashes(backend)
    ):
        raise ValueError('Frozen design, source, or backend changed')


def paired(rows):
    controls = {
        row['case']['pair_id']: row
        for row in rows
        if row['case']['condition'] == 'independent_constant'
    }
    result = {}
    for condition in CONDITIONS:
        part = [row for row in rows if row['case']['condition'] == condition]
        result[condition] = {}
        for method in rows[0]['methods']:
            counts = {'positive_gained': 0, 'positive_lost': 0, 'null_added': 0, 'null_removed': 0}
            for row in part:
                old = controls[row['case']['pair_id']]['methods'][method]
                new = row['methods'][method]
                if row['case']['positive']:
                    counts['positive_gained'] += new and not old
                    counts['positive_lost'] += old and not new
                else:
                    counts['null_added'] += new and not old
                    counts['null_removed'] += old and not new
            result[condition][method] = counts
    return result


def routes(rows):
    result = {}
    for condition in CONDITIONS:
        result[condition] = {}
        for method in ('direct_alone', 'direct_gate'):
            part = [
                row
                for row in rows
                if row['case']['condition'] == condition
                and not row['case']['positive']
                and row['methods'][method]
            ]
            result[condition][method] = {
                route: sum(row['generating_split_route'] == route for row in part)
                for route in ('size_only', 'shape_only', 'both', 'neither')
            }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['freeze', 'run'])
    parser.add_argument('--inherited', type=Path)
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--backend', choices=['native', 'python'], default='native')
    args = parser.parse_args()
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
    t.dump(
        args.output / 'manifest.json',
        {
            'frozen_sha256': hashlib.sha256(args.frozen.read_bytes()).hexdigest(),
            'source_hashes': source_hashes(args.backend),
            'python': platform.python_version(),
            'git_revision': subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], text=True
            ).strip(),
            'design': DESIGN,
            'backend': args.backend,
        },
    )
    rows = []
    with (
        (args.output / 'inputs.jsonl').open('w') as inputs,
        (args.output / 'records.jsonl').open('w') as outputs,
    ):
        for case in cases():
            row = evaluate(case, frozen, args.backend)
            inputs.write(json.dumps(case, allow_nan=False) + '\n')
            outputs.write(json.dumps(r.json_safe(row), allow_nan=False) + '\n')
            inputs.flush()
            outputs.flush()
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
    t.dump(args.output / 'paired.json', paired(rows))
    t.dump(args.output / 'routes.json', routes(rows))


if __name__ == '__main__':
    main()
