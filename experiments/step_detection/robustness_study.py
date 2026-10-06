"""Stress frozen scoring settings without further calibration."""

import argparse
import hashlib
import itertools
import json
import math
import platform
import random
import statistics
import subprocess
from pathlib import Path

from . import ablation_study as a
from . import harness as h
from . import threshold_study as t

HERE = Path(__file__).parent
CONDITIONS = [
    'gaussian',
    'heavy_tails',
    'outliers',
    'missing_random',
    'missing_gap',
    'variance_up',
    'variance_down',
    'outliers_missing',
]
DESIGN = {
    'sizes': [40, 100],
    'changes': [0, 0.04, 0.06, 0.08],
    'noise': [0.005, 0.02, 0.05],
    'correlations': [0, 0.7],
    'locations': ['middle', 'recent'],
    'seeds': list(range(400, 410)),
    'conditions': CONDITIONS,
    'warmup': 256,
}


def make_case(n, change, noise, rho, location, seed, condition):
    if condition not in CONDITIONS:
        raise ValueError('Unknown noise condition')
    rng = random.Random(seed)
    outliers = random.Random(seed + 10000)
    missing = random.Random(seed + 20000)
    position = n // 2 if location == 'middle' else n - max(4, n // 10)
    sigma = 10 * noise
    error, values = 0, []
    for i in range(-DESIGN['warmup'], n):
        innovation = rng.gauss(0, 1)
        if condition == 'heavy_tails':
            innovation /= math.sqrt(sum(rng.gauss(0, 1) ** 2 for _ in range(3)))
        error = rho * error + math.sqrt(1 - rho * rho) * innovation
        if i < 0:
            continue
        multiplier = (
            3
            if (
                (condition == 'variance_up' and i >= n // 2)
                or (condition == 'variance_down' and i < n // 2)
            )
            else 1
        )
        value = 10 * (1 + change * (i >= position)) + sigma * multiplier * error
        if condition in ('outliers', 'outliers_missing') and outliers.random() < 0.03:
            value += 10 * sigma
        if (
            condition in ('missing_random', 'outliers_missing')
            and missing.random() < 0.3
            and 0 < i < n - 1
        ):
            value = None
        values.append(value)
    if condition == 'missing_gap':
        width = max(2, n // 10)
        start = max(1, position - width // 2)
        for i in range(start, min(n - 1, start + width)):
            values[i] = None
    pair_id = f'n{n}-d{change}-s{noise}-r{rho}-{location}-seed{seed}'
    return {
        'id': condition + '-' + pair_id,
        'pair_id': pair_id,
        'condition': condition,
        'n': n,
        'change': change,
        'noise': noise,
        'correlation': rho,
        'location': location,
        'seed': seed,
        'position': position,
        'positive': change > 0.05,
        'values': values,
        'weights': [1] * n,
    }


def cases():
    for condition in CONDITIONS:
        for args in itertools.product(
            DESIGN['sizes'],
            DESIGN['changes'],
            DESIGN['noise'],
            DESIGN['correlations'],
            DESIGN['locations'],
            DESIGN['seeds'],
        ):
            yield make_case(*args, condition)


def assess(module, case, y, indices, fit):
    steps, left = [], 0
    for right, value, cost in zip(fit['right'], fit['values'], fit['costs']):
        steps.append(
            (
                indices[left],
                indices[right - 1] + 1,
                value,
                min(y[left:right]),
                abs(cost / (right - left)),
            )
        )
        left = right
    _, _, alerts = module.detect_regressions(steps, threshold=0.05)
    positions = [r[1] for r in (alerts or [])]
    boundaries = [indices[r] for r in fit['right'][:-1]]
    truth = h.observed_truth([case['position']] if case['change'] else [], indices)
    return {
        'segments': len(steps),
        'boundaries': boundaries,
        'alerts': positions,
        'has_alert': bool(alerts),
        'boundary_metrics': h.match_boundaries(truth, boundaries),
        'alert_metrics': h.match_boundaries(truth if case['positive'] else [], positions),
    }


def evaluate(case, configs, backend='native'):
    module = h.load_detector(backend)
    y, weights, indices = h.prepare(module, case)
    assert weights == [1] * len(y)
    production, pool = t.candidate_pool(module, y, backend)
    scale = statistics.median(abs(b - a) for a, b in zip(y, y[1:]))
    losses = {
        cap: [t.innovations(module, y, fit, cap) for fit in pool]
        for cap in {c['cap'] for c in configs}
    }
    selected = {'production': production}
    for cfg in configs:
        scores = [
            a.candidate_score(fit, loss, len(y), scale, cfg)
            for fit, loss in zip(pool, losses[cfg['cap']])
        ]
        selected[cfg['name']] = pool[min(range(len(pool)), key=scores.__getitem__)]
    return {
        'case': {k: v for k, v in case.items() if k not in ('values', 'weights')},
        'retained': len(y),
        'estimated_difference_scale': scale,
        'candidates': len(pool),
        'methods': {name: assess(module, case, y, indices, fit) for name, fit in selected.items()},
    }


def source_hashes(backend):
    hashes = a.source_hashes(backend)
    for p in (Path(__file__), HERE / 'robustness_protocol.rst'):
        hashes[str(p.relative_to(h.ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return hashes


def freeze(calibration, backend):
    a.check_frozen(calibration, a.source_hashes(backend), backend)
    selections = [
        row for row in calibration['selections'] if row['family'] in ('current-r1', 'shared-r1')
    ]
    if len(selections) != 6 or any(row['config'] is None for row in selections):
        raise ValueError('Expected six inherited feasible configurations')
    configs = {row['config']['name']: row['config'] for row in selections}
    return {
        'source_hashes': source_hashes(backend),
        'design': DESIGN,
        'backend': backend,
        'selections': selections,
        'configurations': list(configs.values()),
    }


def check_frozen(frozen, backend):
    if (
        frozen['source_hashes'] != source_hashes(backend)
        or frozen['design'] != DESIGN
        or frozen['backend'] != backend
    ):
        raise ValueError('Frozen robustness code, protocol, design, or backend changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['freeze', 'run'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--calibration', type=Path)
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--backend', choices=['native', 'python'], default='native')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must not already exist')
    if args.phase == 'freeze':
        if args.calibration is None:
            parser.error('freeze requires --calibration')
        frozen = freeze(json.loads(args.calibration.read_text()), args.backend)
        frozen['calibration_sha256'] = hashlib.sha256(args.calibration.read_bytes()).hexdigest()
        t.dump(args.output, frozen)
        return
    if args.frozen is None:
        parser.error('run requires --frozen')
    frozen = json.loads(args.frozen.read_text())
    check_frozen(frozen, args.backend)
    args.output.mkdir(parents=True)
    t.dump(
        args.output / 'manifest.json',
        {
            'frozen': frozen,
            'frozen_sha256': hashlib.sha256(args.frozen.read_bytes()).hexdigest(),
            'git_revision': subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], cwd=h.ROOT, text=True
            ).strip(),
            'python': platform.python_version(),
            'platform': platform.platform(),
        },
    )
    rows = []
    with (
        (args.output / 'inputs.jsonl').open('w') as inputs,
        (args.output / 'records.jsonl').open('w') as records,
    ):
        for i, case in enumerate(cases()):
            row = evaluate(case, frozen['configurations'], args.backend)
            inputs.write(json.dumps(case, allow_nan=False) + '\n')
            records.write(json.dumps(row, allow_nan=False) + '\n')
            inputs.flush()
            records.flush()
            rows.append(row)
            if (i + 1) % 200 == 0:
                print(f'{i + 1} histories: {case["condition"]}', flush=True)
    t.dump(args.output / 'summary.json', t.summarize(rows))
    breakdown = {
        field: {
            str(value): t.summarize([r for r in rows if r['case'][field] == value])
            for value in sorted({r['case'][field] for r in rows})
        }
        for field in ('condition', 'n', 'change', 'noise', 'correlation', 'location')
    }
    t.dump(args.output / 'breakdown.json', breakdown)
    paired = {}
    controls = {r['case']['pair_id']: r for r in rows if r['case']['condition'] == 'gaussian'}
    for condition in CONDITIONS:
        paired[condition] = {}
        subset = [r for r in rows if r['case']['condition'] == condition]
        for method in rows[0]['methods']:
            totals = {
                'positive_gained': 0,
                'positive_lost': 0,
                'negative_fixed': 0,
                'negative_broken': 0,
            }
            for row in subset:
                old = controls[row['case']['pair_id']]['methods'][method]['has_alert']
                new = row['methods'][method]['has_alert']
                if row['case']['positive']:
                    totals['positive_gained'] += new and not old
                    totals['positive_lost'] += old and not new
                else:
                    totals['negative_fixed'] += old and not new
                    totals['negative_broken'] += new and not old
            paired[condition][method] = totals
    t.dump(args.output / 'paired.json', paired)


if __name__ == '__main__':
    main()
