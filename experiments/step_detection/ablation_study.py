"""Factorial scoring comparison and frozen false-alert-budget calibration."""

import argparse
import hashlib
import itertools
import json
import math
import platform
import statistics
import subprocess
from pathlib import Path

from . import harness as h
from . import threshold_study as t
from .noise_model import selection_score

HERE = Path(__file__).parent
DESIGN = dict(
    t.DESIGN, development_seeds=list(range(200, 210)), heldout_seeds=list(range(300, 320))
)
PENALTIES = [1, 2, 3, 4, 6, 8, 12, 16]
BUDGETS = [0, 0.01, 0.05]


def configurations():
    return [
        {
            'name': f'{floor}-r{cap}-c{c}',
            'family': f'{floor}-r{cap}',
            'floor': floor,
            'cap': cap,
            'penalty': c,
        }
        for floor, cap, c in itertools.product(['current', 'shared'], [0.5, 1], PENALTIES)
    ]


def cases(phase):
    for args in itertools.product(
        DESIGN['sizes'],
        DESIGN['changes'],
        DESIGN['noise'],
        DESIGN['correlations'],
        DESIGN['locations'],
        DESIGN[f'{phase}_seeds'],
    ):
        yield t.make_case(*args)


def candidate_score(fit, error, n, scale, cfg):
    beta = cfg['penalty'] * math.log(n) / n
    if cfg['floor'] == 'shared':
        return selection_score(error, n, len(fit['right']), max(1e-12, 0.5 * scale), beta)
    levels = fit['values']
    if len(levels) > 2:
        floor = 0.1 * min(abs(b - a) for a, b in zip(levels, levels[1:]))
    else:
        floor = 0.001 * min(abs(v) for v in levels)
    return beta * len(levels) + math.log(max(1e-300, floor) + error)


def evaluate(case, configs, backend):
    module = h.load_detector(backend)
    y = case['values']
    production, pool = t.candidate_pool(module, y, backend)
    scale = statistics.median(abs(b - a) for a, b in zip(y, y[1:]))
    losses = {
        cap: [t.innovations(module, y, fit, cap) for fit in pool]
        for cap in {cfg['cap'] for cfg in configs}
    }
    chosen = {'production': production}
    for cfg in configs:
        scores = [
            candidate_score(fit, error, len(y), scale, cfg)
            for fit, error in zip(pool, losses[cfg['cap']])
        ]
        chosen[cfg['name']] = pool[min(range(len(pool)), key=scores.__getitem__)]
    return {
        'case': {k: v for k, v in case.items() if k != 'values'},
        'estimated_difference_scale': scale,
        'candidates': len(pool),
        'methods': {name: t.assess(module, case, fit) for name, fit in chosen.items()},
    }


def calibrate(summary, configs):
    selections = []
    for family, budget in itertools.product(sorted({c['family'] for c in configs}), BUDGETS):
        eligible = [
            c
            for c in configs
            if c['family'] == family and summary[c['name']]['false_rate'] <= budget
        ]

        def key(cfg):
            s = summary[cfg['name']]
            return s['miss_rate'], s['false_rate'], -cfg['penalty']

        selections.append(
            {
                'family': family,
                'budget': budget,
                'config': min(eligible, key=key) if eligible else None,
            }
        )
    return selections


def evaluation_configs(selections):
    chosen = {c['name']: c for c in configurations() if c['penalty'] in (2, 4)}
    for row in selections:
        if row['config'] is not None:
            chosen[row['config']['name']] = row['config']
    return list(chosen.values())


def source_hashes(backend):
    hashes = t.source_hashes(backend)
    for p in (Path(__file__), HERE / 'ablation_protocol.rst'):
        hashes[str(p.relative_to(h.ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return hashes


def check_frozen(frozen, hashes, backend):
    if (
        frozen['source_hashes'] != hashes
        or frozen['design'] != DESIGN
        or frozen['backend'] != backend
        or frozen['penalties'] != PENALTIES
        or frozen['budgets'] != BUDGETS
    ):
        raise ValueError('Frozen code, protocol, design, backend, or calibration grid changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['development', 'heldout'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--backend', choices=['native', 'python'], default='native')
    args = parser.parse_args()
    hashes = source_hashes(args.backend)
    if args.phase == 'heldout':
        if args.frozen is None:
            parser.error('heldout requires --frozen')
        frozen = json.loads(args.frozen.read_text())
        check_frozen(frozen, hashes, args.backend)
        configs = evaluation_configs(frozen['selections'])
    else:
        configs = configurations()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {
        'phase': args.phase,
        'source_hashes': hashes,
        'design': DESIGN,
        'backend': args.backend,
        'configurations': configs,
        'penalties': PENALTIES,
        'budgets': BUDGETS,
        'python': platform.python_version(),
        'platform': platform.platform(),
        'git_revision': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=h.ROOT, text=True
        ).strip(),
        'frozen_sha256': hashlib.sha256(args.frozen.read_bytes()).hexdigest()
        if args.frozen
        else None,
    }
    t.dump(args.output / 'manifest.json', manifest)
    rows = []
    with (
        (args.output / 'inputs.jsonl').open('w') as inputs,
        (args.output / 'records.jsonl').open('w') as records,
    ):
        for i, case in enumerate(cases(args.phase)):
            row = evaluate(case, configs, args.backend)
            inputs.write(json.dumps(case, allow_nan=False) + '\n')
            records.write(json.dumps(row, allow_nan=False) + '\n')
            inputs.flush()
            records.flush()
            rows.append(row)
            if (i + 1) % 100 == 0:
                print(f'{args.phase}: {i + 1} histories', flush=True)
    summary = t.summarize(rows)
    t.dump(args.output / 'summary.json', summary)
    breakdown = {
        field: {
            str(value): t.summarize([r for r in rows if r['case'][field] == value])
            for value in sorted({r['case'][field] for r in rows})
        }
        for field in ('n', 'change', 'noise', 'correlation', 'location')
    }
    t.dump(args.output / 'breakdown.json', breakdown)
    if args.phase == 'development':
        t.dump(
            args.output / 'frozen.json',
            {
                'source_hashes': hashes,
                'design': DESIGN,
                'backend': args.backend,
                'penalties': PENALTIES,
                'budgets': BUDGETS,
                'selections': calibrate(summary, configs),
            },
        )


if __name__ == '__main__':
    main()
