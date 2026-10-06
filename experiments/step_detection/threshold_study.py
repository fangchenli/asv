"""Development/frozen held-out evaluation of near-threshold regression alerts."""

import argparse
import hashlib
import itertools
import json
import math
import platform
import random
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

from . import harness as h
from .exact_reference import fit_independent
from .inspect_comparison import score_candidate
from .noise_model import correlation_cap, selection_score

HERE = Path(__file__).parent
DESIGN = {
    'sizes': [40, 100],
    'changes': [0, 0.04, 0.05, 0.06, 0.08],
    'noise': [0.005, 0.02, 0.05],
    'correlations': [0, 0.7],
    'locations': ['middle', 'recent'],
    'development_seeds': [0, 1, 2],
    'heldout_seeds': [100, 101, 102, 103, 104],
}


def configurations():
    return [
        {
            'name': f'{"independent" if half == 0 else "bounded"}-c{c}-f{floor}-h{half}',
            'family': 'independent' if half == 0 else 'bounded',
            'penalty': c,
            'floor_factor': floor,
            'half_life': half,
        }
        for c, floor, half in itertools.product([2, 4, 8], [0.25, 0.5, 1], [0, 1, 4])
    ]


def make_case(n, change, noise, rho, location, seed):
    rng = random.Random(seed)
    pos = n // 2 if location == 'middle' else n - max(4, n // 10)
    sigma = 10 * noise
    error = rng.gauss(0, sigma)
    values = []
    for i in range(n):
        if i:
            error = rho * error + rng.gauss(0, sigma * math.sqrt(1 - rho * rho))
        values.append(10 * (1 + change * (i >= pos)) + error)
    return {
        'id': f'n{n}-d{change}-s{noise}-r{rho}-{location}-seed{seed}',
        'n': n,
        'change': change,
        'noise': noise,
        'correlation': rho,
        'location': location,
        'seed': seed,
        'position': pos,
        'values': values,
        'positive': None if change == 0.05 else change > 0.05,
    }


def cases(phase):
    for args in itertools.product(
        DESIGN['sizes'],
        DESIGN['changes'],
        DESIGN['noise'],
        DESIGN['correlations'],
        DESIGN['locations'],
        DESIGN[f'{phase}_seeds'],
    ):
        yield make_case(*args)


def candidate_pool(module, y, backend):
    w = [1] * len(y)
    original = module.solve_potts_approx
    collected = []

    def capture(*args, **kwargs):
        right, values, costs = original(*args, **kwargs)
        collected.append({'right': right, 'values': values, 'costs': costs})
        return right, values, costs

    module.solve_potts_approx = capture
    try:
        right, values, costs, _ = module.solve_potts_autogamma(y, w)
    finally:
        module.solve_potts_approx = original
    production = {'right': right, 'values': values, 'costs': costs}
    frontier = fit_independent(y, w, noise_floor=1, beta=0, backend=backend)['frontier']
    unique = {}
    # Production candidates first: retain its tie choice when scores are equal.
    for fit in collected + frontier:
        unique.setdefault((tuple(fit['right']), tuple(fit['values'])), fit)
    return production, list(unique.values())


def innovations(module, y, fit, cap):
    residuals, left = [], 0
    for right, value in zip(fit['right'], fit['values']):
        residuals.extend(x - value for x in y[left:right])
        left = right
    rho = module._fit_ar1(residuals, [1] * len(y), rho_max=cap)
    return abs(residuals[0]) + math.fsum(
        abs(current - rho * previous) for previous, current in zip(residuals, residuals[1:])
    )


def assess(module, case, fit):
    y, left, steps = case['values'], 0, []
    for right, value, cost in zip(fit['right'], fit['values'], fit['costs']):
        steps.append((left, right, value, min(y[left:right]), abs(cost / (right - left))))
        left = right
    _, _, alerts = module.detect_regressions(steps, threshold=0.05)
    positions = [a[1] for a in (alerts or [])]
    truth = [case['position']] if case['change'] else []
    alert_truth = None if case['positive'] is None else truth if case['positive'] else []
    return {
        'segments': len(steps),
        'boundaries': fit['right'][:-1],
        'alerts': positions,
        'has_alert': bool(alerts),
        'boundary_metrics': h.match_boundaries(truth, fit['right'][:-1]),
        'alert_metrics': h.match_boundaries(alert_truth, positions),
    }


def evaluate(case, settings, backend='native'):
    module = h.load_detector(backend)
    y = case['values']
    production, pool = candidate_pool(module, y, backend)
    current_scores = [score_candidate(module, y, [1] * len(y), fit) for fit in pool]
    scale = max(1e-12, statistics.median(abs(b - a) for a, b in zip(y, y[1:])))
    losses = {
        half: [innovations(module, y, fit, correlation_cap(half)) for fit in pool]
        for half in {cfg['half_life'] for cfg in settings}
    }
    selected = {
        'production': production,
        'current_pool': pool[min(range(len(pool)), key=current_scores.__getitem__)],
    }
    for cfg in settings:
        beta = cfg['penalty'] * math.log(len(y)) / len(y)
        floor = max(1e-12, cfg['floor_factor'] * scale)
        scores = [
            selection_score(error, len(y), len(fit['right']), floor, beta)
            for fit, error in zip(pool, losses[cfg['half_life']])
        ]
        selected[cfg['name']] = pool[min(range(len(pool)), key=scores.__getitem__)]
    return {
        'case': {key: value for key, value in case.items() if key != 'values'},
        'estimated_difference_scale': scale,
        'candidates': len(pool),
        'methods': {name: assess(module, case, fit) for name, fit in selected.items()},
    }


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        for method, result in row['methods'].items():
            groups[method].append((row['case'], result))
    summary = {}
    for method, items in groups.items():
        negative = [r for c, r in items if c['positive'] is False]
        positive = [r for c, r in items if c['positive'] is True]
        ambiguous = [r for c, r in items if c['positive'] is None]
        fp = sum(r['has_alert'] for r in negative)
        misses = sum(not r['has_alert'] for r in positive)
        result = {
            'histories': len(items),
            'negative_histories': len(negative),
            'positive_histories': len(positive),
            'false_alert_histories': fp,
            'missed_alert_histories': misses,
            'false_rate': fp / len(negative) if negative else None,
            'miss_rate': misses / len(positive) if positive else None,
            'at_threshold_alerts': sum(r['has_alert'] for r in ambiguous),
            'at_threshold_histories': len(ambiguous),
            'mean_segments': statistics.mean(r['segments'] for _, r in items),
        }
        for field in ('boundary_metrics', 'alert_metrics'):
            ms = [r[field] for _, r in items if r[field] is not None]
            result[field] = {
                key: sum(m[key] for m in ms)
                for key in ('matched', 'missed', 'false', 'location_error_sum')
            }
        summary[method] = result
    return summary


def select_settings(summary, settings):
    selected = []
    for family in ('independent', 'bounded'):

        def key(cfg):
            score = summary[cfg['name']]
            return (
                score['false_rate'] + score['miss_rate'],
                score['false_rate'],
                -cfg['penalty'],
                -cfg['floor_factor'],
                cfg['half_life'],
            )

        selected.append(min((c for c in settings if c['family'] == family), key=key))
    return selected


def source_hashes(backend):
    paths = [
        Path(__file__),
        HERE / 'threshold_protocol.rst',
        HERE / 'harness.py',
        HERE / 'exact_reference.py',
        HERE / 'noise_model.py',
        HERE / 'inspect_comparison.py',
        HERE / 'compare_versions.py',
        h.SOURCE,
    ]
    result = {
        str(p.relative_to(h.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths
    }
    if backend == 'native':
        path = Path(h.load_detector(backend)._rangemedian.__file__)
        result['native_extension'] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def check_frozen(frozen, hashes, backend):
    if (
        frozen['source_hashes'] != hashes
        or frozen['design'] != DESIGN
        or frozen['backend'] != backend
    ):
        raise ValueError('Frozen development protocol, code, design, or backend changed')


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


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
            parser.error('heldout requires --frozen from a completed development run')
        frozen = json.loads(args.frozen.read_text())
        check_frozen(frozen, hashes, args.backend)
        settings = frozen['selected']
    else:
        settings = configurations()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {
        'phase': args.phase,
        'source_hashes': hashes,
        'design': DESIGN,
        'backend': args.backend,
        'settings': settings,
        'python': platform.python_version(),
        'platform': platform.platform(),
        'git_revision': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=h.ROOT, text=True
        ).strip(),
        'frozen_sha256': hashlib.sha256(args.frozen.read_bytes()).hexdigest()
        if args.frozen
        else None,
    }
    dump(args.output / 'manifest.json', manifest)
    rows = []
    with (
        (args.output / 'inputs.jsonl').open('w') as inputs,
        (args.output / 'records.jsonl').open('w') as records,
    ):
        for index, case in enumerate(cases(args.phase)):
            row = evaluate(case, settings, args.backend)
            inputs.write(json.dumps(case, allow_nan=False) + '\n')
            records.write(json.dumps(row, allow_nan=False) + '\n')
            inputs.flush()
            records.flush()
            rows.append(row)
            if (index + 1) % 20 == 0:
                print(f'{args.phase}: {index + 1} histories', flush=True)
    summary = summarize(rows)
    dump(args.output / 'summary.json', summary)
    breakdown = {}
    for field in ('n', 'change', 'noise', 'correlation', 'location'):
        breakdown[field] = {
            str(value): summarize([r for r in rows if r['case'][field] == value])
            for value in sorted({r['case'][field] for r in rows})
        }
    dump(args.output / 'breakdown.json', breakdown)
    if args.phase == 'development':
        dump(
            args.output / 'frozen.json',
            {
                'source_hashes': hashes,
                'design': DESIGN,
                'backend': args.backend,
                'selected': select_settings(summary, settings),
                'selection_rule': 'min false_rate + miss_rate; then false_rate, -penalty, -floor, half_life',
            },
        )


if __name__ == '__main__':
    main()
