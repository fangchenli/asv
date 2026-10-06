"""Freeze sign-only calibration, then evaluate reporting on fresh histories."""

import argparse
import hashlib
import itertools
import json
import math
import platform
import random
import statistics
import subprocess
from fractions import Fraction
from pathlib import Path

from . import ablation_study as a
from . import harness as h
from . import reporting_signs as s
from . import threshold_study as t
from .reporting_reference import single_change_evidence

HERE = Path(__file__).parent
DESIGN = {
    'sizes': [40, 100],
    'changes': [0, 0.04, 0.05, 0.06, 0.08],
    'noise': [0.005, 0.02, 0.05],
    'locations': ['middle', 'recent'],
    'conditions': ['gaussian', 'laplace'],
    'seeds': list(range(600, 620)),
}
CALIBRATION = {'draws': 8191, 'alpha': '1/20', 'total_failure': '1/100', 'seed_base': 550}


def source_hashes(backend):
    hashes = a.source_hashes(backend)
    for name in (
        'reporting_reference.py',
        'reporting_signs.py',
        'reporting_study.py',
        'reporting_protocol.rst',
    ):
        path = HERE / name
        hashes[str(path.relative_to(h.ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def freeze(inherited, backend):
    a.check_frozen(inherited, a.source_hashes(backend), backend)
    cfg = next(
        row['config']
        for row in inherited['selections']
        if row['family'] == 'shared-r1' and row['budget'] == 0.05
    )
    calibrations = {}
    for n in DESIGN['sizes']:
        calibrations[str(n)] = s.calibrate(
            n,
            draws=CALIBRATION['draws'],
            seed=CALIBRATION['seed_base'] + n,
            alpha=Fraction(CALIBRATION['alpha']),
            failure=Fraction(CALIBRATION['total_failure']) / len(DESIGN['sizes']),
        )
        print(f'Calibrated {n} observations', flush=True)
    return {
        'design': DESIGN,
        'calibration_design': CALIBRATION,
        'backend': backend,
        'source_hashes': source_hashes(backend),
        'detector': cfg,
        'calibrations': calibrations,
        'python': platform.python_version(),
    }


def check_frozen(frozen, backend):
    if (
        frozen['design'] != DESIGN
        or frozen['calibration_design'] != CALIBRATION
        or frozen['source_hashes'] != source_hashes(backend)
        or frozen['backend'] != backend
    ):
        raise ValueError('Frozen design, calibration, source, or backend changed')


def make_case(condition, n, change, noise, location, seed):
    identity = f'{condition}-n{n}-d{change}-s{noise}-{location}-seed{seed}'
    # Distinct random streams per complete case; no duplicated flat/location
    # histories and no shared draws across change sizes or noise families.
    rng = random.Random(int.from_bytes(hashlib.sha256(identity.encode()).digest(), 'big'))
    position = n // 2 if location == 'middle' else n - max(4, n // 10)
    values = []
    for i in range(n):
        if condition == 'gaussian':
            error = rng.gauss(0, 1)
        elif condition == 'laplace':
            error = rng.expovariate(math.sqrt(2)) * (1 if rng.getrandbits(1) else -1)
        else:
            raise ValueError('Unknown condition')
        values.append(10 * (1 + change * (i >= position)) + 10 * noise * error)
    return {
        'id': identity,
        'condition': condition,
        'n': n,
        'change': change,
        'noise': noise,
        'location': location,
        'seed': seed,
        'position': position,
        'positive': change > 0.05,
        'values': values,
    }


def cases():
    for args in itertools.product(
        *(DESIGN[key] for key in ('conditions', 'sizes', 'changes', 'noise', 'locations', 'seeds'))
    ):
        yield make_case(*args)


def evaluate(case, frozen, backend):
    module = h.load_detector(backend)
    y = case['values']
    production, pool = t.candidate_pool(module, y, backend)
    scale = statistics.median(abs(b - a) for a, b in zip(y, y[1:]))
    cfg = frozen['detector']
    scores = [
        a.candidate_score(fit, t.innovations(module, y, fit, cfg['cap']), len(y), scale, cfg)
        for fit in pool
    ]
    fit = pool[min(range(len(pool)), key=scores.__getitem__)]
    old = t.assess(module, case, fit)
    reference = single_change_evidence(y, alpha=Fraction(CALIBRATION['alpha']))
    joint = s.evidence(y, frozen['calibrations'][str(len(y))])
    return {
        'case': {k: v for k, v in case.items() if k != 'values'},
        'shared_fit': fit,
        'shared_report': old,
        'reference': reference,
        'joint': joint,
        'methods': {
            'production': t.assess(module, case, production)['has_alert'],
            'shared_existing': old['has_alert'],
            'reference_alone': reference['has_alert'],
            'joint_alone': joint['has_alert'],
            'reference_gate': old['has_alert'] and reference['has_alert'],
            'joint_gate': old['has_alert'] and joint['has_alert'],
        },
    }


def summarize(rows):
    result = {}
    for method in rows[0]['methods']:
        null = [r for r in rows if not r['case']['positive']]
        positive = [r for r in rows if r['case']['positive']]
        equal = [r for r in rows if r['case']['change'] == 0.05]
        false = sum(r['methods'][method] for r in null)
        misses = sum(not r['methods'][method] for r in positive)
        result[method] = {
            'null_histories': len(null),
            'positive_histories': len(positive),
            'false_alerts': false,
            'misses': misses,
            'false_rate': false / len(null) if null else None,
            'miss_rate': misses / len(positive) if positive else None,
            'at_threshold_histories': len(equal),
            'at_threshold_alerts': sum(r['methods'][method] for r in equal),
        }
    return result


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return '-infinity' if value == -math.inf else 'infinity' if value == math.inf else 'NaN'
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    return value


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
            outputs.write(json.dumps(json_safe(row), allow_nan=False) + '\n')
            inputs.flush()
            outputs.flush()
            rows.append(row)
            if len(rows) % 100 == 0:
                print(f'Evaluated {len(rows)} histories', flush=True)
    t.dump(args.output / 'summary.json', summarize(rows))
    t.dump(
        args.output / 'breakdown.json',
        {
            key: {
                str(value): summarize([r for r in rows if r['case'][key] == value])
                for value in sorted({r['case'][key] for r in rows})
            }
            for key in ('condition', 'n', 'change', 'noise', 'location')
        },
    )


if __name__ == '__main__':
    main()
