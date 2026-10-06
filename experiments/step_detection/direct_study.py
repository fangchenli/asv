"""Freeze an analytical direct test, then evaluate fresh independent histories."""

import argparse
import hashlib
import itertools
import json
import platform
import subprocess
from pathlib import Path

import scipy

from . import harness as h
from . import reporting_study as r
from . import threshold_study as t
from .reporting_direct import critical_values, evidence

HERE = Path(__file__).parent
DESIGN = dict(r.DESIGN, seeds=list(range(800, 820)))
CALIBRATION = {'alpha': 0.05, 'size_share': 0.8, 'threshold': 0.05}


def source_hashes(backend):
    hashes = r.source_hashes(backend)
    for name in ('reporting_direct.py', 'direct_study.py', 'direct_protocol.rst'):
        p = HERE / name
        hashes[str(p.relative_to(h.ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return hashes


def freeze(inherited, backend):
    r.check_frozen(inherited, backend)
    return {
        'design': DESIGN,
        'calibration_design': CALIBRATION,
        'backend': backend,
        'source_hashes': source_hashes(backend),
        'inherited': inherited,
        'critical_values': {str(n): critical_values(n, **CALIBRATION) for n in DESIGN['sizes']},
        'python': platform.python_version(),
        'scipy': scipy.__version__,
    }


def check_frozen(frozen, backend):
    r.check_frozen(frozen['inherited'], backend)
    if (
        frozen['design'] != DESIGN
        or frozen['calibration_design'] != CALIBRATION
        or frozen['source_hashes'] != source_hashes(backend)
        or frozen['backend'] != backend
        or frozen['scipy'] != scipy.__version__
    ):
        raise ValueError('Frozen design, calibration, source, backend, or SciPy changed')


def cases():
    for args in itertools.product(
        *(DESIGN[key] for key in ('conditions', 'sizes', 'changes', 'noise', 'locations', 'seeds'))
    ):
        yield r.make_case(*args)


def evaluate(case, frozen, backend):
    # This reuses the complete frozen fit and older evidence pipeline.
    result = r.evaluate(case, frozen['inherited'], backend)
    direct = evidence(case['values'], frozen['critical_values'][str(case['n'])])
    result['direct'] = direct
    result['methods']['direct_alone'] = direct['has_alert']
    result['methods']['direct_gate'] = direct['has_alert'] and result['methods']['shared_existing']
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
            'scipy': scipy.__version__,
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
                print(f'Evaluated {len(rows)} histories', flush=True)
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


if __name__ == '__main__':
    main()
