"""Post-evaluation diagnosis of detections lost relative to indexed-jump.

Numerical tail estimates guide subsequent mathematics; they never change the
frozen decisions and are not rigorous probability or interval certificates.
"""

import argparse
import gzip
import hashlib
import json
import platform
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy

from . import directional_tail as tail


def diagnose(index_path):
    index = json.loads(index_path.read_text())
    inputs, records, sources = {}, {}, {}
    for archive in index['archives']:
        kind = archive['kind']
        if kind not in ('inputs', 'records'):
            continue
        path = index_path.parent / archive['file']
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != archive['sha256']:
            raise ValueError(f'Archive hash mismatch: {path}')
        for line in gzip.decompress(raw).splitlines():
            value = json.loads(line)
            key = value['id'] if kind == 'inputs' else value['case']['id']
            target = inputs if kind == 'inputs' else records
            assert key not in target
            target[key] = value
            sources.setdefault(key, {})[kind] = {
                'file': archive['file'],
                'sha256': archive['sha256'],
            }
    assert inputs.keys() == records.keys() and len(inputs) == 360
    selected = []
    for key, row in records.items():
        if not (
            row['case']['positive']
            and row['methods']['indexed_jump_gate']
            and not row['methods']['directional_tail_gate']
        ):
            continue
        result = row['directional_tail']
        item = {'case': row['case'], 'archives': sources[key], 'status': result['status']}
        if result['status'] == 'surviving_explanation':
            witness = result['witness']
            delta = Fraction(result['confidence']['delta'])
            bounds = witness['confidence_bounds']
            uniform = Fraction(bounds['uniform_p_upper_squared'])
            tail_squared = (
                Fraction(bounds['tail']['p_upper_squared']) if bounds['tail'] else Fraction(1)
            )
            assert not bounds['excluded'] and min(uniform, tail_squared) >= delta**2
            estimated = tail.approximate_tail(
                inputs[key]['values'], witness['split'], float(Fraction(witness['rho']))
            )
            item.update(
                witness={'split': witness['split'], 'rho': witness['rho']},
                confidence_cutoff=float(delta),
                saved_upper_bound_squared=str(min(uniform, tail_squared)),
                numerical_tail=estimated,
                estimated_below_cutoff=estimated['probability'] < float(delta),
            )
        selected.append(item)
    expected = index['comparisons']['all']['directional_tail_vs_indexed_jump_gate'][
        'positive_lost'
    ]
    assert len(selected) == expected
    root = Path(__file__).parent
    return {
        'purpose': 'Post-evaluation development diagnosis; frozen decisions unchanged',
        'numerical_qualification': 'Imhof inversion estimates with heuristic quadrature errors, not rigorous certificates',
        'result_index_sha256': hashlib.sha256(index_path.read_bytes()).hexdigest(),
        'environment': {
            'python': platform.python_version(),
            'numpy': np.__version__,
            'scipy': scipy.__version__,
        },
        'source_hashes': {
            name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in ('directional_losses.py', 'directional_tail.py', 'residual_direction.py')
        },
        'lost_histories': selected,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = diagnose(args.index)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
