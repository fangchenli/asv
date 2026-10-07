"""Verify complete directional reporting results on saved histories and fixtures."""

import argparse
import hashlib
import json
import math
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from pathlib import Path

import numpy as np

from . import directional_tail as tail
from . import rational_polynomial as p
from . import reporting_ar1 as ar1
from . import reporting_covariance as gls
from . import reporting_directional as reporting
from . import residual_direction as direction

HERE = Path(__file__).parent


def verify(values, result):
    """Reconstruct all intervals; independently check surviving pairs with GLS."""
    n, config = len(values), result['calibration']
    assert config == ar1.calibration(
        n,
        alpha=config['total_alpha'],
        confidence_alpha=config['confidence_alpha'],
        threshold=config['threshold'],
    )
    confidence = result['confidence']
    delta = F(confidence['delta'])
    assert delta == F(config['confidence_alpha'])
    assert confidence['tilt'] == str(tail.TILT)
    assert confidence['radius_upper_bits'] == tail.BITS
    assert result['has_alert'] == (result['status'] == 'certified_alert')
    models, states = ar1.Models(values), {}
    for split, saved in confidence['by_split'].items():
        state = direction.state(values, int(split))
        assert saved['ss'] == str(state['ss'])
        assert all(saved[key] == [str(x) for x in state[key]] for key in ('num', 'den', 'b'))
        states[int(split)] = state
    coverage = {split: [] for split in range(1, n)}
    routes = Counter()
    for cell in result['certificate']:
        split, left, right = cell['split'], F(cell['left']), F(cell['right'])
        assert 1 <= split < n and -1 <= left < right <= 1
        coverage[split].append((left, right))
        route = cell['route']
        routes[route] += 1
        if route == 'direction_uniform':
            assert direction.interval_excluded(states[split], left, right, F(1), delta=delta)
        elif route == 'direction_tail':
            assert confidence['tail_enabled'] and left >= 0
            proof = tail.certify_interval(states[split], left, right, delta=delta)
            assert cell['proof'] == proof and proof['status'] == 'certified_excluded'
            assert (
                F(proof['determinant_lower']) * F(proof['density_correction']) ** 2 > 1 / delta**2
            )
        elif route == 'size':
            assert all(p.positive_on(poly, left, right) for poly in models.size(split, config))
        else:
            assert route == 'shape'
            assert p.positive_on(models.shape(split, cell['extra'], config), left, right)
    for intervals in coverage.values():
        intervals.sort()
        assert all(a[1] <= b[0] for a, b in zip(intervals[:-1], intervals[1:], strict=True))
        if result['has_alert']:
            assert intervals and intervals[0][0] == -1 and intervals[-1][1] == 1
            assert all(a[1] == b[0] for a, b in zip(intervals[:-1], intervals[1:], strict=True))
    if result['status'] == 'surviving_explanation':
        witness = result['witness']
        split, rho = witness['split'], F(witness['rho'])
        assert -1 < rho < 1
        bounds = reporting.point_bounds(
            states[split], rho, delta, use_tail=confidence['tail_enabled']
        )
        assert not bounds['excluded'] and bounds == witness['confidence_bounds']
        assert direction.squared_density(states[split], rho) >= delta**2
        distance = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        reference = gls.evidence(values, float(rho) ** distance, config)
        assert reference['status'] == 'ok'
        assert split not in reference['size_rejected_splits']
        assert split not in reference['fit_rejected_splits']
    else:
        assert result['status'] in ('certified_alert', 'unresolved', 'insufficient_variation')
    return {
        'verified': True,
        'certificate_routes': dict(routes),
        'complete_coverage': result['has_alert'],
    }


def examples():
    independent, independent_source = direction.info.load_example()
    positive, positive_source = direction.load_positive_example()
    result = [
        {'id': case['id'], 'values': case['values'], 'archive': source}
        for case, source in ((independent, independent_source), (positive, positive_source))
    ]
    for n, change in ((24, 0.2), (100, 0.2), (24, 0.04)):
        result.append(
            {
                'id': f'deterministic-n{n}-change{change}',
                'description': 'Baseline 10, middle change, deterministic 0.1*sin(1.7*i) variation',
                'values': [
                    10 * (1 + change * (i >= n // 2)) + 0.1 * math.sin(1.7 * i) for i in range(n)
                ],
            }
        )
    return result


def evaluate(job):
    case, use_tail = job
    result = reporting.evidence(case['values'], use_tail=use_tail)
    validation = verify(case['values'], result)
    return {
        'case': case,
        'method': 'directional_tail' if use_tail else 'directional_uniform',
        'result': result,
        'validation': validation,
    }


def diagnose(workers):
    jobs = [(case, use_tail) for case in examples() for use_tail in (False, True)]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for row in pool.map(evaluate, jobs):
            rows.append(row)
            print(
                f"{row['case']['id']} {row['method']}: {row['result']['status']} ({row['result']['cells_visited']} cells)",
                flush=True,
            )
    return {
        'purpose': 'Verified development examples, not a fresh sensitivity evaluation',
        'work': {'max_cells': 4096, 'max_depth': 16},
        'rows': rows,
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in (
                'reporting_directional.py',
                'directional_diagnosis.py',
                'directional_tail.py',
                'residual_direction.py',
                'rational_polynomial.py',
                'reporting_ar1.py',
                'reporting_ar1_full.py',
                'reporting_covariance.py',
                'ar1_information.py',
            )
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    if args.output.exists():
        raise FileExistsError(args.output)
    result = diagnose(args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
