"""Algebra checks and archived-witness diagnostics, not a reporting detector."""

import argparse
import gzip
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np

from . import ar1_information as info
from . import rational_polynomial as p
from . import reporting_ar1 as ar1
from . import reporting_ar1_full as full
from . import reporting_ar1_jump as jump

HERE = Path(__file__).parent
# Fixed before the first archived directional calculation; illustrative only.
COMPONENTS = (F(-9, 10), F(-1, 2), F(0), F(1, 2), F(9, 10))
DELTA = F(1, 100)
POSITIVE_CASE = 'positive-indexed-ar1-study-n100-d0.08-s0.005-early-seed1500'


def residual_basis(n, split):
    """Block Helmert contrasts: orthonormal, perpendicular to both levels."""
    if n < 3 or not 1 <= split < n:
        raise ValueError('Require n>=3 and an interior split')
    basis = np.zeros((n, n - 2))
    col = 0
    for start, length in ((0, split), (split, n - split)):
        for j in range(1, length):
            scale = math.sqrt(j * (j + 1))
            basis[start : start + j, col] = 1 / scale
            basis[start + j, col] = -j / scale
            col += 1
    return basis


def matrix_log_density(values, split, rho):
    """Independent dense evaluation of the projected Gaussian direction."""
    if not -1 < rho < 1:
        raise ValueError('Require stationary correlation')
    y = np.array(values, dtype=float)
    basis = residual_basis(len(y), split)
    y[:split] -= y[:split].mean()
    y[split:] -= y[split:].mean()
    z = basis.T @ y
    if not np.isfinite(z).all() or np.linalg.norm(z) == 0:
        raise ValueError('Require finite, nonzero residuals')
    direction = z / np.linalg.norm(z)
    distance = np.abs(np.arange(len(y))[:, None] - np.arange(len(y))[None, :])
    covariance = basis.T @ (float(rho) ** distance) @ basis
    factor = np.linalg.cholesky(covariance)
    whitened = np.linalg.solve(factor, direction)
    return float(-np.log(np.diag(factor)).sum() - len(z) / 2 * np.log(whitened @ whitened))


def state(values, split):
    """Exact quantities for the scalar identity Q=N/h and density squared."""
    y = [F(x) for x in values]
    n = len(y)
    if n < 3 or not 1 <= split < n:
        raise ValueError('Require n>=3 and an interior split')
    centered = []
    for block in (y[:split], y[split:]):
        center = sum(block) / len(block)
        centered.extend(x - center for x in block)
    ss = sum(x * x for x in centered)
    if ss == 0:
        raise ValueError('Require nonzero residuals')
    numerator, denominator = full.residual_ratio(ar1.Models(centered), split)
    cubic = p.add(
        p.add(p.polynomial([1, 1]), p.scale(p.polynomial([1, -1]), n - 2)),
        p.scale(p.polynomial([1, -3, 3, -1]), (split - 1) * (n - split - 1)),
    )
    return {'n': n, 'split': split, 'ss': ss, 'num': numerator, 'den': denominator, 'b': cubic}


def squared_density(model, rho):
    """Exact squared directional density, also allowing the positive limit."""
    rho = F(rho)
    if not -1 < rho <= 1:
        raise ValueError('Require -1<rho<=1; rho=1 denotes a limit')
    n, split = model['n'], model['split']
    q = p.evaluate(model['num'], rho) / p.evaluate(model['den'], rho)
    return (
        split
        * (n - split)
        * (1 + rho)
        / p.evaluate(model['b'], rho)
        * (model['ss'] / q) ** (n - 2)
    )


def sqrt_lower(value, bits=128):
    """Dyadic lower square root, with absolute error strictly below 2**-bits."""
    value = F(value)
    if value < 0 or not isinstance(bits, int) or bits < 0:
        raise ValueError('Require a nonnegative value and precision')
    scale = 1 << bits
    return F(math.isqrt(value.numerator * scale * scale // value.denominator), scale)


def predictor_lower(model, components=COMPONENTS):
    """Equal-weight proper mixture; bound it without transcendental arithmetic."""
    if not components or any(not -1 < r < 1 for r in components):
        raise ValueError('Require nonempty stationary mixture components')
    return sum(sqrt_lower(squared_density(model, rho)) for rho in components) / len(components)


def interval_excluded(model, left, right, q_lower, delta=DELTA):
    """Check the derived sufficient inequality, without a correlation search."""
    left, right, q_lower, delta = map(F, (left, right, q_lower, delta))
    if not -1 <= left < right <= 1 or q_lower < 0 or not 0 < delta < 1:
        raise ValueError('Invalid interval, predictor bound, or confidence budget')
    lower = min(p.bernstein(model['num'], left, right))
    upper = max(p.bernstein(model['den'], left, right))
    if lower <= 0 or upper <= 0:
        return False
    n, split = model['n'], model['split']
    return (delta * q_lower) ** 2 * p.evaluate(model['b'], right) * lower ** (n - 2) > (
        split * (n - split) * (1 + right) * model['ss'] ** (n - 2) * upper ** (n - 2)
    )


def load_positive_example():
    index_path = HERE / 'data/indexed_v1_results.json'
    index = json.loads(index_path.read_text())
    for item in index['archives']:
        if item['kind'] != 'inputs' or item['condition'] != 'positive' or item['n'] != 100:
            continue
        path = HERE / 'data' / item['file']
        content = path.read_bytes()
        if hashlib.sha256(content).hexdigest() != item['sha256']:
            raise ValueError(f'Archive hash mismatch: {path}')
        for line in gzip.decompress(content).splitlines():
            case = json.loads(line)
            if case['id'] == POSITIVE_CASE:
                return case, {'file': path.name, 'sha256': item['sha256']}
    raise ValueError('Archived positive-correlation example missing')


def log_fraction(value):
    return math.log(value.numerator) - math.log(value.denominator)


def diagnose_case(case, archive, witness):
    values, split = case['values'], case['position']
    model = state(values, split)
    lower = predictor_lower(model)
    upper = lower + F(1, 2**128)
    global_log_q = full.predictive_log_density(values)
    indexed_log_q = jump._indexed_log_density(values, split, global_log_q)
    points = sorted(
        {F(-9, 10), F(0), F(str(case['correlation'])), F(9, 10), F(99, 100), witness, F(1)}
    )
    rows = []
    for rho in points:
        squared = squared_density(model, rho)
        log_g = log_fraction(squared) / 2
        row = {
            'rho': str(rho),
            'positive_endpoint_limit': rho == 1,
            'log_direction_density': log_g,
            'uniform_log_e': -log_g,
            'mixture_lower_log_e': log_fraction(lower) - log_g,
            'uniform_excluded_exact': DELTA**2 > squared,
            'mixture_excluded_exact': (DELTA * lower) ** 2 > squared,
            'mixture_retained_exact': (DELTA * upper) ** 2 <= squared,
        }
        if rho != 1:
            row['dense_log_density'] = matrix_log_density(values, split, float(rho))
            row['indexed_jump_log_e'] = (
                indexed_log_q - info.full_profile(values, split, float(rho))['log_u']
            )
        rows.append(row)
    return {
        'id': case['id'],
        'archive': archive,
        'split': split,
        'witness': str(witness),
        'q_lower': str(lower),
        'q_upper': str(upper),
        'q_lower_log': log_fraction(lower),
        'component_log_densities': {
            str(r): log_fraction(squared_density(model, r)) / 2 for r in COMPONENTS
        },
        'rows': rows,
        'witness_to_one_interval': {
            'left': str(witness),
            'right': '1',
            'uniform_excluded_exact': interval_excluded(model, witness, F(1), F(1)),
            'mixture_excluded_exact': interval_excluded(model, witness, F(1), lower),
        },
    }


def diagnose():
    independent, source = info.load_example()
    positive, positive_source = load_positive_example()
    return {
        'purpose': 'Mathematical development diagnosis; no new histories or detection claims',
        'delta': str(DELTA),
        'mixture_components': [str(r) for r in COMPONENTS],
        'mixture_weights': ['1/5'] * len(COMPONENTS),
        'sqrt_bound_bits': 128,
        'cases': [
            diagnose_case(independent, source, F(4095, 4096)),
            diagnose_case(positive, positive_source, F(16383, 16384)),
        ],
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in (
                'residual_direction.py',
                'rational_polynomial.py',
                'reporting_ar1.py',
                'reporting_ar1_full.py',
                'reporting_ar1_jump.py',
                'ar1_information.py',
            )
        },
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = diagnose()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
