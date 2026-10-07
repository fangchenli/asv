"""Reproduce full-determinant certificates on the two archived lost detections."""

import argparse
import gzip
import hashlib
import json
import math
import platform
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import t as student_t

from . import directional_determinant as determinant
from . import residual_direction as direction

HERE = Path(__file__).parent
HALF_WIDTH = F(1, 4096)


def size_diagnostic(values, split, rho, directional_probability):
    """Floating GLS and independent-p product calculation; not certificates."""
    y = np.array(values, dtype=float)
    n = len(y)
    design = np.column_stack((np.arange(n) < split, np.arange(n) >= split)).astype(float)
    covariance = rho ** np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    precision_design = np.linalg.solve(covariance, design)
    gram_inverse = np.linalg.inv(design.T @ precision_design)
    fitted = gram_inverse @ precision_design.T @ y
    residual = y - design @ fitted
    residual_ss = float(residual @ np.linalg.solve(covariance, residual))
    contrast = np.array([-1.05, 1.0])
    excess = float(contrast @ fitted)
    standard_error = math.sqrt(residual_ss / (n - 2) * float(contrast @ gram_inverse @ contrast))
    statistic = excess / standard_error
    size_p = float(student_t.sf(statistic, n - 2))
    product = directional_probability * size_p
    return {
        'estimated_earlier': float(fitted[0]),
        'estimated_later': float(fitted[1]),
        'estimated_excess_over_five_percent': excess,
        'standard_error': standard_error,
        't_statistic': statistic,
        'size_p': size_p,
        'product_combination_p': product * (1 - math.log(product)) if product else 0.0,
        'joint_budget': 0.042,
    }


def diagnose():
    data = HERE / 'data'
    source = data / 'directional_v1_loss_diagnosis.json'
    saved = json.loads(source.read_text())
    index = data / 'directional_v1_results.json'
    if hashlib.sha256(index.read_bytes()).hexdigest() != saved['result_index_sha256']:
        raise ValueError('Result index hash mismatch')
    results = []
    for item in saved['lost_histories']:
        rows = {}
        for kind in ('inputs', 'records'):
            archive = item['archives'][kind]
            raw = (data / archive['file']).read_bytes()
            if hashlib.sha256(raw).hexdigest() != archive['sha256']:
                raise ValueError(f'Archive hash mismatch: {archive["file"]}')
            entries = [json.loads(line) for line in gzip.decompress(raw).splitlines()]
            rows[kind] = next(
                row
                for row in entries
                if (row['id'] if kind == 'inputs' else row['case']['id']) == item['case']['id']
            )
        witness = rows['records']['directional_tail']['witness']
        split, rho = item['witness']['split'], F(item['witness']['rho'])
        if witness['split'] != split or F(witness['rho']) != rho:
            raise ValueError('Archived witness mismatch')
        values = rows['inputs']['values']
        model = direction.state(values, split)
        results.append(
            {
                'id': item['case']['id'],
                'archives': item['archives'],
                'witness': item['witness'],
                'previous_bound_squared': item['saved_upper_bound_squared'],
                'point_certificate': determinant.certify_interval(model, rho, rho),
                'interval_certificate': determinant.certify_interval(
                    model, rho - HALF_WIDTH, rho + HALF_WIDTH
                ),
                'saved_numerical_tail': item['numerical_tail'],
                'size_diagnostic': size_diagnostic(
                    values, split, float(rho), item['numerical_tail']['probability']
                ),
            }
        )
    return {
        'purpose': 'Mathematical development on two saved witnesses; no new detection counts',
        'numerical_qualification': 'Size and product p-values use floating arithmetic and the saved Imhof estimate; only determinant certificates are rigorous upper bounds',
        'interval_half_width': str(HALF_WIDTH),
        'loss_diagnosis_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in (
                'determinant_diagnosis.py',
                'directional_determinant.py',
                'directional_tail.py',
                'residual_direction.py',
                'rational_polynomial.py',
                'reporting_ar1.py',
                'reporting_ar1_full.py',
            )
        },
        'environment': {
            'python': platform.python_version(),
            'numpy': np.__version__,
            'scipy': scipy.__version__,
        },
        'cases': results,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = diagnose()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
