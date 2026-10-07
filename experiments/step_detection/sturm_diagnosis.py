"""Separate spectral and density losses and certify the saved 51/64 witness."""

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

from . import directional_spectral as spectral
from . import directional_sturm as sturm
from . import residual_direction as direction

HERE = Path(__file__).parent
HALF_WIDTH = F(1, 65536)


def numerical_decomposition(values, split, rho):
    """Floating counterfactuals; these numbers never certify rejection."""
    n = len(values)
    basis = direction.residual_basis(n, split)
    residual = basis.T @ np.asarray(values)
    residual /= np.linalg.norm(residual)
    covariance = (
        basis.T @ (float(rho) ** np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])) @ basis
    )
    radius = float(1 / (residual @ np.linalg.solve(covariance, residual)))
    eigenvalues = np.linalg.eigvalsh(covariance)
    weights = eigenvalues / radius - 1
    per_tilt = []
    for tilt in spectral.TILTS:
        moment = math.exp(-0.5 * float(np.log1p(2 * float(tilt) * weights).sum()))
        stages = {}
        for stage in ('toeplitz', 'full_sturm', 'plateau_sturm', 'projected_eigenvalues'):
            best, selected = min(1.0, moment), None
            for count in spectral.COUNTS:
                if count > n - 2:
                    continue
                eigen_lower = (
                    eigenvalues[-count]
                    if stage == 'projected_eigenvalues'
                    else float(sturm.residual_eigenvalue_bounds(n, split, count, rho, rho)[stage])
                )
                weight = eigen_lower / radius - 1
                coefficient = weight / (1 + 2 * float(tilt) * weight) if weight > 0 else 0.0
                correction = max(
                    1.0, float(tilt) * coefficient / float(spectral.density_constant(count))
                )
                candidate = min(1.0, moment / correction)
                if candidate < best:
                    best, selected = candidate, count
            stages[stage] = {'p_upper_estimate': best, 'count': selected}
        per_tilt.append({'tilt': str(tilt), 'chernoff_estimate': moment, 'stages': stages})
    return {'radius_estimate': radius, 'per_tilt': per_tilt}


def diagnose():
    summary_path = HERE / 'data/spectral_reporting_v1_summary.json'
    summary = json.loads(summary_path.read_text())
    source = summary['archive']
    raw = (HERE / 'data' / source['file']).read_bytes()
    if hashlib.sha256(raw).hexdigest() != source['sha256']:
        raise ValueError('Saved spectral replay hash mismatch')
    original = json.loads(gzip.decompress(raw))
    for name, digest in original['source_hashes'].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Original replay source changed: {name}')
    rows = []
    for row in original['rows']:
        if not row['case']['id'].startswith('negative-'):
            continue
        witness = row['result']['witness']
        split, rho = witness['split'], F(witness['rho'])
        values = row['case']['values']
        model = direction.state(values, split)
        previous = spectral.certify_interval(
            model, rho, rho, delta=F(row['result']['confidence']['delta'])
        )
        assert previous == witness['confidence_bounds']['spectral']
        rows.append(
            {
                'id': row['case']['id'],
                'witness': {'split': split, 'rho': str(rho)},
                'previous_point_bound': previous['p_upper'],
                'point_certificate': sturm.certify_interval(model, rho, rho),
                'interval_certificate': sturm.certify_interval(
                    model, rho - HALF_WIDTH, rho + HALF_WIDTH
                ),
                'numerical_decomposition': numerical_decomposition(values, split, rho),
                'saved_numerical_tail': row['numerical_witness_tail'],
            }
        )
    assert len(rows) == 2
    return {
        'purpose': 'Mathematical diagnosis on saved witnesses; no new complete-search decisions',
        'numerical_qualification': 'Counterfactual spectra and tail estimates use floating arithmetic; only rational certificates determine exclusions',
        'source_archive': source,
        'source_summary_sha256': hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        'half_width': str(HALF_WIDTH),
        'source_hashes': {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in (
                'sturm_diagnosis.py',
                'directional_sturm.py',
                'directional_spectral.py',
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
        'rows': rows,
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
