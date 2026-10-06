"""Decompose prediction cost and isolate prior/location changes on saved inputs."""

import json
import math

import numpy as np

from . import ar1_information as info
from . import reporting_ar1_full as full
from . import reporting_ar1_jump as jump
from . import reporting_covariance as gls


def cost_decomposition(values, split):
    y = np.asarray(values, dtype=float)
    n = len(y)
    rss = sum(
        float((block - block.mean()) @ (block - block.mean())) for block in [y[:split], y[split:]]
    )
    cost, volume = jump.level_terms(y, split, full.PRIOR['mean_precision'] / 2)
    maximum = info.conditional_log_maximum(rss, n)
    noise_terms = [
        jump.noise_log_density(n, rss, 0.5 * 2.0 ** (2 * j)) for j in full.PRIOR['scale_exponents']
    ]
    noise_rss, noise_cost = jump.noise_mixture(n, rss), jump.noise_mixture(n, cost)
    original = full.predictive_log_density(y)
    parts = {
        'variance_integration': maximum - max(noise_terms),
        'scale_mixture': max(noise_terms) - noise_rss,
        'level_shrinkage': noise_rss - noise_cost,
        'level_integration': -volume,
        'location_mixture': noise_cost + volume - original,
    }
    return {
        'rss': rss,
        'prior_adjusted_cost': cost,
        'maximum_log_likelihood': maximum,
        'original_log_density': original,
        'parts': parts,
        'sum': sum(parts.values()),
        'total_gap': maximum - original,
        'jump_components': [
            {
                'jump_scale': scale,
                'log_density_at_location': jump.partition_log_density(y, split, 1 / scale**2),
            }
            for scale in jump.JUMP_SCALES
        ],
    }


def diagnose():
    case, source = info.load_example()
    y, split = case['values'], case['position']
    old_log_q = full.predictive_log_density(y)
    methods = {
        'global_original': full.evidence,
        'global_jump': jump.global_evidence,
        'indexed_original': jump.original_location_evidence,
        'indexed_jump': jump.evidence,
    }
    log_densities = {
        'global_original': old_log_q,
        'global_jump': jump.global_log_density(y),
        'indexed_original': jump._indexed_log_density(y, split, old_log_q, original_levels=True),
        'indexed_jump': jump.predictive_log_density(y, split),
    }
    results = {}
    witness_checks = []
    for name, method in methods.items():
        result = method(y, unit=1)
        results[name] = {
            'log_density_at_generating_location': log_densities[name],
            'result': result,
        }
        if result['status'] == 'surviving_explanation':
            t, rho = result['witness']['split'], result['witness']['rho_float']
            confidence = result['confidence']
            witness_log_q = (
                confidence['by_split'][str(t)]['log_predictive_density']
                if 'by_split' in confidence
                else confidence['log_predictive_density']
            )
            covariance = rho ** np.abs(np.arange(len(y))[:, None] - np.arange(len(y))[None, :])
            report = gls.evidence(y, covariance, result['calibration'])
            witness_checks.append(
                {
                    'method': name,
                    'split': t,
                    'rho': rho,
                    'log_e': witness_log_q - info.full_profile(y, t, rho)['log_u'],
                    'size_t': report['size_t'][t - 1],
                    'shape_f': report['lack_of_fit_f'][t - 1],
                }
            )
    fixtures = []
    for n, change in [(24, 2.0), (100, 2.0), (24, 0.4)]:
        values = [10 + change * (i >= n // 2) + 0.1 * math.sin(i * 1.7) for i in range(n)]
        fixtures.append(
            {
                'n': n,
                'change': change,
                'values': values,
                'indexed_jump': jump.evidence(values, unit=1),
            }
        )
    return {
        'description': 'Development diagnostics; no new random histories; settings fixed within comparison.',
        'case_id': case['id'],
        'source_archive': source,
        'unit': 1,
        'decomposition': cost_decomposition(y, split),
        'methods': results,
        'witness_checks': witness_checks,
        'deterministic_fixtures': fixtures,
    }


if __name__ == '__main__':
    print(json.dumps(diagnose(), indent=2, allow_nan=False))
