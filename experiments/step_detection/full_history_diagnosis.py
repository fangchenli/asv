"""Print diagnostics for one saved case and explicit mathematical fixtures."""

import json
import math

import numpy as np

from . import ar1_information as info
from . import reporting_ar1_full as full


def diagnose():
    case, source = info.load_example()
    y = case['values']
    log_q = full.predictive_log_density(y)
    comparisons = []
    for rho in [0, 0.9, 0.99, 4095 / 4096, 0.9999, 0.99999]:
        profile = info.full_profile(y, case['position'], rho)
        comparisons.append({'rho': rho, 'log_e': log_q - profile['log_u'], **profile})
    residual = np.array(y) - np.where(np.arange(len(y)) < case['position'], 10, 10.8)
    oracle_log_q = -0.5 * (
        len(y) * math.log(2 * math.pi * 0.05**2) + float(residual @ residual) / 0.05**2
    )
    fixtures = []
    for n in [24, 100]:
        values = [10 + 2 * (i >= n // 2) + 0.1 * math.sin(i * 1.7) for i in range(n)]
        fixtures.append(
            {
                'description': '10 to 12 with deterministic 0.1*sin(1.7*i) variation',
                'n': n,
                'values': values,
                'result': full.evidence(values, unit=1),
            }
        )
    return {
        'description': 'Development diagnosis only; no new statistical evaluation.',
        'case_id': case['id'],
        'source_archive': source,
        'unit': 1,
        'log_q': log_q,
        'oracle_log_q_diagnostic_only': oracle_log_q,
        'saved_case_profiles': comparisons,
        'saved_case_result': full.evidence(y, unit=1),
        'deterministic_fixtures': fixtures,
    }


if __name__ == '__main__':
    print(json.dumps(diagnose(), indent=2, allow_nan=False))
