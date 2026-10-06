"""Reproduce a mathematical diagnosis of one archived missed slowdown.

No new histories, fitted tuning parameters, or reporting decisions are generated.
See ar1_information.rst. Profile likelihoods allow unrestricted plateau levels,
so they also upper-bound the model with positive levels.
"""

import gzip
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import nnls
from scipy.special import logsumexp

from . import rational_polynomial as p
from . import reporting_ar1 as a

HERE = Path(__file__).parent
CASE_ID = 'independent_constant-ar1-study-n100-d0.08-s0.005-early-seed1400'


def conditional_sse(values, split, training_count, rho):
    """Exact conditional nuisance fit; rho=1 denotes its profiled limit.

    Write c=(1-rho)*a and d=b-a. The innovation mean is c+d*v,
    where v is 1 at the change and 1-rho afterward. This parameterization
    retains the limiting common drift that a direct design at rho=1 loses.
    """
    y = np.asarray(values, dtype=float)
    indices = np.arange(training_count, len(y))
    z = y[training_count:] - rho * y[training_count - 1 : -1]
    if split < training_count:
        design = np.ones((len(z), 1))
    else:
        v = (indices == split) + (1 - rho) * (indices > split)
        design = np.column_stack((np.ones(len(z)), v))
    residual = z - design @ np.linalg.lstsq(design, z, rcond=None)[0]
    return float(residual @ residual)


def relaxation_gap(values, split, training_count, rho):
    """Additional SSE from linking the three relaxed innovation means."""
    y = np.asarray(values, dtype=float)
    if not training_count < split < len(y) - 1:
        raise ValueError('Both ordinary innovation groups must be nonempty')
    z = y[training_count:] - rho * y[training_count - 1 : -1]
    left, right = z[: split - training_count], z[split - training_count + 1 :]
    crossing = z[split - training_count]
    numerator = ((1 - rho) * crossing + rho * left.mean() - right.mean()) ** 2
    denominator = (1 - rho) ** 2 + rho**2 / len(left) + 1 / len(right)
    return float(numerator / denominator)


def conditional_log_maximum(sse, count):
    if sse == 0:
        return math.inf
    return -count / 2 * (math.log(2 * math.pi) + 1 + math.log(sse / count))


def full_profile(values, split, rho):
    """Stationary full-history likelihood, profiled over levels and variance."""
    if not -1 < rho < 1:
        raise ValueError('Require stationary correlation')
    y = np.asarray(values, dtype=float)
    n = len(y)
    design = np.column_stack((np.ones(n), np.arange(n) >= split))
    first_factor = math.sqrt((1 - rho) * (1 + rho))
    transformed_y = np.r_[first_factor * y[0], y[1:] - rho * y[:-1]]
    transformed_design = np.vstack((first_factor * design[0], design[1:] - rho * design[:-1]))
    beta = np.linalg.lstsq(transformed_design, transformed_y, rcond=None)[0]
    residual = transformed_y - transformed_design @ beta
    sse = float(residual @ residual)
    log_maximum = 0.5 * (math.log1p(-rho) + math.log1p(rho)) + conditional_log_maximum(sse, n)
    return {'sse': sse, 'levels': [float(beta[0]), float(beta.sum())], 'log_u': log_maximum}


def load_example():
    """Read the existing observation, checking its archived source hash."""
    index = json.loads((HERE / 'data/ar1_v1_results.json').read_text())
    for item in index['archives']:
        if item['kind'] != 'inputs' or item['condition'] != 'independent_constant':
            continue
        path = HERE / 'data' / item['file']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != item['sha256']:
            raise ValueError(f'Archive hash mismatch: {path}')
        with gzip.open(path, 'rt') as stream:
            for line in stream:
                case = json.loads(line)
                if case['id'] == CASE_ID:
                    return case, {'file': item['file'], 'sha256': digest}
    raise ValueError('Archived example missing')


def direction_diagnosis(values, split, true_means, true_sigma):
    y = np.asarray(values, dtype=float)
    n = len(y)
    training_count = n // 2
    m = n - training_count
    scale = max(abs(y[:training_count])) or 1.0
    normalized = y / scale
    info = a.confidence_threshold(normalized, 0.005)
    log_q = info['log_predictive_density'] - m * math.log(scale)
    residual = a.confidence_residual(normalized, split, training_count)
    rows = []
    for rho in [0.0, 0.9, 0.99, 4095 / 4096, 1.0]:
        relaxed = float(p.evaluate(residual, Fraction(rho))) * scale**2
        exact = conditional_sse(y, split, training_count, rho)
        rows.append(
            {
                'rho': rho,
                'relaxed_sse': relaxed,
                'exact_sse': exact,
                'relaxed_log_e': log_q - conditional_log_maximum(relaxed, m),
                'exact_log_e': log_q - conditional_log_maximum(exact, m),
            }
        )
        if rho == 4095 / 4096:
            design = np.column_stack((np.arange(n) < split, np.arange(n) >= split)).astype(float)
            conditional = design[training_count:] - rho * design[training_count - 1 : -1]
            z = y[training_count:] - rho * y[training_count - 1 : -1]
            levels, norm = nnls(conditional, z)
            rows[-1]['nonnegative_levels'] = levels.tolist()
            rows[-1]['nonnegative_sse'] = float(norm**2)
            rows[-1]['nonnegative_log_e'] = log_q - conditional_log_maximum(float(norm**2), m)
    # The known generating density is a diagnostic benchmark for this example's
    # prediction loss, not a predictor available to a practical detector.
    true_residual = y[training_count:] - np.asarray(true_means)[training_count:]
    oracle_log_q = -0.5 * (
        m * math.log(2 * math.pi * true_sigma**2)
        + float(true_residual @ true_residual) / true_sigma**2
    )
    gain = conditional_log_maximum(rows[0]['exact_sse'], m) - conditional_log_maximum(
        rows[-1]['exact_sse'], m
    )
    regret = conditional_log_maximum(rows[0]['exact_sse'], m) - log_q
    return {
        'split': split,
        'predictor': info['predictor'],
        'normalization_scale': scale,
        'log_q_original_units': log_q,
        'residual_bound_original_units': info['residual_bound'] * scale**2,
        'exclusion_log_cutoff': math.log(200),
        'exact_near_one_likelihood_gain_from_rho_zero': gain,
        'prediction_regret_against_rho_zero_profile': regret,
        'generating_density_log_q_diagnostic_only': oracle_log_q,
        'generating_density_near_one_log_e_diagnostic_only': (
            oracle_log_q - conditional_log_maximum(rows[-1]['exact_sse'], m)
        ),
        'rows': rows,
    }


def diagnose():
    case, archive = load_example()
    y = np.array(case['values'])
    split, n = case['position'], len(y)
    means = np.where(np.arange(n) < split, 10.0, 10.8)
    forward = direction_diagnosis(y, split, means, 0.05)
    reverse = direction_diagnosis(y[::-1], n - split, means[::-1], 0.05)
    # A single set of nuisance parameters also survives both directions.
    # Choosing a denominator witness from the observations is allowed: this
    # is an existence check, not a data-fitted predictive numerator.
    witness_rho = 4095 / 4096
    costs = []
    for errors in [y - means, (y - means)[::-1]]:
        innovations = errors[n // 2 :] - witness_rho * errors[n // 2 - 1 : -1]
        costs.append(float(innovations @ innovations))
    tau_squared = sum(costs) / n
    shared_logs = [
        direction['log_q_original_units']
        + 0.5 * (n // 2 * math.log(2 * math.pi * tau_squared) + cost / tau_squared)
        for direction, cost in zip([forward, reverse], costs, strict=True)
    ]
    averaged = []
    for left, right in zip(forward['rows'], reverse['rows'], strict=True):
        averaged.append(
            {
                'rho': left['rho'],
                'relaxed_log_e_average': float(
                    logsumexp([left['relaxed_log_e'], right['relaxed_log_e']]) - math.log(2)
                ),
                'exact_log_e_average': float(
                    logsumexp([left['exact_log_e'], right['exact_log_e']]) - math.log(2)
                ),
            }
        )
    differences = np.diff(y)
    limiting_sse = float(np.delete(differences, split - 1) @ np.delete(differences, split - 1))
    return {
        'description': 'Mathematical diagnosis of one archived history; no new evaluation.',
        'case_id': CASE_ID,
        'source_archive': archive,
        'forward': forward,
        'reverse': reverse,
        'averaged_evidence': averaged,
        'average_exclusion_log_cutoff': math.log(100),
        'shared_nuisance_witness': {
            'rho': witness_rho,
            'levels_original_order': [10.0, 10.8],
            'innovation_variance': tau_squared,
            'marginal_variance': tau_squared / (1 - witness_rho**2),
            'forward_log_e': shared_logs[0],
            'reverse_log_e': shared_logs[1],
            'average_log_e': float(logsumexp(shared_logs) - math.log(2)),
        },
        'full_history': {
            'limiting_sse_at_one': limiting_sse,
            'profiles': [
                {'rho': rho, **full_profile(y, split, rho)}
                for rho in [0, 0.9, 0.99, 0.9999, 0.999999, 1 - 1e-10]
            ],
        },
    }


if __name__ == '__main__':
    print(json.dumps(diagnose(), indent=2, allow_nan=False))
