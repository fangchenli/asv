"""Full-history predictor with separate midpoint and jump priors.

The original full-history density receives half the mixture weight. The other
half predicts at the candidate location with separate midpoint and jump priors.
Controls isolate the prior change from candidate-specific prediction.
See baseline_jump_prior.rst. Settings are development choices, not calibrated.
"""

import math

import numpy as np
from scipy.special import gammaln, logsumexp

from . import reporting_ar1_full as full

JUMP_SCALES = (1, 2, 4, 8, 16, 32, 64)
PRIOR = {
    'location_model': 'candidate_specific',
    'level_model': 'midpoint_and_jump',
    'original_density_weight': 0.5,
    'midpoint_precision': 2 * full.PRIOR['mean_precision'],
    'jump_scales_in_noise_sd': list(JUMP_SCALES),
    'original_prior': full.PRIOR,
}


def level_terms(values, split, jump_precision):
    """Ridge residual and determinant from a proper Gaussian level prior."""
    y = np.asarray(values, dtype=float)
    design = np.column_stack((np.ones(len(y)), (np.arange(len(y)) >= split) - 0.5))
    precision = np.diag([PRIOR['midpoint_precision'], jump_precision])
    gram = design.T @ design + precision
    coefficients = np.linalg.solve(gram, design.T @ y)
    residual = y - design @ coefficients
    cost = float(residual @ residual + coefficients @ precision @ coefficients)
    log_volume = 0.5 * (float(np.linalg.slogdet(precision)[1]) - float(np.linalg.slogdet(gram)[1]))
    return cost, log_volume


def noise_log_density(n, cost, beta):
    """Integrated variance factor; useful separately for the cost accounting."""
    alpha = full.PRIOR['shape']
    return float(
        -n / 2 * math.log(2 * math.pi)
        + gammaln(alpha + n / 2)
        - gammaln(alpha)
        + alpha * math.log(beta)
        - (alpha + n / 2) * math.log(beta + cost / 2)
    )


def component_log_density(values, split, beta, jump_precision):
    cost, log_volume = level_terms(values, split, jump_precision)
    return log_volume + noise_log_density(len(values), cost, beta)


def noise_mixture(n, cost):
    terms = [
        noise_log_density(n, cost, 0.5 * 2.0 ** (2 * j)) for j in full.PRIOR['scale_exponents']
    ]
    return float(logsumexp(terms) - math.log(len(terms)))


def partition_log_density(values, split, jump_precision):
    cost, log_volume = level_terms(values, split, jump_precision)
    return log_volume + noise_mixture(len(values), cost)


def global_log_density(values):
    """Average two proper full-history densities; retain all mixture weights."""
    n = len(values)
    original = full.predictive_log_density(values)
    with np.errstate(over='ignore', invalid='ignore'):
        try:
            # The no-change part, variance mixture, and location weights stay
            # exactly as in the original predictor.
            no_change = float(
                logsumexp(
                    [
                        full.component_log_density(values, None, 0.5 * 2.0 ** (2 * j))
                        for j in full.PRIOR['scale_exponents']
                    ]
                )
                - math.log(len(full.PRIOR['scale_exponents']))
            )
            alternatives = [no_change + math.log(0.5)]
            for split in range(1, n):
                jump_mixture = logsumexp(
                    [partition_log_density(values, split, 1 / scale**2) for scale in JUMP_SCALES]
                ) - math.log(len(JUMP_SCALES))
                alternatives.append(float(jump_mixture) + math.log(0.5 / (n - 1)))
            changed = float(logsumexp(alternatives))
            result = float(logsumexp([original, changed]) - math.log(2))
        except (OverflowError, ValueError, np.linalg.LinAlgError):
            return -math.inf
    return result if math.isfinite(result) else -math.inf


def _indexed_log_density(values, split, original, *, original_levels=False):
    with np.errstate(over='ignore', invalid='ignore'):
        try:
            if original_levels:
                changed = partition_log_density(values, split, full.PRIOR['mean_precision'] / 2)
            else:
                changed = float(
                    logsumexp(
                        [
                            partition_log_density(values, split, 1 / scale**2)
                            for scale in JUMP_SCALES
                        ]
                    )
                    - math.log(len(JUMP_SCALES))
                )
            result = float(logsumexp([original, changed]) - math.log(2))
        except (OverflowError, ValueError, np.linalg.LinAlgError):
            return -math.inf
    return result if math.isfinite(result) else -math.inf


def predictive_log_density(values, split):
    """A proper density for each fixed candidate location; no location search."""
    if not 1 <= split < len(values):
        raise ValueError('Require an interior candidate location')
    return _indexed_log_density(values, split, full.predictive_log_density(values))


def global_evidence(values, *, unit, config=None, max_cells=4096, max_depth=16):
    """Development control that retains the original uniform location mixture."""
    return full._evidence(
        values,
        unit=unit,
        config=config,
        max_cells=max_cells,
        max_depth=max_depth,
        predictor=global_log_density,
        prior={**PRIOR, 'location_model': 'uniform_mixture'},
    )


def _indexed_evidence(values, *, unit, config, max_cells, max_depth, original_levels):
    original = None

    def predictor(y, split):
        nonlocal original
        if original is None:
            original = full.predictive_log_density(y)
        return _indexed_log_density(y, split, original, original_levels=original_levels)

    return full._evidence(
        values,
        unit=unit,
        config=config,
        max_cells=max_cells,
        max_depth=max_depth,
        predictor=predictor,
        prior={
            'location_model': 'candidate_specific',
            'level_model': 'original_independent',
            'original_density_weight': 0.5,
            'original_prior': full.PRIOR,
        }
        if original_levels
        else PRIOR,
        indexed=True,
    )


def original_location_evidence(values, *, unit, config=None, max_cells=4096, max_depth=16):
    """Control: indexed prediction with the original independent-level prior."""
    return _indexed_evidence(
        values,
        unit=unit,
        config=config,
        max_cells=max_cells,
        max_depth=max_depth,
        original_levels=True,
    )


def evidence(values, *, unit, config=None, max_cells=4096, max_depth=16):
    return _indexed_evidence(
        values,
        unit=unit,
        config=config,
        max_cells=max_cells,
        max_depth=max_depth,
        original_levels=False,
    )
