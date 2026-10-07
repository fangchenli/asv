"""Complete reporting reference adding the spectral spectral bound.

Frozen reporting sources remain intact; this fork preserves their search order
until a new certificate removes an interval or a saved point explanation.
"""

import math
from fractions import Fraction as F

from . import directional_spectral as spectral
from . import directional_tail as tail
from . import rational_polynomial as p
from . import reporting_ar1 as ar1
from . import reporting_ar1_full as full
from . import residual_direction as direction


def confidence_state(models, split):
    """Reuse exact reporting fits; Q(0) is the within-plateau residual sum."""
    num, den = full.residual_ratio(models, split)
    gram, _, _, _ = models.fit((split,))
    return {
        'n': models.n,
        'split': split,
        'ss': p.evaluate(num, F(0)) / p.evaluate(den, F(0)),
        'num': num,
        'den': den,
        'b': p.remove_stationary_endpoint_factors(gram),
    }


def point_bounds(model, rho, delta, *, use_tail=True, use_spectral=True):
    """All bounds majorize the same p-value; failure to exclude is conservative."""
    uniform_squared = min(F(1), direction.squared_density(model, rho))
    certificate = (
        tail.certify_interval(model, rho, rho, delta=delta) if use_tail and rho >= 0 else None
    )
    tail_squared = F(certificate['p_upper_squared']) if certificate else F(1)
    spectral_proof = (
        spectral.certify_interval(model, rho, rho, delta=delta) if use_spectral else None
    )
    spectral_squared = F(spectral_proof['p_upper_squared']) if spectral_proof else F(1)
    return {
        'spectral': spectral_proof,
        'uniform_p_upper_squared': str(uniform_squared),
        'tail': certificate,
        'excluded': min(uniform_squared, tail_squared, spectral_squared) < delta**2,
    }


def evidence(
    values, config=None, *, max_cells=4096, max_depth=16, use_tail=True, use_spectral=True
):
    """Certify every location and rho in (-1,1), or return a non-alert.

    A surviving explanation belongs to the implemented conservative confidence
    set. It need not survive an ideal, exactly evaluated directional tail test.
    Disabling the spectral reproduces the inherited search. Disabling
    both tail and spectral bounds gives the uniform-direction control.
    """
    y = list(values)
    n = len(y)
    if not 8 <= n <= 200 or any(not math.isfinite(x) for x in y):
        raise ValueError('Use 8 to 200 finite observations')
    if not isinstance(use_tail, bool) or not isinstance(use_spectral, bool):
        raise ValueError('Bound switches must be boolean')
    config = ar1.calibration(n) if config is None else dict(config)
    expected = ar1.calibration(
        n,
        alpha=config['total_alpha'],
        confidence_alpha=config['confidence_alpha'],
        threshold=config['threshold'],
    )
    if config != expected:
        raise ValueError('Use an unmodified matching calibration')
    for count, minimum in ((max_cells, 1), (max_depth, 0)):
        if not isinstance(count, int) or isinstance(count, bool) or count < minimum:
            raise ValueError('Invalid certification budget')
    # Never enlarge the supplied binary floating-point confidence allocation.
    delta = F(config['confidence_alpha'])
    models = ar1.Models(y)
    certificates, states = [], {}
    visited = 0

    def result(status, witness=None):
        return {
            'has_alert': status == 'certified_alert',
            'status': status,
            'witness': witness,
            'cells_visited': visited,
            'certificate': certificates,
            'confidence': {
                'rule': 'directional_full_spectral',
                'spectral_enabled': use_spectral,
                'spectral_bits': spectral.determinant.BITS,
                'spectral_tilts': [str(tilt) for tilt in spectral.TILTS],
                'spectral_counts': list(spectral.COUNTS),
                'tail_enabled': use_tail,
                'delta': str(delta),
                'tilt': str(tail.TILT),
                'radius_upper_bits': tail.BITS,
                'by_split': states,
            },
            'calibration': config,
        }

    try:
        for split in range(1, n):
            model = confidence_state(models, split)
            states[str(split)] = {
                'ss': str(model['ss']),
                **{key: [str(x) for x in model[key]] for key in ('num', 'den', 'b')},
            }
            size = models.size(split, config)
            shapes = {}
            pending = [(F(-1), F(1), 0)]
            while pending:
                left, right, depth = pending.pop()
                if visited >= max_cells:
                    return result(
                        'unresolved',
                        {
                            'split': split,
                            'left': str(left),
                            'right': str(right),
                            'reason': 'max_cells',
                        },
                    )
                visited += 1
                cell = {'split': split, 'left': str(left), 'right': str(right), 'extra': None}
                if direction.interval_excluded(model, left, right, F(1), delta=delta):
                    certificates.append({**cell, 'route': 'direction_uniform'})
                    continue
                if use_tail and left >= 0:
                    proof = tail.certify_interval(model, left, right, delta=delta)
                    if proof['status'] == 'certified_excluded':
                        certificates.append({**cell, 'route': 'direction_tail', 'proof': proof})
                        continue
                midpoint = (left + right) / 2
                candidates = []
                if all(p.evaluate(poly, midpoint) > 0 for poly in size):
                    candidates.append(('size', None, size))
                for extra in range(1, n):
                    if extra == split:
                        continue
                    if extra not in shapes:
                        shapes[extra] = models.shape(split, extra, config)
                    if p.evaluate(shapes[extra], midpoint) > 0:
                        candidates.append(('shape', extra, (shapes[extra],)))
                if not candidates:
                    bounds = point_bounds(
                        model, midpoint, delta, use_tail=use_tail, use_spectral=use_spectral
                    )
                    if not bounds['excluded']:
                        return result(
                            'surviving_explanation',
                            {
                                'split': split,
                                'rho': str(midpoint),
                                'rho_float': float(midpoint),
                                'confidence_bounds': bounds,
                            },
                        )
                for route, extra, polynomials in candidates:
                    if all(p.positive_on(poly, left, right) for poly in polynomials):
                        certificates.append({**cell, 'route': route, 'extra': extra})
                        break
                else:
                    # Keep the inherited inexpensive routes first. Only interior
                    # cells can use the full spectral; endpoints use old routes.
                    if use_spectral and -1 < left <= right < 1:
                        proof = spectral.certify_interval(model, left, right, delta=delta)
                        if proof['status'] == 'certified_excluded':
                            certificates.append(
                                {**cell, 'route': 'direction_spectral', 'proof': proof}
                            )
                            continue
                    if depth >= max_depth:
                        return result(
                            'unresolved',
                            {
                                'split': split,
                                'left': str(left),
                                'right': str(right),
                                'reason': 'max_depth',
                            },
                        )
                    pending.extend(((midpoint, right, depth + 1), (left, midpoint, depth + 1)))
    except ar1._DegenerateResidual:
        return result('insufficient_variation')
    return result('certified_alert')
