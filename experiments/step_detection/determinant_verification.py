"""Verify the new route and reuse the frozen verifier for inherited routes."""

import copy
from collections import Counter
from fractions import Fraction as F

from . import directional_determinant as determinant
from . import directional_diagnosis as inherited
from . import residual_direction as direction


def verify(values, result):
    """Reconstruct certificates from independently centered residual states.

    The interval arithmetic kernel is shared with the detector. Independent
    exact dense determinant checks live in test_directional_determinant.py.
    """
    confidence = result['confidence']
    assert confidence['rule'] == 'directional_full_determinant'
    assert isinstance(confidence['determinant_enabled'], bool)
    assert isinstance(confidence['tail_enabled'], bool)
    assert confidence['determinant_bits'] == determinant.BITS
    assert result['has_alert'] == (result['status'] == 'certified_alert')
    delta = F(confidence['delta'])
    states, routes = {}, Counter()
    coverage = {split: [] for split in range(1, len(values))}

    def state(split):
        if split not in states:
            states[split] = direction.state(values, split)
        return states[split]

    for cell in result['certificate']:
        split, left, right = cell['split'], F(cell['left']), F(cell['right'])
        assert split in coverage and -1 <= left < right <= 1
        assert str(split) in confidence['by_split']
        coverage[split].append((left, right))
        routes[cell['route']] += 1
        if cell['route'] == 'direction_determinant':
            assert confidence['determinant_enabled'] and -1 < left < right < 1
            assert cell['extra'] is None
            proof = determinant.certify_interval(state(split), left, right, delta=delta)
            assert cell['proof'] == proof and proof['status'] == 'certified_excluded'
            assert F(proof['determinant_enclosure'][0]) > 1 / delta**2
    for intervals in coverage.values():
        intervals.sort()
        assert all(a[1] <= b[0] for a, b in zip(intervals[:-1], intervals[1:], strict=True))
        if result['has_alert']:
            assert intervals and intervals[0][0] == -1 and intervals[-1][1] == 1
            assert all(a[1] == b[0] for a, b in zip(intervals[:-1], intervals[1:], strict=True))

    legacy = copy.deepcopy(result)
    legacy['certificate'] = [
        cell for cell in legacy['certificate'] if cell['route'] != 'direction_determinant'
    ]
    if result['status'] == 'surviving_explanation':
        witness = result['witness']
        split, rho = witness['split'], F(witness['rho'])
        assert split in coverage and -1 < rho < 1
        proof = (
            determinant.certify_interval(state(split), rho, rho, delta=delta)
            if confidence['determinant_enabled']
            else None
        )
        assert witness['confidence_bounds']['determinant'] == proof
        assert proof is None or F(proof['p_upper_squared']) >= delta**2
        del legacy['witness']['confidence_bounds']['determinant']
    # The inherited verifier checks its retained intervals as partial coverage.
    # Full coverage, including the new route, was checked above.
    if result['has_alert']:
        legacy['status'], legacy['has_alert'] = 'unresolved', False
    inherited.verify(values, legacy)
    return {
        'verified': True,
        'complete_coverage': result['has_alert'],
        'certificate_routes': dict(routes),
    }
