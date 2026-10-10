"""Check radius consistency, coverage, and conservative partition reporting."""

import hashlib
import json
from fractions import Fraction as F

import numpy as np
import pytest

from experiments.step_detection import interval_polynomial_diagnosis as diagnosis
from experiments.step_detection import residual_direction as direction


@pytest.mark.parametrize('radius_scope', ['parent', 'cell'])
@pytest.mark.parametrize('inconclusive', [False, True])
def test_partition_uses_consistent_radius_and_reports_worst_cell(
    monkeypatch, radius_scope, inconclusive
):
    values = [1, 4, 2, 5, 3, 8]
    split, left, middle, right = 2, F(3, 4), F(13, 16), F(7, 8)
    case = {'id': 'small-partition-test', 'values': values}
    record = {'witness': {'split': split, 'left': str(left), 'right': str(right)}}
    seen = {}

    def polynomial_bound(n, candidate_split, lo, hi, a, b, *, bits):
        assert (n, candidate_split, bits) == (len(values), split, 128)
        tilt = (1 - a) / 2
        seen[lo, hi, tilt] = 2 * tilt / b
        return None if inconclusive and lo == middle else F(4)

    def tail_bound(model, lo, hi, tilt, *, determinant_certifier, rank_groups_on_intervals):
        assert rank_groups_on_intervals
        proof = determinant_certifier(model, lo, hi, tilt=tilt, delta=F(1, 100), bits=192)
        radius = F(proof['radius_upper'])
        assert radius == seen[lo, hi, tilt]
        # An independent dense covariance calculation checks radius containment.
        basis = direction.residual_basis(len(values), split)
        residual = basis.T @ np.array(values)
        residual /= np.linalg.norm(residual)
        distance = np.abs(np.arange(len(values))[:, None] - np.arange(len(values))[None, :])
        for rho in (lo, (lo + hi) / 2, hi):
            covariance = basis.T @ (float(rho) ** distance) @ basis
            exact_radius = 1 / (residual @ np.linalg.solve(covariance, residual))
            assert exact_radius <= float(radius)
        # Tilt choices cross between cells. Picking the global best is unsafe.
        if lo == left:
            p_upper = F(1, 200) if tilt == diagnosis.spectral.TILTS[0] else F(1, 50)
        else:
            p_upper = F(1, 40) if tilt == diagnosis.spectral.TILTS[0] else F(1, 80)
        return {
            'density_correction': '1',
            'selected_count': None,
            'p_upper': str(p_upper),
            'p_upper_squared': str(p_upper**2),
            'status': 'certified_excluded' if p_upper < F(1, 100) else 'not_certified',
        }

    monkeypatch.setattr(diagnosis.polynomial, 'determinant_interval', polynomial_bound)
    monkeypatch.setattr(diagnosis.sturm, 'at_tilt', tail_bound)
    result = diagnosis.diagnose(case, record, bits=128, divisions=2, radius_scope=radius_scope)
    expected = F(1) if inconclusive else F(1, 80)
    assert F(result['p_upper']) == F(result['max_cell_p_upper']) == expected
    assert F(result['p_upper_squared']) == expected**2
    assert result['worst_cell_index'] == 1
    assert result['status'] == 'not_certified'
    cells = result['cells']
    assert [(F(c['left']), F(c['right'])) for c in cells] == [(left, middle), (middle, right)]
    assert cells[0]['best']['tilt'] == str(diagnosis.spectral.TILTS[0])
    if not inconclusive:
        assert cells[1]['best']['tilt'] == str(diagnosis.spectral.TILTS[1])
    parent_radius = F(result['parent_radius_upper'])
    for cell in cells:
        cell_radius = F(cell['radius_upper'])
        assert cell_radius <= parent_radius
        if radius_scope == 'parent':
            assert cell_radius == parent_radius
    if radius_scope == 'cell':
        assert F(cells[1]['radius_upper']) < parent_radius


@pytest.mark.parametrize('divisions', [0, -1, 3, True, 1.5])
def test_invalid_partition_rejected(divisions):
    with pytest.raises(ValueError, match='positive power of two'):
        diagnosis.diagnose({}, {}, bits=128, divisions=divisions)


def test_invalid_radius_scope_rejected():
    with pytest.raises(ValueError, match='radius_scope'):
        diagnosis.diagnose({}, {}, bits=128, radius_scope='invalid')


@pytest.mark.parametrize(
    'location,divisions,expected_certified',
    [('early', 2, False), ('early', 8, True), ('middle', 8, False), ('recent', 8, True)],
)
def test_saved_local_radius_partition_has_complete_consistent_coverage(
    location, divisions, expected_certified
):
    path = diagnosis.HERE / (
        f'data/sturm_fresh_v2_interval_poly_n100_{location}_local_d{divisions}.json'
    )
    proof = json.loads(path.read_text())
    assert proof['schema_version'] == 2
    assert proof['radius_scope'] == 'cell'
    assert (
        proof['input_archive_sha256'] == hashlib.sha256(diagnosis.INPUTS.read_bytes()).hexdigest()
    )
    assert (
        proof['record_archive_sha256']
        == hashlib.sha256(diagnosis.RECORDS.read_bytes()).hexdigest()
    )
    case, record = diagnosis.saved_case(proof['case_id'])
    assert proof['split'] == record['witness']['split']
    assert proof['left'] == record['witness']['left']
    assert proof['right'] == record['witness']['right']
    model = diagnosis.reporting.confidence_state(
        diagnosis.ar1.Models(case['values']), proof['split']
    )
    left, right = F(proof['left']), F(proof['right'])
    assert len(proof['cells']) == divisions
    for index, cell in enumerate(proof['cells']):
        lo, hi = F(cell['left']), F(cell['right'])
        assert lo == left + (right - left) * F(index, divisions)
        assert hi == left + (right - left) * F(index + 1, divisions)
        radius = diagnosis.determinant.certify_interval(model, lo, hi)['radius_upper']
        assert cell['radius_upper'] == radius
        assert F(radius) <= F(proof['parent_radius_upper'])
        assert [F(a['tilt']) for a in cell['attempts']] == list(diagnosis.spectral.TILTS)
        for attempt in cell['attempts']:
            squared = min(
                F(1),
                1 / (F(attempt['determinant_lower']) * F(attempt['density_correction']) ** 2),
            )
            assert squared == F(attempt['p_upper_squared'])
            upper = F(attempt['p_upper'])
            assert (upper - F(1, 2**48)) ** 2 <= squared <= upper**2
        assert cell['best'] == min(cell['attempts'], key=lambda a: F(a['p_upper_squared']))
    worst = max(proof['cells'], key=lambda c: F(c['best']['p_upper_squared']))
    assert proof['worst_cell_index'] == worst['index']
    assert proof['p_upper_squared'] == worst['best']['p_upper_squared']
    assert proof['p_upper'] == proof['max_cell_p_upper'] == worst['best']['p_upper']
    certified = all(F(c['best']['p_upper_squared']) < F(1, 100) ** 2 for c in proof['cells'])
    assert (proof['status'] == 'certified_excluded') == certified == expected_certified
