"""Certify a saved unresolved interval with rounded polynomial coefficients."""

import argparse
import gzip
import hashlib
import json
import sys
from fractions import Fraction as F
from functools import partial
from pathlib import Path

from . import determinant_interval_polynomial as polynomial
from . import directional_determinant as determinant
from . import directional_normalized as normalized
from . import directional_spectral as spectral
from . import directional_sturm as sturm
from . import reporting_ar1 as ar1
from . import reporting_sturm as reporting

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(100_000)

HERE = Path(__file__).parent
INPUTS = HERE / 'data/sturm_fresh_v2_inputs.jsonl.gz'
RECORDS = HERE / 'data/sturm_fresh_v2_records.jsonl.gz'


def _determinant_certifier(
    _model,
    lower,
    upper,
    *,
    tilt,
    delta,
    bits,
    expected_left,
    expected_right,
    expected_tilt,
    radius,
    determinant_lower,
    stationary_normalized=False,
):
    if (F(lower), F(upper), F(tilt)) != (expected_left, expected_right, expected_tilt):
        raise ValueError('Unexpected interval or tilt')
    return {
        'left': str(expected_left),
        'right': str(expected_right),
        'tilt': str(expected_tilt),
        'delta': str(delta),
        'bits': bits,
        'radius_upper': str(radius),
        'stationary_normalized': stationary_normalized,
        'determinant_enclosure': [str(determinant_lower), str(determinant_lower)],
        'status': 'rounded_polynomial_lower_bound',
    }


def saved_case(case_id):
    records = {
        row['case']['id']: row for row in (json.loads(line) for line in gzip.open(RECORDS, 'rt'))
    }
    inputs = {row['id']: row for row in (json.loads(line) for line in gzip.open(INPUTS, 'rt'))}
    if case_id not in records or case_id not in inputs:
        raise ValueError(f'No matching saved case: {case_id}')
    if records[case_id]['methods']['rank_group']['status'] != 'unresolved':
        raise ValueError(f'Saved candidate case is not unresolved: {case_id}')
    return inputs[case_id], records[case_id]['methods']['rank_group']


def _radius_upper(model, left, right, *, stationary_normalized=False):
    """Bound the radius, optionally after cancelling its stationary factor."""
    if not stationary_normalized:
        return F(
            determinant.certify_interval(model, left, right, tilt=spectral.TILTS[0])[
                'radius_upper'
            ]
        )

    return normalized.radius_upper(model, left, right)


def diagnose(case, record, *, bits, divisions=1, radius_scope='cell', stationary_normalized=False):
    if (
        not isinstance(divisions, int)
        or isinstance(divisions, bool)
        or divisions < 1
        or divisions & (divisions - 1)
    ):
        raise ValueError('divisions must be a positive power of two')
    if radius_scope not in ('parent', 'cell'):
        raise ValueError('radius_scope must be parent or cell')
    witness = record['witness']
    left, right = F(witness['left']), F(witness['right'])
    split = int(witness['split'])
    model = reporting.confidence_state(ar1.Models(case['values']), split)
    parent_radius = _radius_upper(model, left, right, stationary_normalized=stationary_normalized)
    cells = []
    width = right - left
    for index in range(divisions):
        cell_left = left + width * index / divisions
        cell_right = left + width * (index + 1) / divisions
        # A parent radius is valid but can stay loose even as cells shrink.
        # Use the same radius in the determinant and density correction.
        radius = parent_radius
        if radius_scope == 'cell' and divisions > 1:
            radius = _radius_upper(
                model, cell_left, cell_right, stationary_normalized=stationary_normalized
            )
        cell_attempts = []
        for tilt in spectral.TILTS:
            matrix_a, matrix_b = 1 - 2 * tilt, 2 * tilt / radius
            determinant_lower = polynomial.determinant_interval(
                model['n'],
                split,
                cell_left,
                cell_right,
                matrix_a,
                matrix_b,
                bits=bits,
                stationary_normalized=stationary_normalized,
            )
            if determinant_lower is None:
                attempt = {
                    'tilt': str(tilt),
                    'determinant_lower': None,
                    'status': 'determinant_inconclusive',
                    'p_upper_squared': '1',
                    'p_upper': '1',
                }
                cell_attempts.append(attempt)
                continue

            attempt_tilt = tilt

            determinant_certifier = partial(
                _determinant_certifier,
                expected_left=cell_left,
                expected_right=cell_right,
                expected_tilt=attempt_tilt,
                radius=radius,
                determinant_lower=determinant_lower,
                stationary_normalized=stationary_normalized,
            )
            proof = sturm.at_tilt(
                model,
                cell_left,
                cell_right,
                tilt,
                determinant_certifier=determinant_certifier,
                rank_groups_on_intervals=True,
            )
            attempt = {
                'tilt': str(tilt),
                'determinant_lower': str(determinant_lower),
                'density_correction': proof['density_correction'],
                'selected_count': proof['selected_count'],
                'p_upper_squared': proof['p_upper_squared'],
                'p_upper': proof['p_upper'],
                'status': proof['status'],
            }
            cell_attempts.append(attempt)
        best_cell = min(cell_attempts, key=lambda attempt: F(attempt['p_upper_squared']))
        cells.append(
            {
                'index': index,
                'left': str(cell_left),
                'right': str(cell_right),
                'radius_upper': str(radius),
                'best': best_cell,
                'attempts': cell_attempts,
            }
        )
    # Each cell may choose its best tilt, but all cells must be covered.
    worst_cell = max(cells, key=lambda cell: F(cell['best']['p_upper_squared']))
    best = worst_cell['best']
    all_certified = all(cell['best']['status'] == 'certified_excluded' for cell in cells)
    worst_p = max(F(cell['best']['p_upper']) for cell in cells)
    return {
        'schema_version': 3,
        'purpose': 'Post-hoc proof diagnostic; frozen outcomes are unchanged',
        'case_id': case['id'],
        'split': split,
        'left': str(left),
        'right': str(right),
        'tilt': best['tilt'],
        'bits': bits,
        'radius_scope': radius_scope,
        'stationary_normalized': stationary_normalized,
        'radius_bits': determinant.BITS,
        'parent_radius_upper': str(parent_radius),
        'radius_upper': worst_cell['radius_upper'],
        'worst_cell_index': worst_cell['index'],
        'determinant_lower': best['determinant_lower'],
        'density_correction': best.get('density_correction'),
        'selected_count': best.get('selected_count'),
        'p_upper_squared': best['p_upper_squared'],
        'p_upper': best['p_upper'],
        'max_cell_p_upper': str(worst_p),
        'status': 'certified_excluded' if all_certified else 'not_certified',
        'divisions': divisions,
        'cells': cells,
        'input_archive_sha256': hashlib.sha256(INPUTS.read_bytes()).hexdigest(),
        'record_archive_sha256': hashlib.sha256(RECORDS.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-id', required=True)
    parser.add_argument('--bits', type=int, default=128)
    parser.add_argument('--divisions', type=int, default=1)
    parser.add_argument('--radius-scope', choices=('parent', 'cell'), default='cell')
    parser.add_argument('--stationary-normalized', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    case, record = saved_case(args.case_id)
    result = diagnose(
        case,
        record,
        bits=args.bits,
        divisions=args.divisions,
        radius_scope=args.radius_scope,
        stationary_normalized=args.stationary_normalized,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(
        json.dumps(
            {
                'case_id': result['case_id'],
                'bits': result['bits'],
                'divisions': result['divisions'],
                'radius_scope': result['radius_scope'],
                'stationary_normalized': result['stationary_normalized'],
                'determinant_lower': result['determinant_lower'],
                'p_upper': result['p_upper'],
                'max_cell_p_upper': result['max_cell_p_upper'],
                'status': result['status'],
                'output': str(args.output),
            },
            indent=2,
        )
    )


if __name__ == '__main__':
    main()
