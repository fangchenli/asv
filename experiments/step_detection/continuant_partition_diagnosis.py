"""Evaluate one saved unresolved interval with a continuant certificate."""

import argparse
import gzip
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

from . import determinant_continuant as continuant
from . import directional_determinant as determinant
from . import directional_sturm as sturm
from . import reporting_ar1 as ar1
from . import reporting_sturm as reporting

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(100_000)

HERE = Path(__file__).parent
INPUTS = HERE / 'data/sturm_fresh_v2_inputs.jsonl.gz'
RECORDS = HERE / 'data/sturm_fresh_v2_records.jsonl.gz'
TILT = F(1, 8)


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


def diagnose(case, original, divisions):
    if divisions < 1 or divisions & (divisions - 1):
        raise ValueError('Require a positive power-of-two partition count')
    witness = original['witness']
    left, right = F(witness['left']), F(witness['right'])
    split = int(witness['split'])
    model = reporting.confidence_state(ar1.Models(case['values']), split)

    parent = determinant.certify_interval(model, left, right, tilt=TILT)
    radius = F(parent['radius_upper'])
    matrix_a = 1 - 2 * TILT
    matrix_b = 2 * TILT / radius
    center, half_width = (left + right) / 2, (right - left) / 2
    numerator, denominator = continuant.determinant_polynomial(
        model['n'], split, matrix_a, matrix_b, center
    )
    numerator_bernstein = continuant._bernstein_coefficients(numerator, -half_width, half_width)
    denominator_bernstein = continuant._bernstein_coefficients(
        denominator, -half_width, half_width
    )
    pieces = [(numerator_bernstein, denominator_bernstein)]
    while len(pieces) < divisions:
        pieces = [
            (num_part, den_part)
            for num, den in pieces
            for num_part, den_part in zip(
                continuant._split_bernstein(num), continuant._split_bernstein(den)
            )
        ]

    rows = []
    for index, (numerator_piece, denominator_piece) in enumerate(pieces):
        cell_left = left + (right - left) * F(index, divisions)
        cell_right = left + (right - left) * F(index + 1, divisions)
        determinant_lower = continuant._bernstein_ratio_lower(numerator_piece, denominator_piece)
        if determinant_lower is None:
            raise ArithmeticError(
                f'Denominator Bernstein coefficients are inconclusive in cell {index}'
            )

        def certifier(
            _model,
            lower,
            upper,
            *,
            tilt,
            delta,
            bits,
            _left=cell_left,
            _right=cell_right,
            _determinant_lower=determinant_lower,
        ):
            if (F(lower), F(upper), F(tilt)) != (_left, _right, TILT):
                raise ValueError('Unexpected interval or tilt in partition certificate')
            return {
                'left': str(_left),
                'right': str(_right),
                'tilt': str(TILT),
                'delta': str(delta),
                'bits': bits,
                'radius_upper': str(radius),
                'determinant_enclosure': [str(_determinant_lower), str(_determinant_lower)],
                'status': 'not_certified',
            }

        proof = sturm.at_tilt(
            model,
            cell_left,
            cell_right,
            TILT,
            determinant_certifier=certifier,
            rank_groups_on_intervals=True,
        )
        rows.append(
            {
                'index': index,
                'left': str(cell_left),
                'right': str(cell_right),
                'determinant_lower': str(determinant_lower),
                'density_correction': proof['density_correction'],
                'p_upper': proof['p_upper'],
                'status': proof['status'],
                'selected_count': proof['selected_count'],
            }
        )

    return {
        'purpose': 'Post-hoc proof diagnostic; original frozen outcomes are unchanged',
        'case_id': case['id'],
        'split': split,
        'tilt': str(TILT),
        'divisions': divisions,
        'shared_radius_upper': str(radius),
        'input_archive_sha256': hashlib.sha256(INPUTS.read_bytes()).hexdigest(),
        'record_archive_sha256': hashlib.sha256(RECORDS.read_bytes()).hexdigest(),
        'all_certified': all(row['status'] == 'certified_excluded' for row in rows),
        'max_p_upper': str(max(F(row['p_upper']) for row in rows)),
        'rows': rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-id', required=True)
    parser.add_argument('--divisions', type=int, default=16)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    case, original = saved_case(args.case_id)
    result = diagnose(case, original, args.divisions)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(
        json.dumps(
            {
                'case_id': result['case_id'],
                'divisions': result['divisions'],
                'all_certified': result['all_certified'],
                'max_p_upper': result['max_p_upper'],
                'output': str(args.output),
            },
            indent=2,
        )
    )


if __name__ == '__main__':
    main()
