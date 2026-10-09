"""Inspect saved unresolved correlation intervals and their midpoint bounds."""

import argparse
import gzip
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

from . import directional_sturm as sturm
from . import directional_tail as tail
from . import reporting_ar1 as ar1
from . import reporting_sturm as reporting

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(100_000)

HERE = Path(__file__).parent
INPUTS = HERE / 'data/sturm_fresh_v2_inputs.jsonl.gz'
RECORDS = HERE / 'data/sturm_fresh_v2_records.jsonl.gz'
DELTA = F(1, 100)


def unresolved_cases():
    records = {
        row['case']['id']: row
        for row in (json.loads(line) for line in gzip.open(RECORDS, 'rt'))
        if row['methods']['rank_group']['status'] == 'unresolved'
    }
    inputs = {row['id']: row for row in (json.loads(line) for line in gzip.open(INPUTS, 'rt'))}
    if records.keys() - inputs.keys():
        raise ValueError('An unresolved record has no matching saved input')
    return [(inputs[key], records[key]['methods']['rank_group']) for key in sorted(records)]


def diagnose_case(case, original):
    stop = original['witness']
    split, left, right = int(stop['split']), F(stop['left']), F(stop['right'])
    midpoint = (left + right) / 2
    model = reporting.confidence_state(ar1.Models(case['values']), split)
    interval = sturm.certify_interval(model, left, right)
    point = sturm.certify_interval(model, midpoint, midpoint)
    halves = [
        sturm.certify_interval(model, lower, upper)
        for lower, upper in ((left, midpoint), (midpoint, right))
    ]
    numeric = tail.approximate_tail(case['values'], split, float(midpoint))

    def proof_summary(proof):
        return {
            'status': proof['status'],
            'p_upper': proof['p_upper'],
            'selected_tilt': proof['selected_tilt'],
            'selected_count': proof['selected_count'],
            'attempts': [
                {
                    'tilt': attempt['tilt'],
                    'p_upper': attempt['p_upper'],
                    'selected_count': attempt['selected_count'],
                    'reason': attempt.get('reason'),
                    'base_reason': attempt.get('base', {}).get('reason'),
                }
                for attempt in proof['attempts']
            ],
        }

    return {
        'id': case['id'],
        'condition': case['condition'],
        'n': case['n'],
        'change': case['change'],
        'noise': case['noise'],
        'location': case['location'],
        'original_stop': {
            'reason': stop['reason'],
            'cells': original['cells_visited'],
            'split': split,
            'left': str(left),
            'right': str(right),
            'width': str(right - left),
        },
        'interval_certificate': proof_summary(interval),
        'midpoint_certificate': {'rho': str(midpoint), **proof_summary(point)},
        'one_level_subdivision': [proof_summary(proof) for proof in halves],
        'numerical_tail_at_midpoint': numeric,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    rows = [diagnose_case(case, original) for case, original in unresolved_cases()]
    result = {
        'purpose': 'Post-hoc interval diagnosis only; no decisions or settings were changed',
        'confidence_cutoff': str(DELTA),
        'input_archive_sha256': hashlib.sha256(INPUTS.read_bytes()).hexdigest(),
        'record_archive_sha256': hashlib.sha256(RECORDS.read_bytes()).hexdigest(),
        'case_count': len(rows),
        'rows': rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'cases': len(rows), 'output': str(args.output)}, indent=2))


if __name__ == '__main__':
    main()
