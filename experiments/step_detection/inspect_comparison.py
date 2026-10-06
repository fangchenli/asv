"""Explain selected small histories using exact fits and production scoring."""

import argparse
import json

from . import compare_versions as c
from . import harness as h
from .exact_reference import fit_independent


def score_candidate(module, y, w, fit):
    """Use the revision's actual scoring closure on one externally supplied fit."""
    original_fit, original_search = module.solve_potts_approx, module.golden_search
    scores = []

    def fixed_fit(*args, **kwargs):
        return fit['right'], fit['values'], fit['costs']

    def search(f, a, b, **kwargs):
        if kwargs.get('ftol') != 0:
            return original_search(f, a, b, **kwargs)
        scores.append(f(0.0))
        return 0.0

    module.solve_potts_approx, module.golden_search = fixed_fit, search
    try:
        module.solve_potts_autogamma(y, w)
    finally:
        module.solve_potts_approx, module.golden_search = original_fit, original_search
    assert len(scores) == 1
    return scores[0]


def inspect(run, case_ids, backend='native'):
    cases = {case['id']: case for case in json.loads((run / 'inputs.json').read_text())}
    rows = []
    for case_id in case_ids:
        case = cases[case_id]
        row = {'case_id': case_id, 'backend': backend, 'versions': {}}
        for version in ('before', 'after'):
            module = c.load_source(run / f'{version}.py', backend)
            y, w, indices = h.prepare(module, case)
            diagnostic = c.diagnose(module, case, backend)
            frontier = fit_independent(y, w, noise_floor=1, beta=0, backend=backend)['frontier']
            scored = [
                dict(fit, production_score=score_candidate(module, y, w, fit)) for fit in frontier
            ]
            best = min(scored, key=lambda fit: fit['production_score'])
            diagnostic['best_frontier_candidate'] = best
            diagnostic['selected_minus_best_frontier_score'] = (
                diagnostic['selected']['score'] - best['production_score']
            )
            diagnostic['frontier_candidate_boundary_metrics'] = h.match_boundaries(
                h.observed_truth(case['true_boundaries'], indices),
                [indices[r] for r in best['right'][:-1]],
            )
            diagnostic['frontier_candidates'] = scored
            row['versions'][version] = diagnostic
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=c.Path)
    parser.add_argument(
        '--cases',
        nargs='+',
        default=[
            'weak-n100-seed0',
            'weak-n100-seed3',
            'weak-n100-seed5',
            'weak-n100-seed7',
            'weak-n100-seed9',
            'dip-n100-seed0',
            'outlier-n100-seed0',
        ],
    )
    parser.add_argument('--backend', choices=['native', 'python'], default='native')
    args = parser.parse_args()
    c.dump(
        args.run / f'inspection-{args.backend}.json', inspect(args.run, args.cases, args.backend)
    )


if __name__ == '__main__':
    main()
