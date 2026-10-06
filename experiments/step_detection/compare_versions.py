"""Paired, uninstrumented comparison of saved detector revisions.

Run as ``python -m experiments.step_detection.compare_versions --output ...``.
Every timing sample uses a fresh subprocess; diagnostics run separately.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

from . import harness as h
from .exact_reference import fit_independent

BEFORE = 'd526c2e'
AFTER = '123405a'
SNAPSHOT = Path(__file__).with_name('data') / 'numpy_snapshot.json'


def load_source(path, backend):
    spec = importlib.util.spec_from_file_location('_comparison_detector', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module._rangemedian = h.load_detector(backend)._rangemedian
    return module


def evaluate(module, case):
    """Time the public detector and alert pipeline, including input preparation."""
    started = time.perf_counter()
    steps = module.detect_steps(case['values'], case['weights'])
    revisions = case['revisions']
    graph_steps = [(revisions[l], revisions[r - 1] + 1, v, m, e) for l, r, v, m, e in steps]
    latest, best, alerts = module.detect_regressions(graph_steps, threshold=0.05)
    seconds = time.perf_counter() - started
    rss = h.peak_rss_bytes()
    indices = [
        i
        for i, (y, w) in enumerate(zip(case['values'], case['weights']))
        if y is not None and math.isfinite(y) and (w is None or w > 0)
    ]
    boundaries = [s[0] for s in steps[1:]]
    positions = {rev: i for i, rev in enumerate(revisions)}
    alert_positions = [positions[a[1]] for a in (alerts or [])]
    return {
        'seconds': seconds,
        'peak_rss_bytes': rss,
        'steps': steps,
        'boundaries': boundaries,
        'alerts': alerts or [],
        'alert_positions': alert_positions,
        'latest': latest,
        'best': best,
        'boundary_metrics': h.match_boundaries(
            h.observed_truth(case['true_boundaries'], indices), boundaries
        ),
        'alert_metrics': h.match_boundaries(
            h.observed_truth(case['true_alerts'], indices), alert_positions
        ),
    }


def diagnose(module, case, backend):
    """Record actual outer-search candidates without changing its decisions."""
    y, w, _ = h.prepare(module, case)
    original_search = module.golden_search
    original_fit = module.solve_potts_approx
    candidates = []
    pending = {}

    def fit(*args, **kwargs):
        result = original_fit(*args, **kwargs)
        pending.update(gamma=kwargs['gamma'], right=result[0], values=result[1], costs=result[2])
        return result

    def search(f, a, b, **kwargs):
        # Only the outer gamma search sets ftol=0. Preserve old inner rho search.
        if kwargs.get('ftol') != 0:
            return original_search(f, a, b, **kwargs)

        def record(x):
            score = f(x)
            candidates.append(dict(pending, score=score))
            return score

        return original_search(record, a, b, **kwargs)

    module.golden_search, module.solve_potts_approx = search, fit
    try:
        result = module.solve_potts_autogamma(y, w)
    finally:
        module.golden_search, module.solve_potts_approx = original_search, original_fit
    selected = min(candidates, key=lambda c: c['score'])
    assert selected['right'] == result[0] and selected['gamma'] == result[3]
    diagnostic = {'selected': selected, 'candidates': candidates}
    if len(y) <= 200:
        exact = module.solve_potts(y, w, gamma=result[3])
        actual_objective = h.check_solution(y, w, result[3], result[:3])
        exact_objective = h.check_solution(y, w, result[3], exact)
        # Only D_K is used here. Floor/beta select no model for this diagnostic.
        frontier = fit_independent(y, w, noise_floor=1, beta=0, backend=backend)['frontier']
        same_count = next(f for f in frontier if f['segments'] == len(result[0]))
        diagnostic.update(
            fixed_gamma_gap=actual_objective - exact_objective,
            fixed_gamma_exact_right=exact[0],
            same_count_error_gap=math.fsum(result[2]) - same_count['error_sum'],
            same_count_exact_right=same_count['right'],
        )
    return diagnostic


def make_cases(seeds, include_real=True, include_scaling=True):
    cases = []
    for scenario in h.SCENARIOS:
        for seed in seeds:
            case = h.make_case(scenario, 100, seed)
            case['corpus'] = 'baseline' if seed < 4 else 'new_seeds'
            cases.append(case)
    if include_real:
        cases.extend(json.loads(SNAPSHOT.read_text())['cases'])
    if include_scaling:
        for size in (1000, 10000):
            for scenario in ('flat', 'step'):
                case = h.make_case(scenario, size, 0)
                case['corpus'] = 'scaling'
                cases.append(case)
    return cases


def summarize(records):
    groups = defaultdict(list)
    for row in records:
        groups[row['case_id'], row['backend'], row['version']].append(row)
    reduced = {}
    for key, rows in groups.items():
        # Detect any nondeterminism before aggregating accuracy or timing.
        for row in rows[1:]:
            for field in ('steps', 'alerts', 'boundary_metrics', 'alert_metrics'):
                if row[field] != rows[0][field]:
                    raise ValueError(f'Nondeterministic {field}: {key}')
        reduced[key] = dict(
            rows[0],
            seconds=statistics.median(r['seconds'] for r in rows),
            peak_rss_bytes=statistics.median(r['peak_rss_bytes'] for r in rows),
        )
    pairs = []
    for (case_id, backend, version), before in reduced.items():
        if version != 'before':
            continue
        after = reduced[case_id, backend, 'after']
        pairs.append(
            {
                'case_id': case_id,
                'backend': backend,
                'corpus': before['corpus'],
                'before': before,
                'after': after,
                'runtime_ratio': after['seconds'] / before['seconds'],
                'rss_delta_bytes': after['peak_rss_bytes'] - before['peak_rss_bytes'],
                'fit_changed': after['steps'] != before['steps'],
                'alert_locations_changed': after['alert_positions'] != before['alert_positions'],
                'alert_records_changed': after['alerts'] != before['alerts'],
            }
        )
    summaries = []
    grouped_pairs = defaultdict(list)
    for pair in pairs:
        grouped_pairs[pair['corpus'], pair['backend']].append(pair)
    for (corpus, backend), rows in sorted(grouped_pairs.items()):
        summary = {
            'corpus': corpus,
            'backend': backend,
            'histories': len(rows),
            'changed_fits': sum(r['fit_changed'] for r in rows),
            'changed_alert_locations': sum(r['alert_locations_changed'] for r in rows),
            'changed_alert_records': sum(r['alert_records_changed'] for r in rows),
            'median_paired_runtime_ratio': statistics.median(r['runtime_ratio'] for r in rows),
            'median_rss_delta_bytes': statistics.median(r['rss_delta_bytes'] for r in rows),
        }
        for version in ('before', 'after'):
            totals = {'median_seconds': statistics.median(r[version]['seconds'] for r in rows)}
            for field in ('boundary_metrics', 'alert_metrics'):
                metrics = [r[version][field] for r in rows if r[version][field] is not None]
                totals[field] = {
                    key: sum(m[key] for m in metrics)
                    for key in ('matched', 'missed', 'false', 'location_error_sum')
                }
                totals[field]['scorable_histories'] = len(metrics)
                totals[field]['histories_with_false'] = sum(m['false'] > 0 for m in metrics)
            summary[version] = totals
        summaries.append(summary)
    return {'groups': summaries, 'pairs': pairs}


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--before', default=BEFORE)
    parser.add_argument('--after', default=AFTER)
    parser.add_argument('--seeds', nargs='+', type=int, default=list(range(12)))
    parser.add_argument(
        '--backends', nargs='+', choices=('native', 'python'), default=['native', 'python']
    )
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--no-real', action='store_true')
    parser.add_argument('--no-scaling', action='store_true')
    args = parser.parse_args()
    if args.worker:
        request = json.load(sys.stdin)
        module = load_source(request['source'], request['backend'])
        result = evaluate(module, request['case'])
        print(json.dumps(result, allow_nan=False))
        return
    if args.output is None or args.repeats < 1:
        parser.error('--output and positive --repeats are required')
    args.output.mkdir(parents=True, exist_ok=False)
    sources, revisions = {}, {}
    for version, revision in [('before', args.before), ('after', args.after)]:
        revisions[version] = subprocess.check_output(
            ['git', 'rev-parse', revision], cwd=h.ROOT, text=True
        ).strip()
        source = subprocess.check_output(
            ['git', 'show', f'{revisions[version]}:asv/step_detect.py'], cwd=h.ROOT
        )
        sources[version] = args.output.resolve() / f'{version}.py'
        sources[version].write_bytes(source)
    cases = make_cases(args.seeds, not args.no_real, not args.no_scaling)
    dump(args.output / 'inputs.json', cases)
    native = h.load_detector('native')._rangemedian if 'native' in args.backends else None
    manifest = {
        'revisions': revisions,
        'python': sys.version,
        'platform': platform.platform(),
        'command': sys.argv,
        'repeats': args.repeats,
        'backends': args.backends,
        'source_sha256': {
            v: hashlib.sha256(p.read_bytes()).hexdigest() for v, p in sources.items()
        },
        'native_sha256': hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest()
        if native
        else None,
        'inputs_sha256': hashlib.sha256((args.output / 'inputs.json').read_bytes()).hexdigest(),
        'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'timing': 'fresh process per sample; public detector plus alert mapping; imports excluded',
        'order': 'alternate version order by case index plus repeat index',
        'threshold': 0.05,
        'boundary_tolerance': 2,
    }
    dump(args.output / 'manifest.json', manifest)
    records = []
    with (args.output / 'records.jsonl').open('w') as stream:
        for i, case in enumerate(cases):
            for backend in args.backends:
                for repeat in range(args.repeats):
                    versions = (
                        ['before', 'after'] if (i + repeat) % 2 == 0 else ['after', 'before']
                    )
                    for version in versions:
                        payload = {
                            'source': str(sources[version]),
                            'backend': backend,
                            'case': case,
                        }
                        worker = subprocess.run(
                            [
                                sys.executable,
                                '-m',
                                'experiments.step_detection.compare_versions',
                                '--worker',
                            ],
                            input=json.dumps(payload),
                            text=True,
                            capture_output=True,
                            cwd=h.ROOT,
                            timeout=120,
                            check=True,
                        )
                        row = dict(
                            json.loads(worker.stdout),
                            case_id=case['id'],
                            corpus=case['corpus'],
                            backend=backend,
                            version=version,
                            repeat=repeat,
                        )
                        records.append(row)
                        stream.write(json.dumps(row, allow_nan=False) + '\n')
                        stream.flush()
            print(f'{i + 1}/{len(cases)} {case["id"]}', flush=True)
    summary = summarize(records)
    dump(args.output / 'summary.json', summary)
    # Diagnostics are separate from every time/RSS measurement.
    indexed = {c['id']: c for c in cases}
    diagnostics = []
    for pair in summary['pairs']:
        if not pair['fit_changed']:
            continue
        case = indexed[pair['case_id']]
        row = {'case_id': case['id'], 'backend': pair['backend']}
        for version in sources:
            module = load_source(sources[version], pair['backend'])
            row[version] = diagnose(module, case, pair['backend'])
        diagnostics.append(row)
    dump(args.output / 'diagnostics.json', diagnostics)


if __name__ == '__main__':
    main()
