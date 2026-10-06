"""Reproducible step-detection experiments using ASV's production scoring code.

Run this file directly; the Python backend needs only the standard library.
Each measured method runs in a fresh subprocess to isolate caches and peak RSS.
"""

import argparse
import hashlib
import importlib.machinery
import importlib.util
import itertools
import json
import math
import platform
import random
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'asv' / 'step_detect.py'
SCENARIOS = (
    'flat',
    'step',
    'weak',
    'recent',
    'recovered',
    'partial',
    'dip',
    'multiple',
    'outlier',
    'correlated',
    'weighted',
    'missing',
    'sparse',
    'quantized',
    'drift',
)
METHODS = ('exact', 'approximate', 'current', 'grid', 'hybrid')


def load_detector(backend='python'):
    """Load an isolated module so instrumentation never changes imported ASV."""
    spec = importlib.util.spec_from_file_location('_step_detection_experiment', SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if backend == 'native':
        candidates = [
            SOURCE.with_name('_rangemedian' + suffix)
            for suffix in importlib.machinery.EXTENSION_SUFFIXES
        ]
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            raise RuntimeError('Build the local extension first: pip install -e .')
        native_spec = importlib.util.spec_from_file_location('_rangemedian', path)
        native = importlib.util.module_from_spec(native_spec)
        native_spec.loader.exec_module(native)
        module._rangemedian = native
    elif backend != 'python':
        raise ValueError(f'Unknown backend: {backend}')
    return module


def make_case(scenario, size, seed):
    """Generate data and explicit truth; retain arrays rather than only seeds."""
    if scenario not in SCENARIOS or size < 12:
        raise ValueError('Use a known scenario and size >= 12')
    rng = random.Random(seed)
    levels = [10.0] * size
    changes = []
    alerts = []
    transitions = []
    mid = size // 2
    if scenario in ('step', 'weighted', 'missing', 'sparse', 'quantized'):
        transitions = [(mid, 12)]
        alerts = [mid]
    elif scenario == 'weak':
        transitions = [(mid, 10.2)]
    elif scenario == 'recent':
        transitions = [(size - 4, 12)]
        alerts = [size - 4]
    elif scenario == 'recovered':
        transitions = [(size // 3, 12), (2 * size // 3, 10)]
    elif scenario == 'partial':
        transitions = [(size // 3, 14), (2 * size // 3, 12)]
        alerts = [size // 3]
    elif scenario == 'dip':
        transitions = [(mid, 9), (mid + 1, 10)]
    elif scenario == 'multiple':
        transitions = [(size // 3, 12), (2 * size // 3, 14)]
        alerts = [p for p, _ in transitions]
    for pos, level in transitions:
        levels[pos:] = [float(level)] * (size - pos)
        changes.append(pos)
    if scenario == 'drift':
        levels = [10 + 2 * i / (size - 1) for i in range(size)]
        changes = alerts = None  # No arbitrary discrete truth for a continuous drift.
    weights = [1.0] * size
    values = []
    residual = 0.0
    for i, level in enumerate(levels):
        sigma = 0.12
        if scenario == 'weighted':
            sigma = 0.05 if i % 2 == 0 else 0.8
            weights[i] = 1 / sigma
        if scenario == 'correlated':
            residual = 0.85 * residual + rng.gauss(0, 0.2)
        else:
            residual = rng.gauss(0, sigma)
        value = level + residual
        if scenario == 'quantized':
            value = round(value, 1)
        values.append(value)
    if scenario == 'outlier':
        values[mid] += 8
    if scenario == 'missing':
        for i in range(size):
            if abs(i - mid) <= 2 or rng.random() < 0.1:
                values[i] = None
            if i % 11 == 0:
                weights[i] = None
        weights[2] = 0
    revisions = [100 + (3 if scenario == 'sparse' else 1) * i for i in range(size)]
    return {
        'id': f'{scenario}-n{size}-seed{seed}',
        'scenario': scenario,
        'seed': seed,
        'values': values,
        'weights': weights,
        'revisions': revisions,
        'true_levels': levels,
        'true_boundaries': changes,
        'true_alerts': alerts,
        'split': 'development' if seed % 2 == 0 else 'held_out',
    }


def prepare(module, case):
    """Use detect_steps itself for filtering and normalization, capturing its input."""
    captured = {}
    original = module.solve_potts_autogamma

    def capture(y, w):
        captured.update(values=y, weights=w)
        return [], [], [], None

    module.solve_potts_autogamma = capture
    try:
        module.detect_steps(case['values'], case['weights'])
    finally:
        module.solve_potts_autogamma = original
    indices = [
        i
        for i, (y, w) in enumerate(zip(case['values'], case['weights']))
        if y is not None and not math.isnan(y) and (w is None or not w <= 0)
    ]
    assert len(indices) == len(captured['values'])
    return captured['values'], captured['weights'], indices


def direct_interval(y, w):
    """Independent oracle: try each observed value, without a median algorithm."""
    if not y:
        raise ValueError('An interval must contain an observation')
    candidates = [(sum(ww * abs(yy - level) for yy, ww in zip(y, w)), level) for level in set(y)]
    cost, level = min(candidates)
    return level, cost


def exhaustive_fit(y, w, gamma, min_size=1, max_size=None, min_pos=0, max_pos=None):
    """Enumerate every valid partition of a tiny series as an independent oracle."""
    end = len(y) if max_pos is None else max_pos
    max_size = len(y) if max_size is None else max_size
    if not 0 <= min_pos < end <= len(y) or not 1 <= min_size <= max_size:
        raise ValueError('Invalid oracle constraints')
    if end - min_pos > 12:
        raise ValueError('Exhaustive oracle is limited to 12 observations')
    best = None
    for mask in range(1 << (end - min_pos - 1)):
        right = [i for i in range(min_pos + 1, end) if mask & (1 << (i - min_pos - 1))]
        right.append(end)
        left = min_pos
        values, costs = [], []
        for r in right:
            if not min_size <= r - left <= max_size:
                break
            value, cost = direct_interval(y[left:r], w[left:r])
            values.append(value)
            costs.append(cost)
            left = r
        else:
            objective = sum(costs) + gamma * (len(right) - 1)
            if best is None or objective < best['objective']:
                best = {'right': right, 'values': values, 'costs': costs, 'objective': objective}
    if best is None:
        raise ValueError('No feasible partition')
    return best


def check_solution(y, w, gamma, solution, min_pos=0, max_pos=None, min_size=1, max_size=None):
    """Recompute returned costs and enforce coverage, bounds, and segment lengths."""
    right, values, costs = solution
    end = len(y) if max_pos is None else max_pos
    max_size = len(y) if max_size is None else max_size
    if not (len(right) == len(values) == len(costs) and right and right[-1] == end):
        raise AssertionError('Invalid output lengths or endpoint')
    left = min_pos
    recomputed = []
    for r, value, reported in zip(right, values, costs):
        if not left < r <= end or not min_size <= r - left <= max_size:
            raise AssertionError('Invalid interval bounds or length')
        actual = sum(ww * abs(yy - value) for yy, ww in zip(y[left:r], w[left:r]))
        if not math.isclose(actual, reported, rel_tol=1e-10, abs_tol=1e-10):
            raise AssertionError('Reported interval cost differs from direct calculation')
        recomputed.append(actual)
        left = r
    return sum(recomputed) + gamma * (len(right) - 1)


def match_boundaries(expected, detected, tolerance=2):
    """Ordered one-to-one matching: maximize matches, then minimize total distance."""
    if expected is None:
        return None
    if tolerance < 0:
        raise ValueError('Tolerance must be nonnegative')
    expected, detected = sorted(expected), sorted(detected)
    # Each state is (number matched, total error, matched position pairs).
    row = [(0, 0, []) for _ in range(len(detected) + 1)]
    for target in expected:
        new = [(0, 0, [])]
        for j, found in enumerate(detected, 1):
            options = [row[j], new[j - 1]]
            if abs(target - found) <= tolerance:
                n, error, pairs = row[j - 1]
                options.append((n + 1, error + abs(target - found), pairs + [[target, found]]))
            new.append(max(options, key=lambda item: (item[0], -item[1])))
        row = new
    matched, error, pairs = row[-1]
    return {
        'matched': matched,
        'missed': len(expected) - matched,
        'false': len(detected) - matched,
        'location_error_sum': error,
        'pairs': pairs,
    }


def expanded_bounds(a, b):
    """The full bracket used by production golden_search(expand_bounds=True)."""
    ratio = 2 / (1 + math.sqrt(5))
    return (
        (ratio * a - (1 - ratio) * b) / (2 * ratio - 1),
        (ratio * b - (1 - ratio) * a) / (2 * ratio - 1),
    )


class SearchExperiment:
    """Intercept only the search around the production objective closure.

    The inner fitter, noise floor, correlation residual formula and complexity
    penalty remain those in solve_potts_autogamma. All patches affect an isolated
    module and are restored after each experiment.
    """

    def __init__(self, module, method, budget=12):
        self.module = module
        self.method = method
        self.budget = budget
        self.trace = []
        self.counts = defaultdict(int)
        self.seconds = defaultdict(float)
        self.outer_active = False
        self.last_fit = None
        self.rho_cache = {}
        self.originals = {}
        self.bounds = None

    def __enter__(self):
        for name in ('golden_search', 'solve_potts_approx', 'solve_potts', 'merge_pieces'):
            self.originals[name] = getattr(self.module, name)
        self.module.golden_search = self.search
        self.module.solve_potts_approx = self.fit
        for name, phase in [('solve_potts', 'dynamic_program'), ('merge_pieces', 'merge')]:
            original = self.originals[name]

            def timed(*args, _original=original, _phase=phase, **kwargs):
                start = time.perf_counter()
                try:
                    return _original(*args, **kwargs)
                finally:
                    self.seconds[_phase] += time.perf_counter() - start

            setattr(self.module, name, timed)
        return self

    def __exit__(self, *args):
        for name, original in self.originals.items():
            setattr(self.module, name, original)

    def fit(self, *args, **kwargs):
        start = time.perf_counter()
        solution = self.originals['solve_potts_approx'](*args, **kwargs)
        self.seconds['candidate_fit'] += time.perf_counter() - start
        self.counts['solver_calls'] += 1
        self.last_fit = {
            'gamma': kwargs['gamma'],
            'right': solution[0],
            'values': solution[1],
            'costs': solution[2],
        }
        return solution

    def search(self, f, a, b, **kwargs):
        native_search = self.originals['golden_search']
        if self.outer_active:
            # Nested call: preserve production rho search, optionally reuse rho for an identical fit.
            key = (tuple(self.last_fit['right']), tuple(self.last_fit['values']))
            if self.method != 'current' and key in self.rho_cache:
                self.counts['rho_cache_hits'] += 1
                return self.rho_cache[key]
            start = time.perf_counter()

            def counted(rho):
                self.counts['rho_objective_evaluations'] += 1
                return f(rho)

            rho = native_search(counted, a, b, **kwargs)
            self.seconds['rho_search'] += time.perf_counter() - start
            self.rho_cache[key] = rho
            return rho

        self.outer_active = True
        try:
            self.bounds = expanded_bounds(a, b)

            def evaluate(x):
                start = time.perf_counter()
                score = f(x)  # Production closure also updates its selected solution.
                self.trace.append(
                    {
                        **self.last_fit,
                        'log_gamma_ratio': x,
                        'score': score,
                        'seconds': time.perf_counter() - start,
                    }
                )
                return score

            if self.method == 'current':
                return native_search(evaluate, a, b, **kwargs)
            if self.method == 'hybrid':
                native_search(evaluate, a, b, **kwargs)
            if self.method not in ('grid', 'hybrid'):
                raise ValueError(f'Unknown search method: {self.method}')
            if self.budget < 4:
                raise ValueError('Grid budget must be at least 4')
            initial = len(self.trace)
            target_calls = initial + self.budget
            lo, hi = self.bounds
            # Spend half the budget on coverage, including production's initial points.
            points = {lo, hi, a, b}
            coarse_count = max(4, self.budget // 2)
            for i in range(1, coarse_count - 1):
                if len(points) >= coarse_count:
                    break
                points.add(lo + (hi - lo) * i / (coarse_count - 1))
            seen = {r['log_gamma_ratio'] for r in self.trace}
            for x in sorted(points):
                if len(self.trace) >= target_calls:
                    break
                if x not in seen:
                    evaluate(x)
                    seen.add(x)
            while len(self.trace) < target_calls:
                samples = sorted(self.trace, key=lambda row: row['log_gamma_ratio'])
                intervals = []
                for left, right in zip(samples, samples[1:]):
                    x = (left['log_gamma_ratio'] + right['log_gamma_ratio']) / 2
                    if x in seen:
                        continue
                    changed = (left['right'], left['values']) != (right['right'], right['values'])
                    width = right['log_gamma_ratio'] - left['log_gamma_ratio']
                    # Refine partition transitions, favoring promising scores; otherwise cover gaps.
                    priority = (
                        changed,
                        -min(left['score'], right['score']) if changed else 0,
                        width,
                    )
                    intervals.append((priority, x))
                if not intervals:
                    break
                _, x = max(intervals)
                evaluate(x)
                seen.add(x)
            return min(self.trace, key=lambda row: row['score'])['log_gamma_ratio']
        finally:
            self.outer_active = False


def peak_rss_bytes():
    try:
        import resource
    except ImportError:
        return None
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == 'darwin' else peak * 1024)


def observed_truth(positions, indices):
    """A change inside a missing-data gap is located at the next usable observation."""
    if positions is None:
        return None
    result = []
    for pos in positions:
        after = next((i for i in indices if i >= pos), None)
        if after is not None and any(i < pos for i in indices):
            result.append(after)
    return sorted(set(result))


def run_case(case, backend, method, budget=12, tolerance=2):
    module = load_detector(backend)
    y, w, indices = prepare(module, case)
    if not y:
        raise ValueError('Measured cases require at least one retained observation')
    scale = module.get_mu_dist(y, w).dist(0, len(y) - 1)
    gamma = case.get('gamma', 3 * scale * math.log(len(y)) / len(y))
    started = time.perf_counter()
    with SearchExperiment(module, method, budget) as experiment:
        if method in ('exact', 'approximate'):
            solver = module.solve_potts if method == 'exact' else module.solve_potts_approx
            result = solver(y, w, gamma=gamma)
            if method == 'exact':
                experiment.counts['solver_calls'] = 1
            selected_score = None
        else:
            right, values, costs, gamma = module.solve_potts_autogamma(y, w)
            result = (right, values, costs)
            selected_score = min(row['score'] for row in experiment.trace)
        seconds = time.perf_counter() - started
    rss = peak_rss_bytes()
    objective = check_solution(y, w, gamma, result)
    right, values, costs = result
    steps, left = [], 0
    for r, value, cost in zip(right, values, costs):
        steps.append(
            (indices[left], indices[r - 1] + 1, value, min(y[left:r]), abs(cost / (r - left)))
        )
        left = r
    # Reproduce graph._compute_graph_steps's coordinate conversion.
    revisions = case['revisions']
    graph_steps = [(revisions[l], revisions[r - 1] + 1, v, m, e) for l, r, v, m, e in steps]
    latest, best, alerts = module.detect_regressions(graph_steps, threshold=0.05)
    revision_positions = {rev: i for i, rev in enumerate(revisions)}
    alert_positions = [revision_positions[a[1]] for a in (alerts or [])]
    detected = [indices[r] for r in right[:-1]]
    truth = observed_truth(case['true_boundaries'], indices)
    truth_alerts = observed_truth(case['true_alerts'], indices)
    return {
        'case_id': case['id'],
        'scenario': case['scenario'],
        'split': case['split'],
        'backend': backend,
        'method': method,
        'size': len(case['values']),
        'retained_size': len(y),
        'seconds': seconds,
        'process_peak_rss_bytes': rss,
        'gamma': gamma,
        'fixed_objective': objective,
        'selection_score': selected_score,
        'right': right,
        'values': values,
        'costs': costs,
        'steps': steps,
        'boundaries': detected,
        'observable_truth': truth,
        'boundary_metrics': match_boundaries(truth, detected, tolerance),
        'alerts': alerts or [],
        'latest': latest,
        'best': best,
        'alert_metrics': match_boundaries(truth_alerts, alert_positions, tolerance),
        'trace': experiment.trace,
        'counts': dict(experiment.counts),
        'phase_seconds': dict(experiment.seconds),
        'search_bounds': experiment.bounds,
        'additional_grid_budget': budget if method == 'hybrid' else None,
    }


def verify_oracle(backend, count=40):
    module = load_detector(backend)
    rng = random.Random(1007)
    failures, checked = [], 0
    for trial in range(count):
        size = rng.randrange(3, 9)
        y = [rng.choice([0.0, 1.0, 3.0]) for _ in range(size)]
        w = [rng.choice([1.0, 2.0, 4.0]) for _ in y]
        gamma = rng.choice([0.0, 0.5, 2.0, 10.0])
        end = size - 1
        configurations = [
            {},
            {'max_size': 2},
            {'min_size': 2},
            {'min_pos': 1, 'max_pos': end},
            {'min_pos': 1, 'max_pos': end, 'min_size': end - 1},
        ]
        for kwargs in configurations:
            try:
                expected = exhaustive_fit(y, w, gamma, **kwargs)
            except ValueError:
                continue
            checked += 1
            try:
                result = module.solve_potts(y, w, gamma, **kwargs)
                actual = check_solution(y, w, gamma, result, **kwargs)
                if not math.isclose(actual, expected['objective'], abs_tol=1e-9, rel_tol=1e-9):
                    raise AssertionError(f'Objective {actual} != {expected["objective"]}')
            except (AssertionError, IndexError, ValueError) as exc:
                failures.append(
                    {
                        'trial': trial,
                        'values': y,
                        'weights': w,
                        'gamma': gamma,
                        'constraints': kwargs,
                        'error': str(exc),
                    }
                )
    return {'backend': backend, 'checked': checked, 'failures': failures}


def wilson(successes, total):
    if not total:
        return None
    z = 1.96
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    radius = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return [center - radius, center + radius]


def summarize(records):
    groups = defaultdict(list)
    for row in records:
        groups[row['backend'], row['method'], row['split']].append(row)
    summary = []
    for (backend, method, split), rows in sorted(groups.items()):
        entry = {
            'backend': backend,
            'method': method,
            'split': split,
            'runs': len(rows),
            'median_seconds': statistics.median(r['seconds'] for r in rows),
            'max_peak_rss_bytes': max((r['process_peak_rss_bytes'] or 0) for r in rows),
            'mean_solver_calls': statistics.mean(r['counts'].get('solver_calls', 0) for r in rows),
        }
        for field in ('boundary_metrics', 'alert_metrics'):
            metrics = [r[field] for r in rows if r[field] is not None]
            total = {
                k: sum(m[k] for m in metrics)
                for k in ('matched', 'missed', 'false', 'location_error_sum')
            }
            false_histories = sum(m['false'] > 0 for m in metrics)
            total.update(
                histories=len(metrics),
                histories_with_false=false_histories,
                false_history_rate=false_histories / len(metrics) if metrics else None,
                false_history_rate_ci95=wilson(false_histories, len(metrics)),
            )
            entry[field] = total
        summary.append(entry)
    comparisons = []
    indexed = {(r['case_id'], r['backend'], r['method']): r for r in records}
    for row in records:
        other_method = 'exact' if row['method'] == 'approximate' else 'current'
        reference = indexed.get((row['case_id'], row['backend'], other_method))
        if reference is None or row['method'] not in ('approximate', 'grid', 'hybrid'):
            continue
        field = 'fixed_objective' if row['method'] == 'approximate' else 'selection_score'
        comparisons.append(
            {
                'case_id': row['case_id'],
                'backend': row['backend'],
                'method': row['method'],
                'reference': other_method,
                'score_delta': row[field] - reference[field],
                'seconds_ratio': row['seconds'] / reference['seconds'],
                'boundaries_changed': row['boundaries'] != reference['boundaries'],
            }
        )
    return {'groups': summary, 'comparisons': comparisons}


def write_report(output, manifest, records):
    summary = summarize(records)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    lines = [
        'Step detection experiment results',
        '=================================',
        '',
        f'Source revision: ``{manifest["git_revision"]}``.',
        '',
        f'Detector SHA256: ``{manifest["detector_sha256"]}``.',
        '',
        f'Command: ``{manifest["command"]}``.',
        '',
        'Methods and interpretation',
        '--------------------------',
        '',
        'Exact and approximate use the same fixed gamma. Current runs the production',
        'automatic search. Grid uses the production score with the same number of',
        'candidate fits as current for that case. Hybrid retains current and adds',
        'that many grid trials. Grid and hybrid reuse correlation searches for',
        'identical fits. Lower selection scores do not necessarily improve detection.',
        '',
        'Even seeds are development data; odd seeds are held out. This first report',
        'evaluates the prespecified methods without tuning them on either split.',
        '',
        'Timing includes search instrumentation but excludes imports and input preparation.',
        'Peak RSS is the whole fresh worker process, including native allocations and',
        'interpreter overhead. Phase timings overlap: candidate fitting contains dynamic',
        'programming and merging. C++ internal cache hits are not instrumented.',
        '',
        'Boundary positions use input observation indices with a fixed matching tolerance.',
        'Changes inside missing-data gaps map to the next retained observation. Alerts',
        'are evaluated after revision mapping and the existing 5 percent reporting rule.',
        'Drift cases have no discrete truth and are excluded from accuracy counts.',
        '',
        'Aggregate results',
        '-----------------',
        '',
        '.. list-table::',
        '   :header-rows: 1',
        '',
        '   * - Backend / method / split',
        '     - Runs',
        '     - Boundary missed / false',
        '     - Alert missed / false',
        '     - Median seconds',
    ]
    for group in summary['groups']:
        bm, am = group['boundary_metrics'], group['alert_metrics']
        lines += [
            f'   * - {group["backend"]} / {group["method"]} / {group["split"]}',
            f'     - {group["runs"]}',
            f'     - {bm["missed"]} / {bm["false"]}',
            f'     - {am["missed"]} / {am["false"]}',
            f'     - {group["median_seconds"]:.6f}',
        ]
    lines += ['', 'Paired comparisons', '------------------', '']
    for backend, method in sorted({(c['backend'], c['method']) for c in summary['comparisons']}):
        pairs = [
            c for c in summary['comparisons'] if (c['backend'], c['method']) == (backend, method)
        ]
        improved = sum(c['score_delta'] < -1e-9 for c in pairs)
        worse = sum(c['score_delta'] > 1e-9 for c in pairs)
        lines += [
            f'* {backend} {method}: {len(pairs)} pairs, {improved} lower scores,',
            f'  {worse} higher scores, {sum(c["boundaries_changed"] for c in pairs)} changed boundaries.',
        ]
    lines += [
        '',
        'Artifacts',
        '---------',
        '',
        '``inputs.json`` contains the complete generated inputs and truth.',
        '``records.jsonl`` contains every result and evaluated penalty with its fit',
        'and production score. ``summary.json`` includes paired comparisons,',
        'aggregate counts, and Wilson 95 percent intervals for false-history rates.',
        '``oracle.json`` records independent small-input correctness checks.',
        '',
    ]
    (output / 'report.rst').write_text('\n'.join(lines))


def call_worker(payload, timeout):
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), '--worker'],
        input=json.dumps(payload),
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
        cwd=ROOT,
    )
    return json.loads(completed.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument(
        '--output', type=Path, default=Path('experiments/step_detection/runs/baseline')
    )
    parser.add_argument(
        '--backends', nargs='+', choices=['python', 'native'], default=['python', 'native']
    )
    parser.add_argument('--scenarios', nargs='+', choices=SCENARIOS, default=list(SCENARIOS))
    parser.add_argument('--sizes', nargs='+', type=int, default=[100])
    parser.add_argument('--seeds', nargs='+', type=int, default=[0, 1, 2, 3])
    parser.add_argument('--methods', nargs='+', choices=METHODS, default=list(METHODS))
    parser.add_argument('--max-exact-size', type=int, default=200)
    parser.add_argument('--tolerance', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=120)
    args = parser.parse_args()
    if args.worker:
        payload = json.load(sys.stdin)
        result = (
            verify_oracle(payload['backend'])
            if payload.get('oracle')
            else run_case(
                payload['case'],
                payload['backend'],
                payload['method'],
                payload['budget'],
                payload['tolerance'],
            )
        )
        print(json.dumps(result, allow_nan=False))
        return
    if args.tolerance < 0 or min(args.sizes) < 12:
        parser.error('Use sizes >= 12 and nonnegative tolerance')
    if any(m in args.methods for m in ('grid', 'hybrid')) and 'current' not in args.methods:
        parser.error('Grid and hybrid require current for a matched per-case budget')
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('Output directory must be empty; use a new run directory')
    manifest = {
        'schema_version': 1,
        'command': ' '.join(sys.argv),
        'python': sys.version,
        'platform': platform.platform(),
        'settings': {**vars(args), 'output': str(args.output)},
        'git_revision': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True
        ).strip(),
        'git_status': subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
        'detector_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'harness_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    checks = [call_worker({'oracle': True, 'backend': b}, args.timeout) for b in args.backends]
    (args.output / 'oracle.json').write_text(json.dumps(checks, indent=2) + '\n')
    if any(check['failures'] for check in checks):
        raise SystemExit(f'Oracle disagreements: inspect {args.output / "oracle.json"}')
    cases = [
        make_case(s, n, seed)
        for s, n, seed in itertools.product(args.scenarios, args.sizes, args.seeds)
    ]
    (args.output / 'inputs.json').write_text(json.dumps(cases, indent=2, allow_nan=False) + '\n')
    records = []
    with (args.output / 'records.jsonl').open('w') as stream:
        for case in cases:
            for backend in args.backends:
                budget = 12
                for method in METHODS:
                    if method not in args.methods:
                        continue
                    if method == 'exact' and len(case['values']) > args.max_exact_size:
                        continue
                    result = call_worker(
                        {
                            'case': case,
                            'backend': backend,
                            'method': method,
                            'budget': budget,
                            'tolerance': args.tolerance,
                        },
                        args.timeout,
                    )
                    if method == 'current':
                        budget = len(result['trace'])
                    stream.write(json.dumps(result, allow_nan=False) + '\n')
                    stream.flush()
                    records.append(result)
                print(f'{case["id"]} / {backend}: complete', flush=True)
    write_report(args.output, manifest, records)
    print(f'Report: {args.output / "report.rst"}')


if __name__ == '__main__':
    main()
