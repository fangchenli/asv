"""Save six public NumPy histories, chosen by operation before comparing fits."""

import datetime
import hashlib
import json
import re
import urllib.parse
import urllib.request

from .compare_versions import SNAPSHOT, dump

BASE = 'https://pv.github.io/numpy-bench/'
BENCHMARKS = (
    'bench_app.MaxesOfDots.time_it',
    'bench_core.Core.time_arange_100',
    'bench_core.Core.time_array_1',
    'bench_function_base.Sort.time_sort_worst',
    'bench_linalg.Eindot.time_dot_a_b',
    'bench_random.Shuffle.time_100000',
)


def download(url):
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read()


def main():
    if SNAPSHOT.exists():
        raise FileExistsError('Keep the saved snapshot; move it explicitly to refresh')
    raw_index = download(BASE + 'index.json')
    index = json.loads(raw_index)
    params = next(
        p
        for p in index['graph_param_list']
        if p['machine'] == 'i7' and p['python'] == '3.7' and p['Cython'] == '0.29.21'
    )
    parts = []
    for key, value in params.items():
        part = key + ('-null' if value is None else '-' + value if value else '')
        parts.append(re.sub('[<>:"/\\\\^|?*\x00-\x1f]', '_', part))
    prefix = BASE + 'graphs/' + '/'.join(urllib.parse.quote(p) for p in sorted(parts)) + '/'
    cases = []
    for name in BENCHMARKS:
        url = prefix + name + '.json'
        raw = download(url)
        graph = json.loads(raw)
        revisions, values = zip(*graph)
        if list(revisions) != sorted(set(revisions)):
            raise ValueError('Expected sorted distinct revisions')
        if not all(v is None or isinstance(v, (int, float)) for v in values):
            raise ValueError('Expected scalar benchmark')
        cases.append(
            {
                'id': 'numpy-' + name,
                'scenario': name,
                'corpus': 'numpy',
                'split': 'unlabeled',
                'values': values,
                'weights': [1] * len(values),
                'revisions': revisions,
                'true_levels': None,
                'true_boundaries': None,
                'true_alerts': None,
                'provenance': {
                    'url': url,
                    'sha256': hashlib.sha256(raw).hexdigest(),
                    'unit': index['benchmarks'][name]['unit'],
                    'commit_hashes': [index['revision_to_hash'][str(r)] for r in revisions],
                    'dates_ms': [index['revision_to_date'][str(r)] for r in revisions],
                },
            }
        )
        print(name, len(values), flush=True)
    snapshot = {
        'source': BASE,
        'index_url': BASE + 'index.json',
        'index_sha256': hashlib.sha256(raw_index).hexdigest(),
        'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'params': params,
        'selection': 'Six named scalar operations on i7/Python 3.7, all published observations',
        'weights': 'Unit weights: published graphs contain no uncertainty weights',
        'cases': cases,
    }
    SNAPSHOT.parent.mkdir(exist_ok=True)
    dump(SNAPSHOT, snapshot)


if __name__ == '__main__':
    main()
