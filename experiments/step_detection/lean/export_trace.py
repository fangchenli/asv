"""Export the archived determinant's interval operations as compositional Lean proofs."""

# Lean type notation in the generated source is intentional.
# ruff: noqa: RUF001

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent
STUDY = ROOT.parent
sys.path.insert(0, str(STUDY.parent.parent))

from experiments.step_detection import directional_determinant as determinant  # noqa: E402

CASE_ID = 'negative-directional-ar1-study-n100-d0.08-s0.02-early-seed1601'
CHUNK_SIZE = 64


class Node(tuple):
    def __new__(cls, bounds, index):
        obj = super().__new__(cls, bounds)
        obj.index = index
        return obj


class Recorder:
    def __init__(self, bits):
        self.arithmetic = determinant.Arithmetic(bits)
        self.nodes = []
        self.cache = {}
        self.calls = 0

    def record(self, op, args, bounds, value=None):
        self.calls += 1
        key = (op, tuple(a.index for a in args), value)
        if key in self.cache:
            node = self.cache[key]
            if tuple(node) != bounds:
                raise ValueError('Repeated expression has inconsistent bounds')
            return node
        node = Node(bounds, len(self.nodes))
        self.nodes.append({'op': op, 'args': key[1], 'value': value, 'bounds': tuple(bounds)})
        self.cache[key] = node
        return node

    def interval(self, lower, upper=None):
        lower = F(lower)
        upper = lower if upper is None else F(upper)
        op = 'constant' if lower == upper else 'input'
        return self.record(op, (), self.arithmetic.interval(lower, upper), (lower, upper))

    def add(self, x, y):
        return self.record('add', (x, y), self.arithmetic.add(x, y))

    def sub(self, x, y):
        return self.record('sub', (x, y), self.arithmetic.sub(x, y))

    def mul(self, x, y):
        return self.record('mul', (x, y), self.arithmetic.mul(x, y))

    def div(self, x, y):
        if y[0] <= 0 <= y[1]:
            raise ArithmeticError('Denominator enclosure contains zero')
        inverse = self.record('inv', (y,), self.arithmetic.interval(1 / y[1], 1 / y[0]))
        result = self.mul(x, inverse)
        if tuple(result) != self.arithmetic.div(x, y):
            raise ValueError('Division decomposition differs from the original operation')
        return result

    def square(self, x):
        return self.record('square', (x,), self.arithmetic.square(x))


def fraction(q):
    return f'({q.numerator} / {q.denominator} : ℝ)'


def collect(*, capture=None):
    archive = STUDY / 'data' / 'sturm_v1_diagnosis.json'
    raw = archive.read_bytes()
    diagnosis = json.loads(raw)
    for name, digest in diagnosis['source_hashes'].items():
        if hashlib.sha256((STUDY / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Frozen research source changed: {name}')
    row = next(r for r in diagnosis['rows'] if r['id'] == CASE_ID)
    cert = row['interval_certificate']
    attempt = next(a for a in cert['attempts'] if a['tilt'] == cert['selected_tilt'])
    source = diagnosis['source_archive']
    source_bytes = (STUDY / 'data' / source['file']).read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != source['sha256']:
        raise ValueError('Saved history archive hash mismatch')
    histories = json.loads(gzip.decompress(source_bytes))
    history = next(r for r in histories['rows'] if r['case']['id'] == CASE_ID)
    n, split = len(history['case']['values']), row['witness']['split']
    left, right = F(cert['left']), F(cert['right'])
    tilt, radius = F(cert['selected_tilt']), F(attempt['base']['radius_upper'])
    a, b = 1 - 2 * tilt, 2 * tilt / radius
    recorder = Recorder(cert['bits'])
    expected = tuple(map(F, attempt['base']['determinant_enclosure']))
    # Run the unchanged routine once normally and once with a recording arithmetic object.
    original = determinant.determinant_interval(n, split, left, right, a, b, bits=cert['bits'])
    with patch.object(determinant, 'Arithmetic', return_value=recorder):
        if capture is None:
            result = determinant.determinant_interval(
                n, split, left, right, a, b, bits=cert['bits']
            )
        else:
            previous_profile = sys.getprofile()

            def profile(frame, event, _arg):
                if event != 'return' or frame.f_code.co_filename != determinant.__file__:
                    return
                if frame.f_code.co_name == 'determinant_interval':
                    capture.update(frame.f_locals)
                elif frame.f_code.co_name == 'solve':
                    capture.setdefault('solves', []).append(dict(frame.f_locals))

            try:
                sys.setprofile(profile)
                result = determinant.determinant_interval(
                    n, split, left, right, a, b, bits=cert['bits']
                )
            finally:
                sys.setprofile(previous_profile)
    if tuple(result) != expected or original != expected:
        raise ValueError('Recomputed determinant differs from the saved enclosure')
    inputs = [node for node in recorder.nodes if node['op'] == 'input']
    if len(inputs) != 1 or inputs[0]['value'] != (left, right):
        raise ValueError('Expected exactly one correlation input')
    metadata = {
        'archive': archive.name,
        'sha256': hashlib.sha256(raw).hexdigest(),
        'case': CASE_ID,
        'n': n,
        'split': split,
        'bits': cert['bits'],
        'left': str(left),
        'right': str(right),
        'a': str(a),
        'b': str(b),
        'result_node': result.index,
        'recorded_calls': recorder.calls,
        'unique_nodes': len(recorder.nodes),
        'operations': dict(sorted(Counter(node['op'] for node in recorder.nodes).items())),
        'determinant_enclosure': list(map(str, expected)),
    }
    return recorder.nodes, metadata


def node_source(index, node):
    op, args = node['op'], node['args']
    lo, hi = node['bounds']
    definition = f'noncomputable def b{index} : Interval := ⟨{fraction(lo)}, {fraction(hi)}⟩'
    if op == 'input':
        value = 'rho'
        proof = '  exact hrho'
    elif op == 'constant':
        value = fraction(node['value'][0])
        proof = f'  norm_num [b{index}, v{index}, Interval.Contains]'
    else:
        x = args[0]
        if op in ('add', 'sub', 'mul'):
            y = args[1]
            symbol = {'add': '+', 'sub': '-', 'mul': '*'}[op]
            value = f'v{x} rho {symbol} v{y} rho'
            sound = f'Interval.{op}_sound (h{x} hrho) (h{y} hrho)'
            unfolded = f'b{index}, b{x}, b{y}, Interval.{op}, Interval.Enlarges'
        else:
            value = f'1 / v{x} rho' if op == 'inv' else f'v{x} rho ^ 2'
            sound = f'Interval.{op}_sound (h{x} hrho)'
            if op == 'inv':
                sound += f' (by norm_num [b{x}])'
            unfolded = f'b{index}, b{x}, Interval.{op}, Interval.Enlarges'
        proof = f'  exact Interval.enlarge_sound ({sound})\n    (by norm_num [{unfolded}])'
    parameter = '_rho' if op == 'constant' else 'rho'
    hypothesis = '_hrho' if op == 'constant' else 'hrho'
    return f'''{definition}
noncomputable def v{index} ({parameter} : ℝ) : ℝ := {value}
theorem h{index} {{rho : ℝ}} ({hypothesis} : inputBounds.Contains rho) :
    b{index}.Contains (v{index} rho) := by
{proof}
'''


def render(nodes, metadata):
    header = (
        '-- Generated by export_trace.py; do not edit by hand.\n'
        f'-- Archive SHA256: {metadata["sha256"]}\n'
    )
    files = {}
    left, right = F(metadata['left']), F(metadata['right'])
    files['TraceInputs.lean'] = (
        header
        + f'''import StepDetection.Trace

namespace StepDetection.DeterminantTrace
noncomputable def inputBounds : Interval := ⟨{fraction(left)}, {fraction(right)}⟩
end StepDetection.DeterminantTrace
'''
    )
    previous = 'TraceInputs'
    for start in range(0, len(nodes), CHUNK_SIZE):
        name = f'TraceChunk{start // CHUNK_SIZE:03}'
        body = '\n'.join(
            node_source(i, nodes[i]) for i in range(start, min(start + CHUNK_SIZE, len(nodes)))
        )
        files[f'{name}.lean'] = (
            header + f'import StepDetection.DeterminantTrace.{previous}\n\n'
            'namespace StepDetection.DeterminantTrace\n\n'
            + body
            + '\nend StepDetection.DeterminantTrace\n'
        )
        previous = name
    result = metadata['result_node']
    lo, hi = map(F, metadata['determinant_enclosure'])
    files['Result.lean'] = (
        header
        + f'''import StepDetection.DeterminantTrace.{previous}
import StepDetection.SavedCertificate

namespace StepDetection.DeterminantTrace

-- This expression is the recorded arithmetic circuit; the matrix identity is separate.
noncomputable def determinantExpression (rho : ℝ) : ℝ := v{result} rho

theorem saved_enclosure {{rho : ℝ}} (hrho : inputBounds.Contains rho) :
    {fraction(lo)} ≤ determinantExpression rho ∧
    determinantExpression rho ≤ {fraction(hi)} := h{result} hrho

theorem saved_lower {{rho : ℝ}} (hrho : inputBounds.Contains rho) :
    (Saved.detLower : ℝ) ≤ determinantExpression rho := by
  have h := (saved_enclosure hrho).1
  norm_num [Saved.detLower] at h ⊢
  exact h

theorem cutoff_with_trace {{rho p correction : ℝ}}
    (hrho : inputBounds.Contains rho)
    (hcorr : (Saved.correctionLower : ℝ) ≤ correction)
    (hprob : p ^ 2 ≤ 1 / (determinantExpression rho * correction ^ 2)) :
    p < (1 : ℝ) / 100 :=
  Saved.saved_model_cutoff (saved_lower hrho) hcorr hprob

end StepDetection.DeterminantTrace
'''
    )
    files['manifest.json'] = json.dumps(metadata, indent=2) + '\n'
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    nodes, metadata = collect()
    files = render(nodes, metadata)
    directory = ROOT / 'StepDetection' / 'DeterminantTrace'
    if args.check:
        existing = {p.name for p in directory.iterdir()} if directory.exists() else set()
        if existing != set(files) or any(
            (directory / name).read_text() != text for name, text in files.items()
        ):
            parser.exit(1, 'Determinant trace is stale; run export_trace.py\n')
    else:
        directory.mkdir(exist_ok=True)
        extra = {p.name for p in directory.iterdir()} - set(files)
        if extra:
            parser.exit(
                1, f'Refusing to overwrite a directory with unexpected files: {sorted(extra)}\n'
            )
        for name, text in files.items():
            (directory / name).write_text(text)
    print(f'Determinant trace: {len(nodes)} unique nodes; saved enclosure matches.')


if __name__ == '__main__':
    main()
