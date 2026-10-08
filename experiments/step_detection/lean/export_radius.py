"""Check the saved radius proof data against its frozen exact inputs."""

# Lean's real-number type notation is intentional.
# ruff: noqa: RUF001

import gzip
import hashlib
import json
import re
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STUDY = ROOT.parent
sys.path.insert(0, str(STUDY.parent.parent))

from experiments.step_detection import rational_polynomial as poly  # noqa: E402
from experiments.step_detection import reporting_ar1, reporting_sturm  # noqa: E402

CASE_ID = 'negative-directional-ar1-study-n100-d0.08-s0.02-early-seed1601'
DIAGNOSIS = STUDY / 'data' / 'sturm_v1_diagnosis.json'
LEAN_DATA = ROOT / 'StepDetection' / 'SavedRadius.lean'


def rational_present(source, value):
    value = F(value)
    if value.denominator == 1:
        sign = r'-\s*' if value.numerator < 0 else ''
        return re.search(rf'{sign}(?<!\d){abs(value.numerator)}(?!\d)', source) is not None
    sign = r'-\s*' if value.numerator < 0 else ''
    pattern = rf'{sign}\(\s*{abs(value.numerator)}\s*:\s*ℝ\s*\)\s*/\s*{value.denominator}'
    return re.search(pattern, source) is not None


def check():
    diagnosis_bytes = DIAGNOSIS.read_bytes()
    diagnosis = json.loads(diagnosis_bytes)
    for name, digest in diagnosis['source_hashes'].items():
        if hashlib.sha256((STUDY / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Frozen source changed: {name}')

    archived = diagnosis['source_archive']
    archive_path = STUDY / 'data' / archived['file']
    archive_bytes = archive_path.read_bytes()
    if hashlib.sha256(archive_bytes).hexdigest() != archived['sha256']:
        raise ValueError('Frozen history archive hash mismatch')
    histories = json.loads(gzip.decompress(archive_bytes))
    history = next(row for row in histories['rows'] if row['case']['id'] == CASE_ID)
    row = next(row for row in diagnosis['rows'] if row['id'] == CASE_ID)

    split = row['witness']['split']
    model = reporting_sturm.confidence_state(
        reporting_ar1.Models(history['case']['values']), split
    )
    certificate = row['interval_certificate']
    left, right = F(certificate['left']), F(certificate['right'])
    attempt = next(
        item for item in certificate['attempts'] if item['tilt'] == certificate['selected_tilt']
    )
    radius = F(attempt['base']['radius_upper'])
    num_bernstein = poly.bernstein(model['num'], left, right)
    den_bernstein = poly.bernstein(model['den'], left, right)

    if len(model['num']) != 4 or len(model['den']) != 2:
        raise ValueError('Saved example polynomial degrees changed')
    if model['ss'] != poly.evaluate(model['num'], F(0)) / poly.evaluate(model['den'], F(0)):
        raise ValueError('Saved residual scale does not match the residual quadratic at rho=0')
    source = LEAN_DATA.read_text()
    required = [
        left,
        right,
        model['ss'],
        *model['num'],
        *model['den'],
        *num_bernstein,
        *den_bernstein,
        min(num_bernstein),
        max(den_bernstein),
        radius,
    ]
    missing = [str(value) for value in required if not rational_present(source, value)]
    if missing:
        raise ValueError(f'Lean radius data no longer matches the archived calculation: {missing}')
    if hashlib.sha256(diagnosis_bytes).hexdigest() not in source:
        raise ValueError('Lean radius module lacks the diagnosis provenance hash')
    if archived['sha256'] not in source:
        raise ValueError('Lean radius module lacks the source archive provenance hash')
    print('Saved radius polynomial data matches the frozen calculation.')


if __name__ == '__main__':
    check()
