"""Export exact sufficient statistics for the saved one-split GLS fit."""

# Lean's polynomial notation is intentional.
# ruff: noqa: RUF001

import gzip
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STUDY = ROOT.parent
sys.path.insert(0, str(STUDY.parent.parent))

from experiments.step_detection import reporting_ar1 as ar1  # noqa: E402
from experiments.step_detection import reporting_ar1_full as full  # noqa: E402

DIAGNOSIS = STUDY / 'data' / 'sturm_v1_diagnosis.json'
OUTPUT = ROOT / 'StepDetection' / 'SavedGLS.lean'
CASE_ID = 'negative-directional-ar1-study-n100-d0.08-s0.02-early-seed1601'


def qreal(value):
    value = F(value)
    if value.denominator == 1:
        return f'({value.numerator} : ℝ)'
    if value.numerator < 0:
        return f'-(({abs(value.numerator)} : ℝ) / {value.denominator})'
    return f'({value.numerator} : ℝ) / {value.denominator}'


def iexpr(coefficients):
    terms = []
    for degree, coefficient in enumerate(coefficients):
        coefficient = F(coefficient)
        if coefficient.denominator != 1:
            raise ValueError('Scaled sufficient-statistic coefficient is not integral')
        coefficient = coefficient.numerator
        if not coefficient:
            continue
        term = f'{abs(coefficient)}'
        if degree:
            term += f' * rho ^ {degree}'
        terms.append(('-' if coefficient < 0 else '+', term))
    if not terms:
        return '0'
    sign, first = terms[0]
    expression = ('-' if sign == '-' else '') + first
    for sign, term in terms[1:]:
        expression += f' {sign} {term}'
    return expression


def scalar_expression(coefficients, variable='rho'):
    terms = []
    for degree, coefficient in enumerate(coefficients):
        coefficient = F(coefficient)
        if not coefficient:
            continue
        term = qreal(abs(coefficient))
        if degree:
            term += f' * {variable} ^ {degree}'
        terms.append(('-' if coefficient < 0 else '+', term))
    if not terms:
        return '0'
    else:
        sign, first = terms[0]
        expression = ('-' if sign == '-' else '') + first
        for sign, term in terms[1:]:
            expression += f' {sign} {term}'
    return expression


def render():
    diagnosis_bytes = DIAGNOSIS.read_bytes()
    diagnosis = json.loads(diagnosis_bytes)
    for name, digest in diagnosis['source_hashes'].items():
        if hashlib.sha256((STUDY / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Frozen source changed: {name}')
    archive = diagnosis['source_archive']
    archive_bytes = (STUDY / 'data' / archive['file']).read_bytes()
    if hashlib.sha256(archive_bytes).hexdigest() != archive['sha256']:
        raise ValueError('Frozen history archive hash mismatch')
    histories = json.loads(gzip.decompress(archive_bytes))
    history = next(
        row for row in histories['rows'] if row['case']['id'] == CASE_ID
    )
    row = next(row for row in diagnosis['rows'] if row['id'] == CASE_ID)
    split = row['witness']['split']
    models = ar1.Models(history['case']['values'])
    columns = [
        tuple(int(i == 0) for i in range(models.n)),
        tuple(int(i >= split) for i in range(models.n)),
    ]
    gram = [[ar1.precision_product(x, y) for y in columns] for x in columns]
    linear = [ar1.precision_product(x, models.values) for x in columns]
    saved_num, saved_den = full.residual_ratio(models, split)
    scale = 2**49
    saved_values = [int(F(value) * scale) for value in history['case']['values']]
    saved_values_lean = ', '.join(str(value) for value in saved_values)
    rational_stats = [
        ('total', models.total),
        ('gram00', gram[0][0]),
        ('gram01', gram[0][1]),
        ('gram11', gram[1][1]),
        ('linear0', linear[0]),
        ('linear1', linear[1]),
    ]
    rational_defs = '\n'.join(
        f'def saved{name.title()}Scaled (rho : ℤ) : ℤ := '
        f'{iexpr([F(coefficient) * scale**2 for coefficient in coefficients])}'
        for name, coefficients in rational_stats
    )

    definitions = [
        f'noncomputable def fitTotal (rho : ℝ) : ℝ := {scalar_expression(models.total)}',
        f'noncomputable def fitGram00 (_rho : ℝ) : ℝ := {scalar_expression(gram[0][0])}',
        f'noncomputable def fitGram01 (rho : ℝ) : ℝ := {scalar_expression(gram[0][1])}',
        f'noncomputable def fitGram11 (rho : ℝ) : ℝ := {scalar_expression(gram[1][1])}',
        f'noncomputable def fitLinear0 (rho : ℝ) : ℝ := {scalar_expression(linear[0])}',
        f'noncomputable def fitLinear1 (rho : ℝ) : ℝ := {scalar_expression(linear[1])}',
        f'noncomputable def fitSavedNumerator (rho : ℝ) : ℝ := {scalar_expression(saved_num)}',
        f'noncomputable def fitSavedDenominator (rho : ℝ) : ℝ := {scalar_expression(saved_den)}',
        'noncomputable def fitGramDet (rho : ℝ) : ℝ := fitGram00 rho * fitGram11 rho - fitGram01 rho ^ 2',
        'noncomputable def fitProfiledNumerator (rho : ℝ) : ℝ := fitTotal rho * fitGramDet rho -\n'
        '  (fitGram11 rho * fitLinear0 rho ^ 2 - 2 * fitGram01 rho *\n'
        '    fitLinear0 rho * fitLinear1 rho + fitGram00 rho * fitLinear1 rho ^ 2)',
    ]
    return f'''-- Generated by export_gls.py; do not edit by hand.
-- Case: {CASE_ID}; split: {split}
-- Diagnosis SHA256: {hashlib.sha256(diagnosis_bytes).hexdigest()}
-- Source archive SHA256: {archive['sha256']}
import StepDetection.SavedRadius
import Mathlib.Algebra.Polynomial.Eval.Algebra
import Mathlib.Tactic

set_option maxRecDepth 4096
set_option maxHeartbeats 5000000

namespace StepDetection.DeterminantTrace
open Polynomial

open scoped BigOperators

def observationScale : ℤ := {scale}
def savedObservationsI : Fin 100 → ℤ := ![{saved_values_lean}]
def savedEntryI (i : ℕ) : ℤ := if h : i < 100 then savedObservationsI ⟨i, h⟩ else 0
def savedPlateau0I (i : ℕ) : ℤ := if i = 0 then observationScale else 0
def savedPlateau1I (i : ℕ) : ℤ := if i = 0 then 0 else observationScale
def precisionProductI (rho : ℤ) (a b : ℕ → ℤ) : ℤ :=
  (∑ i ∈ Finset.range 100, a i * b i) -
    rho * (∑ i ∈ Finset.range 99, (a i * b (i + 1) + a (i + 1) * b i)) +
    rho ^ 2 * (∑ i ∈ Finset.range 98, a (i + 1) * b (i + 1))

{rational_defs}

theorem saved_total_from_observations (rho : ℤ) :
    precisionProductI rho savedEntryI savedEntryI = savedTotalScaled rho := by
  norm_num [precisionProductI, savedTotalScaled, savedEntryI, savedObservationsI,
    Finset.sum_range_succ] <;> ring

theorem saved_gram00_from_observations (rho : ℤ) :
    precisionProductI rho savedPlateau0I savedPlateau0I = savedGram00Scaled rho := by
  norm_num [precisionProductI, savedGram00Scaled, savedPlateau0I, observationScale,
    Finset.sum_range_succ] <;> ring

theorem saved_gram01_from_observations (rho : ℤ) :
    precisionProductI rho savedPlateau0I savedPlateau1I = savedGram01Scaled rho := by
  norm_num [precisionProductI, savedGram01Scaled, savedPlateau0I, savedPlateau1I,
    observationScale,
    Finset.sum_range_succ] <;> ring

theorem saved_gram11_from_observations (rho : ℤ) :
    precisionProductI rho savedPlateau1I savedPlateau1I = savedGram11Scaled rho := by
  norm_num [precisionProductI, savedGram11Scaled, savedPlateau1I, observationScale,
    Finset.sum_range_succ] <;> ring

theorem saved_linear0_from_observations (rho : ℤ) :
    precisionProductI rho savedPlateau0I savedEntryI = savedLinear0Scaled rho := by
  norm_num [precisionProductI, savedLinear0Scaled, savedPlateau0I, savedEntryI,
    savedObservationsI, observationScale, Finset.sum_range_succ] <;> ring

theorem saved_linear1_from_observations (rho : ℤ) :
    precisionProductI rho savedPlateau1I savedEntryI = savedLinear1Scaled rho := by
  norm_num [precisionProductI, savedLinear1Scaled, savedPlateau1I, savedEntryI,
    savedObservationsI, observationScale, Finset.sum_range_succ] <;> ring

{chr(10).join(definitions)}

theorem fit_gram_determinant_identity (rho : ℝ) :
    fitGramDet rho = (1 - rho) * SavedRadius.denominator rho := by
  norm_num [fitGramDet, fitGram00, fitGram01, fitGram11, SavedRadius.denominator]
  <;> ring

theorem fit_profiled_numerator_identity (rho : ℝ) :
    fitProfiledNumerator rho = (1 - rho) * SavedRadius.numerator rho := by
  norm_num [fitProfiledNumerator, fitGramDet, fitTotal, fitGram00, fitGram01,
    fitGram11, fitLinear0, fitLinear1, fitSavedNumerator, SavedRadius.numerator]
  <;> ring

theorem saved_numerator_eval (rho : ℝ) :
    fitSavedNumerator rho = SavedRadius.numerator rho := by
  norm_num [fitSavedNumerator, SavedRadius.numerator]

theorem saved_denominator_eval (rho : ℝ) :
    fitSavedDenominator rho = SavedRadius.denominator rho := by
  norm_num [fitSavedDenominator, SavedRadius.denominator]

/-- The archived polynomial ratio is the profiled GLS residual quadratic
    computed from the exact saved sufficient statistics. -/
theorem saved_profiled_quadratic_eq_ratio {{rho : ℝ}}
    (hrho : rho ≠ 1) (hden : SavedRadius.denominator rho ≠ 0) :
    fitProfiledNumerator rho / fitGramDet rho =
      SavedRadius.numerator rho / SavedRadius.denominator rho := by
  rw [fit_profiled_numerator_identity, fit_gram_determinant_identity]
  field_simp [hrho, hden]

end StepDetection.DeterminantTrace
'''


def main():
    expected = render()
    if len(sys.argv) > 1 and sys.argv[1] == '--check':
        if not OUTPUT.exists() or OUTPUT.read_text() != expected:
            raise SystemExit('SavedGLS.lean is stale; run export_gls.py')
        print('Saved GLS sufficient statistics match the frozen history.')
    else:
        OUTPUT.write_text(expected)
        print(f'Wrote {OUTPUT.name}')


if __name__ == '__main__':
    main()
