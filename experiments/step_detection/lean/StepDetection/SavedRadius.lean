import StepDetection.Bernstein
import StepDetection.Radius
import Mathlib.Tactic

/- Generated radius data for the frozen negative-directional n=100 example.
   Diagnosis SHA256: ed44a41ad3108fddc7f8066b2f62362cba222e1a6c93500ea637907f96d896b6
   Source archive SHA256: ae567451c2cce383cbdea0427ae093256ddf7256fd5a30758d9d3fc5ddae2d3a -/

namespace StepDetection.SavedRadius

noncomputable section

def left : ℝ := (52223 : ℝ) / 65536
def right : ℝ := (52225 : ℝ) / 65536
def residualSS : ℝ :=
  (4337245619167860129096951655 : ℝ) / 239367312283696576591822848
def numerator (rho : ℝ) : ℝ :=
  (4337245619167860129096951655 : ℝ) / 2417851639229258349412352
    - (577924052090121653902176239771249 : ℝ) / 158456325028528675187087900672 * rho
    + (568172494753268793708650004473907 : ℝ) / 158456325028528675187087900672 * rho ^ 2
    - (269249086060844397951836192221989 : ℝ) / 158456325028528675187087900672 * rho ^ 3
def denominator (rho : ℝ) : ℝ := 99 - 97 * rho
def parameter (t : ℝ) : ℝ := left + (right - left) * t

def numeratorBernstein (i : ℕ) : ℝ :=
  if i = 0 then
    (13585128148614303003326318082147249248056522533 : ℝ) /
      44601490397061246283071436545296723011960832
  else if i = 1 then
    (40753792522813802634153913575391570910164447889 : ℝ) /
      133804471191183738849214309635890169035882496
  else if i = 2 then
    (40752200579995766466664043631436843252060026223 : ℝ) /
      133804471191183738849214309635890169035882496
  else if i = 3 then
    (13583536205794112843147961383008906900414324955 : ℝ) /
      44601490397061246283071436545296723011960832
  else 0

def denominatorBernstein (i : ℕ) : ℝ :=
  if i = 0 then (1422433 : ℝ) / 65536
  else if i = 1 then (1422239 : ℝ) / 65536
  else 0

def numeratorLower : ℝ :=
  (13583536205794112843147961383008906900414324955 : ℝ) /
    44601490397061246283071436545296723011960832
def numeratorUpper : ℝ :=
  (13585128148614303003326318082147249248056522533 : ℝ) /
    44601490397061246283071436545296723011960832
def denominatorUpper : ℝ := (1422433 : ℝ) / 65536
def radiusUpper : ℝ :=
  (2958740906340992790304509298794373175588860052530009324199 : ℝ) /
    6277101735386680763835789423207666416102355444464034512896

theorem numerator_bernstein_identity (t : ℝ) :
    numerator (parameter t) = (bernsteinExpansion 3 numeratorBernstein).eval t := by
  norm_num [numerator, parameter, left, right, bernsteinExpansion,
    numeratorBernstein, bernsteinPolynomial, Nat.choose, Polynomial.eval_finsetSum,
    Polynomial.eval_mul, Polynomial.eval_C, Finset.sum_range_succ]
  ring_nf

theorem denominator_bernstein_identity (t : ℝ) :
    denominator (parameter t) = (bernsteinExpansion 1 denominatorBernstein).eval t := by
  norm_num [denominator, parameter, left, right, bernsteinExpansion,
    denominatorBernstein, bernsteinPolynomial, Nat.choose, Polynomial.eval_finsetSum,
    Polynomial.eval_mul, Polynomial.eval_C, Finset.sum_range_succ]
  ring_nf

theorem numerator_lower_bound {t : ℝ} (ht : 0 ≤ t) (htone : t ≤ 1) :
    numeratorLower ≤ numerator (parameter t) := by
  rw [numerator_bernstein_identity]
  have hcoeff : ∀ i ∈ Finset.range 4,
      numeratorLower ≤ numeratorBernstein i ∧ numeratorBernstein i ≤ numeratorUpper := by
    intro i hi
    simp only [Finset.mem_range] at hi
    interval_cases i <;> norm_num [numeratorLower, numeratorUpper, numeratorBernstein]
  exact (bernsteinExpansion_bounds (lower := numeratorLower) (upper := numeratorUpper)
    ht htone hcoeff).1

/-- The saved ratio `numerator / denominator` is the precision-weighted
    residual quadratic; the radius uses its reciprocal. -/
theorem saved_radius_as_precision_ratio {t : ℝ} (ht : 0 ≤ t) (htone : t ≤ 1) :
    residualRadius (parameter t) residualSS (numerator (parameter t))
        (denominator (parameter t)) =
      (1 - (parameter t) ^ 2) * residualSS /
        (numerator (parameter t) / denominator (parameter t)) := by
  have hrho : parameter t ≤ right := by
    unfold parameter left right
    norm_num
    nlinarith
  have hrho_one : right < 1 := by norm_num [right]
  have hden : 0 < denominator (parameter t) := by
    rw [denominator]
    nlinarith
  have hnum : 0 < numerator (parameter t) := by
    exact lt_of_lt_of_le (by norm_num [numeratorLower])
      (numerator_lower_bound ht htone)
  unfold residualRadius
  field_simp [ne_of_gt hden, ne_of_gt hnum]

/-- At zero correlation, the saved scale is the ordinary residual sum of
    squares, expressed as the same precision numerator/denominator ratio. -/
theorem saved_residualSS_eq_zero_precision_ratio :
    residualSS = numerator 0 / denominator 0 := by
  norm_num [residualSS, numerator, denominator]

theorem denominator_upper_bound {t : ℝ} (ht : 0 ≤ t) (htone : t ≤ 1) :
    denominator (parameter t) ≤ denominatorUpper := by
  rw [denominator_bernstein_identity]
  have hcoeff : ∀ i ∈ Finset.range 2,
      0 ≤ denominatorBernstein i ∧ denominatorBernstein i ≤ denominatorUpper := by
    intro i hi
    simp only [Finset.mem_range] at hi
    interval_cases i <;> norm_num [denominatorUpper, denominatorBernstein]
  exact (bernsteinExpansion_bounds (lower := 0) (upper := denominatorUpper)
    ht htone hcoeff).2

theorem denominator_lower_bound {t : ℝ} (ht : 0 ≤ t) (htone : t ≤ 1) :
    0 ≤ denominator (parameter t) := by
  rw [denominator_bernstein_identity]
  have hcoeff : ∀ i ∈ Finset.range 2,
      0 ≤ denominatorBernstein i ∧ denominatorBernstein i ≤ denominatorUpper := by
    intro i hi
    simp only [Finset.mem_range] at hi
    interval_cases i <;> norm_num [denominatorUpper, denominatorBernstein]
  exact (bernsteinExpansion_bounds (lower := 0) (upper := denominatorUpper)
    ht htone hcoeff).1

theorem saved_radius_bound {t : ℝ} (ht : 0 ≤ t) (htone : t ≤ 1) :
    residualRadius (parameter t) residualSS (numerator (parameter t))
      (denominator (parameter t)) ≤ radiusUpper := by
  have hleft_nonneg : 0 ≤ left := by norm_num [left]
  have hrho_left : left ≤ parameter t := by
    unfold parameter left right
    norm_num
    nlinarith
  have hrho_right : parameter t ≤ right := by
    unfold parameter left right
    norm_num
    nlinarith
  have hrho_le : right ≤ 1 := by norm_num [right]
  have hrho_le_one : parameter t ≤ 1 := hrho_right.trans hrho_le
  have hnumlo : 0 < numeratorLower := by norm_num [numeratorLower]
  have hnum : numeratorLower ≤ numerator (parameter t) := numerator_lower_bound ht htone
  have hdenlo : 0 ≤ denominator (parameter t) := denominator_lower_bound ht htone
  have hden : denominator (parameter t) ≤ denominatorUpper :=
    denominator_upper_bound ht htone
  have hdenupper : 0 ≤ denominatorUpper := by norm_num [denominatorUpper]
  have hss : 0 ≤ residualSS := by norm_num [residualSS]
  have hgeneric := residualRadius_le_of_interval_ratio_bounds hleft_nonneg hrho_left
    hrho_le_one
    hss hnumlo hnum hdenlo hden hdenupper
  calc
    residualRadius (parameter t) residualSS (numerator (parameter t))
        (denominator (parameter t)) ≤
      (1 - left ^ 2) * residualSS * denominatorUpper / numeratorLower := hgeneric
    _ ≤ radiusUpper := by norm_num [right, residualSS, denominatorUpper,
      numeratorLower, radiusUpper, left]

end
end StepDetection.SavedRadius
