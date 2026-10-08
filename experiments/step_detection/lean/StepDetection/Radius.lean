import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith

namespace StepDetection

/-- The residual radius written using the numerator and denominator of the
    AR(1) precision quadratic form. -/
noncomputable def residualRadius (rho ss numerator denominator : ℝ) : ℝ :=
  (1 - rho ^ 2) * ss * denominator / numerator

/-- A polynomial enclosure for the precision quadratic gives an upper bound
    on the directional radius. The factor `1 + right` is valid because the
    correlation interval is nonnegative and `1 - rho ≤ 1` there. -/
theorem residualRadius_le_of_ratio_bounds
    {rho right ss numerator denominator numeratorLower denominatorUpper : ℝ}
    (hrho : 0 ≤ rho) (hright : rho ≤ right) (hright_le : right ≤ 1)
    (hss : 0 ≤ ss) (hnum_lower : 0 < numeratorLower)
    (hnum : numeratorLower ≤ numerator)
    (hden_nonneg : 0 ≤ denominator) (hden : denominator ≤ denominatorUpper)
    (hden_upper : 0 ≤ denominatorUpper) :
    residualRadius rho ss numerator denominator ≤
      (1 + right) * ss * denominatorUpper / numeratorLower := by
  have hfactor_nonneg : 0 ≤ 1 - rho ^ 2 := by nlinarith
  have hfactor_upper : 1 - rho ^ 2 ≤ 1 + right := by
    nlinarith [mul_nonneg (show 0 ≤ 1 - rho by linarith)
      (show 0 ≤ 1 + rho by linarith)]
  have hcoeff_nonneg : 0 ≤ (1 + right) * ss :=
    mul_nonneg (by linarith) hss
  have hproduct :
      (1 - rho ^ 2) * ss * denominator ≤
        (1 + right) * ss * denominatorUpper := by
    calc
      (1 - rho ^ 2) * ss * denominator ≤
          (1 + right) * ss * denominator := by
            exact mul_le_mul_of_nonneg_right
              (mul_le_mul_of_nonneg_right hfactor_upper hss) hden_nonneg
      _ ≤ (1 + right) * ss * denominatorUpper := by
            exact mul_le_mul_of_nonneg_left hden hcoeff_nonneg
  have hnum_pos : 0 < numerator := lt_of_lt_of_le hnum_lower hnum
  have hupper_nonneg : 0 ≤ (1 + right) * ss * denominatorUpper :=
    mul_nonneg hcoeff_nonneg hden_upper
  unfold residualRadius
  exact (div_le_div_of_nonneg_right hproduct (le_of_lt hnum_pos)).trans
    (div_le_div_of_nonneg_left hupper_nonneg hnum_lower hnum)

end StepDetection
