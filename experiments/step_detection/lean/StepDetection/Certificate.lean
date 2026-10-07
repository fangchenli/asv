import StepDetection.Interval

namespace StepDetection

/-! The probability inequality is an explicit hypothesis, not a new axiom. -/
theorem cutoff_of_squared_bound {p bound upper delta : ℝ}
    (hprob : p ^ 2 ≤ bound) (henclosure : bound ≤ upper ^ 2)
    (hupper : 0 ≤ upper) (hcutoff : upper < delta) : p < delta := by
  have : p ≤ upper := by nlinarith
  exact this.trans_lt hcutoff

theorem reciprocal_square_bound {det correction detLower correctionLower : ℝ}
    (hd : 0 < detLower) (hc : 0 < correctionLower)
    (hdet : detLower ≤ det) (hcorr : correctionLower ≤ correction) :
    1 / (det * correction ^ 2) ≤ 1 / (detLower * correctionLower ^ 2) := by
  apply one_div_le_one_div_of_le (mul_pos hd (sq_pos_of_pos hc))
  apply mul_le_mul hdet
  · nlinarith
  · positivity
  · linarith

end StepDetection
