import StepDetection.MatrixIdentity

namespace StepDetection
open Matrix

variable {n k : Type*} [Fintype n] [DecidableEq n]

/-- Moving a covariance shift to the precision matrix preserves this factorization. -/
theorem precision_shift_factor (A T : Matrix n n ℝ) (a c : ℝ) [Invertible A]
    (hT : T = a • A + c • (1 : Matrix n n ℝ)) :
    a • (1 : Matrix n n ℝ) + c • A⁻¹ = T * A⁻¹ := by
  rw [hT, Matrix.add_mul, Matrix.smul_mul, Matrix.smul_mul,
    Matrix.mul_inv_of_invertible, Matrix.one_mul]

theorem precision_shift_det (A T : Matrix n n ℝ) (a c : ℝ) [Invertible A]
    (hT : T = a • A + c • (1 : Matrix n n ℝ)) :
    (a • (1 : Matrix n n ℝ) + c • A⁻¹).det = T.det / A.det := by
  rw [precision_shift_factor A T a c hT, Matrix.det_mul, Matrix.det_nonsing_inv,
    Ring.inverse_eq_inv, div_eq_mul_inv]

theorem precision_shift_inverse (A T : Matrix n n ℝ) (a c : ℝ)
    [Invertible A] [Invertible T] (ha : a ≠ 0)
    (hT : T = a • A + c • (1 : Matrix n n ℝ)) :
    (a • (1 : Matrix n n ℝ) + c • A⁻¹)⁻¹ =
      a⁻¹ • ((1 : Matrix n n ℝ) - c • T⁻¹) := by
  have hinv : (a • (1 : Matrix n n ℝ) + c • A⁻¹)⁻¹ = A * T⁻¹ := by
    apply Matrix.inv_eq_right_inv
    rw [precision_shift_factor A T a c hT, Matrix.mul_assoc,
      Matrix.inv_mul_cancel_left_of_invertible, Matrix.mul_inv_of_invertible]
  have h := congrArg (fun B : Matrix n n ℝ => B * T⁻¹) hT
  simp only [Matrix.mul_inv_of_invertible, Matrix.add_mul, Matrix.smul_mul,
    Matrix.one_mul] at h
  have hdiff : (1 : Matrix n n ℝ) - c • T⁻¹ = a • (A * T⁻¹) := by
    rw [h, add_sub_cancel_right]
  rw [hinv, hdiff, smul_smul, inv_mul_cancel₀ ha, one_smul]

theorem precision_shift_solve (A T : Matrix n n ℝ) (a c : ℝ)
    (X Y : Matrix n k ℝ) [Invertible A] [Invertible T] (ha : a ≠ 0)
    (hT : T = a • A + c • (1 : Matrix n n ℝ)) (hY : T * Y = X) :
    (a • (1 : Matrix n n ℝ) + c • A⁻¹)⁻¹ * X = a⁻¹ • (X - c • Y) := by
  rw [precision_shift_inverse A T a c ha hT, Matrix.smul_mul, Matrix.sub_mul,
    Matrix.one_mul, Matrix.smul_mul, ← hY, Matrix.inv_mul_cancel_left_of_invertible]

end StepDetection
