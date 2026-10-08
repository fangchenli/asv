import StepDetection.ProfiledPrecision

set_option linter.style.haveILetI false
set_option linter.unusedSectionVars false

namespace StepDetection
open Matrix

variable {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]

/-- Completing the square in the first block leaves the Schur-complement
    quadratic in the residual coordinates. The fitted coefficient is the
    one that makes the first block of the precision-weighted residual zero. -/
theorem block_quadratic_at_profiled_coefficients
    (A : Matrix m m ℝ) (B : Matrix m n ℝ)
    (C : Matrix n m ℝ) (D : Matrix n n ℝ) [Invertible A]
    (z : n → ℝ) :
    let beta := -((A⁻¹ * B) *ᵥ z)
    Sum.elim beta z ⬝ᵥ
      (fromBlocks A B C D *ᵥ Sum.elim beta z) =
        z ⬝ᵥ ((D - C * A⁻¹ * B) *ᵥ z) := by
  let beta := -((A⁻¹ * B) *ᵥ z)
  have hAcomp : A *ᵥ ((A⁻¹ * B) *ᵥ z) = B *ᵥ z := by
    calc
      A *ᵥ ((A⁻¹ * B) *ᵥ z) = (A * (A⁻¹ * B)) *ᵥ z :=
        Matrix.mulVec_mulVec z A (A⁻¹ * B)
      _ = B *ᵥ z := by
        rw [← Matrix.mul_assoc, Matrix.mul_inv_of_invertible, Matrix.one_mul]
  have hCcomp : C *ᵥ ((A⁻¹ * B) *ᵥ z) = (C * A⁻¹ * B) *ᵥ z := by
    calc
      C *ᵥ ((A⁻¹ * B) *ᵥ z) = (C * (A⁻¹ * B)) *ᵥ z :=
        Matrix.mulVec_mulVec z C (A⁻¹ * B)
      _ = (C * A⁻¹ * B) *ᵥ z := by rw [Matrix.mul_assoc]
  have htop : A *ᵥ beta + B *ᵥ z = 0 := by
    dsimp [beta]
    rw [Matrix.mulVec_neg, hAcomp, neg_add_cancel]
  have hbottom : C *ᵥ beta + D *ᵥ z = (D - C * A⁻¹ * B) *ᵥ z := by
    dsimp [beta]
    rw [Matrix.mulVec_neg, hCcomp, Matrix.sub_mulVec]
    ext i
    simp only [Pi.add_apply, Pi.neg_apply, Pi.sub_apply]
    ring
  change Sum.elim beta z ⬝ᵥ
      (fromBlocks A B C D *ᵥ Sum.elim beta z) = _
  simp [Matrix.fromBlocks_mulVec, sumElim_dotProduct_sumElim, htop, hbottom]

/-- The precision-weighted residual after profiling the plateau coordinates
    is the inverse-covariance quadratic on the orthonormal residual space. -/
theorem profiled_quadratic_eq_projected_covariance_inverse
    (M : Matrix (m ⊕ n) (m ⊕ n) ℝ)
    (X : Matrix (m ⊕ n) m ℝ) (U : Matrix (m ⊕ n) n ℝ)
    (G : Matrix m m ℝ)
    [Invertible M] [Invertible (fromCols X U)] [Invertible G]
    [Invertible ((fromCols X U)ᵀ * M * (fromCols X U)).toBlocks₁₁]
    [Invertible (Uᵀ * M⁻¹ * U)]
    (hgram : (fromCols X U)ᵀ * fromCols X U =
      fromBlocks G 0 0 (1 : Matrix n n ℝ))
    (z : n → ℝ) :
    let H := (fromCols X U)ᵀ * M * (fromCols X U)
    let beta := -((H.toBlocks₁₁⁻¹ * H.toBlocks₁₂) *ᵥ z)
    Sum.elim beta z ⬝ᵥ (H *ᵥ Sum.elim beta z) =
      z ⬝ᵥ ((Uᵀ * M⁻¹ * U)⁻¹ *ᵥ z) := by
  let H := (fromCols X U)ᵀ * M * (fromCols X U)
  let beta := -((H.toBlocks₁₁⁻¹ * H.toBlocks₁₂) *ᵥ z)
  have hprofile := compressed_covariance_inverse_eq_profiled_precision
    M X U G hgram
  have hquad := block_quadratic_at_profiled_coefficients H.toBlocks₁₁ H.toBlocks₁₂
    H.toBlocks₂₁ H.toBlocks₂₂ z
  have hprof : (Uᵀ * M⁻¹ * U)⁻¹ =
      H.toBlocks₂₂ - H.toBlocks₂₁ * H.toBlocks₁₁⁻¹ * H.toBlocks₁₂ := by
    simpa [H] using hprofile
  calc
    Sum.elim beta z ⬝ᵥ (H *ᵥ Sum.elim beta z) =
        z ⬝ᵥ ((H.toBlocks₂₂ - H.toBlocks₂₁ * H.toBlocks₁₁⁻¹ *
          H.toBlocks₁₂) *ᵥ z) := by
            rw [show H = fromBlocks H.toBlocks₁₁ H.toBlocks₁₂ H.toBlocks₂₁ H.toBlocks₂₂ by
              exact (Matrix.fromBlocks_toBlocks H).symm]
            change Sum.elim beta z ⬝ᵥ
              (fromBlocks H.toBlocks₁₁ H.toBlocks₁₂ H.toBlocks₂₁ H.toBlocks₂₂ *ᵥ
                Sum.elim beta z) = _
            simpa [beta] using hquad
    _ = z ⬝ᵥ ((Uᵀ * M⁻¹ * U)⁻¹ *ᵥ z) := by rw [← hprof]

end StepDetection
