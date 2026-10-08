import StepDetection.SchurInverse

namespace StepDetection
open Matrix

variable {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]

/-- Orthogonal residual coordinates complete a full-rank design to a basis.
    The determinant of the combined Gram matrix is the design Gram
    determinant, so the combined square matrix has nonzero determinant. -/
theorem residual_coordinates_det_ne_zero
    (X : Matrix (m ⊕ n) m ℝ) (U : Matrix (m ⊕ n) n ℝ)
    (hUU : Uᵀ * U = 1) (hXU : Xᵀ * U = 0)
    (hXX : (Xᵀ * X).det ≠ 0) : (fromCols X U).det ≠ 0 := by
  have hUX : Uᵀ * X = 0 := by
    have h := congrArg Matrix.transpose hXU
    simpa only [Matrix.transpose_mul, Matrix.transpose_transpose,
      Matrix.transpose_zero] using h
  have hgram : (fromCols X U)ᵀ * fromCols X U =
      fromBlocks (Xᵀ * X) 0 0 (1 : Matrix n n ℝ) := by
    rw [Matrix.transpose_fromCols, Matrix.fromRows_mul_fromCols, hUU, hXU, hUX]
  have hdet := congrArg Matrix.det hgram
  simp only [Matrix.det_mul, Matrix.det_transpose, Matrix.det_fromBlocks_zero₁₂,
    Matrix.det_one, mul_one] at hdet
  have hp : (fromCols X U).det ≠ 0 := by
    intro hz
    rw [hz, zero_mul] at hdet
    exact hXX hdet.symm
  exact hp

/-- Profiling out the first block of a precision matrix gives the inverse of
    the covariance compressed to the orthonormal complementary directions.
    The first block's Gram matrix is retained throughout the proof. -/
theorem compressed_covariance_inverse_eq_profiled_precision
    (M : Matrix (m ⊕ n) (m ⊕ n) ℝ)
    (X : Matrix (m ⊕ n) m ℝ) (U : Matrix (m ⊕ n) n ℝ)
    (G : Matrix m m ℝ)
    [Invertible M] [Invertible (fromCols X U)] [Invertible G]
    [Invertible ((fromCols X U)ᵀ * M * (fromCols X U)).toBlocks₁₁]
    [Invertible (Uᵀ * M⁻¹ * U)]
    (hgram : (fromCols X U)ᵀ * fromCols X U =
      fromBlocks G 0 0 (1 : Matrix n n ℝ)) :
    (Uᵀ * M⁻¹ * U)⁻¹ =
      ((fromCols X U)ᵀ * M * (fromCols X U)).toBlocks₂₂ -
        ((fromCols X U)ᵀ * M * (fromCols X U)).toBlocks₂₁ *
          (((fromCols X U)ᵀ * M * (fromCols X U)).toBlocks₁₁)⁻¹ *
            ((fromCols X U)ᵀ * M * (fromCols X U)).toBlocks₁₂ := by
  let P : Matrix (m ⊕ n) (m ⊕ n) ℝ := fromCols X U
  let H : Matrix (m ⊕ n) (m ⊕ n) ℝ := Pᵀ * M * P
  let C : Matrix (m ⊕ n) (m ⊕ n) ℝ := Pᵀ * M⁻¹ * P
  let K : Matrix (m ⊕ n) (m ⊕ n) ℝ := fromBlocks G 0 0 (1 : Matrix n n ℝ)
  let Ki : Matrix (m ⊕ n) (m ⊕ n) ℝ := fromBlocks G⁻¹ 0 0 (1 : Matrix n n ℝ)
  have hPK : Pᵀ * P = K := by simpa [P, K] using hgram
  letI : Invertible (Pᵀ * M) := invertibleMul _ _
  letI : Invertible H := invertibleMul _ _
  have hKschur : (1 : Matrix n n ℝ) - (0 : Matrix n m ℝ) * G⁻¹ *
      (0 : Matrix m n ℝ) = 1 := by simp
  letI : Invertible (1 - (0 : Matrix n m ℝ) * G⁻¹ * (0 : Matrix m n ℝ)) :=
    Invertible.copy invertibleOne _ hKschur
  letI : Invertible (1 - (0 : Matrix n m ℝ) * ⅟G * (0 : Matrix m n ℝ)) := by
    letI : Invertible (1 : Matrix n n ℝ) := invertibleOne
    apply Invertible.copy (inferInstance : Invertible (1 : Matrix n n ℝ))
    simp
  letI : Invertible K := Matrix.fromBlocks₁₁Invertible G 0 0 1
  have hC : C = K * H⁻¹ * K := by
    have hinv : Pᵀ * M⁻¹ * P = (Pᵀ * P) * H⁻¹ * (Pᵀ * P) := by
      dsimp [H]
      simp only [Matrix.mul_inv_rev, Matrix.mul_assoc,
        Matrix.mul_inv_cancel_left_of_invertible,
        Matrix.inv_mul_cancel_left_of_invertible]
    calc
      C = (Pᵀ * P) * H⁻¹ * (Pᵀ * P) := by simpa [C] using hinv
      _ = K * H⁻¹ * K := by rw [hPK]
  have hKKi : K * Ki = 1 := by
    dsimp [K, Ki]
    simp only [Matrix.fromBlocks_multiply, Matrix.mul_zero, Matrix.zero_mul,
      Matrix.one_mul, Matrix.mul_one, add_zero, zero_add]
    simp [Matrix.mul_inv_cancel_left_of_invertible]
  have hKinv : K⁻¹ = Ki := Matrix.inv_eq_right_inv hKKi
  have hCinv : C⁻¹ = Ki * H * Ki := by
    rw [hC, Matrix.mul_inv_rev, Matrix.mul_inv_rev]
    simp only [Matrix.inv_inv_of_invertible, hKinv, Matrix.mul_assoc]
  have hHblocks : Ki * H * Ki = fromBlocks
      (G⁻¹ * H.toBlocks₁₁ * G⁻¹) (G⁻¹ * H.toBlocks₁₂)
      (H.toBlocks₂₁ * G⁻¹) H.toBlocks₂₂ := by
    rw [← Matrix.fromBlocks_toBlocks H]
    dsimp [Ki]
    simp only [Matrix.fromBlocks_multiply, Matrix.mul_zero, Matrix.zero_mul,
      Matrix.one_mul, Matrix.mul_one, add_zero, zero_add]
  have hCinv₁₁ : C⁻¹.toBlocks₁₁ = G⁻¹ * H.toBlocks₁₁ * G⁻¹ := by
    rw [hCinv, hHblocks]
    rfl
  have hCinv₁₂ : C⁻¹.toBlocks₁₂ = G⁻¹ * H.toBlocks₁₂ := by
    rw [hCinv, hHblocks]
    rfl
  have hCinv₂₁ : C⁻¹.toBlocks₂₁ = H.toBlocks₂₁ * G⁻¹ := by
    rw [hCinv, hHblocks]
    rfl
  have hCinv₂₂ : C⁻¹.toBlocks₂₂ = H.toBlocks₂₂ := by
    rw [hCinv, hHblocks]
    rfl
  have hblocks (N : Matrix (m ⊕ n) (m ⊕ n) ℝ) :
      Pᵀ * N * P = fromBlocks (Xᵀ * N * X) (Xᵀ * N * U)
        (Uᵀ * N * X) (Uᵀ * N * U) := by
    dsimp [P]
    rw [Matrix.transpose_fromCols, Matrix.fromRows_mul, Matrix.fromRows_mul_fromCols]
  have hC₂₂ : C.toBlocks₂₂ = Uᵀ * M⁻¹ * U := by
    rw [show C = Pᵀ * M⁻¹ * P by rfl, hblocks]
    rfl
  letI : Invertible C := by
    letI : Invertible M⁻¹ := {
      invOf := M
      invOf_mul_self := by simp [Matrix.invOf_eq_nonsing_inv]
      mul_invOf_self := by simp [Matrix.invOf_eq_nonsing_inv]
    }
    letI : Invertible (Pᵀ * M⁻¹) := invertibleMul _ _
    exact Invertible.copy (invertibleMul (Pᵀ * M⁻¹) P) C (by rfl)
  letI : Invertible C.toBlocks₂₂ :=
    Invertible.copy (inferInstance : Invertible (Uᵀ * M⁻¹ * U)) _ hC₂₂
  letI : Invertible C⁻¹.toBlocks₁₁ := by
    rw [hCinv₁₁]
    letI : Invertible (G⁻¹ * H.toBlocks₁₁) := invertibleMul _ _
    exact invertibleMul _ _
  have hschur := inverse_bottomRight_eq_schur C
  rw [hC₂₂, hCinv₁₁, hCinv₁₂, hCinv₂₁, hCinv₂₂] at hschur
  have hnorm : (G⁻¹ * H.toBlocks₁₁ * G⁻¹)⁻¹ =
      G * H.toBlocks₁₁⁻¹ * G := by
    letI : Invertible (G⁻¹ * H.toBlocks₁₁) := invertibleMul _ _
    letI : Invertible (G⁻¹ * H.toBlocks₁₁ * G⁻¹) := invertibleMul _ _
    apply Matrix.inv_eq_right_inv
    simp only [Matrix.mul_assoc]
    simp [Matrix.mul_inv_cancel_left_of_invertible,
      Matrix.inv_mul_cancel_left_of_invertible]
  have hcancel : (H.toBlocks₂₁ * G⁻¹) *
      (G⁻¹ * H.toBlocks₁₁ * G⁻¹)⁻¹ * (G⁻¹ * H.toBlocks₁₂) =
      H.toBlocks₂₁ * H.toBlocks₁₁⁻¹ * H.toBlocks₁₂ := by
    rw [hnorm]
    simp only [Matrix.mul_assoc]
    simp [Matrix.mul_inv_cancel_left_of_invertible,
      Matrix.inv_mul_cancel_left_of_invertible]
  rw [hcancel] at hschur
  have hHblocks' := hblocks M
  change (Pᵀ * M * P) = _ at hHblocks'
  simpa only [H, P, hHblocks'] using hschur

end StepDetection
