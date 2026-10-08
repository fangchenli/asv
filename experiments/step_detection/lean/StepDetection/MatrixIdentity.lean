import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Data.Real.Basic
import Mathlib.Data.Matrix.ColumnRowPartitioned
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

namespace StepDetection
open Matrix

variable {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]

/-- The complementary block of the inverse recovers a block determinant. -/
theorem complementary_block_det
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (C : Matrix n m ℝ)
    (D : Matrix n n ℝ) [Invertible D] [Invertible (fromBlocks A B C D)] :
    D.det = (fromBlocks A B C D).det * ((fromBlocks A B C D)⁻¹).toBlocks₁₁.det := by
  let _ := Matrix.invertibleOfFromBlocks₂₂Invertible A B C D
  have hinv := Matrix.invOf_fromBlocks₂₂_eq A B C D
  have hblock := congrArg Matrix.toBlocks₁₁ hinv
  simp only [Matrix.toBlocks_fromBlocks₁₁, Matrix.invOf_eq_nonsing_inv] at hblock
  rw [hblock, Matrix.det_fromBlocks₂₂]
  simp only [Matrix.invOf_eq_nonsing_inv, Matrix.det_nonsing_inv]
  have hunit := Matrix.isUnit_det_of_invertible (A - B * ⅟D * C)
  simp only [Matrix.invOf_eq_nonsing_inv, isUnit_iff_ne_zero] at hunit
  rw [Ring.inverse_eq_inv]
  field_simp

/-- A basis with Gram matrix `diag(G, I)` accounts for the non-normalized columns. -/
theorem compression_det_of_gram
    (M P : Matrix (m ⊕ n) (m ⊕ n) ℝ) (G : Matrix m m ℝ)
    [Invertible M] [Invertible P]
    [Invertible (Pᵀ * M * P).toBlocks₂₂]
    (hgram : Pᵀ * P = fromBlocks G 0 0 (1 : Matrix n n ℝ)) :
    (Pᵀ * M * P).toBlocks₂₂.det =
      M.det * (Pᵀ * M⁻¹ * P).toBlocks₁₁.det / G.det := by
  let H := Pᵀ * M * P
  let _ : Invertible (Pᵀ * M) := invertibleMul _ _
  let _ : Invertible H := invertibleMul _ _
  have hinv : Pᵀ * M⁻¹ * P = (Pᵀ * P) * H⁻¹ * (Pᵀ * P) := by
    dsimp [H]
    simp only [Matrix.mul_inv_rev, Matrix.mul_assoc,
      Matrix.mul_inv_cancel_left_of_invertible, Matrix.inv_mul_cancel_left_of_invertible]
  have hb : (Pᵀ * M⁻¹ * P).toBlocks₁₁ = G * H⁻¹.toBlocks₁₁ * G := by
    rw [hinv, hgram, ← Matrix.fromBlocks_toBlocks H⁻¹]
    simp only [Matrix.fromBlocks_multiply, Matrix.mul_zero, Matrix.zero_mul,
      add_zero, zero_add, Matrix.one_mul, Matrix.mul_one, Matrix.toBlocks_fromBlocks₁₁]
  have hG : P.det * P.det = G.det := by
    have := congrArg Matrix.det hgram
    simpa only [Matrix.det_mul, Matrix.det_transpose, Matrix.det_fromBlocks_zero₁₂,
      Matrix.det_one, mul_one] using this
  have hG0 : G.det ≠ 0 := by
    rw [← hG]
    exact mul_ne_zero (Matrix.isUnit_det_of_invertible P).ne_zero
      (Matrix.isUnit_det_of_invertible P).ne_zero
  have hH : H.det = G.det * M.det := by
    dsimp [H]
    rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose, ← hG]
    ring
  have hc : H.toBlocks₂₂.det = H.det * H⁻¹.toBlocks₁₁.det := by
    have hh := Matrix.fromBlocks_toBlocks H
    let _ : Invertible (fromBlocks H.toBlocks₁₁ H.toBlocks₁₂ H.toBlocks₂₁ H.toBlocks₂₂) :=
      Invertible.copy (inferInstance : Invertible H) _ hh
    simpa only [hh] using complementary_block_det H.toBlocks₁₁ H.toBlocks₁₂
      H.toBlocks₂₁ H.toBlocks₂₂
  change H.toBlocks₂₂.det = _
  rw [hc, hb, hH, Matrix.det_mul, Matrix.det_mul]
  field_simp [hG0]

/-- The projection determinant identity, retaining the original plateau normalization. -/
theorem residual_det_identity
    (M : Matrix (m ⊕ n) (m ⊕ n) ℝ)
    (X : Matrix (m ⊕ n) m ℝ) (U : Matrix (m ⊕ n) n ℝ)
    [Invertible M] [Invertible (fromCols X U)] [Invertible (Uᵀ * M * U)]
    (hUU : Uᵀ * U = 1) (hXU : Xᵀ * U = 0) :
    (Uᵀ * M * U).det = M.det * (Xᵀ * M⁻¹ * X).det / (Xᵀ * X).det := by
  have hUX : Uᵀ * X = 0 := by
    have := congrArg Matrix.transpose hXU
    simpa only [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.transpose_zero] using this
  have hblocks (N : Matrix (m ⊕ n) (m ⊕ n) ℝ) :
      (fromCols X U)ᵀ * N * fromCols X U =
        fromBlocks (Xᵀ * N * X) (Xᵀ * N * U) (Uᵀ * N * X) (Uᵀ * N * U) := by
    rw [Matrix.transpose_fromCols, Matrix.fromRows_mul, Matrix.fromRows_mul_fromCols]
  let _ : Invertible ((fromCols X U)ᵀ * M * fromCols X U).toBlocks₂₂ :=
    Invertible.copy (inferInstance : Invertible (Uᵀ * M * U)) _ (by rw [hblocks]; rfl)
  have hgram : (fromCols X U)ᵀ * fromCols X U = fromBlocks (Xᵀ * X) 0 0 1 := by
    rw [Matrix.transpose_fromCols, Matrix.fromRows_mul_fromCols, hUU, hXU, hUX]
  have h := compression_det_of_gram M (fromCols X U) (Xᵀ * X) hgram
  simpa only [hblocks, Matrix.toBlocks_fromBlocks₁₁, Matrix.toBlocks_fromBlocks₂₂] using h

end StepDetection
