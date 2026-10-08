import StepDetection.MatrixIdentity

namespace StepDetection
open Matrix

variable {r m n : Type*} [Fintype r] [Fintype m] [Fintype n]
  [DecidableEq r] [DecidableEq m] [DecidableEq n]

theorem residual_det_reindex (e : m ⊕ n ≃ r) (M : Matrix r r ℝ)
    (X : Matrix r m ℝ) (U : Matrix r n ℝ)
    [Invertible M] [Invertible (Uᵀ * M * U)]
    (hUU : Uᵀ * U = 1) (hXU : Xᵀ * U = 0) (hX : (Xᵀ * X).det ≠ 0) :
    (Uᵀ * M * U).det = M.det * (Xᵀ * M⁻¹ * X).det / (Xᵀ * X).det := by
  let V := U.submatrix e id
  let W := X.submatrix e id
  let N := M.submatrix e e
  let P := fromCols W V
  have hprod (A : Matrix r r ℝ) (Y : Matrix r m ℝ) :
      (Y.submatrix e id)ᵀ * A.submatrix e e * Y.submatrix e id = Yᵀ * A * Y := by
    simp only [Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv, Matrix.submatrix_id_id]
  have hVprod (A : Matrix r r ℝ) : Vᵀ * A.submatrix e e * V = Uᵀ * A * U := by
    dsimp [V]
    simp only [Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv, Matrix.submatrix_id_id]
  have hVV : Vᵀ * V = 1 := by
    simpa only [V, Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv,
      Matrix.submatrix_id_id] using hUU
  have hWV : Wᵀ * V = 0 := by
    simpa only [W, V, Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv,
      Matrix.submatrix_id_id] using hXU
  have hVW : Vᵀ * W = 0 := by
    simpa only [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.transpose_zero] using
      congrArg Matrix.transpose hWV
  have hWW : Wᵀ * W = Xᵀ * X := by
    simp only [W, Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv, Matrix.submatrix_id_id]
  have hgram : Pᵀ * P = fromBlocks (Xᵀ * X) 0 0 (1 : Matrix n n ℝ) := by
    dsimp [P]
    rw [Matrix.transpose_fromCols, Matrix.fromRows_mul_fromCols, hVV, hWV, hVW, hWW]
  have hp : P.det ≠ 0 := by
    have h := congrArg Matrix.det hgram
    simp only [Matrix.det_mul, Matrix.det_transpose, Matrix.det_fromBlocks_zero₁₂,
      Matrix.det_one, mul_one] at h
    intro hz
    rw [hz, zero_mul] at h
    exact hX h.symm
  let _ := Matrix.invertibleOfIsUnitDet P (isUnit_iff_ne_zero.mpr hp)
  let _ : Invertible N := Matrix.submatrixEquivInvertible M e e
  let _ : Invertible (Vᵀ * N * V) :=
    Invertible.copy (inferInstance : Invertible (Uᵀ * M * U)) _ (hVprod M)
  have h := residual_det_identity N W V hVV hWV
  simpa only [N, W, hVprod, hWW, Matrix.det_submatrix_equiv_self,
    Matrix.inv_submatrix_equiv, hprod] using h

end StepDetection
