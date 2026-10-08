import StepDetection.SavedMatrix.Gram
import StepDetection.AR1
import StepDetection.SavedPositive
import StepDetection.MatrixReindex
import StepDetection.ProfiledPrecision

set_option maxRecDepth 2048

namespace StepDetection.DeterminantTrace
open Matrix

/-- The saved two-plateau design and any orthonormal residual basis form a
    complete coordinate system; the plateau Gram determinant is 99. -/
theorem saved_residual_coordinates_det_ne_zero
    (U : Matrix (Fin 100) (Fin 98) ℝ)
    (hU : Uᵀ * U = 1) (hXU : plateauColumnsᵀ * U = 0) :
    (fromCols (plateauColumns.submatrix
      (finSumFinEquiv : Fin 2 ⊕ Fin 98 ≃ Fin 100) id)
      (U.submatrix (finSumFinEquiv : Fin 2 ⊕ Fin 98 ≃ Fin 100) id)).det ≠ 0 := by
  let e : Fin 2 ⊕ Fin 98 ≃ Fin 100 := finSumFinEquiv
  let X := plateauColumns.submatrix e id
  let V := U.submatrix e id
  have hVV : Vᵀ * V = 1 := by
    simpa only [V, Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv,
      Matrix.submatrix_id_id] using hU
  have hXV : Xᵀ * V = 0 := by
    simpa only [X, V, Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv,
      Matrix.submatrix_id_id] using hXU
  have hXX : (Xᵀ * X).det ≠ 0 := by
    have hgram : Xᵀ * X = plateauColumnsᵀ * plateauColumns := by
      simp only [X, Matrix.transpose_submatrix, Matrix.submatrix_mul_equiv,
        Matrix.submatrix_id_id]
    rw [hgram, plateau_gram_det]
    norm_num
  exact residual_coordinates_det_ne_zero X V hVV hXV hXX

theorem recorded_gram_is_inverse {rho : ℝ} (hrho : inputBounds.Contains rho) :
    recordedGram rho = plateauColumnsᵀ * (shiftedCovariance rho)⁻¹ * plateauColumns := by
  let _ := (ar_posDef hrho).isUnit.invertible
  let _ := Matrix.invertibleOfIsUnitDet (precisionShift rho)
    (isUnit_iff_ne_zero.mpr (ne_of_gt (saved_det_positive hrho)))
  have ha : v2 rho ≠ 0 := by norm_num [v2]
  have hs : (shiftedCovariance rho)⁻¹ * plateauColumns =
      (v2 rho)⁻¹ • (plateauColumns - v7 rho • solutionMatrix rho) := by
    rw [shifted_covariance_eq]
    exact precision_shift_solve _ _ _ _ _ _ ha (shift_matrix rho) (saved_solve hrho)
  rw [recorded_gram, Matrix.mul_assoc, hs, Matrix.mul_smul, Matrix.mul_sub, Matrix.mul_smul]

theorem shifted_determinant {rho : ℝ} (hrho : inputBounds.Contains rho) :
    (shiftedCovariance rho).det = (precisionShift rho).det / (1 - rho ^ 2) := by
  let _ := (ar_posDef hrho).isUnit.invertible
  rw [shifted_covariance_eq, precision_shift_det _ _ _ _ (shift_matrix rho), ar_determinant]

theorem recorded_expression_identity {rho : ℝ} (hrho : inputBounds.Contains rho) :
    determinantExpression rho = (shiftedCovariance rho).det *
      (plateauColumnsᵀ * (shiftedCovariance rho)⁻¹ * plateauColumns).det / 99 := by
  rw [recorded_final_formula, ← saved_determinant hrho, recorded_gram_is_inverse hrho,
    shifted_determinant hrho]
  simp only [v6, v1, v5, v4, div_one]
  field_simp [ne_of_gt (stationary_positive hrho)]

/-- The complete recorded expression equals the determinant on the residual subspace. -/
theorem residual_determinant_eq_trace {rho : ℝ} (hrho : inputBounds.Contains rho)
    (U : Matrix (Fin 100) (Fin 98) ℝ)
    (hU : Uᵀ * U = 1) (hXU : plateauColumnsᵀ * U = 0) :
    (v2 rho • (1 : Matrix (Fin 98) (Fin 98) ℝ) +
      v3 rho • (Uᵀ * arMatrix rho * U)).det = determinantExpression rho := by
  let _ := (shifted_posDef hrho).isUnit.invertible
  let _ := (residual_posDef hrho U hU).isUnit.invertible
  have h := residual_det_reindex (finSumFinEquiv : Fin 2 ⊕ Fin 98 ≃ Fin 100)
    (shiftedCovariance rho) plateauColumns U hU hXU (by rw [plateau_gram_det]; norm_num)
  rw [plateau_gram_det, ← recorded_expression_identity hrho] at h
  have hc : Uᵀ * shiftedCovariance rho * U =
      v2 rho • (1 : Matrix (Fin 98) (Fin 98) ℝ) + v3 rho • (Uᵀ * arMatrix rho * U) := by
    rw [shiftedCovariance, covariance_entries hrho, Matrix.mul_add, Matrix.mul_smul,
      Matrix.mul_smul, Matrix.mul_one, Matrix.add_mul, Matrix.smul_mul, Matrix.smul_mul, hU]
  rwa [hc] at h

theorem saved_residual_enclosure {rho : ℝ} (hrho : inputBounds.Contains rho)
    (U : Matrix (Fin 100) (Fin 98) ℝ)
    (hU : Uᵀ * U = 1) (hXU : plateauColumnsᵀ * U = 0) :
    b1722.Contains ((v2 rho • (1 : Matrix (Fin 98) (Fin 98) ℝ) +
      v3 rho • (Uᵀ * arMatrix rho * U)).det) := by
  rw [residual_determinant_eq_trace hrho U hU hXU]
  exact h1722 hrho

theorem cutoff_with_residual_matrix {rho p correction : ℝ}
    (hrho : inputBounds.Contains rho) (U : Matrix (Fin 100) (Fin 98) ℝ)
    (hU : Uᵀ * U = 1) (hXU : plateauColumnsᵀ * U = 0)
    (hcorr : (Saved.correctionLower : ℝ) ≤ correction)
    (hprob : p ^ 2 ≤ 1 / ((v2 rho • (1 : Matrix (Fin 98) (Fin 98) ℝ) +
      v3 rho • (Uᵀ * arMatrix rho * U)).det * correction ^ 2)) : p < (1 : ℝ) / 100 := by
  rw [residual_determinant_eq_trace hrho U hU hXU] at hprob
  exact cutoff_with_trace hrho hcorr hprob

end StepDetection.DeterminantTrace
