import StepDetection.SavedMatrix.Data
import StepDetection.PrecisionIdentity

namespace StepDetection.DeterminantTrace
open Matrix

noncomputable def arDiagonal (rho : ℝ) (i : Fin 100) : ℝ :=
  if i.val = 0 ∨ i.val = 99 then 1 else 1 + rho ^ 2

noncomputable def arPrecision (rho : ℝ) : Matrix (Fin 100) (Fin 100) ℝ :=
  Tridiagonal.matrix (arDiagonal rho) (fun _ => -rho)

noncomputable def arPivots (rho : ℝ) (i : Fin 100) : ℝ :=
  if i = 99 then 1 - rho ^ 2 else 1

theorem ar_pivot_equations (rho : ℝ) :
    (∀ i, arDiagonal rho i = arPivots rho i + if i.val = 0 then 0 else
      (-rho) * arPivots rho (Tridiagonal.previous i) * (-rho)) ∧
    (∀ i : Fin 100, i.val ≠ 0 → -rho = (-rho) * arPivots rho (Tridiagonal.previous i)) := by
  have hdiag (i : Fin 100) : arDiagonal rho i = arPivots rho i +
      if i.val = 0 then 0 else (-rho) * arPivots rho (Tridiagonal.previous i) * (-rho) := by
    have hprev : Tridiagonal.previous i ≠ 99 := by
      intro h
      have hval := congrArg Fin.val h
      have := i.isLt
      dsimp [Tridiagonal.previous] at hval
      omega
    have hlast : i = 99 ↔ i.val = 99 := Fin.ext_iff
    simp only [arDiagonal, arPivots, hprev, if_false, hlast]
    split_ifs <;> first | omega | ring
  have hoff (i : Fin 100) (_hi : i.val ≠ 0) : -rho = (-rho) * arPivots rho (Tridiagonal.previous i) := by
    have hprev : Tridiagonal.previous i ≠ 99 := by
      intro h
      have hval := congrArg Fin.val h
      have := i.isLt
      dsimp [Tridiagonal.previous] at hval
      omega
    simp [arPivots, hprev]
  exact ⟨hdiag, hoff⟩

theorem ar_factorization (rho : ℝ) :
    arPrecision rho = Tridiagonal.lower (fun _ => -rho) * Matrix.diagonal (arPivots rho) *
      (Tridiagonal.lower (fun _ : Fin 100 => -rho))ᵀ :=
  Tridiagonal.factorization _ _ _ _ (ar_pivot_equations rho).1 (ar_pivot_equations rho).2

theorem ar_determinant (rho : ℝ) : (arPrecision rho).det = 1 - rho ^ 2 := by
  rw [arPrecision, Tridiagonal.determinant _ _ (arPivots rho) (fun _ => -rho)
    (ar_pivot_equations rho).1 (ar_pivot_equations rho).2]
  simp [arPivots]

theorem shift_matrix (rho : ℝ) :
    precisionShift rho = v2 rho • arPrecision rho + v7 rho • (1 : Matrix (Fin 100) (Fin 100) ℝ) := by
  ext i j
  simp only [precisionShift, Tridiagonal.matrix, Matrix.add_apply,
    Matrix.smul_apply, smul_eq_mul, arPrecision, Matrix.one_apply, arDiagonal,
    diagonal, offDiagonal]
  split_ifs <;>
    simp_all [v8, v11, v10, v9, v13, v12, v0, v1, v4, v5, v6, v7]

theorem stationary_positive {rho : ℝ} (hrho : inputBounds.Contains rho) :
    0 < 1 - rho ^ 2 := by
  rcases hrho with ⟨hlo, hhi⟩
  dsimp [inputBounds] at hlo hhi
  have hl : 0 < 1 + rho := by linarith
  have hh : 0 < 1 - rho := by linarith
  nlinarith [mul_pos hl hh]

noncomputable def ar_covariance (rho : ℝ) : Matrix (Fin 100) (Fin 100) ℝ :=
  (1 - rho ^ 2) • (arPrecision rho)⁻¹

noncomputable def shiftedCovariance (rho : ℝ) : Matrix (Fin 100) (Fin 100) ℝ :=
  v2 rho • (1 : Matrix (Fin 100) (Fin 100) ℝ) + v3 rho • ar_covariance rho

theorem shifted_covariance_eq (rho : ℝ) :
    shiftedCovariance rho = v2 rho • (1 : Matrix (Fin 100) (Fin 100) ℝ) +
      v7 rho • (arPrecision rho)⁻¹ := by
  simp [shiftedCovariance, ar_covariance, smul_smul, v7, v6, v5, v4, v1]

end StepDetection.DeterminantTrace
