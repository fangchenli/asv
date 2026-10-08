import StepDetection.SavedPrecision
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Algebra.Order.Star.Real

namespace StepDetection.DeterminantTrace
open Matrix

theorem ar_posDef {rho : ℝ} (hrho : inputBounds.Contains rho) : (arPrecision rho).PosDef := by
  have hd : (Matrix.diagonal (arPivots rho)).PosDef := Matrix.PosDef.diagonal (by
    intro i
    dsimp [arPivots]
    split_ifs
    · exact stationary_positive hrho
    · norm_num)
  let L := Tridiagonal.lower (fun _ : Fin 100 => -rho)
  have hL : IsUnit L := by
    apply (Matrix.isUnit_iff_isUnit_det L).mpr
    rw [Tridiagonal.det_lower]
    exact isUnit_one
  have h := hd.mul_mul_conjTranspose_same (B := L) (Matrix.vecMul_injective_of_isUnit hL)
  rw [ar_factorization]
  simpa only [Matrix.conjTranspose_eq_transpose_of_trivial] using h

theorem shifted_posDef {rho : ℝ} (hrho : inputBounds.Contains rho) :
    (shiftedCovariance rho).PosDef := by
  rw [shifted_covariance_eq]
  have ha : 0 < v2 rho := by norm_num [v2]
  have hc : 0 < v7 rho := by
    have hb : 0 < v3 rho := by norm_num [v3]
    simpa only [v7, v6, v1, v5, v4, div_one] using mul_pos hb (stationary_positive hrho)
  exact (Matrix.PosDef.one.smul ha).add ((ar_posDef hrho).inv.smul hc)

theorem residual_posDef {rho : ℝ} (hrho : inputBounds.Contains rho)
    (U : Matrix (Fin 100) (Fin 98) ℝ) (hU : Uᵀ * U = 1) :
    (Uᵀ * shiftedCovariance rho * U).PosDef := by
  have hinj : Function.Injective U.mulVec := by
    intro x y h
    have h' := congrArg (fun v => Uᵀ *ᵥ v) h
    simpa only [Matrix.mulVec_mulVec, hU, Matrix.one_mulVec] using h'
  simpa only [Matrix.conjTranspose_eq_transpose_of_trivial] using
    (shifted_posDef hrho).conjTranspose_mul_mul_same hinj

end StepDetection.DeterminantTrace
