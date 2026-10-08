import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.Data.Matrix.ColumnRowPartitioned
import Mathlib.Data.Real.Basic

set_option linter.style.haveILetI false

namespace StepDetection
open Matrix

variable {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]

/-- The inverse of the lower-right covariance block is the Schur complement
    of the upper-left precision block. -/
theorem inverse_bottomRight_eq_schur
    (C : Matrix (m ⊕ n) (m ⊕ n) ℝ)
    [Invertible C] [Invertible C.toBlocks₂₂]
    [Invertible (C⁻¹).toBlocks₁₁] :
    (C.toBlocks₂₂)⁻¹ = (C⁻¹).toBlocks₂₂ -
      (C⁻¹).toBlocks₂₁ * ((C⁻¹).toBlocks₁₁)⁻¹ * (C⁻¹).toBlocks₁₂ := by
  let A := (C⁻¹).toBlocks₁₁
  let B := (C⁻¹).toBlocks₁₂
  let D := (C⁻¹).toBlocks₂₁
  let E := (C⁻¹).toBlocks₂₂
  let H : Matrix (m ⊕ n) (m ⊕ n) ℝ := fromBlocks A B D E
  have hH : H = C⁻¹ := by
    dsimp [H, A, B, D, E]
    exact Matrix.fromBlocks_toBlocks (C⁻¹)
  letI : Invertible H := Invertible.copy (inferInstance : Invertible (C⁻¹)) H hH
  letI : Invertible (E - D * ⅟A * B) :=
    Matrix.invertibleOfFromBlocks₁₁Invertible A B D E
  let S := E - D * A⁻¹ * B
  letI : Invertible S := Invertible.copy
    (inferInstance : Invertible (E - D * ⅟A * B)) S (by simp [S])
  have hblock := Matrix.invOf_fromBlocks₁₁_eq A B D E
  have hHC : H * C = 1 := by rw [hH]; simp
  have hinv : H⁻¹ = C := Matrix.inv_eq_right_inv hHC
  have hbottom : C.toBlocks₂₂ = S⁻¹ := by
    have hparts := congrArg Matrix.toBlocks₂₂ hblock
    change (⅟H).toBlocks₂₂ = _ at hparts
    rw [Matrix.invOf_eq_nonsing_inv, hinv] at hparts
    simpa only [Matrix.toBlocks_fromBlocks₂₂, Matrix.invOf_eq_nonsing_inv] using hparts
  have hright : C.toBlocks₂₂ * (E - D * A⁻¹ * B) = 1 := by
    rw [hbottom]
    exact Matrix.inv_mul_of_invertible S
  have hinv := Matrix.inv_eq_right_inv hright
  simpa only [A, B, D, E] using hinv

end StepDetection
