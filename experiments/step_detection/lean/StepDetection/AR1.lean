import StepDetection.SavedPrecision
import StepDetection.TridiagonalEntries
import Mathlib.Data.Nat.Dist

namespace StepDetection.DeterminantTrace
open Matrix

theorem ar_left (rho : ℝ) (j : ℕ) :
    rho ^ j - rho * rho ^ Nat.dist 1 j = if j = 0 then 1 - rho ^ 2 else 0 := by
  by_cases hj : j = 0
  · subst j
    norm_num [Nat.dist]
    ring
  · have hdist : Nat.dist 1 j = j - 1 := Nat.dist_eq_sub_of_le (by omega)
    have hexp : j = (j - 1) + 1 := by omega
    rw [if_neg hj, hdist]
    conv_lhs => lhs; rw [hexp, pow_succ]
    ring

theorem ar_interior (rho : ℝ) (i j : ℕ) (hi : 0 < i) :
    (1 + rho ^ 2) * rho ^ Nat.dist i j - rho * rho ^ Nat.dist (i - 1) j -
      rho * rho ^ Nat.dist (i + 1) j = if i = j then 1 - rho ^ 2 else 0 := by
  rcases lt_trichotomy i j with hij | hij | hij
  · have h0 : Nat.dist i j = (j - i - 1) + 1 := by simp only [Nat.dist]; omega
    have h1 : Nat.dist (i - 1) j = (j - i - 1) + 2 := by simp only [Nat.dist]; omega
    have h2 : Nat.dist (i + 1) j = j - i - 1 := by simp only [Nat.dist]; omega
    rw [if_neg (ne_of_lt hij), h0, h1, h2]
    simp only [pow_succ]
    ring
  · subst j
    have h1 : Nat.dist (i - 1) i = 1 := by simp only [Nat.dist]; omega
    have h2 : Nat.dist (i + 1) i = 1 := by simp only [Nat.dist]; omega
    simp only [if_true, Nat.dist_self, pow_zero, h1, h2, pow_one]
    ring
  · have h0 : Nat.dist i j = (i - j - 1) + 1 := by simp only [Nat.dist]; omega
    have h1 : Nat.dist (i - 1) j = i - j - 1 := by simp only [Nat.dist]; omega
    have h2 : Nat.dist (i + 1) j = (i - j - 1) + 2 := by simp only [Nat.dist]; omega
    rw [if_neg (ne_of_gt hij), h0, h1, h2]
    simp only [pow_succ]
    ring

theorem ar_right (rho : ℝ) (i j : ℕ) (hi : 0 < i) (hj : j ≤ i) :
    rho ^ Nat.dist i j - rho * rho ^ Nat.dist (i - 1) j =
      if i = j then 1 - rho ^ 2 else 0 := by
  by_cases hij : i = j
  · subst j
    have h1 : Nat.dist (i - 1) i = 1 := by simp only [Nat.dist]; omega
    simp only [if_true, Nat.dist_self, h1, pow_zero, pow_one]
    ring
  · have h0 : Nat.dist i j = (i - j - 1) + 1 := by simp only [Nat.dist]; omega
    have h1 : Nat.dist (i - 1) j = i - j - 1 := by simp only [Nat.dist]; omega
    rw [if_neg hij, h0, h1, pow_succ]
    ring

noncomputable def arMatrix (rho : ℝ) : Matrix (Fin 100) (Fin 100) ℝ :=
  fun i j => rho ^ Nat.dist i.val j.val

theorem precision_covariance_product (rho : ℝ) :
    arPrecision rho * arMatrix rho = (1 - rho ^ 2) • (1 : Matrix (Fin 100) (Fin 100) ℝ) := by
  ext i j
  rw [arPrecision, Tridiagonal.matrix_mul_apply]
  simp only [arMatrix, Matrix.smul_apply, smul_eq_mul, Matrix.one_apply, Fin.ext_iff]
  by_cases hi : i.val = 0
  · have hn : i.val + 1 < 100 := by omega
    have h := ar_left rho j.val
    simpa [arDiagonal, hi, hn, Tridiagonal.next, Nat.dist, sub_eq_add_neg,
      eq_comm, mul_ite] using h
  · by_cases hlast : i.val = 99
    · have hn : ¬i.val + 1 < 100 := by omega
      have h := ar_right rho i.val j.val (by omega) (by have := j.isLt; omega)
      simpa [arDiagonal, hi, hlast, hn, Tridiagonal.previous, sub_eq_add_neg,
        mul_ite] using h
    · have hn : i.val + 1 < 100 := by have := i.isLt; omega
      have h := ar_interior rho i.val j.val (by omega)
      simpa [arDiagonal, hi, hlast, hn, Tridiagonal.previous, Tridiagonal.next,
        sub_eq_add_neg, mul_ite, add_assoc] using h

theorem covariance_entries {rho : ℝ} (hrho : inputBounds.Contains rho) :
    ar_covariance rho = arMatrix rho := by
  let _ := Matrix.invertibleOfIsUnitDet (arPrecision rho)
    (isUnit_iff_ne_zero.mpr (by rw [ar_determinant]; exact ne_of_gt (stationary_positive hrho)))
  have h := congrArg (fun B : Matrix (Fin 100) (Fin 100) ℝ => (arPrecision rho)⁻¹ * B)
    (precision_covariance_product rho)
  simpa only [Matrix.inv_mul_cancel_left_of_invertible, Matrix.mul_smul,
    Matrix.mul_one, ar_covariance] using h.symm

end StepDetection.DeterminantTrace
