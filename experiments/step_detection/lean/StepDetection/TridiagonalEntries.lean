import StepDetection.Tridiagonal

namespace StepDetection.Tridiagonal
open Matrix

variable {n : ℕ}

theorem matrix_as_sum (d off : Fin n → ℝ) :
    matrix d off = Matrix.diagonal d + shift off + (shift off)ᵀ := by
  ext i j
  simp only [matrix, shift, Matrix.add_apply, Matrix.diagonal_apply, Matrix.transpose_apply]
  have hval : i = j ↔ i.val = j.val := Fin.ext_iff
  split_ifs <;> first | omega | ring

theorem matrix_mul_apply {k : Type*} (d off : Fin n → ℝ) (A : Matrix (Fin n) k ℝ)
    (i : Fin n) (j : k) :
    (matrix d off * A) i j = d i * A i j +
      (if i.val = 0 then 0 else off i * A (previous i) j) +
      (if h : i.val + 1 < n then off (next i h) * A (next i h) j else 0) := by
  have hl : (shift off * A) i j = if i.val = 0 then 0 else off i * A (previous i) j := by
    have h := lower_mul_apply off A i j
    simp only [lower, Matrix.add_mul, Matrix.one_mul, Matrix.add_apply] at h
    exact add_left_cancel h
  have hu : ((shift off)ᵀ * A) i j =
      if h : i.val + 1 < n then off (next i h) * A (next i h) j else 0 := by
    have h := upper_mul_apply off A i j
    simp only [lower, Matrix.transpose_add, Matrix.transpose_one, Matrix.add_mul,
      Matrix.one_mul, Matrix.add_apply] at h
    exact add_left_cancel h
  rw [matrix_as_sum, Matrix.add_mul, Matrix.add_mul]
  simp only [Matrix.add_apply, Matrix.diagonal_mul, hl, hu]

end StepDetection.Tridiagonal
