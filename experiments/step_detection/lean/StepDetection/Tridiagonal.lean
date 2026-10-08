import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.FieldSimp

namespace StepDetection.Tridiagonal
open Matrix

variable {n : ℕ}

def previous (i : Fin n) : Fin n := ⟨i.val - 1, lt_of_le_of_lt (Nat.sub_le _ _) i.isLt⟩

def shift (m : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => if j.val + 1 = i.val then m i else 0

def lower (m : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ := 1 + shift m

theorem lower_mul_apply {k : Type*} (m : Fin n → ℝ) (A : Matrix (Fin n) k ℝ)
    (i : Fin n) (j : k) :
    (lower m * A) i j = A i j + if i.val = 0 then 0 else m i * A (previous i) j := by
  rw [lower, Matrix.add_mul, Matrix.one_mul, Matrix.add_apply]
  congr 1
  rw [Matrix.mul_apply]
  by_cases hi : i.val = 0
  · simp [shift, hi]
  · have hindex (k : Fin n) : k.val + 1 = i.val ↔ k = previous i := by
      constructor
      · intro h
        apply Fin.ext
        dsimp [previous]
        omega
      · intro h
        subst k
        dsimp [previous]
        omega
    simp [shift, hindex, hi]

theorem mul_lower_transpose_apply {k : Type*} (m : Fin n → ℝ) (A : Matrix k (Fin n) ℝ)
    (i : k) (j : Fin n) :
    (A * (lower m)ᵀ) i j = A i j + if j.val = 0 then 0 else m j * A i (previous j) := by
  have h := lower_mul_apply m Aᵀ j i
  calc
    (A * (lower m)ᵀ) i j = (lower m * Aᵀ) j i := by
      rw [← Matrix.transpose_apply (A * (lower m)ᵀ), Matrix.transpose_mul,
        Matrix.transpose_transpose]
    _ = _ := by simpa only [Matrix.transpose_apply] using h

theorem det_lower (m : Fin n → ℝ) : (lower m).det = 1 := by
  have htri : (lower m).IsLowerTriangular := by
    intro i j hij
    have hne : i ≠ j := ne_of_lt hij
    have hnot : ¬j.val + 1 = i.val := by
      have := hij
      change i.val < j.val at this
      omega
    simp [lower, shift, hne, hnot]
  rw [Matrix.det_of_isLowerTriangular _ htri]
  simp [lower, shift]

def matrix (d off : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ := fun i j =>
  if i = j then d i else if j.val + 1 = i.val then off i
  else if i.val + 1 = j.val then off j else 0

/-- The local pivot equations imply the complete matrix factorization. -/
theorem factorization (d off p m : Fin n → ℝ)
    (hdiag : ∀ i, d i = p i + if i.val = 0 then 0 else m i * p (previous i) * m i)
    (hoff : ∀ i, i.val ≠ 0 → off i = m i * p (previous i)) :
    matrix d off = lower m * Matrix.diagonal p * (lower m)ᵀ := by
  ext i j
  rw [mul_lower_transpose_apply, lower_mul_apply, lower_mul_apply]
  simp only [Matrix.diagonal_apply, matrix]
  by_cases hij : i = j
  · subst j
    by_cases hi : i.val = 0
    · simp [hi, hdiag]
    · have hprev : previous i ≠ i := by
        intro h
        have := congrArg Fin.val h
        dsimp [previous] at this
        omega
      simp [hi, hprev, Ne.symm hprev, hdiag]
      ring
  · by_cases hi : i.val = 0 <;> by_cases hj : j.val = 0
    all_goals
      have hprev_i : previous i = j ↔ j.val + 1 = i.val := by
        dsimp [previous]
        constructor <;> intro h
        · have := congrArg Fin.val h
          simp only at this
          omega
        · apply Fin.ext
          simp only
          omega
      have hprev_j : i = previous j ↔ i.val + 1 = j.val := by
        dsimp [previous]
        constructor <;> intro h
        · have := congrArg Fin.val h
          simp only at this
          omega
        · apply Fin.ext
          simp only
          omega
      simp only [hij, if_false, hi, hj, hprev_i, hprev_j, ite_true, add_zero, zero_add]
      try
        have hprev_ne : previous i ≠ previous j := by
          intro h
          have := congrArg Fin.val h
          dsimp [previous] at this
          apply hij
          apply Fin.ext
          omega
        simp only [hprev_ne, if_false, mul_zero, add_zero]
      split_ifs <;> simp_all <;> try omega

theorem determinant (d off p m : Fin n → ℝ)
    (hdiag : ∀ i, d i = p i + if i.val = 0 then 0 else m i * p (previous i) * m i)
    (hoff : ∀ i, i.val ≠ 0 → off i = m i * p (previous i)) :
    (matrix d off).det = ∏ i, p i := by
  rw [factorization d off p m hdiag hoff, Matrix.det_mul, Matrix.det_mul,
    Matrix.det_transpose, det_lower, Matrix.det_diagonal, one_mul, mul_one]

def next (i : Fin n) (h : i.val + 1 < n) : Fin n := ⟨i.val + 1, h⟩

theorem upper_mul_apply {k : Type*} (m : Fin n → ℝ) (A : Matrix (Fin n) k ℝ)
    (i : Fin n) (j : k) :
    ((lower m)ᵀ * A) i j = A i j +
      if h : i.val + 1 < n then m (next i h) * A (next i h) j else 0 := by
  rw [lower, Matrix.transpose_add, Matrix.transpose_one, Matrix.add_mul,
    Matrix.one_mul, Matrix.add_apply]
  congr 1
  rw [Matrix.mul_apply]
  by_cases hi : i.val + 1 < n
  · have hindex (k : Fin n) : i.val + 1 = k.val ↔ k = next i hi := by
      constructor <;> intro h
      · apply Fin.ext
        exact h.symm
      · subst k
        rfl
    simp [shift, Matrix.transpose_apply, hindex, hi]
  · have hindex (k : Fin n) : ¬i.val + 1 = k.val := by
      have := k.isLt
      omega
    simp [shift, Matrix.transpose_apply, hindex, hi]

/-- The forward and backward recurrences solve the original tridiagonal system. -/
theorem solve {k : Type*} (d off p m : Fin n → ℝ) (X F Y : Matrix (Fin n) k ℝ)
    (hdiag : ∀ i, d i = p i + if i.val = 0 then 0 else m i * p (previous i) * m i)
    (hoff : ∀ i, i.val ≠ 0 → off i = m i * p (previous i))
    (hp : ∀ i, p i ≠ 0)
    (hf : ∀ i j, F i j = X i j - if i.val = 0 then 0 else m i * F (previous i) j)
    (hy : ∀ i j, Y i j = F i j / p i -
      if h : i.val + 1 < n then m (next i h) * Y (next i h) j else 0) :
    matrix d off * Y = X := by
  have hforward : lower m * F = X := by
    ext i j
    rw [lower_mul_apply, hf i j]
    exact sub_add_cancel _ _
  have hbackward : Matrix.diagonal p * ((lower m)ᵀ * Y) = F := by
    ext i j
    rw [Matrix.diagonal_mul, upper_mul_apply, hy i j]
    simp only [sub_add_cancel]
    exact mul_div_cancel₀ _ (hp i)
  rw [factorization d off p m hdiag hoff, Matrix.mul_assoc, Matrix.mul_assoc,
    hbackward, hforward]

theorem inverse_solve {k : Type*} (T : Matrix (Fin n) (Fin n) ℝ)
    (X Y : Matrix (Fin n) k ℝ) [Invertible T] (h : T * Y = X) : T⁻¹ * X = Y := by
  rw [← h, Matrix.inv_mul_cancel_left_of_invertible]

end StepDetection.Tridiagonal
