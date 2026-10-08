import Mathlib.RingTheory.Polynomial.Bernstein
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Positivity

namespace StepDetection

open Polynomial

/-- A polynomial represented in the Bernstein basis on the unit interval. -/
noncomputable def bernsteinExpansion (degree : ℕ) (coeff : ℕ → ℝ) : ℝ[X] :=
  ∑ i ∈ Finset.range (degree + 1), Polynomial.C (coeff i) *
    bernsteinPolynomial ℝ degree i

private theorem bernsteinBasis_nonneg {degree i : ℕ} {t : ℝ}
    (ht : 0 ≤ t) (htone : t ≤ 1) :
    0 ≤ (bernsteinPolynomial ℝ degree i).eval t := by
  rw [bernsteinPolynomial]
  simp only [eval_mul, eval_natCast, eval_X, eval_pow, eval_sub, eval_one]
  positivity

/-- Bernstein coefficients bound every value on `[0,1]` because the basis
    functions are nonnegative and sum to one. -/
theorem bernsteinExpansion_bounds {degree : ℕ} {coeff : ℕ → ℝ}
    {lower upper t : ℝ} (ht : 0 ≤ t) (htone : t ≤ 1)
    (hcoeff : ∀ i ∈ Finset.range (degree + 1), lower ≤ coeff i ∧ coeff i ≤ upper) :
    lower ≤ (bernsteinExpansion degree coeff).eval t ∧
      (bernsteinExpansion degree coeff).eval t ≤ upper := by
  have hweights : ∀ i ∈ Finset.range (degree + 1),
      0 ≤ (bernsteinPolynomial ℝ degree i).eval t := by
    intro i hi
    exact bernsteinBasis_nonneg ht htone
  have hsum : ∑ i ∈ Finset.range (degree + 1),
      (bernsteinPolynomial ℝ degree i).eval t = 1 := by
    calc
      ∑ i ∈ Finset.range (degree + 1),
          (bernsteinPolynomial ℝ degree i).eval t =
        (∑ i ∈ Finset.range (degree + 1), bernsteinPolynomial ℝ degree i).eval t := by
          rw [eval_finsetSum]
      _ = (1 : ℝ[X]).eval t := congrArg (fun p : ℝ[X] => p.eval t)
        (bernsteinPolynomial.sum (R := ℝ) degree)
      _ = 1 := by simp
  constructor
  · simp only [bernsteinExpansion, eval_finsetSum, eval_mul, eval_C]
    change lower ≤ ∑ i ∈ Finset.range (degree + 1),
      coeff i * (bernsteinPolynomial ℝ degree i).eval t
    calc
      lower = lower * (∑ i ∈ Finset.range (degree + 1),
          (bernsteinPolynomial ℝ degree i).eval t) := by rw [hsum, mul_one]
      _ = ∑ i ∈ Finset.range (degree + 1),
          lower * (bernsteinPolynomial ℝ degree i).eval t := by rw [Finset.mul_sum]
      _ ≤ ∑ i ∈ Finset.range (degree + 1),
          coeff i * (bernsteinPolynomial ℝ degree i).eval t := by
            apply Finset.sum_le_sum
            intro i hi
            exact mul_le_mul_of_nonneg_right (hcoeff i hi).1 (hweights i hi)
  · simp only [bernsteinExpansion, eval_finsetSum, eval_mul, eval_C]
    change ∑ i ∈ Finset.range (degree + 1),
      coeff i * (bernsteinPolynomial ℝ degree i).eval t ≤ upper
    calc
      ∑ i ∈ Finset.range (degree + 1),
          coeff i * (bernsteinPolynomial ℝ degree i).eval t ≤
        ∑ i ∈ Finset.range (degree + 1),
          upper * (bernsteinPolynomial ℝ degree i).eval t := by
            apply Finset.sum_le_sum
            intro i hi
            exact mul_le_mul_of_nonneg_right (hcoeff i hi).2 (hweights i hi)
      _ = upper * (∑ i ∈ Finset.range (degree + 1),
          (bernsteinPolynomial ℝ degree i).eval t) := by rw [← Finset.mul_sum]
      _ = upper := by rw [hsum, mul_one]

end StepDetection
