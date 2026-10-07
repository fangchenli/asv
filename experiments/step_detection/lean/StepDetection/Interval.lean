import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Algebra.Order.Floor.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity

/-! Real enclosures for the mathematical operations used by the Python reference.
The endpoints may be rational, but the enclosed value need not be. -/
namespace StepDetection

noncomputable section

structure Interval where
  lo : ℝ
  hi : ℝ

def Interval.Contains (a : Interval) (x : ℝ) : Prop := a.lo ≤ x ∧ x ≤ a.hi

def Interval.add (a b : Interval) : Interval := ⟨a.lo + b.lo, a.hi + b.hi⟩
def Interval.sub (a b : Interval) : Interval := ⟨a.lo - b.hi, a.hi - b.lo⟩
def Interval.mul (a b : Interval) : Interval :=
  ⟨min (min (a.lo * b.lo) (a.lo * b.hi)) (min (a.hi * b.lo) (a.hi * b.hi)),
   max (max (a.lo * b.lo) (a.lo * b.hi)) (max (a.hi * b.lo) (a.hi * b.hi))⟩
def Interval.inv (a : Interval) : Interval := ⟨1 / a.hi, 1 / a.lo⟩
def Interval.div (a b : Interval) : Interval := a.mul b.inv
def Interval.square (a : Interval) : Interval :=
  ⟨if a.lo ≤ 0 ∧ 0 ≤ a.hi then 0 else min (a.lo ^ 2) (a.hi ^ 2),
   max (a.lo ^ 2) (a.hi ^ 2)⟩

theorem Interval.add_sound {a b : Interval} {x y : ℝ}
    (hx : a.Contains x) (hy : b.Contains y) : (a.add b).Contains (x + y) :=
  ⟨add_le_add hx.1 hy.1, add_le_add hx.2 hy.2⟩

theorem Interval.sub_sound {a b : Interval} {x y : ℝ}
    (hx : a.Contains x) (hy : b.Contains y) : (a.sub b).Contains (x - y) :=
  ⟨sub_le_sub hx.1 hy.2, sub_le_sub hx.2 hy.1⟩

private theorem mul_bounds {a b x y : ℝ} (ha : a ≤ x) (hb : x ≤ b) :
    min (a * y) (b * y) ≤ x * y ∧ x * y ≤ max (a * y) (b * y) := by
  by_cases h : 0 ≤ y
  · exact ⟨(min_le_left _ _).trans (mul_le_mul_of_nonneg_right ha h),
      (mul_le_mul_of_nonneg_right hb h).trans (le_max_right _ _)⟩
  · have h' : y ≤ 0 := le_of_not_ge h
    exact ⟨(min_le_right _ _).trans (mul_le_mul_of_nonpos_right hb h'),
      (mul_le_mul_of_nonpos_right ha h').trans (le_max_left _ _)⟩

theorem Interval.mul_sound {a b : Interval} {x y : ℝ}
    (hx : a.Contains x) (hy : b.Contains y) : (a.mul b).Contains (x * y) := by
  have hlo := mul_bounds (y := a.lo) hy.1 hy.2
  have hhi := mul_bounds (y := a.hi) hy.1 hy.2
  have hxy := mul_bounds (y := y) hx.1 hx.2
  simp only [mul_comm y, mul_comm b.lo, mul_comm b.hi] at hlo hhi
  exact ⟨(min_le_min hlo.1 hhi.1).trans hxy.1,
    hxy.2.trans (max_le_max hlo.2 hhi.2)⟩

theorem Interval.inv_sound {a : Interval} {x : ℝ}
    (hx : a.Contains x) (hzero : 0 < a.lo ∨ a.hi < 0) :
    a.inv.Contains (1 / x) := by
  rcases hzero with hpos | hneg
  · exact ⟨one_div_le_one_div_of_le (hpos.trans_le hx.1) hx.2,
      one_div_le_one_div_of_le hpos hx.1⟩
  · have h1 := one_div_le_one_div_of_le (neg_pos.mpr hneg) (neg_le_neg hx.2)
    have h2 := one_div_le_one_div_of_le (neg_pos.mpr (hx.2.trans_lt hneg))
      (neg_le_neg hx.1)
    simp only [one_div_neg_eq_neg_one_div] at h1 h2
    constructor <;> dsimp [Interval.inv] <;> linarith

theorem Interval.div_sound {a b : Interval} {x y : ℝ}
    (hx : a.Contains x) (hy : b.Contains y) (hzero : 0 < b.lo ∨ b.hi < 0) :
    (a.div b).Contains (x / y) := by
  simpa only [Interval.div, div_eq_mul_inv, one_mul] using
    a.mul_sound hx (b.inv_sound hy hzero)

theorem Interval.square_sound {a : Interval} {x : ℝ} (hx : a.Contains x) :
    a.square.Contains (x ^ 2) := by
  rcases hx with ⟨hxlo, hxhi⟩
  constructor
  · dsimp [Interval.square]
    split_ifs with h
    · exact sq_nonneg x
    · by_cases hlo : a.lo ≤ 0
      · have hhi : a.hi < 0 := lt_of_not_ge (fun hhi => h ⟨hlo, hhi⟩)
        exact (min_le_right _ _).trans (by nlinarith)
      · have hlo' : 0 < a.lo := lt_of_not_ge hlo
        exact (min_le_left _ _).trans (by nlinarith)
  · dsimp [Interval.square]
    by_cases h : 0 ≤ x
    · exact (show x ^ 2 ≤ a.hi ^ 2 by nlinarith).trans (le_max_right _ _)
    · have h' : x ≤ 0 := le_of_not_ge h
      exact (show x ^ 2 ≤ a.lo ^ 2 by nlinarith).trans (le_max_left _ _)

def roundDown (scale x : ℝ) : ℝ := (⌊x * scale⌋ : ℝ) / scale
def roundUp (scale x : ℝ) : ℝ := (⌈x * scale⌉ : ℝ) / scale
def Interval.outward (scale : ℝ) (a : Interval) : Interval :=
  ⟨roundDown scale a.lo, roundUp scale a.hi⟩

theorem roundDown_le {scale x : ℝ} (h : 0 < scale) : roundDown scale x ≤ x := by
  exact (div_le_iff₀ h).2 (Int.floor_le (x * scale))

theorem le_roundUp {scale x : ℝ} (h : 0 < scale) : x ≤ roundUp scale x := by
  exact (le_div_iff₀ h).2 (Int.le_ceil (x * scale))

theorem Interval.outward_sound {a : Interval} {x scale : ℝ}
    (hx : a.Contains x) (h : 0 < scale) : (a.outward scale).Contains x :=
  ⟨(roundDown_le h).trans hx.1, hx.2.trans (le_roundUp h)⟩

theorem Interval.dyadic_sound {a : Interval} {x : ℝ} (bits : ℕ)
    (hx : a.Contains x) : (a.outward (2 ^ bits)).Contains x :=
  a.outward_sound hx (pow_pos (by norm_num) bits)

end
end StepDetection
