import StepDetection.Interval

namespace StepDetection

def Interval.Enlarges (outer inner : Interval) : Prop :=
  outer.lo ≤ inner.lo ∧ inner.hi ≤ outer.hi

theorem Interval.enlarge_sound {outer inner : Interval} {x : ℝ}
    (hx : inner.Contains x) (h : outer.Enlarges inner) : outer.Contains x :=
  ⟨h.1.trans hx.1, hx.2.trans h.2⟩

end StepDetection
