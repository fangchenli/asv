Current status of the step-detection experiment
================================================

This page gathers the current result from the step-detection research. It
separates the experimental changes from ASV's production detector and records
what the latest saved-case replay does and does not establish.

What the detector is deciding
------------------------------

ASV measures how a benchmark's runtime changes across software revisions.
The step detector looks for stretches of similar runtimes separated by a
lasting change. A timing history can also be explained by ordinary noise, so
the detector checks whether each proposed step and noise model is convincing
enough to report.

In this experiment, each possible step location and correlation is treated
as a candidate explanation. A mathematical upper bound is computed for the
probability of seeing data this unusual under that explanation. If the bound
falls below the 0.01 confidence cutoff, that candidate is ruled out. A history
is called an alert only when the search rules out every candidate. One
surviving candidate, or an interval the search cannot resolve, is enough to
withhold an alert.

What changed in the experiment
------------------------------

The experimental certificate now uses three refinements:

* **More informative covariance bounds.** Exact Sturm counts give lower
  bounds on covariance eigenvalues. The proposed step split helps identify a
  long plateau whose directions retain more information.
* **A sharper density bound.** The Fourier calculation keeps separate
  strengths for successive groups of four directions and integrates their
  product directly with partial fractions. This avoids replacing all groups
  with the weakest one. The inputs and final bound are rational; equal group
  strengths are separated by a tiny downward adjustment so the exact formula
  remains valid.
* **Faster exact sign counting.** The Sturm recurrence uses integer-scaled
  minors instead of repeatedly simplifying fractions. The signs, and
  therefore the eigenvalue counts, are unchanged.

The experimental reporting search also uses depth 17 instead of 16. One
additional split was enough to certify the interval that previously stopped
the search as unresolved.

What the saved replay found
---------------------------

At correlation 101/128, the previous probability upper bound was about
0.012103. The direct Fourier integral lowers it to **0.0087203** at the point
and **0.0087751** over a small interval. Both are below 0.01, so this
particular explanation is rejected.

The search then finds another candidate: a step after the second reading,
with correlation 99/128. Its probability upper bound is **0.012290**, above
the cutoff, so it remains a plausible explanation. The history still does
not alert.

Across the seven saved histories, four alerts remain alerts and the other
three remain non-alerts. The replay verifies every result. The new bound
changes which explanations survive, but it does not change the alert
decisions in this saved set.

What this contributes to ASV
----------------------------

These changes live under ``experiments/step_detection/``. They are a
mathematical prototype for a possible future improvement. **ASV's production
step detector has not been changed.** The latest work also has not been
evaluated on a fresh benchmark set, so it does not establish a better
detection rate or false-alert rate in real use.

The Lean pilot checks some exact interval and determinant arithmetic. It does
not yet formalize the new Fourier probability inequality or the complete
statistical argument.

What remains
------------

The next mathematical target is the 99/128 candidate, whose bound remains
above 0.01. After improving the certificate, it should be replayed on saved
cases and then evaluated on fresh histories with false-alert rates measured
alongside detection rates. Only that broader evidence can tell us whether
the prototype belongs in ASV's production implementation.

Further detail
--------------

The `plain-language introduction <README.rst>`_ explains how ASV fits steps
to benchmark histories. The `mathematical derivation <sturm_refinement.rst>`_
develops the eigenvalue and Fourier bounds. The `replay report
<sturm_reporting.rst>`_ describes the seven saved cases, and the
`verified replay summary <data/sturm_reporting_v3_depth17_summary.json>`_
records their outcomes. The initial `experiment guide <running.rst>`_ gives
reproduction commands.
