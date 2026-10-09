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

The latest point-certificate refinement bounds individual covariance ranks.
Previously, each group of four directions used the weakest rank bound for all
four. The new point calculation bounds each rank separately, then combines
four lower bounds through their geometric mean in the Fourier integral. The
interval search still uses its existing interval certificate, so this addition
affects candidate points without changing interval coverage.

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
with correlation 99/128. Its earlier probability upper bound was **0.012290**.
The individual-rank point refinement lowers it to **0.010850**, an 11.7%
reduction, but it is still above 0.01. The numerical tail diagnostic is about
0.00773; that floating estimate is not a certificate, so the candidate
remains and the history still does not alert.

Across the seven saved histories, four alerts remain alerts and the other
three remain non-alerts. The replay verifies every result. The new bound
changes which explanations survive, but it does not change the alert
decisions in this saved set.

What this contributes to ASV
----------------------------

These changes live under ``experiments/step_detection/``. They are a
mathematical prototype for a possible future improvement. **ASV's production
step detector has not been changed.** On the independent 360-history
comparison, the Sturm-enabled reference detected 53/144 positive histories,
versus 52/144 for the control: one paired gain and no losses. Both had zero
alerts among 216 null histories. The candidate also left nine positive
searches unresolved, versus one for the control, so this is a small
descriptive gain with more abstentions, not evidence of a broad improvement.
See the `fresh evaluation report <sturm_fresh_report.rst>`_.

The Lean pilot checks some exact interval and determinant arithmetic. It does
not yet formalize the new Fourier probability inequality or the complete
statistical argument.

What remains
------------

The rank-group bound still misses the 0.01 cutoff at the saved 99/128
candidate. The fresh study's single gained alert does not resolve that
mathematical gap or justify a production change. Any next step should target
the remaining looseness in the tail bound, then use a new frozen study before
considering integration. The `unresolved-search diagnosis
<sturm_unresolved_diagnosis.rst>`_ shows that one extra interval split
certifies five of nine saved witnesses; four still have at least one
uncertified half. Those four fail at the determinant enclosure near
``rho = 1``; increasing precision to 512 bits and extensive subdivision did
not fix a representative case. A continuant/Bernstein prototype now
certifies 13 of 16 equal pieces in one 40-reading witness interval; three
pieces near ``rho = 1`` remain just above 0.01. Exact rational setup takes
tens of seconds, and the 32-piece follow-up exceeded two minutes. The next
mathematical task is to reduce coefficient cost and refine only the endpoint
pieces before another benchmark comparison.

Further detail
--------------

The `plain-language introduction <README.rst>`_ explains how ASV fits steps
to benchmark histories. The `mathematical derivation <sturm_refinement.rst>`_
develops the eigenvalue and Fourier bounds. The `replay report
<sturm_reporting.rst>`_ describes the seven saved cases, and the
`verified replay summary <data/sturm_reporting_v4_rank_groups_summary.json>`_
records their latest outcomes. The `fresh evaluation report
<sturm_fresh_report.rst>`_ describes the independent comparison. The
`experiment guide <running.rst>`_ gives reproduction commands.
