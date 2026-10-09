Fresh evaluation of the Sturm reporting reference
=================================================

This study compares the complete reporting search with its split-aware Sturm
confidence route disabled and enabled. The enabled route includes the new
rank-group point bound. The control and candidate use the same generated
history, confidence allocation, reporting tests, and work limits. This is an
experimental Python reference; ASV's production C++ detector is unchanged.

Result
------

The candidate produced **53 alerts among 144 positive histories**, compared
with **52 of 144** for the control. It gained one alert and lost none. Neither
method alerted on any of the **216 null histories**.

The candidate also left **9 positive searches unresolved**, compared with
**1** for the control. Seven stopped at the depth limit and two at the cell
limit. An unresolved search is a non-alert. The extra alert therefore comes
with more abstentions elsewhere; it is not a broad increase in detection.

.. list-table:: Overall decisions
   :header-rows: 1

   * - Method
     - Positive alerts
     - Null alerts
     - Positive unresolved
   * - Sturm disabled (control)
     - 52 / 144
     - 0 / 216
     - 1
   * - Sturm enabled (candidate)
     - 53 / 144
     - 0 / 216
     - 9

The one gained alert occurs for a 6% slowdown at the middle split, with
positive AR(1) correlation and the higher noise level. The candidate's
remaining unresolved searches are concentrated in 100-reading histories
with 6% or 8% changes, especially under positive correlation. Full counts
by correlation, length, change, noise, and location are in the committed
summary JSON.

.. list-table:: Positive histories by generating correlation
   :header-rows: 1

   * - Correlation condition
     - Positive alerts, control
     - Positive alerts, candidate
     - Unresolved, control
     - Unresolved, candidate
   * - Negative AR(1)
     - 24 / 48
     - 24 / 48
     - 1
     - 2
   * - Independent
     - 22 / 48
     - 22 / 48
     - 0
     - 1
   * - Positive AR(1)
     - 6 / 48
     - 7 / 48
     - 0
     - 6

The paired result is small: one of 144 positive histories changed from a
non-alert to an alert, and none changed in the opposite direction. The 360
rows represent 120 base histories, each evaluated at three correlations;
those three versions share their random stream and are dependent. These
counts are descriptive and do not establish a reliable improvement in
detection probability.

The Python reference's median elapsed time over all histories was 0.45 s for
the control and 0.57 s for the candidate. Positive histories took longer
because more of their searches continued. These measurements describe this
research implementation and are not ASV production timing results.

Reproduction and artifacts
--------------------------

The source, settings, and protocol were frozen before evaluation in
`sturm_fresh_v2_frozen.json <data/sturm_fresh_v2_frozen.json>`_. Reproduce the
run using the `protocol <sturm_fresh_protocol.rst>`_. The archived
`manifest <data/sturm_fresh_v2_manifest.json>`_ records the environment and
source hashes. Compressed `inputs <data/sturm_fresh_v2_inputs.jsonl.gz>`_ and
`full per-history records <data/sturm_fresh_v2_records.jsonl.gz>`_ preserve
the generated cases and both method outcomes. The `summary
<data/sturm_fresh_v2_summary.json>`_ contains all subgroup counts, statuses,
paired alert changes, and integrity hashes.

The inputs were regenerated from the frozen generator and match the saved
records. The subgroup summaries were recomputed from those records, and both
compressed archive hashes and round trips were checked. The study does not
formalize the complete probability argument in Lean and does not justify a
production change. Its result is a small gain in saved-study sensitivity,
more unresolved searches, and no observed null alerts in this sample.
