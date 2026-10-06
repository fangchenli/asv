Does the shared-floor improvement survive messier data?
=======================================================

Mostly, but the original false-alert budgets do not hold uniformly. The
shared-floor setting calibrated at the 5% budget still has fewer misses
and fewer below-threshold alerts overall than the calibrated current-floor
setting. With missing observations, however, its extra sensitivity also
produces more false alerts. Positive outliers and falling noise variance
push its below-threshold alert rate above 5%.

Keep the shared floor as the leading experimental formulation. These results
do not establish a production false-alert guarantee or justify replacing the
current conservative default.

Frozen settings, new conditions
-------------------------------

We reused the full-correlation-cap settings selected in the
`component study <ablation_report.rst>`_, without recalibration. The protocol,
runner, and inherited settings were committed as ``7980fb5`` before generating
the new evaluation inputs. Seeds 400--409 were unused in earlier studies.

The `protocol <robustness_protocol.rst>`_ defines a Gaussian control plus
heavy-tailed innovations, positive outliers, random missing readings, a gap
at the change, increasing noise amplitude, decreasing noise amplitude, and
outliers combined with missing readings. There are 960 histories per condition:
480 above the 5% reporting threshold and 480 flat or with a 4% slowdown.
The total is 7680 histories. Exactly 5% is omitted from this study.

All methods keep ASV's actual filtering and regression-reporting function.
Missing readings are removed before fitting, and fits are mapped back to
original revision coordinates before reporting. Evaluation moves a true
boundary inside a missing gap to the first retained observation afterward.
This measures the observable boundary, not exact revision recovery in the gap.

Overall results
---------------

Each error column below has 3840 eligible histories. Calibration budgets
refer to the earlier development data; percentages here are measured on
the new control and stress cases.

.. list-table:: Aggregate error rates
   :header-rows: 1

   * - Method / inherited budget
     - Missed above-threshold histories
     - Below-threshold histories alerted
   * - Current production default
     - 2242 (58.4%)
     - 7 (0.18%)
   * - Current floor / 5%
     - 1696 (44.2%)
     - 208 (5.42%)
   * - Shared floor / 5%
     - 1347 (35.1%)
     - 163 (4.24%)
   * - Shared floor / 1%
     - 1881 (49.0%)
     - 39 (1.02%)
   * - Shared floor / 0%
     - 2798 (72.9%)
     - 0

Against the calibrated current-floor setting at the same 5% development
budget, the shared floor reduces misses by 349 and below-threshold alert
histories by 45. Against the unchanged production default, it gains 896
detected positive histories and loses one, while adding 156 below-threshold
alert histories. Those are different comparisons: it improves the calibrated
tradeoff, but produces substantially more false alerts than the default.

The stricter shared-floor setting preserves lower false-alert rates at the
cost of more misses. Its 1% setting gains 361 detected positives over
production and adds 32 below-threshold alert histories. The zero-budget
setting has no observed below-threshold alerts, but misses nearly three
quarters of the true above-threshold slowdowns. The corresponding current-floor
zero-budget setting misses 2806, so the two are very similar at that end.

Behavior differs by condition
-----------------------------

This table compares the two settings inherited from the 5% development
budget. Each miss count is out of 480 positive histories, and each
false-alert count is out of 480 negative histories.

.. list-table:: Condition-specific counts
   :header-rows: 1

   * - Condition
     - Current-floor misses
     - Shared-floor misses
     - Current-floor false alerts
     - Shared-floor false alerts
   * - Gaussian control
     - 199
     - 152
     - 30
     - 13
   * - Heavy tails
     - 136
     - 92
     - 23
     - 21
   * - Positive outliers
     - 177
     - 118
     - 32
     - 27
   * - Random missing
     - 229
     - 162
     - 5
     - 15
   * - Gap at change
     - 218
     - 189
     - 29
     - 13
   * - Noise amplitude increases
     - 286
     - 264
     - 32
     - 22
   * - Noise amplitude decreases
     - 244
     - 229
     - 47
     - 28
   * - Outliers plus missing
     - 207
     - 141
     - 10
     - 24

The shared-floor 5% setting reaches 5.63% below-threshold alerts with positive
outliers and 5.83% when noise amplitude falls. Combined outliers and missing
readings reach 5.00%. Its 1% setting reaches 2.08% with heavy tails, 1.46%
when noise grows, and 1.67% when noise falls. The inherited calibration is
therefore not a stable per-condition bound.

Changing noise amplitude hurts sensitivity markedly. The shared-floor 5%
setting misses 31.7% of positive Gaussian-control histories, versus 55.0%
when noise increases and 47.7% when it decreases. The noise deviations are
tripled in one half of those histories; the underlying performance level
does not change unless the case includes a genuine slowdown.

Heavy tails are not uniformly harder in this construction. The innovations
have the same variance as the Gaussian control, but differ in their central
spread and rare extremes. The observed median adjacent-difference scale is
about 0.667 times the paired control's scale, and more changes are detected.
This result applies to the declared Student-t construction; it does not
establish generic robustness to arbitrary heavy-tailed measurements.

Random missingness retains about 50 of the original 70 observations on
average across the two lengths. The gap condition retains 63 on average.
Fitting and penalty selection use these retained lengths. The scale estimator
also uses adjacent retained observations, without correcting correlation for
the number of skipped revisions.

What the remaining false alerts look like
-----------------------------------------

Of the shared-floor 5% setting's 163 below-threshold alert histories,
12 are truly flat and 151 contain a genuine 4% slowdown. That does not mean
all 151 alerts locate the true change: only 65 have an alert within two original
observation indices of the observable boundary. Incorrect boundary placement
and overstated change size both need attention.

All 39 false-alert histories from the stricter shared-floor 1% setting
contain a genuine 4% slowdown and use two fitted segments. Nineteen locate
the real change within the tolerance. Fewer segments alone do not eliminate
uncertainty about the boundary or the size of the change.

For a concrete correctly located example, consider
``missing_random-n40-d0.04-s0.02-r0-recent-seed401``. The true level rises
from 10 to 10.4 at observation 36. After missing readings are removed:

.. list-table:: Fitted levels in that example
   :header-rows: 1

   * - Segment
     - Usable readings
     - Fitted level
     - Mean absolute residual
   * - Before the change
     - 28
     - 9.904
     - 0.145
   * - After the change
     - 3
     - 10.543
     - 0.094

The fitted increase is 6.45%, even though the true increase is 4%. The
reporter compares the fitted level difference, about 0.639, with the larger
of its noise allowances and 5% of the earlier level, about 0.495. The
difference exceeds that threshold, so it emits an alert. Both shared-floor
penalties 2 and 3 select this fit. The current production default selects
one level and emits none.

This example identifies a reporting question: how much evidence do three
post-change readings provide that the increase exceeds 5%? A point estimate
and a within-segment residual scale do not directly answer that question.
Any confidence rule also needs to account for the boundary being selected
from the same data; treating it as known would miss part of the uncertainty.

Next step
---------

Keep the settings fixed and treat this as the end of their current evaluation,
not another opportunity to adjust them to the observed failures. The shared
floor remains promising, but a default needs a clearer statistical meaning
for an alert near the threshold.

Investigate a reporting rule that evaluates evidence for an above-threshold
change, accounting for usable observation counts, serial correlation, and
boundary selection. Derive its assumptions before implementing it, then
calibrate and test it on new histories. Continue to assess the noise-scale
estimator separately, particularly under changing variance and missing gaps.
Production defaults remain unchanged.

Evidence and limitations
------------------------

The study uses ten new random seeds and paired configurations. Histories
share random draws, and some flat cases repeat under different location
labels. The percentages describe this synthetic design, not independent
binomial trials or measured real-project error rates. The control uses a
256-observation warmup; it is not a byte-for-byte replay of the previous
Gaussian sample. All comparison methods receive identical inputs here.

The `running instructions <running.rst>`_ describe reproduction.
``data/robustness_v1_frozen.json`` pins the inherited settings and source
hashes. ``data/robustness_v1_results.json`` contains summaries, condition
breakdowns, paired changes against the Gaussian control and production,
scale diagnostics, and the example above. Thirty-two condition-specific
``robustness_v1_*.jsonl.gz`` files preserve all inputs and per-history results.

The runner verifies frozen source/protocol/extension hashes before execution.
The robustness, component, threshold-pipeline, and exact-reference tests pass
88 checks, including missing-gap coordinate mapping and agreement with the
actual public production reporting pipeline.
