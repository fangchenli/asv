Which scoring change actually helps?
====================================

The shared noise-floor formulation is the strongest finding. Simply reducing
the segment penalty under the current floor formula causes extreme overfitting.
The shared floor makes lower penalties usable. The correlation cap improves
sensitivity at a fixed penalty, but its advantage largely disappears after
calibrating false-alert rates.

This updates the `previous study's recommendation <threshold_report.rst>`_.
Prioritize the shared floor and penalty calibration; keep the cap as an
experimental option. Production defaults remain unchanged.

What was frozen and measured
----------------------------

The `protocol <ablation_protocol.rst>`_ specifies eight combinations of:

* The current or shared-floor scoring formulation.
* A residual correlation cap of 1 or 0.5.
* A segment penalty of ``4*log(n)/n`` or ``2*log(n)/n``.

All scores use the same candidate pool. This includes production's visited
fits and the exact independent-error optimum at every segment count. The
independent reference does not establish global optimality for correlated
scores. The current detector also runs unchanged as a separate control.

For calibration, each floor/cap family searches eight penalty coefficients
on 1200 development histories. At false-alert budgets of 0%, 1%, and 5%, it
chooses the most sensitive qualifying setting, using predetermined tie rules.
Those choices, the source, and protocol were committed as ``1ae1b3e`` before
evaluating 2400 new histories. No held-out result changes a selected setting.

The new seeds are disjoint from the preceding study. Each phase uses the same
noise families and lengths as before. The held-out set contains 960
above-threshold histories, 960 flat or below-threshold histories, and 480
exactly at the 5% threshold. The last group does not count toward either
error rate. Counts below refer to histories, not individual alert records.

Separate the three components
-----------------------------

This table is the fixed factorial comparison. Each error column has 960
eligible histories. The penalty column is the coefficient multiplying
``log(n)/n``.

.. list-table:: Held-out factorial results
   :header-rows: 1

   * - Floor formulation
     - Correlation cap
     - Penalty
     - Missed slowdowns
     - Below-threshold alerts
   * - Current
     - 1
     - 4
     - 499
     - 2
   * - Current
     - 0.5
     - 4
     - 391
     - 9
   * - Shared
     - 1
     - 4
     - 506
     - 2
   * - Shared
     - 0.5
     - 4
     - 398
     - 9
   * - Current
     - 1
     - 2
     - 14
     - 456
   * - Current
     - 0.5
     - 2
     - 14
     - 456
   * - Shared
     - 1
     - 2
     - 298
     - 24
   * - Shared
     - 0.5
     - 2
     - 198
     - 56

Holding the floor and penalty at their current values, the cap reduces
misses from 499 to 391 but raises below-threshold alerts from 2 to 9.
With the shared floor and penalty 2, it reduces misses from 298 to 198
while raising below-threshold alerts from 24 to 56. It consistently changes
the sensitivity/false-alert tradeoff in these comparisons.

Changing only the floor at penalty 4 has little benefit here. At penalty 2,
however, the floor formulation is decisive: the current formula produces
below-threshold alerts on 47.5% of eligible histories. The shared formula
reduces that to 2.5% with cap 1, or 5.8% with cap 0.5. Floor and penalty
interact; their effects cannot be summarized as independent improvements.

The unmodified production search misses 506 histories and alerts on two
below-threshold histories. Expanding its candidate pool reduces misses by
seven. That small search effect is separated from the larger scoring effects
in the table above.

Why lowering the penalty breaks the current floor
-------------------------------------------------

With penalty 2 and the current floor, 2370 of 2400 histories get one segment
per observation. Every reading becomes its own level, including measurement
noise. The correlation cap cannot help these fits: their residuals are all
zero, so changing residual correlation changes nothing.

The current formula sets its floor using the smallest difference between
adjacent fitted levels when there are more than two levels. With one level
per reading, that is the smallest adjacent measurement difference. Two
readings happening to be close make this floor tiny. Taking its logarithm
then gives the exact fit a very favorable score, enough to overcome the
reduced complexity penalty.

For example, ``n100-d0-s0.02-r0-middle-seed301`` is a flat history near 10
with independent noise of standard deviation 0.2. At penalty 2 and cap 1:

.. list-table:: Lower wins within each score column
   :header-rows: 1

   * - Fit
     - Current candidate-dependent floor
     - Current score
     - Shared-floor score
   * - One level
     - 0.0100275
     - 2.98810
     - 1.50414
   * - 100 levels
     - 0.0001794
     - 0.58447
     - 9.21034

The shared floor is about 0.119821 for both candidates. Its profile loss
has the same exact-fit baseline for every candidate; fitting each observation
exactly leaves the full 100-segment penalty. The shared model selects one
level and reports no alert. The current formula selects 100 levels and
reports an alert in this genuinely flat history.

Scores have different normalizations across the two columns. Their numerical
values only rank candidates within a column. Also, the shared-floor factor
includes the profiled loss formula, not merely substitution of a constant
into the old logarithm.

Compare calibrated false-alert budgets
--------------------------------------

Matching development budgets does not guarantee equal held-out rates. These
tables show the rates actually observed, with each configuration selected
before held-out evaluation.

.. list-table:: Settings chosen under the 1% development budget
   :header-rows: 1

   * - Floor / cap
     - Selected penalty
     - Held-out miss rate
     - Held-out below-threshold alert rate
   * - Current / 1
     - 4
     - 52.0% (499)
     - 0.21% (2)
   * - Current / 0.5
     - 6
     - 53.3% (512)
     - 0.10% (1)
   * - Shared / 1
     - 3
     - 45.0% (432)
     - 1.15% (11)
   * - Shared / 0.5
     - 6
     - 53.3% (512)
     - 0.10% (1)

These realized false-alert rates are too different to call the sensitivity
comparison an equal-rate improvement. The shared/cap-1 setting slightly
exceeds its 1% development budget on new data. The strict cap needs a larger
penalty during calibration, giving back the sensitivity seen at a fixed
penalty.

.. list-table:: Settings chosen under the 5% development budget
   :header-rows: 1

   * - Floor / cap
     - Selected penalty
     - Held-out miss rate
     - Held-out below-threshold alert rate
   * - Current / 1
     - 3
     - 39.9% (383)
     - 5.00% (48)
   * - Current / 0.5
     - 3
     - 28.9% (277)
     - 6.04% (58)
   * - Shared / 1
     - 2
     - 31.0% (298)
     - 2.50% (24)
   * - Shared / 0.5
     - 3
     - 31.9% (306)
     - 2.40% (23)

The shared-floor settings have nearly matched realized rates. Cap 1 finds
eight more affected histories, with one extra below-threshold history alerted.
This is not evidence of a substantial advantage from cap 0.5. Relative to
the calibrated current-floor/cap-1 family, shared-floor/cap-1 reduces both
misses (383 to 298) and below-threshold alerts (48 to 24).

At the zero-false-alert development budget, all selected settings also have
zero observed below-threshold alerts on the new data. Cap-1 settings miss
655 histories with the current floor and 652 with the shared floor, using
penalty 6. Cap-0.5 settings need penalty 16 and miss 767 with either floor.
An empirical zero here is not a guarantee for future histories.

What this means for the next implementation
-------------------------------------------

Keep the shared-floor/cap-1 family as the leading experimental candidate.
Its penalty-2 setting detects 208 more affected histories than unchanged
production, with 22 additional below-threshold histories alerted. No
previously detected positive history loses its alert. Of its 24
below-threshold alerts, two occur on flat histories and 22 on 4% slowdowns.
This is a useful tradeoff, but it does not improve both metrics relative
to the current production default.

Location accuracy also improves: this candidate reports an alert within two
observations of the true change in 629 of 960 positive histories, compared
with 450 for production. The cap-0.5 shared-floor setting calibrated at the
5% budget locates 619. Full localization and noise-condition breakdowns are
saved with the evidence.

Before adopting a default, stress-test the frozen shared-floor settings on
outliers, missing observations, non-Gaussian noise, and changing noise levels.
The adjacent-difference floor estimator is still a heuristic. A robust floor
formulation and a well-calibrated floor estimator are separate requirements.
Keep the reporting threshold fixed while checking how reliably the chosen
false-alert budget carries over. Runtime is not the priority for this work.

Scope and reproducibility
-------------------------

These are paired factorial simulations using 10 development seeds and 20
fresh evaluation seeds. Configurations share random draws, and flat histories
are duplicated under the location labels. The 2400 histories are not 2400
independent trials, and the percentages do not estimate error rates across
real benchmark projects. The simulated noise families remain Gaussian AR(1).

The `running instructions <running.rst>`_ describe reproduction. The frozen
settings are in ``data/ablation_v1_frozen.json``. ``data/ablation_v1_results.json``
contains summaries, breakdowns, paired wins/losses against production, hashes,
and the example above. Eighteen ``ablation_v1_*.jsonl.gz`` archives preserve
all inputs and per-history results in ordered shards of 400 histories.

The source/protocol/extension hashes were checked before evaluation and again
after archiving. The component, threshold-pipeline, and exact-reference tests
passed 74 checks. The copied current score is tested directly against the
production scoring closure across every candidate in representative histories.
