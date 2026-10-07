Evaluating directional confidence on fresh histories
====================================================

The combined directional reference detects **58 of 144** true slowdowns,
compared with **27** for the previous indexed-jump reference. Both produce
**zero false alerts among 216 null histories**. Detection rises from 18.8%
to 40.3%, a gain of **21.5 percentage points** on the same fresh inputs.

The simpler uniform-direction rule detects 47. Adding the tail bound
recovers **11 more histories**, all with 100 readings and positive noise
correlation. It preserves every uniform-only detection. Relative to
indexed-jump, the combined method gains 33 detections and loses two.

This is a substantial improvement under the declared model, with a large
remaining gap: the matching known-correlation oracle detects 108. Of the
combined method's 86 misses, 84 retain a verified explanation and two have
unfinished searches.

Frozen comparison
------------------

The `protocol <directional_protocol.rst>`_ and complete evaluation settings
were frozen in commit ``2278f6d`` before generating observations. The study
uses 360 versions of 120 fresh base histories, with seed labels 1600 and
1601 and the new ``directional-ar1-study`` random streams.

The design combines correlations -0.5, 0, and 0.7; lengths 40 and 100;
noise standard deviations 0.05 and 0.2; changes of 0%, 4%, 5%, 6%, and 8%
from baseline 10; and early, middle, and recent locations. There are 144
above-threshold histories and 216 null histories, including 72 exactly at
5%. The reporting claim is a slowdown greater than 5%, so exactly 5% counts
as null. Each correlation condition contains 48 positives and 72 nulls.

The correlation versions share their Gaussian innovations. Their outcomes
are dependent, and the small subgroup counts describe this sample rather
than precise population error rates. Every history follows the declared
stationary Gaussian AR(1), common-variance, at-most-one-change model.

What the comparisons mean
--------------------------

The previous indexed-jump reference predicts complete measurement histories,
integrating over possible plateau levels and noise amplitudes. The new
directional reference removes those quantities before checking the noise
pattern. This avoids spending predictive probability on levels and scales
that the data can eliminate directly.

The uniform-only variant uses a simple upper bound on the probability of an
unusual residual pattern. The combined variant also uses the tighter tail
bound near strong positive correlation. Both bounds concern the same
probability and share the same confidence budget. Their comparison isolates
what the tail calculation contributes beyond removing levels and scale.

Every unknown-correlation method keeps the 5% slowdown threshold, 1%
confidence budget, 4% reporting budget, and work limits of 4096 cells and
depth 16. The known-correlation oracle receives the true correlation and
uses the matching 4% reporting budget. Its count measures sensitivity with
that extra information; it is not a sample-by-sample ceiling on other rules.

The main tables use gated decisions: both the evidence rule and the
existing shared-floor reporter must accept the history. Standalone evidence
is also saved. Unresolved searches count as non-alerts and as misses on
above-threshold histories.

Overall results
----------------

.. list-table:: Gated detections and false alerts
   :header-rows: 1

   * - Method
     - Detected / 144
     - False alerts / 216
   * - ASV production reporting
     - 113
     - 25
   * - Shared-floor fit with existing reporting
     - 129
     - 37
   * - Conditional confidence
     - 15
     - 0
   * - Global full-history confidence
     - 23
     - 0
   * - Indexed-original
     - 24
     - 0
   * - Indexed-jump
     - 27
     - 0
   * - Uniform directional confidence
     - 47
     - 0
   * - Combined directional confidence
     - 58
     - 0
   * - Known-correlation oracle, reporting alpha=0.04
     - 108
     - 4

The production and shared-floor rows use their existing reporting rules;
they have no matching per-history 5% error guarantee. Their higher detection
counts come with more alerts on histories at or below the threshold.

All six unknown-correlation methods have identical standalone and gated
decisions in this study. The matching oracle has 111 standalone detections
and five false alerts; the shared-floor gate reduces these to 108 and four.
All four gated oracle false alerts occur on changes exactly at 5%.

The earlier indexed study's 28/144 result used different observations.
The 27/144 comparator here was rerun on these fresh histories, making the
31-detection net gain a paired comparison.

Where the gain comes from
-------------------------

.. list-table:: Detections by generating noise correlation
   :header-rows: 1

   * - Correlation
     - Indexed-jump / 48
     - Uniform / 48
     - Combined / 48
     - Oracle / 48
   * - -0.5
     - 17
     - 24
     - 24
     - 41
   * - 0
     - 10
     - 23
     - 23
     - 36
   * - 0.7
     - 0
     - 0
     - 11
     - 31

Removing levels and scale improves the negative-correlation and independent
groups. The tail bound supplies the entire improvement in the positively
correlated group. Its 11 additional detections include five 6% changes and
six 8% changes, with seven in quieter noise and four in louder noise.

History length matters. At 40 readings, indexed-jump detects 3/72 and both
directional variants detect 23/72. At 100 readings, indexed-jump and
uniform-direction confidence each detect 24/72; the combined method detects
35/72. Neither directional variant detects any of the 24 positively
correlated, 40-reading slowdowns, while the oracle detects 15.

Most sensitivity remains in the quieter histories: combined confidence
detects 53/72 at noise standard deviation 0.05, versus 24/72 for indexed-jump.
At standard deviation 0.2, the counts are only 5/72 versus 3/72. The oracle
detects 72/72 and 36/72, respectively. The improvement therefore leaves
substantial difficulty with noisier measurements.

Across early, middle, and recent locations, combined confidence detects
20/48, 19/48, and 19/48; indexed-jump detects 10/48, 10/48, and 7/48.
The result index retains all subgroup counts and paired comparisons.

What still prevents an alert
----------------------------

Combined confidence produces 58 certified alerts, 299 surviving
explanations, and three unresolved searches across all 360 histories.
Two unresolved histories are positive and one is null. Two hit the
4096-cell limit; one hits subdivision depth 16. All remain non-alerts.
There are no insufficient-variation outcomes.

Among its 86 missed true slowdowns, 84 have a concrete surviving location
and correlation. Raising the work limits addresses only the two unfinished
searches. The oracle detects 50 positive histories that combined confidence
misses; combined confidence gains none over the oracle on this sample.

The confidence check excludes the generating pair in one null history,
a positively correlated 4% change with 40 readings. Another explanation
survives, so the complete search returns no alert. This illustrates why
confidence exclusion at one pair is not itself a regression report.

The two detections lost relative to indexed-jump are both 8% changes in
100-reading, negatively correlated histories with noise standard deviation
0.2 and seed 1601: one early change and one middle change. Both retain a
candidate change after the first reading with correlation 7/8, or 0.875.

The `post-evaluation loss diagnosis <data/directional_v1_loss_diagnosis.json>`_
estimates the directional tail probability at those saved witnesses:

.. list-table:: Numerical diagnosis of the two lost detections
   :header-rows: 1

   * - Actual change location
     - Implemented probability upper bound
     - Numerical tail estimate
     - Confidence cutoff
   * - Early
     - 1
     - 0.0003212
     - 0.01
   * - Middle
     - 1
     - 0.01682
     - 0.01

The first estimate suggests a loose bound is retaining the witness. The
second suggests that even an exact evaluation of this same tail probability
would retain its witness. These are floating-point diagnostics, with
quadrature error estimates about 3.3e-9 and 2.0e-10; they are not interval
certificates and do not change the study decisions.

The next mathematical task is to derive useful certificates at moderate
correlations such as 0.875. The present matrix bound is designed for values
very close to 1. Then examine whether combining size and noise-pattern
evidence more efficiently can address the remaining retained explanations.
The middle example suggests that tighter probability bounds alone will not
recover every lost detection. These saved histories can now guide development;
any revised reporting rule needs a new freeze and fresh evaluation.

Evidence and validation
------------------------

The `result index <data/directional_v1_results.json>`_ contains all 24 method
summaries, paired comparisons, subgroup results, diagnostics, and hashes for
78 compressed archives. Those archives preserve all observations, method
outputs, and six covariance matrices, using about 9 MB in total.

Validation regenerated all 360 inputs exactly, reproduced every summary,
and checked the complete frozen source and environment chain. It verified
2160 method/history outcomes across the six unknown-correlation references:
194 alerts, 1955 surviving explanations, and 11 unresolved outcomes.
Alert checks reconstruct full interval coverage; witness checks independently
recompute reporting evidence with GLS. The directional verifier also checks
partial certificates. The combined alert certificates contain 15 intervals
proved by the new tail route.

Every compressed archive passed a round-trip check. The recovered run's
first 91 inputs and results remain byte-identical to the saved prefix, the
original manifest is intact, and every uniform-only alert is preserved by
the combined rule. Ten new study tests and five recovery tests pass, in
addition to the mathematical tests completed before this evaluation.

Execution recovery
-------------------

The first attempt stopped after saving 91 histories because an exact
fraction exceeded Python's default 4300-digit integer-to-text limit.
Recovery set ``PYTHONINTMAXSTRDIGITS=0`` for this trusted-data calculation,
preserved the saved prefix, and evaluated the remaining histories with
the same frozen sources, observations, calibration, and work limits.
Each completed case was then saved immediately before final assembly in
the declared order.

This permits longer integer strings while preserving numerical values.
The `recovery record <data/directional_v1_recovery.json>`_ preserves hashes
of the initial input and result files, original manifest, frozen settings,
and recovery wrapper. The original manifest remains intact. The
`running guide <running.rst>`_ includes replay, recovery, and verification
commands with the required serialization setting.

Scope
------

These are exact-arithmetic research references. Per-method elapsed times
include contention from six evaluation workers and concurrent verification;
they are not isolated benchmarks of ASV's production detector.

The interval inequalities use exact rational arithmetic for the supplied
observations. The inherited t/F critical values are computed numerically,
so reporting certificates remain conditional on those cutoffs. The model
does not cover changing variance, outliers, missing readings, or several
genuine changes. ASV's production detector is unchanged.
