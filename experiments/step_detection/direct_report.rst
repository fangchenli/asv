Does a direct threshold test recover useful sensitivity?
========================================================

Yes, under an explicit noise model. The direct test, used to gate the fixed
shared-floor detector, detects 530 of 960 above-threshold slowdowns versus
119 for the joint-sign reference. It produces 10 false alerts among 1440
null histories, versus zero for that reference. Recent-change detections
rise from zero to 225 of 480.

This is a substantial sensitivity gain over the mathematical references.
It remains a tradeoff against existing reporting: the unchanged production
default detects 605 true slowdowns with 147 false alerts on the same data.
The direct gate is promising for further study, with production unchanged.

What the new rule does
----------------------

For every possible change location, the test checks whether the increase
exceeds 5% and whether the proposed two plateaus leave an unmodeled step.
It alerts only if every location is rejected by at least one of those checks.
A surviving below-threshold explanation prevents an alert.

The `protocol and derivation <direct_protocol.rst>`_ allocate 4% error
probability to a Student t test of ``later - 1.05*earlier`` and 1% to an
extra-boundary F scan. At any true below-threshold location, their combined
rejection probability is at most 5%. A global false alert must also reject
that true location, so the same bound holds while accounting for the unknown
boundary. All locations are checked; none is treated as known after fitting.

This proof assumes independent Gaussian noise with one common unknown
variance, positive plateau levels, and at most one true change. Sharing the
variance lets the longer plateau help estimate noise for a short later
plateau. Gaussian means and medians coincide, so the population target is
consistent with ASV's median-based fitting. The evidence calculation uses
means; ASV's fit stays unchanged.

The Gaussian assumption is stronger than the earlier sign references.
Laplace results below are a stress check, without the formal Gaussian
guarantee. Correlation, changing variance, and multiple true changes remain
outside this model.

This comparison changes both the testing construction and the noise
assumptions. Any gain belongs to that combination; the experiment does not
isolate how much comes from testing the target directly versus borrowing
precision through the Gaussian common-variance model.

Frozen design
--------------

Code, protocol, inherited detector, and analytical critical values were
committed as ``861abe3`` before generating evaluation histories. No cutoff
was adjusted afterward. The Student t cutoffs are approximately 1.7988 at
length 40 and 1.7690 at length 100; the corresponding maximum-F cutoffs are
16.2796 and 16.4271. Unlike the sign calibration, these critical values do
not have Monte Carlo calibration uncertainty.

The study uses 2400 fresh histories with seed labels 800--819: 1440 null
histories at 0%/4%/5%, and 960 positive histories at 6%/8%. Half are Gaussian
and half Laplace, with equal variance. Lengths, noise levels, and middle/recent
locations match the preceding study. Complete case identifiers produce
distinct random streams, so previous histories are not reused.

False alerts are counted against the claim that the underlying increase
exceeds 5%, including when a sample estimate happens to exceed 5% on a true
5% history. The existing reporter uses that sample estimate, so its decisions
answer a weaker question. These rates also cannot be compared directly with
older studies that omitted exactly-5% histories from their null cases.

The shared-floor detector retains ``shared-r1-c2``. All existing fits and
reporting decisions are reused unchanged. Each evidence gate retains an
existing alert only if its history-level test passes. Standalone evidence
is also measured to identify any limit imposed by the existing detector.
The certified claim concerns the initial-to-final change size, not the
detector's reported revision or individual steps.

Overall comparison
-------------------

True-detection counts have denominator 960. False-alert counts have
denominator 1440 and include true 5% changes.

.. list-table:: Results on the same fresh histories
   :header-rows: 1

   * - Method
     - True slowdowns detected
     - False alerts
   * - Production default
     - 605 (63.0%)
     - 147 (10.21%)
   * - Shared floor, existing reporter
     - 801 (83.4%)
     - 232 (16.11%)
   * - Shared floor, original reference gate
     - 70 (7.3%)
     - 0
   * - Shared floor, joint-sign gate
     - 119 (12.4%)
     - 0
   * - Shared floor, direct-test gate
     - 530 (55.2%)
     - 10 (0.69%)

The direct gate gains 411 true detections over the joint-sign gate and loses
none, while adding 10 false alerts. Compared with the existing shared-floor
reporter, it removes 222 false alerts and loses 271 true detections. Compared
with production, it gains 25 true detections, loses 100, and removes 137
false alerts without adding any. It therefore improves sensitivity over
the references, but does not dominate the existing reporters.

The direct test alone detects 531 positive histories and alerts on 14 null
histories. Gating removes four false alerts, all on true 4% changes, and loses
one true 6% detection. Thus the existing reporter has a small but observable
effect on the combined decision. Both older references have identical
standalone and gated decisions on these histories.

Gaussian results and Laplace stress results
-------------------------------------------

Each family has 480 positive and 720 null histories.

.. list-table:: Model-specific comparison
   :header-rows: 1

   * - Method
     - Gaussian detections
     - Gaussian false alerts
     - Laplace detections
     - Laplace false alerts
   * - Existing shared-floor reporter
     - 382
     - 111
     - 419
     - 121
   * - Joint-sign gate
     - 55
     - 0
     - 64
     - 0
   * - Direct-test gate
     - 259
     - 1
     - 271
     - 9

Under Gaussian noise, the direct gate detects 54.0% of positive histories
and alerts on 0.14% of null histories. Under Laplace noise, those figures
are 56.5% and 1.25%. All ten gated false alerts are on true 5% changes:
one of 240 Gaussian threshold histories and nine of 240 Laplace threshold
histories. There are no gated false alerts on flat or 4% histories.

The Laplace aggregate hides a less favorable small subgroup. Among its
80 low-noise, exactly-5% histories, eight alert: 10%. The corresponding
Gaussian count is one of 80. This is an exploratory subgroup description,
not a calibrated estimate of a uniform error bound. It warrants further
investigation and does not extend the Gaussian theorem to Laplace noise.

Where sensitivity returns
--------------------------

Each row below contains 240 positive histories across both noise families.

.. list-table:: Detections by length and change location
   :header-rows: 1

   * - History
     - Existing shared-floor reporter
     - Joint-sign gate
     - Direct-test gate
   * - 40 readings, middle change
     - 218
     - 1
     - 133
   * - 40 readings, recent change
     - 157
     - 0
     - 101
   * - 100 readings, middle change
     - 230
     - 118
     - 172
   * - 100 readings, recent change
     - 196
     - 0
     - 124

The recent-change cases no longer have the sign reference's short-segment
barrier. Across lengths, the direct gate detects 225 of 480 recent changes
(46.9%), and 305 of 480 middle changes (63.5%).

Remaining misses concentrate in noisier histories. At standard deviations
0.05, 0.2, and 0.5, it detects 318, 170, and 42 of 320 positive histories,
respectively. For true 6% changes it detects 197 of 480; for true 8% changes,
333 of 480. The rule still pays a sensitivity cost to establish that an
increase exceeds 5%, especially when the excess is small relative to noise.

A recent change with four later readings
-----------------------------------------

Consider ``gaussian-n40-d0.06-s0.005-recent-seed800``. Its true timing rises
from 10 to 10.6 at observation 36, leaving four later readings. The fixed
detector finds that split and fits medians 10.00094 and 10.60878. The joint-sign
reference cannot certify the increase, while the direct test can.

At length 40 the size statistic must exceed 1.7988, or the shape statistic
must exceed 16.2796. Three of the 39 checked locations illustrate the logic:

.. list-table:: How alternative locations are rejected
   :header-rows: 1

   * - Proposed location
     - Size t statistic
     - Maximum shape F statistic
     - Rejection route
   * - 35
     - -0.774
     - 140.62
     - An extra step is needed
   * - 36
     - 3.412
     - 7.69
     - The increase exceeds 5%
   * - 37
     - 1.405
     - 95.42
     - An extra step is needed

At locations 35 and 37, adding boundary 36 greatly improves the fit, so
their two-plateau explanations are rejected. At 36, the estimated increase
is large relative to its uncertainty. All other locations are rejected too,
leaving no below-threshold explanation within the model. This example
illustrates the method; its frequency and errors are measured across the
full evaluation rather than inferred from this one success.

Why retain the existing alert gate
----------------------------------

A fixed-location null can also be rejected through its shape test, even
when its size statistic is negative. Under the Gaussian model, that error
is covered by the allocated 1% shape allowance. It helps explain why the
standalone test sometimes alerts on a below-threshold history.

For example, ``gaussian-n100-d0.04-s0.005-recent-seed813`` has a true 4%
increase. At its true boundary the size statistic is -8.207, but the shape
statistic is 16.824, just above its cutoff 16.427. The standalone direct test
rejects all locations. The existing reporter sees fitted levels about
9.9973 and 10.3733, below its 5% threshold, so the combined gate emits no
alert. The archives preserve all rejection routes rather than presenting
every direct rejection as a positive size statistic.

Next decision
--------------

Keep the direct gate as the leading candidate for further reporting work
under its declared model. Its improvement is large enough to justify testing
the assumptions that made it possible: one common noise variance and
independence. Investigate varying variance and serial correlation with the
cutoffs held fixed before considering a production default. Derive an
appropriate extension before claiming coverage under either condition.

Laplace behavior also needs more scrutiny, especially near exactly 5% with
short plateaus. Do not adjust the completed study's cutoffs to its observed
false alerts. Multiple changes, recoveries, weights, and repeated publication
remain separate extensions. The current evidence supports continued
development, not a general-purpose false-alert guarantee for ASV histories.

Evidence and reproduction
--------------------------

The `running instructions <running.rst>`_ give the commands.
``data/direct_v1_frozen.json`` preserves the inherited calibration, analytical
critical values, source/protocol hashes, design, and SciPy version.
``data/direct_v1_results.json`` records the aggregate results, breakdowns,
paired comparisons, examples, and numerical archive checksums. Compressed
``direct_v1_*_inputs_*.jsonl.gz`` and ``direct_v1_*_records_*.jsonl.gz`` files
preserve all observations and per-history results.

The new direct-test and pipeline checks pass 25 tests. Together with the
earlier reporting and study checks, 157 tests pass. Independent regressions
verify every t and F statistic on small examples; rational calculations
verify interval moments. Tests also check short post-change plateaus, the
all-locations decision, exact-threshold fixtures, timing-unit invariance,
and the unchanged detector pipeline. Production behavior is unchanged.
