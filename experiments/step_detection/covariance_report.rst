What does knowing the noise structure recover?
==============================================

The covariance-aware gate reduces false alerts from 85 to 8 among 1440 null
histories, while true detections fall from 585 to 499 among 960 positive
histories. It removes much of the overconfidence seen in the stress study,
with a measurable sensitivity cost. Independent noise that becomes quieter
is an exception: there, it improves both detection and false-alert counts.

This is useful evidence for pursuing covariance estimation. It is not a
uniform improvement in sensitivity, and it relies on supplying the true
noise structure. Production remains unchanged.

Frozen comparison
-----------------

The `oracle protocol <covariance_protocol.rst>`_, runner, and covariance
matrices were frozen in commit ``decb689`` before generating observations.
The evaluation uses 2400 fresh versions of 400 base histories, with the
same six conditions as the previous stress study. Each condition has 160
positive histories at 6%/8% and 240 null histories at 0%/4%/5%, including
80 exactly at 5%. Histories have 40 or 100 readings, with changes in the
middle or leaving four/ten later readings.

The shared-floor fit and all previous reporting rules are unchanged. The
new test receives the true covariance shape, which describes relative
variances and correlations, but estimates the overall noise scale from
the history. It uses the same t/F cutoffs and the same 5% reporting threshold.
``Oracle alone`` means this new evidence test by itself. ``Oracle gate``
also requires an existing shared-floor regression report.

All six conditions satisfy the oracle's Gaussian one-change model. Only
the independent constant-variance control satisfies the old direct test's
model. The `derivation <known_covariance.rst>`_ gives the error argument;
the simulation measures sensitivity and observed false alerts. A finite
sample's rate is not itself the mathematical guarantee.

Results by noise condition
--------------------------

Each cell is ``true detections / false alerts``, with denominators 160 and
240 respectively. The direct and oracle gates use the same shared-floor fit
and existing reporter; only their evidence calculations differ.

.. list-table:: Paired comparison on fresh histories
   :header-rows: 1

   * - Noise
     - Old direct gate
     - Oracle gate
     - Oracle alone
   * - Independent, constant
     - 123 / 2
     - 123 / 2
     - 123 / 2
   * - Independent, rises
     - 103 / 7
     - 72 / 1
     - 72 / 1
   * - Independent, falls
     - 68 / 5
     - 89 / 1
     - 90 / 1
   * - Correlated, constant
     - 118 / 19
     - 112 / 2
     - 120 / 4
   * - Correlated, rises
     - 96 / 35
     - 48 / 1
     - 48 / 1
   * - Correlated, falls
     - 77 / 17
     - 55 / 1
     - 59 / 1

Every decision and the numerical statistics agree in the independent
constant-variance control. There were no numerical abstentions in any
condition. The low false-alert counts therefore do not come from failed
numerical calculations suppressing alerts.

Across the paired histories, the oracle gate gains 32 true detections and
loses 118 relative to the old gate. It removes 77 false alerts and adds none.
The standalone comparison is different: the oracle removes 204 false alerts
and adds two, reducing the total from 212 to 10. Existing reporting filters
both newly added standalone false alerts, which occur at true 4% changes.
The oracle is not simply a subset of the old evidence decisions.

All eight remaining gated false alerts occur at exactly 5%: 8/480 in that
subset. They are distributed as two in the control, two under correlated
constant variance, and one in each remaining condition. Exactly 5% remains
a null case even if its sample estimate happens to exceed the threshold.

For context, unchanged production detects 508 positive histories with 135
false alerts on these same observations. Existing shared-floor reporting
detects 738 with 246 false alerts, and the joint-sign gate detects 112 with
five. The oracle's 499 detections and eight false alerts retain much more
sensitivity than the sign reference while substantially reducing false alerts
relative to existing reporting. The sign reference's guarantee still excludes
the correlated conditions.

Where sensitivity improves and where it falls
---------------------------------------------

The strongest sensitivity gain occurs for a short quiet plateau following
a long noisy one. Under independent falling variance, recent-change detections
rise from 23/80 to 46/80, with zero false alerts in either method's 120 recent
null histories. The old common-variance estimate lets earlier noise dominate
the uncertainty assigned to the quiet recent readings. Supplying the correct
relative variances removes that penalty.

The reverse pattern costs sensitivity. Under independent rising variance,
recent-change detections fall from 57/80 to 27/80, while false alerts fall
from 6/120 to 1/120. The later readings really are noisier than the old pooled
estimate suggested. Accounting for that uncertainty prevents some correct
sample decisions as well as incorrect ones.

With both rising variance and correlation, detections fall from 96/160 to
48/160 while false alerts fall from 35/240 to 1/240. Known covariance does
not make four noisy, related readings as informative as many independent
ones. This is the largest sensitivity loss and needs to remain visible when
judging the benefit of the correction.

Across all conditions, the oracle detects 231/480 recent changes and 268/480
middle changes. It detects 165/480 true 6% increases and 334/480 true 8%
increases. Knowing the noise structure supports useful sensitivity, but does
not reliably distinguish every small excess over the 5% threshold.

The failures change in the way the mathematics predicts
-------------------------------------------------------

At the generating split, the old gate's 85 false alerts reject through size
alone in 26 cases, shape alone in 44, and both in 15. For the oracle's eight
false alerts, those counts are seven, one, and zero. Both routes improve.
These counts describe rejection of a valid explanation at that split; they
do not uniquely attribute the global alert across all searched locations.

Three saved examples make the effects concrete. At length 40, the size
cutoff is approximately 1.799 and the shape cutoff is 16.280.

* ``independent_up-direct-stress-n40-d0.05-s0.005-middle-seed1209`` has a
  true 5% increase at split 20. The old size statistic is only 0.800, but
  its shape statistic is 26.628: noisy later readings make an extra boundary
  look compelling. With the correct covariance, the shape statistic falls
  to 9.081 and the size statistic is 0.751. Split 20 survives, preventing
  the false alert.
* ``independent_down-direct-stress-n40-d0.06-s0.005-recent-seed1200`` has
  a true 6% increase at split 36. Its size statistic rises from 1.494 to
  3.203 when the calculation recognizes that the last four readings are
  quieter. The old test retains split 36 as a below-threshold explanation;
  the oracle rejects all explanations and recovers the true detection.
* ``correlated_up-direct-stress-n40-d0.06-s0.005-middle-seed1200`` also
  has a true 6% increase. The size statistic falls from 2.289 to 1.393
  when noise structure is accounted for. The corrected shape statistic
  is only 3.247. The true split survives, and this positive history becomes
  a miss. Correct uncertainty can reduce sensitivity on an individual case.

What the oracle is allowed to know
----------------------------------

The supplied covariance matrix is fixed throughout the search over mean
locations. The true mean split is never supplied as a search constraint.
However, the variance transition occurs at the nominal mean-change location
in this design, so knowing the covariance reveals the noise-transition
location. That is additional information. These results measure performance
with that information available; they do not establish the same performance
when covariance must be estimated from observations.

The six versions of a base history share Gaussian innovations. Pairing
supports direct comparisons of noise transformations, but makes the combined
2400 versions dependent. Aggregate counts are descriptive. The fresh seeds
also mean that counts from this run should be compared within this run,
rather than subtracted from the preceding stress study's counts.

Next: uncertainty about covariance itself
-----------------------------------------

The oracle results justify continuing: correcting covariance removes most
observed false alerts while retaining useful detection, and it can recover
changes hidden by an excessively large pooled variance estimate. The next
mathematical problem is how much of this survives without knowing R.

Start with a restricted family: constant marginal variance and an unknown
AR(1) correlation. This isolates one uncertain parameter before adding
unknown variance ratios and noise-transition locations. Derive a valid
set of possible correlations in the presence of an unknown mean step, and
require evidence to survive every covariance in that set. Charge the set's
coverage failure to the total error budget as described in the
`derivation <known_covariance.rst>`_. A finite grid or a fitted residual
correlation alone does not establish that guarantee.

Keep this completed evaluation as a diagnostic set. A practical extension
should be frozen before comparison on new histories, with both the oracle
and existing reporter retained so losses from estimating covariance can
be separated from the limitations of the reporting construction.

Artifacts and validation
------------------------

The `frozen configuration <data/covariance_v1_frozen.json>`_ contains the
inherited settings, source hashes, environment versions, and matrix hashes.
The `result index <data/covariance_v1_results.json>`_ contains all summaries,
paired comparisons, rejection routes, example statistics, and archive hashes.
Its 49 compressed archives preserve all 2400 numerical inputs and outputs
plus the 24 covariance matrices.

All 230 focused tests passed. Archive verification regenerated the inputs
exactly, matched every input to its covariance, reproduced the summaries,
and checked gate subsets and generating-split rejections. Frozen hashes,
identity-control agreement, and document parsing were checked. See
`reproduction instructions <running.rst>`_ for commands and output fields.
