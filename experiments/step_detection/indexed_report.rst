Evaluating candidate-specific prediction on fresh histories
===========================================================

The new construction improves detection, but it is still too conservative
for a practical replacement. Indexed-jump detects **28 of 144** true slowdowns,
compared with **20** for the conditional reference, **25** for global
full-history prediction, and **27** for indexed-original. All four produce
zero false alerts among 216 null histories. The matching known-correlation
oracle detects **114 of 144**.

The gain over the conditional reference is eight detections, or **5.6
percentage points** (13.9% to 19.4%). Most remaining misses have a concrete
surviving noise explanation. Increasing the search budget would not remove
those explanations.

Frozen comparison
-----------------

The evaluator, priors, external unit, reporting cutoffs, and search limits
were frozen in commit ``7d251fb`` before generating these observations.
The `protocol <indexed_protocol.rst>`_ specifies 360 versions of 120 fresh
paired histories, using the new ``indexed-ar1-study`` random streams and
seed labels 1500 and 1501. No development history is reused.

The design combines correlations -0.5, 0, and 0.7; lengths 40 and 100;
noise standard deviations 0.05 and 0.2; changes of 0%, 4%, 5%, 6%, and 8%
from baseline 10; and early, middle, and recent locations. All versions
satisfy the declared stationary Gaussian AR(1), constant-variance,
single-change model.

There are 144 positive histories with slowdowns above 5% and 216 null
histories, including 72 exactly at 5%. Each correlation condition has 48
positives and 72 nulls. The three versions of a base history share Gaussian
innovations; their outcomes are dependent. Small subgroup counts describe
this sample and do not establish precise population error rates.

The comparison retains the existing shared-floor fit and reporter. A gated
alert requires both the evidence rule and that reporter to accept the history.
Standalone evidence decisions are preserved separately in the artifacts.

The conditional reference is the earlier training/validation construction.
The global full-history reference uses a proper density mixed over all
locations. The two indexed references give each candidate location its own
proper density, either with the original independent-level prior or with
the midpoint/jump prior. Both retain half the original global density.
The oracle is supplied the generating correlation.

All four unknown-correlation references allocate 0.01 to confidence
construction and 0.04 to reporting. The matching oracle uses reporting
alpha=0.04. The new full-history references use external unit=1. Every search
keeps the frozen limits of 4096 cells and depth 16. Unresolved searches
count as non-alerts and as misses on positive histories.

What the comparison answers
---------------------------

A non-alert can mean that the data still admit an explanation with a slowdown
of at most 5%. For example, strong noise persistence may explain much of an
apparent step. Each reference searches for such an explanation across change
locations and correlations. It alerts only after ruling out every candidate.

The predictor helps decide which candidates remain plausible. The global
predictor spreads probability across all change locations. When checking one
particular location, much of that probability may describe changes elsewhere.
The indexed construction supplies a predictor for the location being checked,
while retaining half the global predictor. The original-prior indexed control
isolates this change. Comparing it with indexed-jump then measures what the
new jump prior adds.

This keeps the same error-budget argument: a false alert must eliminate the
true explanation along with all the others. At that true location, its own
predictor supplies the required confidence bound. We do not need to divide
the confidence budget by the number of locations. This argument depends on
using a proper predictor for each location, as derived in the
`baseline/jump analysis <baseline_jump_prior.rst>`_.

The oracle knows the true correlation, so it can skip this uncertainty search.
Its detection count shows how much sensitivity is available with that extra
information. It is a useful reference, not a guaranteed upper bound on the
other methods' decisions on every sample.

What changed in the alerts
--------------------------

All tables below use gated decisions: both the evidence check and the
shared-floor reporter must accept the history. The four unknown-correlation
methods have identical standalone and gated decisions on all 360 histories.
The reporter removes one standalone oracle detection (115 becomes 114).

.. list-table:: Overall results
   :header-rows: 1

   * - Method
     - Detected / 144
     - False alerts / 216
   * - Current production
     - 110
     - 33
   * - Shared-floor fit with existing reporting
     - 135
     - 39
   * - Direct independent-noise test, alpha=0.05
     - 115
     - 10
   * - Conditional
     - 20
     - 0
   * - Global full-history
     - 25
     - 0
   * - Indexed-original
     - 27
     - 0
   * - Indexed-jump
     - 28
     - 0
   * - Oracle, alpha=0.04
     - 114
     - 0

The direct test ignores correlation, explaining why its larger detection
count comes with false alerts in this design. It also uses reporting
alpha=0.05; the oracle at 0.04 is the main matched-budget comparison.
Zero observed false alerts does not establish a zero population error rate.

.. list-table:: Paired changes on the same positive histories
   :header-rows: 1

   * - Comparison
     - New detections
     - Lost detections
   * - Global full-history vs conditional
     - 9
     - 4
   * - Indexed-original vs global full-history
     - 2
     - 0
   * - Indexed-jump vs indexed-original
     - 1
     - 0
   * - Indexed-jump vs conditional
     - 10
     - 2

Restoring the full-history information supplies most of the net gain.
Location-specific prediction adds two detections beyond that, and changing
the jump prior adds one more. The small extra gain does not establish a
broad advantage for the jump prior. Indexed-jump still loses two detections
that the conditional construction found.

.. list-table:: Detections by generating correlation (48 positives per row)
   :header-rows: 1

   * - Correlation
     - Conditional
     - Global full-history
     - Indexed-original
     - Indexed-jump
     - Oracle, alpha=0.04
   * - -0.5
     - 12
     - 15
     - 16
     - 17
     - 42
   * - 0
     - 7
     - 9
     - 10
     - 10
     - 37
   * - 0.7
     - 1
     - 1
     - 1
     - 1
     - 35

All methods in this table have zero false alerts among the 72 null histories
in each row. Positive correlation remains especially difficult: each
unknown-correlation method finds only one of 48 true changes, versus 35 for
the oracle.

.. list-table:: Where detections occur
   :header-rows: 1

   * - Group
     - Positives
     - Conditional
     - Global full-history
     - Indexed-original
     - Indexed-jump
     - Oracle, alpha=0.04
   * - 40 readings
     - 72
     - 2
     - 1
     - 2
     - 3
     - 54
   * - 100 readings
     - 72
     - 18
     - 24
     - 25
     - 25
     - 60
   * - 6% change
     - 72
     - 8
     - 9
     - 9
     - 9
     - 47
   * - 8% change
     - 72
     - 12
     - 16
     - 18
     - 19
     - 67
   * - Noise SD 0.05
     - 72
     - 17
     - 21
     - 23
     - 24
     - 71
   * - Noise SD 0.2
     - 72
     - 3
     - 4
     - 4
     - 4
     - 43
   * - Early change
     - 48
     - 1
     - 8
     - 8
     - 8
     - 38
   * - Middle change
     - 48
     - 11
     - 11
     - 12
     - 12
     - 44
   * - Recent change
     - 48
     - 8
     - 6
     - 7
     - 8
     - 32

These rows overlap: each history contributes to its length, change size,
noise, and location groups. All four unknown-correlation references and the
matched oracle have zero false alerts within every subgroup, including all
72 histories exactly at the threshold. The result index preserves all
breakdowns, including the 0% and 4% null groups.

Why most misses remain
----------------------

.. list-table:: Search outcomes over all 360 histories
   :header-rows: 1

   * - Method
     - Certified alert
     - Surviving explanation
     - Unresolved
   * - Conditional
     - 20
     - 340
     - 0
   * - Global full-history
     - 25
     - 334
     - 1
   * - Indexed-original
     - 27
     - 331
     - 2
   * - Indexed-jump
     - 28
     - 330
     - 2

Every method retains all 360 generating location/correlation pairs in its
confidence set. No case has insufficient variation, and none of the three
full-history methods disables confidence exclusion at the true location
because of numerical fallback. No false-alert rejection routes occur.

Both unresolved indexed searches reach the 4096-cell limit. One is a null
5% history, and one is a true 6% slowdown; both have negative correlation,
100 readings, noise SD 0.2, and a recent change. The global full-history
reference is unresolved only on the null history. Thus only **one of the
116 indexed-jump misses** could potentially change by increasing the work
limit. The other 115 retain a verified explanation.

For a concrete example, consider
``positive-indexed-ar1-study-n100-d0.08-s0.005-early-seed1500``.
The true levels are 10 and 10.8, the change is after reading 25, the noise
standard deviation is 0.05, and the generating correlation is 0.7. The oracle
alerts. Indexed-jump retains location 25 with correlation
``16383/16384``, approximately 0.999939. At that pair, neither reporting test
rejects the below-threshold explanation. This is a surviving mathematical
explanation under the current rule, not an unfinished search. The other
seed with the same design has the same witness.

Next mathematical question
--------------------------

The immediate task is to understand why these highly persistent explanations
survive. Keep this evaluation fixed and use its saved witnesses for a new
development analysis. Decompose the evidence at those witnesses into what
the observations say about correlation and what is lost by integrating over
unknown plateau levels and noise scale.

A specific candidate is to remove those nuisance quantities before building
the correlation confidence set. At a fixed candidate location, subtract each
plateau's arithmetic mean, then normalize the resulting residual vector.
If that location is correct, the direction no longer depends on the two true
levels or the overall noise scale; its distribution still depends on
correlation. Derive that distribution
and a valid confidence construction, including its behavior near correlations
-1 and 1, before implementing or evaluating it. The residuals are correlated,
so an ordinary independent-residual test would not suffice.

This asks whether the remaining prediction cost can be reduced without
choosing priors from the evaluation data. Keep both existing indexed variants
as controls. Another parameter sweep or a larger search budget would not
answer this question.

The follow-up `residual-direction derivation <residual_direction.rst>`_ now
establishes that nuisance-free distribution and its endpoint limits. It also
checks the saved correlated witness and derives a route to direct tail
calibration without changing this evaluation.

Evidence and reproduction
-------------------------

The `result index <data/indexed_v1_results.json>`_ contains all 20 method
summaries, paired comparisons, subgroup breakdowns, confidence diagnostics,
and SHA-256 hashes for 74 compressed archives. Those archives preserve every
input, result, and the six covariance matrices. See the `running instructions
<running.rst>`_ for evaluation and verification commands.

Before generating observations, the focused analytical and harness suite
passed 130 tests. After evaluation, verification regenerated all 360 inputs,
reproduced every summary, checked the frozen source/environment chain, and
checked all six covariance matrices. It reconstructed all 100 alert
certificates across the four references with exact interval bounds and
complete coverage. It checked all 1335 surviving witnesses for exact
confidence membership and independently recomputed their reporting decisions
with GLS. The five unresolved method/history outcomes remain non-alerts.
All 120 independent-control direct/oracle decisions agree at alpha=0.05.
Every archive was decompressed and checked against its uncompressed content.

.. list-table:: Reference elapsed times in seconds
   :header-rows: 1

   * - Method
     - Median, n=40
     - Median, n=100
     - Maximum
   * - Conditional
     - 0.281
     - 13.191
     - 86.635
   * - Global full-history
     - 0.181
     - 13.641
     - 70.055
   * - Indexed-original
     - 0.200
     - 13.011
     - 75.894
   * - Indexed-jump
     - 0.196
     - 12.638
     - 108.592

These exact-arithmetic research references take seconds rather than the
milliseconds of the ordinary detector. Times were recorded during six-worker
execution and include contention; they are not isolated speed benchmarks.
The current priority is sensitivity, since faster execution cannot reject
a valid witness under the same evidence rule.

Scope of the result
-------------------

This comparison measures evidence for a slowdown under a single-change,
constant-variance Gaussian AR(1) model. It does not establish behavior under
changing variance, outliers, gaps, or several genuine changes. The priors and
external unit were fixed before evaluation; choosing them from these outcomes
would require another fresh comparison.

The certificate calculations use exact rational interval bounds. Predictive
densities and analytical cutoffs still originate in floating-point numerical
calculations, so the reference does not supply an end-to-end exact-arithmetic
coverage guarantee. See the numerical qualifications in the
`full-history derivation <full_history_confidence.rst>`_. Production defaults
remain unchanged.
