How much sensitivity survives unknown correlation?
==================================================

The first unknown-correlation reference is too conservative to be a useful
replacement. It detects **16 of 144** above-threshold slowdowns, compared with
**103 of 144** when the true correlation is supplied at the same reporting
budget. Both have zero false alerts among 216 null histories in this sample.
The search completes every case. The loss comes from allowing too many noise
explanations, rather than from running out of search time.

Frozen comparison
-----------------

Code, protocol, analytical cutoffs, and certification budgets were frozen in
commit ``601bed9`` before generating observations. The `protocol
<ar1_protocol.rst>`_ uses 360 versions of 120 fresh paired histories:
correlation -0.5, 0, or 0.7; lengths 40 and 100; noise standard deviations
0.05 and 0.2; changes of 0%, 4%, 5%, 6%, or 8%; and early, middle, or
recent locations. Early changes occur inside the training prefix.

Each correlation condition has 48 positive and 72 null histories, including
24 exactly at 5%. Across all versions there are 144 positives and 216 nulls.
Exactly 5% remains null. The three versions of a base history share Gaussian
innovations, so aggregate counts are descriptive rather than independent
binomial trials. Small subdivisions have limited precision.

All cases satisfy the declared stationary Gaussian AR(1), constant-variance,
one-change model. Only the oracle receives the true correlation. The unknown-
correlation rule scans the full stationary range using its joint confidence
set and interval certificates. Work limits remain 4096 cells and depth 16;
unresolved searches count as non-alerts.

The old direct test and the original oracle both use reporting alpha=0.05.
A second oracle uses alpha=0.04, matching the unknown rule's reporting budget.
The unknown rule spends the remaining 0.01 on constructing its confidence set.
Comparing these two oracle budgets separates an allocation cost from the
cost of uncertain covariance. All gated methods use the same shared-floor
fit and existing reporter.

What the alerts show
--------------------

An alert is correct here when the generating increase exceeds 5%. A gated
method reports only when both its evidence check and the existing reporter
accept the history. All counts below are gated:

.. list-table:: Overall comparison
   :header-rows: 1
   :widths: 52 24 24

   * - Method
     - Detected / 144
     - False alerts / 216
   * - Shared-floor fit with existing reporting
     - 132
     - 33
   * - Direct independent-noise test, alpha=0.05
     - 106
     - 4
   * - Known-correlation oracle, alpha=0.05
     - 105
     - 2
   * - Known-correlation oracle, alpha=0.04
     - 103
     - 0
   * - Unknown-correlation reference
     - 16
     - 0

Reducing the oracle's reporting budget from 0.05 to 0.04 costs two detections.
Allowing uncertainty in correlation then costs another 87, with no gained
detections. Thus the budget split explains little of the sensitivity loss.
Standalone and gated unknown-correlation decisions agree on all 360 histories;
the existing reporter removes none of its alerts.

.. list-table:: By generating correlation
   :header-rows: 1
   :widths: 25 25 25 25

   * - Correlation
     - Direct test
     - Oracle, alpha=0.04
     - Unknown correlation
   * - -0.5
     - 36 / 0
     - 39 / 0
     - 12 / 0
   * - 0
     - 36 / 0
     - 35 / 0
     - 4 / 0
   * - 0.7
     - 34 / 4
     - 29 / 0
     - 0 / 0

Each cell is detections / false alerts, with 48 positives and 72 nulls per
row. The unknown rule detects one of 72 positives at length 40 and 15 of 72
at length 100. It detects four of 72 true 6% changes and 12 of 72 true 8%
changes. None of the positively correlated slowdowns survives its evidence
requirement, including the larger changes.

Why the confidence set is the bottleneck
----------------------------------------

The reference keeps pairs of possible change locations and correlations.
To alert, it must reject every retained pair as an explanation of a slowdown
no larger than 5%. A broad set gives noise more ways to explain the data.

All 360 generating pairs remain in their confidence sets. But the median
width of the union over locations is 2, the entire stationary range from -1
to 1. Even at the generating location alone, the median width is 1.913.
In 273 histories, including 116 of 144 positives, some retained location
allows correlations arbitrarily close to 1. This is much more uncertainty
than the three generating values, -0.5, 0, and 0.7, might suggest.

There are 16 certified alerts and 344 surviving explanations, with zero
unresolved searches and zero insufficient-variation cases. The largest
search visits 586 cells, below the frozen limit of 4096. Each non-alert has
an explicit retained pair that fails both rejection tests. More search time
cannot remove such a pair under the current rule.

A concrete missed slowdown
--------------------------

Consider ``independent_constant-ar1-study-n100-d0.08-s0.005-early-seed1400``:
100 readings, independent noise with standard deviation 0.05, and an 8%
increase from a baseline of 10 after reading 25.

At the true correlation of zero, the size statistic is 25.23, far above
its cutoff of 1.873. The oracle alerts. The unknown rule retains the full
correlation range at location 25 and finds a surviving explanation at
rho=4095/4096, approximately 0.999756.

At such strong correlation, consecutive readings can move together while
providing little information about the common baseline. The overall noise
scale is also unknown. This matters for a percentage claim: to establish
``later > 1.05 * earlier``, the baseline must be known well enough, even if
the absolute jump is easy to see.

Holding the location at 25, an independent GLS calculation gives:

.. list-table:: The same readings under increasingly persistent noise
   :header-rows: 1

   * - Assumed correlation
     - Size statistic
     - Extra-step statistic
   * - 0.99
     - 3.821
     - 6.340
   * - 0.9999
     - 1.112
     - 6.332
   * - 0.999999
     - 0.115
     - 6.332

The size statistic drops below 1.873; the extra-step statistic stays below
its cutoff of 16.924. Neither test can reject the persistent-noise explanation.
The predictor trained on the first half also sees the early mean change and
chooses correlation 0.95, its allowed prediction limit. That limit applies
only to prediction; it does not exclude larger candidate correlations from
the confidence set.

The mathematical obstruction near correlation one
-------------------------------------------------

This is a follow-up diagnosis made after evaluation. It changes no decisions
or settings. It examines the generating location of the 87 detections lost
relative to the oracle at alpha=0.04.

Let ``t`` be a candidate boundary and define consecutive differences
``Delta_i = y_i - y_(i-1)``. As rho approaches 1, the rescaled AR(1)
precision matrix approaches the first-difference precision matrix. A
two-level fit explains the difference at ``t``. An extra boundary can explain
one additional difference. Define::

    Q = sum(Delta_i**2 for i != t)
    M = max(Delta_i**2 for i != t)
    limiting extra-step statistic = (n - 3) * M / (Q - M)

The remaining differences determine whether a second step is needed. They
need not look unusual just because correlation approaches 1.

Meanwhile, for a positive percentage threshold and nonzero residual
variation, the size statistic tends to zero. The two-level estimates remain
finite, but uncertainty about the common level grows without bound. The
contrast ``later - 1.05*earlier`` contains that common level with coefficient
-0.05, so its uncertainty grows too.

The confidence condition is ``S_t(rho) <= K``, where ``S_t`` is the quadratic
conditional residual cost from the `derivation <unknown_correlation.rst>`_.
If ``S_t(1) < K`` and the limiting extra-step statistic is below its cutoff,
continuity gives surviving explanations arbitrarily close to 1. The endpoint
itself is excluded from the stationary model; nearby values suffice.

This sufficient obstruction holds at the generating location for **72 of
the 87** lost oracle detections: 19 under negative correlation, 24 under
independent noise, and 29 under positive correlation. The other 15 losses
have surviving explanations but are not explained by this particular test.
For the example above, ``S_t(1)=0.002257 < K=0.004781`` and the limiting
extra-step statistic is 6.332, well below 16.924.

The next mathematical step
--------------------------

Improve the confidence construction before optimizing its runtime. A concrete
candidate is to use both ends of the history:

1. Construct the current conditional confidence set in forward order with
   failure budget 0.005.
2. Construct it on the reversed history with failure budget 0.005, mapping
   each reversed location back to the original location.
3. Retain only pairs accepted in both directions. Keep reporting alpha=0.04.

Stationary Gaussian AR(1) noise is reversible. A single mean change remains
a single change after reversal, although its direction reverses. Each
conditional confidence argument therefore applies. The probability that
either set excludes the true pair is at most 0.005 + 0.005 = 0.01; no
independence between the two sets is needed.

This may help when a change contaminates one training prefix: the opposite
end can provide a cleaner predictor. Each direction still contributes a
quadratic confidence inequality, so interval certification can reject a
region whenever either inequality excludes it. The stricter budget in each
direction could also reduce sensitivity. This proposal needs a derivation,
analytical checks, and a new frozen evaluation before claiming improvement.

Evidence and practical limits
-----------------------------

The `result index <data/ar1_v1_results.json>`_ contains all method summaries,
paired comparisons, confidence diagnostics, the follow-up limit calculations,
and SHA-256 hashes for 37 compressed archives. The archives preserve all
inputs, records, and six covariance matrices. `Running instructions
<running.rst>`_ describe how to reproduce the frozen study.

Validation regenerated all 360 inputs, checked every covariance matrix and
summary, and checked that gated alerts are subsets of the existing reporter's
alerts. Every certified alert was reconstructed with exact polynomial
arithmetic, including complete interval coverage. All 344 surviving pairs
were checked for exact confidence membership and independently recomputed
with GLS. All 120 independent-control direct/oracle decisions agree at
alpha=0.05. The focused test suite passed 295 tests.

The exact-arithmetic reference's median elapsed time was 0.226 seconds for
length 40 and 10.87 seconds for length 100, with a maximum of 34.66 seconds.
These are elapsed times during six-worker execution and include contention;
they are not production timing benchmarks. The immediate obstacle is the
large sensitivity loss.

Zero false alerts in this small paired sample does not establish a zero
population error rate. The mathematical error bound depends on stationary
Gaussian AR(1) noise with constant marginal variance and at most one mean
change. Changing variance, heavy tails, multiple changes, and recovery remain
outside this evaluation. Production defaults remain unchanged.
