Joint sign calibration for reporting a slowdown above 5%
========================================================

This protocol is frozen before generating the new benchmark histories. The
question is whether joint calibration of the median constraints improves
sensitivity over the conservative reference while retaining a stated
false-alert bound. Production code remains unchanged.

The reporting target and model
------------------------------

Use the target from `the derivation <reporting_uncertainty.rst>`_:
``D = final_median - 1.05*initial_median``. An alert claims ``D > 0``.
The null includes exactly 5% as well as smaller changes. Assume independent
observations, positive plateau medians, unit weights, and at most one true
change. The true boundary may be anywhere. The guarantee concerns one history,
not the exact change location, repeated publication, or multiple benchmarks.

The original reference divides its error allowance among all intervals. Many
of those intervals overlap and their failures occur together. The new method
calibrates their joint maximum, using the common distribution of signs at
the true medians. It changes the interval bounds, while retaining all possible
boundaries and the same confidence-set inversion.

Derivation of the joint statistic
---------------------------------

For an interval of length ``m`` with ``c`` positive signs, define::

    p(m,c) = min(1, 2*P(Binomial(m, 1/2) <= min(c,m-c)))

Smaller values mean a more extreme imbalance. Rank the distinct exact rational
probabilities across lengths 1 through n in descending order. The largest
rank is the most extreme. Let ``T`` be the largest rank over every contiguous
interval in a sequence of n signs. Integer ranks preserve ties without
floating-point comparisons of very small probabilities.

For independent continuous observations, signs relative to their respective
true plateau medians are independent fair coins, even when the noise spreads
or shapes differ. The same is true across the boundary: each sign is centered
on its own true median. Thus the joint distribution of ``T`` can be calibrated
without knowing the levels, noise scales, or boundary.

For each proposed constant interval and level, require its rank to be no
larger than a calibrated cutoff. Inverting the count constraints gives an
order-statistic interval for that level. Intersect the constraints over all
subintervals of each proposed plateau, keep every feasible split and the
constant model, then minimize D over those explanations, exactly as in the
reference. An empty set produces ``incompatible_model`` and no alert.

Why this covers an unknown boundary: if the true median sign sequence has
``T <= cutoff``, every interval wholly inside a true plateau accepts the
true level. The true step function therefore belongs to the confidence set.
Any false alert under ``D <= 0`` requires ``T > cutoff``. Intervals crossing
the true boundary need no common-median coverage claim. This proof handles
the boundary search without estimating one boundary and treating it as fixed.
Ties at the median can be handled conservatively by coupling to independent
randomized median signs; the evaluation below uses continuous noise.

Calibration uncertainty is included
-----------------------------------

For each retained length n in {40, 100}, generate B = 8191 independent fair
sign sequences, using ``random.Random(550+n)``. Save every maximum T. The
per-history error allowance is alpha = 1/20. The total probability that
calibration fails its stated bound is at most 1/100, allocated equally to
the two lengths: epsilon = 1/200 for each.

Choose the smallest integer k such that::

    P(Binomial(B, 1-alpha) >= k) <= epsilon

The cutoff is the kth sorted calibration maximum. The rank and its binomial
failure bound are computed with integer arithmetic and rational comparisons.
If no finite rank qualifies, accept all sign sequences rather than claim
unsupported evidence.

To see why this works, let q be a population (1-alpha) quantile of T. A cutoff
whose true exceedance probability is greater than alpha lies below q. That
requires at least k calibration draws below q. Their probability is bounded
by the binomial tail above. Ties only make the bound conservative.

Consequently, with probability at least 99% over the two calibrations, both
frozen cutoffs have per-history false-alert probability at most 5% for every
history satisfying the stated model at their respective lengths. This is a
calibration confidence statement. It is not an unconditional 5% guarantee
with no calibration uncertainty, nor a claim that each empirical evaluation
cell must contain at most 5% alerts.

No benchmark outcomes select the statistic, cutoff, or detector. There is no
development sweep. Freeze the source hashes, protocol, design, inherited
detector setting, calibration maxima, ranks, and cutoffs in
``data/reporting_v1_frozen.json`` and commit them before evaluation.

Fixed detector and reporting comparisons
----------------------------------------

Inherit ``shared-r1-c2`` from the earlier component study's 5% development
budget: shared floor, correlation cap 1, penalty coefficient 2, and floor
factor 0.5 times the median adjacent absolute difference. This fitting choice
stays fixed even though the new evaluation observations are independent.
Reuse the original candidate pool and scoring implementation. Verify their
frozen hashes before calibration and evaluation.

Compare six history-level decisions:

* Unchanged production fitting and reporting.
* The fixed shared-floor fit with the existing reporter.
* The original union-bound reference, on its own.
* The jointly calibrated reference, on its own.
* Existing shared-floor alerts gated by the original reference.
* Existing shared-floor alerts gated by the joint reference.

The last two methods require both the existing alert and the new evidence.
Gating can only remove existing alerts, so it preserves the evidence rule's
false-alert bound regardless of the detector's selection. It certifies an
above-threshold initial-to-final change, not the locations or individual
steps returned by the existing reporter. Record those existing positions as
detector output; do not call them confidence-certified locations. Standalone
decisions reveal cases where fitting or the existing reporter limits the gate.

Fresh evaluation histories
---------------------------

Generate the Cartesian product of:

* Lengths 40 and 100.
* Relative increases 0%, 4%, 5%, 6%, and 8% from a baseline of 10.
* Noise standard deviations 0.05, 0.2, and 0.5.
* A middle change or a recent change leaving max(4, n/10) later readings.
* Independent Gaussian noise or independent Laplace noise of equal variance.
* Seed labels 600 through 619.

There are 2400 histories: 1440 null and 960 positive. Exactly 5% contributes
480 null histories and is also reported separately. Laplace scale is
``sigma/sqrt(2)``. Seed the random generator with the SHA-256 integer of the
complete case identifier, including distribution, length, change, noise,
location, and seed label. Cases do not reuse random draws across conditions
or duplicate flat histories under different location labels. Methods receive
the same history. Save complete numerical inputs, not just seed labels.

Report false-alert and miss counts overall and by noise family, length,
change size, noise magnitude, and location. Also report confidence-set
incompatibility and paired gains/losses versus the original reference.
Do not choose another cutoff after inspecting these results. Gaussian and
Laplace observations satisfy the declared independence model; this study
does not validate correlated, weighted, missing, or multiple-change histories.

Acceptance and next decision
----------------------------

Analytical checks must verify exact binomial tolerance ranks, enumerate tiny
sign sequences to check the joint statistic and its inversion, and reproduce
the original reference when given its cutoff. Pipeline checks must verify
that the shared-floor fit and original reporting are unchanged and that
gating cannot create a new alert.

The practical question is whether the joint method retains materially more
6%/8% detections than the original reference at the same nominal error
allowance. Report the answer even if the remaining power is poor. Successful
calibration does not establish a useful production default. Correlation,
multiple changes, and recovery semantics still require separate extensions.
