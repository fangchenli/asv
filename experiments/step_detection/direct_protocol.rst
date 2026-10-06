Testing the below-threshold explanations directly
=================================================

The last `reporting study <reporting_report.rst>`_ found that a confidence
set for the whole step function was too conservative. This protocol instead
tests whether any plausible change location still supports an increase of
at most 5%. Derive the rule, freeze its critical values and code, then
evaluate new histories. Production fitting and reporting stay unchanged.

The idea in ordinary terms
--------------------------

Suppose a history looks like timings near 10 followed by timings near 10.6.
One explanation places the change in the middle. Another places it after
the first reading, so its long later plateau mixes readings near both 10
and 10.6.

For each possible location, ask two questions:

* Are the estimated levels precise enough to establish an increase above 5%?
* Does the proposed pair of plateaus leave a further step in its residuals?

The first question can rule out a below-threshold explanation at the correct
location. The second can discard a badly placed location whose apparent
plateau actually contains two levels. An alert requires every possible
location to be rejected by at least one of the two checks. A surviving
below-threshold explanation is enough to prevent an alert.

This avoids bounding all possible plateau levels at once. It also avoids
treating the detector's selected location as known.

The stronger model assumption
------------------------------

Assume ``y_i = mu_i + error_i``, where the errors are independent Gaussian
variables with the same unknown variance ``sigma**2 > 0``. The positive
level ``mu_i`` is constant or changes once, at an unknown location. The
initial and final levels are ``a`` and ``b``. Gaussian means and medians
coincide, so the claim remains::

    D = b - 1.05*a > 0

Sharing one variance lets the long plateau help estimate the noise affecting
a short plateau. This is useful when only four later readings are available.
It is also a real assumption: the procedure does not have its stated
guarantee if noise variance changes, observations are correlated, or the
underlying history contains multiple changes. The evaluation includes
Laplace noise as a stress test outside the Gaussian guarantee.

The shared-floor detector continues to fit medians. Means are used only in
this evidence calculation. For non-Gaussian asymmetric distributions, a mean
contrast would also change the target, so this rule cannot be presented as
a general test about medians.

Test the size at one specified location
---------------------------------------

Let t be a proposed split, leaving t earlier observations and n-t later
observations. Fit their means ``a_hat`` and ``b_hat``. Let ``Q_t`` be the sum
of squared residuals around those two means. Then::

    s_t**2 = Q_t / (n-2)

    T_t = (b_hat - 1.05*a_hat)
          / [s_t * sqrt(1/(n-t) + 1.05**2/t)]

Both levels contribute uncertainty; that is why the baseline factor is
squared in the denominator. Under the fixed-location null ``D = 0``, the
numerator is centered Gaussian and independent of the residual sum, while
``Q_t/sigma**2`` has a chi-squared distribution with n-2 degrees of freedom.
Their ratio has a Student t distribution with n-2 degrees of freedom.
For ``D < 0``, its upper-tail rejection probability is smaller.

This uses the pooled-variance principle described in the
`NIST t-test reference <https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/t_test.htm>`_,
with the contrast changed to ``b - 1.05*a``. Allocate 0.04 error probability
to this check and reject when ``T_t > t_quantile(0.96, n-2)``.

Test whether that location gives two constant plateaus
------------------------------------------------------

For each extra boundary j distinct from t, fit three means on the segments
defined by the sorted pair (t, j). Let ``Q_tj`` be their residual sum of
squares. The improvement from the extra boundary gives::

    F_tj = (Q_t - Q_tj) / [Q_tj/(n-3)]

The two-level model is nested in this three-level model. Under a true
two-level model at t, the additional parameter explains only Gaussian noise.
Its squared projection divided by sigma squared is chi-squared with one
degree of freedom, independently of the remaining residual sum with n-3
degrees of freedom. Therefore ``F_tj`` has an F(1, n-3) distribution. This is
the residual-sum-of-squares logic of an
`ANOVA F test <https://itl.nist.gov/div898/handbook/eda/section3/eda354.htm>`_.

There are n-2 choices of extra boundary. Allocate a total error probability
0.01 to their scan and reject the two-level shape when any F statistic
exceeds ``F_quantile(1 - 0.01/(n-2), 1, n-3)``. Adding the individual tail
bounds controls this scan at 0.01 without assuming its tests are independent.
Computationally, the largest F comes from the smallest three-segment residual
sum, so save that extra boundary and F statistic for each t.

Why checking every location controls false alerts
-------------------------------------------------

At a specified true location t with ``D <= 0``, either rejection route can
be wrong. Their combined probability is at most ``0.04 + 0.01 = 0.05``.

The full null is the union of these fixed-location nulls. Reject it only
when every fixed-location null is rejected. If there is a true location t,
a global false alert necessarily rejects that particular true null. Its
probability is therefore at most 0.05. Checking more locations makes it
harder to alert because all must be rejected, so there is no additional
multiple-testing factor across t. This is an intersection-union test.

A constant positive history is included by taking a=b at any split; its
contrast is ``-0.05*a <= 0``. No separate constant-history test is required.
All n-1 splits are considered, including one-observation plateaus. The common
variance assumption still gives residual degrees of freedom n-2, with noise
estimated from the other observations. Require n >= 4 so the extra-boundary
test has positive residual degrees of freedom.

The proof concerns one history under the declared model and does not depend
on which fit ASV selects. Requiring an existing detector alert as an additional
gate can only reduce the false-alert probability. Surviving ``null_splits``
are explanations the test has not rejected, not a confidence interval for
the true boundary. A rejection certifies neither an exact revision nor
ASV's individual multiple-change reports. If the one-change assumption is
wrong, rejection may reflect that model failure instead of an excessive
initial-to-final increase.

Calibration is analytical and frozen
------------------------------------

Use alpha=0.05, with 80% assigned to the size test and 20% to the shape test.
This allocation is a design choice, not an optimized value. It reserves most
of the error allowance for the target contrast while permitting rejection
of clearly misplaced boundaries. No development sweep is run.

Compute Student t and F critical values with SciPy and save them for lengths
40 and 100. The Gaussian model supplies finite-sample distributions, so no
Monte Carlo calibration or extra calibration-failure allowance is needed.
Floating-point quantiles and statistics are checked numerically; mathematical
coverage refers to the exact distributions and real-arithmetic tests.

Freeze the source/protocol hashes, SciPy version, design, critical values,
inherited shared-floor setting, and old joint-sign calibration in
``data/direct_v1_frozen.json``. Commit the freeze before generating evaluation
histories. The new source must leave the previous frozen modules unchanged.

Numerical implementation uses incremental interval means and squared errors,
after scaling observations by their largest absolute value to avoid square
overflow. Zero residual sums in noiseless fixtures use limiting ratios:
a positive numerator over zero is infinity, and zero over zero is zero.
The statistical model assumes positive variance. These conventions make
deterministic examples reviewable without claiming a separate degenerate
noise theorem.

Evaluation and comparisons
---------------------------

Reuse the preceding study's generator with fresh seed labels 800--819.
Each full case identifier seeds a distinct SHA-256-derived random stream.
The design is the Cartesian product of lengths 40/100, changes
0%/4%/5%/6%/8%, noise standard deviations 0.05/0.2/0.5, middle/recent
locations, and independent Gaussian/Laplace errors of equal variance.
There are 2400 histories: 1440 null and 960 positive, including 480 exactly
at the 5% threshold. No earlier evaluation histories are reused.

Gaussian cases evaluate the stated model. Laplace cases test sensitivity to
the distributional assumption; their measured rates cannot establish the
Gaussian guarantee for Laplace noise. Report the two families separately.
Recent changes leave four later readings at length 40 and ten at length 100.

Hold the inherited shared-floor detector fixed at ``shared-r1-c2``. Reuse
the complete previous reporting pipeline and save its fit and reports.
Compare production, the existing shared-floor reporter, both earlier
references alone and gated, and the new direct test alone and gated.
Gating keeps an existing alert only when the history-level evidence test
passes; it does not certify the detector's reported positions.

Record per-history t statistics, maximum F statistics, selected extra
boundaries, surviving null splits, and each rejection route. Preserve
complete numerical inputs and outputs. Report false alerts, misses, exactly-
5% alerts, paired gains/losses, and breakdowns by family, length, location,
noise, and change size. Keep the cutoffs unchanged after evaluation.

Validation must compare interval moments with exact rational calculations,
all test statistics with independent least-squares regressions, analytical
tail probabilities with the allocated error budget, and the reused detector
pipeline with its existing implementation. Check the all-locations decision,
unit invariance, constant histories, exact-threshold examples, and short
post-change plateaus before freezing.

The next decision is whether useful sensitivity returns under the explicit
Gaussian model, and how much survives the Laplace stress test. A favorable
result would motivate separate work on unequal variance and correlation;
it would not establish a production-ready default.
