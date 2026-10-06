Stress the direct test's variance and independence assumptions
==============================================================

The `direct threshold test <direct_protocol.rst>`_ restored useful sensitivity
under independent Gaussian noise with one common variance. This study keeps
that test and the shared-floor detector fixed, then separates violations of
the common-variance and independence assumptions. Derive the expected failure
mechanisms first, freeze the study, and generate fresh paired histories.

What unequal variance changes
-----------------------------

At a specified split with A earlier and B later readings, the estimated
threshold excess is ``D_hat = mean_after - 1.05*mean_before``. With independent
Gaussian noise of variances va and vb in the two plateaus, its actual variance
is::

    V = 1.05**2 * va/A + vb/B

The frozen test pools the two residual sums to estimate one common variance.
Its reported variance for D_hat has expectation::

    V_reported_mean = [(A-1)*va + (B-1)*vb]/(A+B-2)
                      * [1.05**2/A + 1/B]

These match when va=vb. If only four of forty readings occur after the change,
the long earlier plateau dominates the pooled noise estimate. Tripling the
later noise amplitude sets vb=9*va. The ratio ``V/V_reported_mean`` is then
about 4.981: the contrast varies much more than its reported variance suggests.
Reversing the variance change gives a ratio about 0.224, which suggests a
loss of sensitivity instead. With equal plateau lengths the ratio is close
to one, though the pooled t distribution is still not exact.

This ratio compares a true variance with an expected estimate. It is not
the variance of the random t statistic and does not directly predict a
false-alert probability. Even for independent Gaussian observations, unequal
plateau variances make the pooled residual sum a weighted sum of chi-squared
variables, rather than the common-scale chi-squared variable used by the
frozen calibration.

What correlation changes
------------------------

Write the contrast as ``D_hat = sum(c_i*y_i)``, with ``c_i = -1.05/A`` before
the split and ``c_i = 1/B`` after it. For a covariance matrix Sigma::

    V = c' * Sigma * c

Let P project observations onto the two constant plateaus, and M=I-P. At
the true split the mean lies in that two-level model, so its squared-residual
sum Q satisfies::

    E(Q) = trace(M * Sigma)

    V_reported_mean = trace(M * Sigma)/(n-2) * (c' * c)

The first equality follows by expanding the quadratic form of the centered
noise and summing its covariances. These formulas include covariance across
the boundary. Ignoring that covariance would omit part of the uncertainty.

Positive correlation allows noise excursions to persist. It changes both
the variance of the contrast and the residual variance estimate. It can
also destroy independence between the contrast and residuals. The frozen
Student t calculation therefore needs more than a replacement sample count.

The extra-boundary F scan has its own assumption problem. With general
Sigma, the fitted improvement and remaining residual are Gaussian quadratic
forms with different scales and potentially dependent projections. Their
ratio need not follow the calibrated F distribution. Fixing only the size
standard error would leave this rejection route uncorrected.

A declared covariance model for the stress cases
------------------------------------------------

Use a stationary Gaussian AR(1) process z with marginal variance one and
``Cov(z_i,z_j) = rho**abs(i-j)``. Start z_0 as N(0,1), then use
``z_i = rho*z_(i-1) + sqrt(1-rho**2)*innovation_i`` with independent standard
Gaussian innovations. The stationary initial draw avoids a startup transient.

Observed noise is ``sigma*s_i*z_i``, where s_i is one or three. Therefore::

    Sigma_ij = sigma**2 * s_i * s_j * rho**abs(i-j)

Multiply noise deviations, not the underlying timing level. A variance
change does not count as a true slowdown. The amplitude changes at the
nominal mean-change position; this deliberately includes a noisy short
later plateau. For a flat history that position marks only a possible
noise-amplitude change.

The following ratios are deterministic calculations from that covariance,
before generating any evaluation data. Sigma's absolute scale cancels.

.. list-table:: True contrast variance / mean reported variance, length 40
   :header-rows: 1

   * - Noise condition
     - Middle split, 20/20
     - Recent split, 36/4
   * - Independent, constant amplitude
     - 1.000
     - 1.000
   * - Independent, amplitude rises threefold
     - 0.961
     - 4.981
   * - Independent, amplitude falls threefold
     - 1.039
     - 0.224
   * - Correlation 0.7, constant amplitude
     - 5.661
     - 3.177
   * - Correlation 0.7, amplitude rises
     - 5.616
     - 19.188
   * - Correlation 0.7, amplitude falls
     - 6.095
     - 0.906

At length 100 with a 90/10 split, the rise-plus-correlation ratio reaches
about 23.903. The study saves all 24 length/location/condition calculations
in its freeze. They identify mechanisms worth checking, without changing
the reporting rule or claiming corrected error control.

Frozen experimental design
--------------------------

Use six conditions: rho=0 or 0.7 crossed with constant, rising, or falling
noise amplitude. Reuse the same independent Gaussian initial draw and
innovation sequence across the six versions of each base history. Pairing
makes differences attributable to those transformations instead of a new
random draw. Within a condition, distinct base histories have separate
SHA-256-seeded random streams.

Base histories form the Cartesian product of:

* Lengths 40 and 100.
* True increases 0%, 4%, 5%, 6%, and 8% from baseline 10.
* Base noise standard deviations 0.05 and 0.2.
* Middle and recent positions, the latter leaving max(4,n/10) later readings.
* Fresh seed labels 1000 through 1009.

There are 400 base histories and 2400 condition versions. Each condition has
240 null histories, including 80 exactly at 5%, and 160 positive histories.
Totals are 1440 null and 960 positive versions, but versions sharing a base
history are dependent. Aggregate counts are descriptive, not independent
binomial trials. The largest base noise from the preceding study is omitted
to keep this study focused; multiplying 0.2 by three already gives a 0.6
standard deviation in the noisier plateau.

The pair identifier includes length, true change, noise, location, and seed,
with a study-specific prefix. It excludes the noise condition. This avoids
duplicating flat histories under two position labels while preserving
pairing across conditions. Save numerical observations as well as seeds.

Reuse the entire frozen direct-study pipeline: production, existing
shared-floor reporting, both sign references alone and gated, and the direct
test alone and gated. Inherit its critical values and ``shared-r1-c2`` fit
settings without recalibration. Freeze source/protocol/build hashes,
the inherited configuration, design, and covariance diagnostics in
``data/direct_stress_v1_frozen.json`` and commit before generating histories.

The direct test's Gaussian common-variance guarantee applies to the
independent constant-amplitude control only. The sign references retain
their independent-observation guarantee under the two independent
variance-change conditions, because a common median within a plateau does
not require common variance. Positive correlation removes that independence
guarantee for both families of tests.

Measurements and interpretation
-------------------------------

Report false alerts, misses, exactly-5% alerts, and breakdowns by condition,
length, location, noise, and true change. Save paired gains/losses relative
to the independent constant-amplitude control for every method. Gated
decisions must remain a subset of existing shared-floor alerts.

For each null history alerted by the direct test, inspect its rejection at
the generating split. Record size-only, shape-only, or both. The generating
split is known to the evaluator and is not supplied to the reporting test.
For a flat history, use the nominal split as a designated witness. A global
direct rejection must reject this witness; a ``neither`` route on any such
history indicates an implementation error. These route counts describe
which safeguard fails at that witness, not a unique causal attribution
among all searched locations.

Validate covariance formulas against an independent Gaussian innovation
matrix, variance ratios against the independent unequal-variance formula,
and noise transformations against the paired AR recurrence. Check that
evaluation reuses all frozen fitting and reporting outputs without changes.
Do not generate evaluation histories in pre-freeze tests.

After evaluation, report failures without adjusting cutoffs. A ratio-based
rescaling alone is not an established repair: both size and shape tests
need a valid covariance treatment, including uncertainty in any estimated
variance or correlation. The results should determine which mathematical
extension to pursue next. Production behavior stays unchanged.
