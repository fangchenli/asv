Allowing uncertainty in AR(1) correlation
=========================================

The `covariance oracle study <covariance_report.rst>`_ reduced false alerts
when the true noise structure was supplied. This reference removes that
knowledge for one restricted model: stationary Gaussian AR(1) noise with
constant marginal variance, unknown correlation -1<rho<1, and at most one
change in the positive benchmark level. Variance changes remain outside
this model.

We keep possible correlations and change locations together. A slowdown is
reported only when every remaining pair is rejected. This avoids treating
an estimated residual correlation as exact, and avoids choosing a change
location before accounting for uncertainty in that choice.

The result below is a derivation and an experimental reference. Analytical
checks establish its calculations; no new simulation study measures its
sensitivity yet.

Why correlation and location belong together
--------------------------------------------

A proposed step explains some variation that would otherwise be residual
noise. Moving that step changes how much persistence remains in the noise.
A correlation estimate obtained after selecting a step therefore depends
on an uncertain modeling choice.

Write C for a set of pairs (t,rho) containing the true pair with probability
at least 1-delta. For every pair in C, apply the fixed-covariance size and
shape tests at that pair's own location. Require all pairs to be rejected.
This is less restrictive than first projecting C onto a set of correlations
and then restoring every possible location at each correlation.

With total budget alpha=0.05, allocate delta=0.01 to construction of C and
alpha-delta=0.04 to reporting. Retain the previous 80/20 size/shape allocation
within reporting: 0.032 for size and 0.008 for the extra-boundary scan.
These are declared allocations, not settings selected by simulation.

Constructing C without assuming independent halves
--------------------------------------------------

Use the first floor(n/2) readings as a training prefix. On the remaining m
readings, choose a predictive density q conditional on that prefix. Its
parameters and mixture weights must be set from training data alone. It
must be a proper density over the entire remaining sequence, not a density
fitted afterward to the same values it scores.

For a candidate pair (t,rho), let p be the conditional AR(1) likelihood of
the remaining observations. Maximize it over nuisance levels and innovation
variance, or use a computable upper bound U on that maximum. Define::

    E(t,rho) = q / U(t,rho)
    C = {(t,rho): E(t,rho) <= 1/delta}

At the true pair and nuisance parameters, U is at least the true conditional
density p_true. Conditional on the training prefix::

    E[E(t_true,rho_true) | training]
        <= integral (q/p_true)*p_true = integral q = 1

Markov's inequality then bounds the probability of excluding the true pair
by delta. The halves need not be independent: both densities condition on
the observed prefix. Averaging the conditional bound gives the same
unconditional coverage.

This follows the likelihood-ratio principle of
`Universal Inference <https://arxiv.org/abs/1912.11436>`_
(Wasserman, Ramdas, and Balakrishnan). The conditional AR(1) construction,
nuisance enlargement, and continuous rejection checks below are derived
here for this step-detection problem.

A simple upper bound on the null likelihood
-------------------------------------------

For validation index i, form ``z_i(rho) = y_i - rho*y_(i-1)``. In a stationary
AR(1) history, the innovations have independent common variance
``tau**2 = sigma**2*(1-rho**2)``. Within each plateau the conditional mean
of z_i is constant. At i=t, it is ``b-rho*a`` instead of either ordinary
plateau intercept.

Enlarge the conditional mean model: fit an independent intercept before
the candidate split, an independent intercept after it, and a free mean
for the single transition innovation if it is in validation. The enlarged
model contains every true conditional mean. Its maximum likelihood is
therefore an upper bound U, which makes E smaller and coverage conservative.

The transition observation is perfectly fitted in this enlarged model.
It still counts among the m density factors; it is not deleted from the
likelihood or from the estimated-variance divisor. Empty groups contribute
zero residual error. Let S_t(rho) be the sum of centered squared z_i in
the two remaining groups. Then::

    log U(t,rho) = -m/2 * [log(2*pi) + 1 + log(S_t(rho)/m)]

When S=0, U is unbounded and the pair remains in C. Each z_i is affine in
rho, so centering within fixed groups gives a quadratic polynomial::

    S_t(rho) = A_t - 2*B_t*rho + D_t*rho**2

The membership rule becomes ``S_t(rho) <= K``, where the same K applies
to every t::

    log K = log(m) - log(2*pi) - 1
            + 2/m * [log(1/delta) - log q]

Thus every location has an explicitly defined quadratic confidence region
over correlation. Large uncertainty leaves a broad region; it does not
force a single estimate.

A predictive density that allows a future step
----------------------------------------------

Fit a baseline mean, an AR coefficient, and a positive innovation variance
using training data. The coefficient is clipped to [-0.95,0.95] for the
predictor only; the confidence set still covers the full stationary range.
After scaling by the training prefix's largest absolute value, the predictor
variance has a 1e-12 floor. These choices affect efficiency, not the coverage
argument, because any proper training-chosen density q is valid.

The predictor mixes equally over no future mean change and a change at each
validation index. A change has a zero-mean Gaussian prior on its size, with
variance chosen from the training level and variation. Conditional on a
change at j, the innovation-mean shift has coefficient 1 at j and
1-rho_hat at subsequent indices. Integrating out the change size gives a
Gaussian density with covariance ``tau_hat**2*I + v_jump*v_j*v_j'``.
Here v_jump is the prior variance and v_j is that coefficient vector.

Each mixture component is a proper conditional density. The map from future
observations to conditional innovations is triangular with determinant one,
so the mixture is proper in observation coordinates too. Future values
evaluate the density but never choose its fitted parameters or prior weights.
Allowing a future step prevents one real jump from automatically ruining
the predictive likelihood and making the confidence region uninformative.

Testing a fixed pair
--------------------

For stationary AR(1), the inverse correlation matrix is
``P(rho)/(1-rho**2)``, where P is tridiagonal: endpoint diagonals are 1,
interior diagonals are 1+rho squared, and neighboring off-diagonals are -rho.
The common factor cancels from both reporting statistics. Use P directly.

For a two- or three-plateau design X, define::

    G = X' * P * X
    u = X' * P * y
    v = y' * P * y
    d = determinant(G)
    N = v*d - u'*adjugate(G)*u
    residual_sum = N/d

These are low-degree polynomials in rho. At a two-level split, with contrast
``h = [-(1+threshold), 1]``, also define::

    B = h' * adjugate(G) * u
    C = h' * adjugate(G) * h

For -1<rho<1, G is positive definite and its determinant is positive.
The size statistic exceeds its positive cutoff c_t exactly when::

    B > 0
    (n-2)*B**2 - c_t**2*N*C > 0

For an extra boundary with three-level quantities N3,d3, its F statistic
exceeds cutoff c_f exactly when::

    (n-3)*N*d3 - (n-3+c_f)*N3*d > 0

The reference abstains on exactly zero residual variation in a candidate
model. The probability model assumes positive noise variance; such fixtures
have probability zero. These polynomial conditions otherwise reproduce the
known-covariance t and F calculations.

Covering the continuum
----------------------

For every split, cover the full open interval (-1,1) with subintervals.
A subinterval is eliminated if one of these statements is proved throughout:

* Its confidence quadratic exceeds K.
* Both size-test polynomial inequalities hold.
* One extra-boundary F polynomial is positive.

The implementation represents polynomial coefficients as exact rational
numbers after reading the binary floating-point inputs and cutoffs. It
converts a polynomial to the Bernstein basis on a rational interval. Since
Bernstein basis functions are nonnegative and sum to one, positive
coefficients prove positivity throughout that interval. If a coarse interval
cannot be certified, bisect it and try again. Factors 1-rho and 1+rho can
be removed exactly because they are positive in the stationary domain.

Midpoint evaluations can find a surviving explanation. They cannot certify
that none exists between grid points. If subdivision reaches its work or
depth limit, return an unresolved non-alert. An alert requires certificates
covering every interval for every split. Store the intervals and their
rejection routes so the result is reviewable.

Strict inequalities matter. Equality at an internal point prevents that
route from certifying the whole interval. The endpoints -1 and 1 are excluded
because they are not stationary AR(1) models. Numerical density evaluation
and distribution quantiles remain floating-point calculations; exact
polynomial certificates address the continuous search, not a claim of exact
floating-point coverage for all possible inputs.

The final error bound
---------------------

Suppose the true increase is at most the reporting threshold. Either the true
pair is excluded from C, an event with probability at most 0.01, or it remains.
On the latter event, a global alert must reject the reporting null at that
fixed true pair. The known-covariance tests bound that event unconditionally
by 0.032+0.008=0.04. Adding the bounds gives at most 0.05.

This does not assume independence between confidence construction and
reporting, and does not condition the reporting error bound on selection
into C. Both use the same observations. There is no multiplicity charge
over candidate main splits or correlations because all must be eliminated.
The n-2 extra-boundary alternatives still need the existing Bonferroni charge
because any one can reject a fixed shape.

What this establishes
---------------------

The reference provides a mathematically specified way to handle one unknown
covariance parameter and a continuous nuisance search. Tests check the
predictive mixture against direct multivariate Gaussian densities, the
conditional residuals against independent least-squares fits, the reporting
polynomials against GLS, and interval certificates against exact examples.

The API is ``reporting_ar1.evidence(values, config=None, max_cells=4096,
max_depth=16)`` for 8 to 200 observations. ``calibration(n)`` supplies the
declared error budget. A modified cutoff dictionary is rejected. The result
includes all confidence quadratics as exact fraction strings, their shared
residual bound K, the training normalization, and interval certificates.
The coefficients describe S_t in normalized units; compare it with K in
those same units. Every candidate location gets a quadratic even if the
reporting search stops early.

Results distinguish four outcomes:

* ``certified_alert``: every pair has been eliminated with interval certificates.
* ``surviving_explanation``: a saved rational correlation and location remain
  inside C and fail both reporting rejection routes.
* ``unresolved``: work or depth limits prevent a complete certificate.
* ``insufficient_variation``: an exactly fitted candidate makes the ordinary
  positive-variance statistics degenerate.

Only the first outcome alerts. Partial certificates in other outcomes do
not constitute a full-domain rejection. An unrepresentable predictive density
retains the entire correlation range rather than excluding any pair.

A deterministic example alternates between 10.1 and 9.9 for six readings,
then between 10.9 and 10.7 for six. The reference certifies an alert without
choosing one correlation as the truth. With sixteen readings and only a 4%
increase halfway through, it instead retains split 8 at rho=-0.75. These
examples verify execution and certificate coverage; they are not estimates
of sensitivity or false-alert rates. See `reproduction commands <running.rst>`_.

Sensitivity, confidence-region width, numerical abstentions, and the work
needed for certification still require a frozen evaluation on fresh data.
Unknown variance changes, non-Gaussian noise, multiple mean changes, and
repeated publication are separate extensions. A plug-in estimate can be
included in that evaluation as a comparison, with its lack of this guarantee
made explicit.
