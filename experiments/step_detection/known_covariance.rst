Testing steps when the noise structure is known
===============================================

The `stress study <direct_stress_report.rst>`_ found two failures: changing
variance makes a short noisy plateau look too precise, and correlation makes
noise excursions look like extra steps. We can correct both mathematically
if we know the relative variances and correlations of the observations.
The overall amount of noise may remain unknown.

This reference establishes the test and its error bound under that assumption.
It does not yet estimate the noise structure from a benchmark history or
measure sensitivity on fresh simulations. Those are the next two questions.

What information the model needs
--------------------------------

Write the observations as ``y = X_t * beta + error`` at a proposed split t.
The two columns of X_t indicate the earlier and later plateaus, so
``beta = [a, b]`` contains their timing levels. Assume::

    error ~ Normal(0, sigma**2 * R)

R is a supplied symmetric positive-definite n-by-n matrix. Its diagonal
describes relative noise variances; its off-diagonal entries describe how
errors move together. Sigma is an unknown positive scale. For example,
``R = diag(1, ..., 1, 9, ..., 9)`` says the later readings have three times
the noise amplitude. It does not specify whether the earlier standard
deviation is 0.05 or 0.2.

R must describe the actual noise and be fixed independently of these observed
errors. The same R is used at every proposed mean split and every extra
boundary. Choosing a different covariance to improve each fit would be a
different procedure. In particular, knowing a variance-change location does
not imply that a mean change occurs there; we still test every mean split.

As before, true timing levels are positive and have at most one change.
For a reporting threshold theta, set k=1+theta. The target is::

    D = b - k*a > 0

Gaussian means and medians agree. This makes D compatible with the population
target of ASV's median fit, although this evidence calculation fits means.
The supplied matrix describes noise in retained-observation order, including
correlation across the proposed boundary.

Why the fit must change with the uncertainty calculation
--------------------------------------------------------

With independent equal-variance noise, squared Euclidean distance is a
suitable residual measure. With covariance R, use::

    Q_t = (y - X_t*beta_hat)' * inverse(R) * (y - X_t*beta_hat)

The inverse covariance discounts directions where noise naturally varies
more. A long correlated excursion therefore gets a different weight from
the same number of independent deviations. Minimizing Q_t gives generalized
least squares (GLS)::

    G_t = inverse(X_t' * inverse(R) * X_t)
    beta_hat = G_t * X_t' * inverse(R) * y

This changes both the estimated levels and their uncertainty. Ordinary
plateau means still work for diagonal R that is constant within each plateau.
For general correlation or varying diagonal entries, they need not be the
GLS estimates. Fitting each segment independently also omits cross-boundary
covariances.

An equivalent construction makes the probability argument easier. Factor
``R = L*L'`` and solve with L to obtain::

    W = inverse(L)
    z = W*y
    B_t = W*X_t
    W*R*W' = I

The transformed errors are independent Gaussian variables with common
variance sigma squared. Ordinary least squares of z on B_t is exactly GLS
in the original coordinates. This is called whitening. Transform both
the observations and the design: B_t's columns are generally not ordinary
plateau indicators.

For constant-amplitude AR(1) noise with correlation rho, the transformation
leaves the first error unchanged and replaces subsequent errors by
``(error_i - rho*error_(i-1))/sqrt(1-rho**2)``. Applied to a plateau indicator,
it also creates a special entry at a boundary. Fitting ordinary constant
segments to transformed observations would therefore test the wrong model.

The implementation uses Cholesky solves and QR least squares rather than
explicit inverses. SciPy documents the
`Cholesky factorization <https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cholesky.html>`_
and `triangular solves <https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_triangular.html>`_.
The inverses above specify the mathematical quantities.

The corrected size test
-----------------------

Let ``h = [-k, 1]``. At a fixed split whose two-level model is correct::

    D_hat = h' * beta_hat
    Var(D_hat) = sigma**2 * h' * G_t * h
    sigma_hat**2 = Q_t / (n-2)

    T_t = D_hat / sqrt(sigma_hat**2 * h' * G_t * h)

Here is why the Student t distribution still applies. In whitened coordinates
the coefficient error is a projection onto the two-dimensional column space
of B_t. The residual is a projection onto its orthogonal complement. Gaussian
projections onto orthogonal spaces are independent: their cross-covariance
is zero, and their joint distribution is Gaussian. Consequently::

    (D_hat - D) / (sigma * sqrt(h' * G_t * h)) ~ Normal(0, 1)
    Q_t / sigma**2 ~ chi_squared(n-2)

The two variables are independent. At D=0 their ratio is Student t with
n-2 degrees of freedom. For D<0, the numerator has a negative shift; for
each positive denominator this can only reduce upper-tail rejection.
Reject the size null at the existing upper-tail cutoff with probability
budget ``alpha_size = 0.04``.

For the 36/4 example with relative variances 1 and 9, the variance factor
``h' * G_t * h`` at the true split is::

    1.05**2 / 36 + 9 / 4 = 2.280625

The standardized residual sum Q_t has expectation ``38*sigma**2``. Thus
the expected estimated contrast variance is exactly its true variance.
The old pooled calculation understated that variance by a factor of about
4.981. The new formula accounts for the noisy later plateau explicitly.

The corrected shape test
------------------------

For every extra boundary j distinct from t, build the three-plateau design
X_tj and transform it by the same W. Fit it and call its weighted residual
sum Q_tj. The three-level space contains the original two-level space.

Let P_t and P_tj be their orthogonal projectors in whitened coordinates.
Because the spaces are nested, ``P_tj - P_t`` is an orthogonal projector
of rank one. ``I - P_tj`` has rank n-3, and the two projectors multiply to
zero. Under a true two-level mean at t, both remove that mean. Therefore::

    (Q_t - Q_tj) / sigma**2 ~ chi_squared(1)
    Q_tj / sigma**2 ~ chi_squared(n-3)

These quadratic forms are independent, so::

    F_tj = (Q_t - Q_tj) / (Q_tj/(n-3)) ~ F(1, n-3)

There are n-2 extra boundaries. Keep the existing Bonferroni cutoff at upper
tail probability ``alpha_shape/(n-2)``, with ``alpha_shape = 0.01``.
The probability that any extra boundary rejects a correct two-plateau shape
is at most 0.01. The F statistics need not be independent of each other or
of the size statistic. The largest F corresponds to the smallest Q_tj.

Why the complete rule still has a 5% false-alert bound
------------------------------------------------------

At a true below-threshold split t, the probability of rejection by either
route is at most ``alpha_size + alpha_shape = 0.05``. The full rule alerts
only if every proposed split is rejected by at least one route. A global
false alert must therefore reject this particular true split, giving the
same bound. There is no additional multiplicity factor over the main split:
all main splits must be rejected, whereas any extra boundary may reject a
fixed shape. That difference explains the different handling of the scans.

A constant positive history has a=b and satisfies the null at every split.
All splits from 1 to n-1 remain eligible, including one-reading plateaus.
We need n>=4 for the three-level residual degrees of freedom. Requiring
an existing ASV alert as an additional gate can only lower the probability
of a false alert.

The guarantee is per history under the stated Gaussian one-change model.
It does not certify the reported revision or individual steps in a
multiple-change history. A surviving split is an explanation not rejected
by these tests, not a boundary confidence interval.

Reduction to the existing test
------------------------------

When R=I, GLS levels become the ordinary plateau means, Q_t becomes their
ordinary squared-residual sum, and::

    h' * G_t * h = k**2/t + 1/(n-t)

The size and shape statistics are exactly the previous direct test in real
arithmetic. Multiplying R by any positive scalar also leaves both statistics
unchanged: the fitted levels stay the same, Q is divided by that scalar,
and G is multiplied by it. Changing timing units multiplies D_hat by the
unit factor and Q by its square, again preserving both tests for positive
unit factors. These identities give direct implementation checks.

Reference implementation and analytical checks
----------------------------------------------

``reporting_covariance.evidence(values, covariance, calibration)`` accepts
4 to 200 finite observations, one supplied R, and the existing direct-test
calibration. Outputs retain the old split-indexed statistics and rejection
lists, with an additional numerical status. This helper fits levels for
evidence only; ASV's median fit and reporting pipeline are unchanged.
It takes R as input and does not estimate it.

The tests use deterministic fixtures and matrix identities, not a new
statistical evaluation. They compare against independently computed GLS,
verify the Gaussian projection ranks and orthogonality behind the proof,
check the diagonal-variance example, and recover the old statistics at R=I.
They also check timing-unit and covariance-scale invariance and invalid
covariance handling.

Numerical inputs are normalized to avoid overflow. After dividing R by its
largest absolute entry, asymmetry above 64 machine epsilons is rejected;
smaller differences are symmetrized. Non-positive-definite covariance and
condition numbers above 10**12 are rejected. A residual norm at most
``64*machine_epsilon*n*norm(z)`` causes ``status='insufficient_precision'``
and ``has_alert=False``, with diagnostic arrays set to None. This is an
abstention, not a finding that all splits are plausible. Exactly noiseless
fixtures may differ from the old helper's limiting-ratio convention; the
probability model assumes sigma>0. The real-arithmetic proof does not claim
exact floating-point coverage at every extreme input or cutoff equality.

What remains before using this on real histories
------------------------------------------------

First, freeze the reference and evaluate fresh data with the actual generating
R supplied to it. This measures sensitivity when covariance is known. If R
contains a variance change at the true mean boundary, its known structure
itself supplies information about the noise transition; label that an oracle
comparison and keep the distinction visible. Use a single matrix per history
throughout the location scan.

Second, decide how R is learned. Estimating R from these residuals makes
whitening data-dependent. The projection and independence argument above
then no longer proves Student t and F calibration. Independent external
estimates still have estimation error unless the model conditions on them
with the correct conditional covariance. Their independence alone does
not make a plug-in matrix exact.

One possible mathematical extension is a confidence set for covariance
shapes. If it covers the true shape with probability at least 1-delta,
require an alert for every shape in the set and use test budget alpha-delta.
On coverage, a global alert implies rejection by the test at the fixed true
shape; outside coverage, the error probability is at most delta. A union
bound then gives at most alpha without requiring independence between the
set and the tests. Constructing a valid set and checking every matrix in
it are substantial remaining problems; a finite grid without a coverage
argument is insufficient.

Unknown covariance, non-Gaussian tails, multiple changes, and repeated
publication remain outside this implementation. This derivation supplies
a reference against which those extensions can be assessed.
