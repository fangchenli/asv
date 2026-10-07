Using more residual directions to tighten the probability bound
==============================================================

The saved explanation at correlation 0.8125 now has a rigorous probability
bound of **0.006590**, down from 0.032029. The bound stays below **0.007357**
throughout a small interval around that correlation, crossing the 0.01 cutoff.

Two changes make this possible: using eight residual directions in the
density bound, and evaluating a second exponential tilt. Neither change
spends additional false-alert budget, because every calculation bounds the
same directional p-value.

The complete search still does not recover the missed detection. It finds
another explanation at correlation **0.796875**, whose bound is **0.010687**.
All seven saved-case outcomes verify: four alerts remain alerts, and three
non-alerts remain non-alerts. The fresh-study detection count remains 58/144.

Why another approximation was needed
------------------------------------

The `full-determinant calculation <directional_determinant.rst>`_ retains
the covariance information that an earlier bound discarded. But it still
uses a probability inequality that can overestimate how often an unfavorable
residual pattern occurs. Accurate matrix arithmetic does not remove that
statistical looseness.

At the proposed split after reading 1 and correlation 13/16, the previous
bound was 0.032029, while numerical integration estimated a probability
of 0.0026718. We need a proved upper bound below 0.01 before the search may
reject the explanation. This checkpoint strengthens the inequality rather
than treating the numerical estimate as a certificate.

What the exponential tilt does
------------------------------

Keep the earlier notation: C is the covariance of the residual coordinates
after removing both plateau levels, r is the observed directional threshold,
and the a_j are the eigenvalues of C/r-I. The exact directional probability
is the Gaussian-square tail::

    p = P(W <= 0),  W = sum_j a_j * xi_j**2.

For 0 < tau < 1/2, exponential reweighting gives::

    M = det(I + 2*tau*(C/r-I))**(-1/2)
    p = M * E_tilted[exp(tau*W) * indicator(W <= 0)].

The ordinary Chernoff bound replaces the expectation by 1. The tilt changes
which part of the distribution receives the most weight; choosing a more
useful tau can reduce M. Under that tilted law, the independent coefficients
become::

    b_j = a_j / (1 + 2*tau*a_j).

The density refinement also bounds the expectation, instead of dropping it.
If the tilted distribution cannot concentrate much probability near zero,
the expectation must be appreciably below 1.

More directions give a stronger density bound
---------------------------------------------

Suppose k=4m positive tilted coefficients are each at least b0>0. The
magnitude of the tilted characteristic function is then bounded by::

    |phi(u)| <= (1 + 4*b0**2*u**2)**(-m).

Each Gaussian-square coordinate contributes a factor with exponent -1/4;
k such coordinates give exponent -m. The remaining coordinates have
characteristic-function magnitude at most 1, so ignoring them is conservative.
Integrating this bound in Fourier inversion gives a bound on the density::

    sup_w f_tilted(w) <= c_m / b0
    c_m = binomial(2*m-2, m-1) / 4**m.

These constants follow from the elementary integral recurrence for
(1+x**2)**(-m). Four, eight, twelve, and sixteen directions give constants
1/4, 1/8, 3/32, and 5/64 respectively. Thus::

    p <= M * min(1, c_m / (tau*b0)).

Using more directions lowers c_m, but the weakest included direction may
also lower b0. More is therefore not always better. The implementation
checks k in {4, 8, 12, 16}, skipping counts larger than n-2, and uses the
strongest valid correction.

Bounding the residual eigenvalues without numerical eigensolvers
-----------------------------------------------------------------

The density refinement needs a lower bound on the k-th largest residual
covariance eigenvalue. The former near-endpoint certificate supplied one
using a coarse matrix inequality that failed at moderate correlations.
Here we use the structure of the full AR(1) matrix.

Write v=|rho|. Positive and negative correlations of the same magnitude
have the same covariance eigenvalues: an alternating diagonal sign change
transforms one covariance matrix into the other. Let A_v be the tridiagonal
precision numerator, so R_v=(1-v**2)*A_v**(-1).

Replace A_v's two endpoint diagonal entries 1 by 1+v**2 to obtain a
tridiagonal Toeplitz matrix T_v. This increases the matrix in positive
semidefinite order, hence::

    A_v <= T_v
    R_v >= (1-v**2) * inverse(T_v).

T_v has known eigenvalues. Numbering its eigenvalues from smallest to
largest, the j-th is::

    (1-v)**2 + 4*v*sin(j*pi/(2*(n+1)))**2.

The residual space removes two dimensions. Eigenvalue interlacing therefore
gives the following bound, where residual and covariance eigenvalues are
numbered from largest to smallest::

    lambda_k(C) >= lambda_(k+2)(R_v)
                >= (1-v**2) /
                   ((1-v)**2 + v*(22*(k+2)/(7*(n+1)))**2).

The last step uses sin(x)<=x and pi<22/7. The result is rational and applies
to every split. It is deliberately conservative about the two removed
directions; no numerical eigenvalue estimate enters the certificate.

For a correlation interval, let v_low and v_high enclose |rho|. Bound the
numerator below by 1-v_high**2. The denominator is a convex quadratic in v,
so its maximum occurs at one of those two endpoints. This gives a uniform
lower bound L_k over the entire correlation interval, including intervals
that cross zero.

From eigenvalues to the final certificate
-----------------------------------------

The existing determinant routine supplies a radius upper bound r_upper
and a determinant lower bound D_lower for each chosen tilt. Therefore::

    a_k >= L_k / r_upper - 1 = a_lower.

When a_lower>0, the map a/(1+2*tau*a) is increasing, so the k-th largest
positive tilted coefficient is at least::

    b0 = a_lower / (1 + 2*tau*a_lower).

Combining that with the density bound gives the exact rational comparison::

    correction = max(1, tau*b0/c_m)
    p**2 <= min(1, 1/(D_lower * correction**2)).

The correction must be applied before clipping the Chernoff bound at 1.
Clipping first and then dividing could manufacture an invalid probability
bound when the original exponential moment is very large. An analytical
test explicitly covers that case.

The implementation evaluates tau=1/16 and tau=1/8. Taking the minimum of
valid bounds is safe because both bound the same p-value pointwise and
uniformly over the declared interval. It is not a union of separate tests.
The direction counts and tilts are development choices made using saved
cases, not settings evaluated on fresh histories.

What each ingredient contributes
--------------------------------

At the old witness, split 1 and rho=13/16, the point bounds are:

.. list-table:: Bounds rounded upward
   :header-rows: 1

   * - Calculation
     - Probability upper bound
   * - Previous Chernoff bound, tau=1/16
     - 0.032029
   * - Eight-direction refinement, tau=1/16
     - 0.018915
   * - Full determinant alone, tau=1/8
     - 0.015681
   * - Eight-direction refinement, tau=1/8
     - 0.006590

Neither ingredient alone crosses 0.01 here. Together they do. On the full
interval [3327/4096, 3329/4096], the combined upper bound is 0.007357.
This provides a valid interval certificate as well as a point exclusion.

The complete replay
--------------------

The separate `reporting reference <reporting_spectral.py>`_ keeps the existing
uniform and endpoint tail routes and adds the new spectral certificate.
It preserves the confidence allocation, reporting cutoffs, and work limits
of 4096 cells and depth 16. Reporting still requires every split and every
stationary correlation to be excluded; incomplete searches produce non-alerts.

The first lost detection now reaches a survivor at split 1, rho=51/64,
after 11 cells instead of 8. The new bound there is 0.010687; a floating
quadrature estimate is 0.0043470, with estimated integration error 4.4e-9
after rounding upward. The estimate is not used to make the decision.

The second lost detection retains rho=7/8. Its refined bound is 0.055928,
while the earlier numerical tail estimate remains about 0.01682. All four
previous alerts keep their complete certificates, and the 4% deterministic
example remains a non-alert. None of the seven runs hits a work limit.

As in the preceding replay, these seven completed searches do not record
a whole-interval certificate through the new route: the changed early case
finds another survivor first, and inherited routes cover the existing alerts.
The saved calculations separately certify the nonzero interval around
13/16, and tests validate that interval's proof record and reject tampering.

Verification and next work
---------------------------

The 37 mathematical tests check closed-form Gaussian-square tail cases,
the eigenvalue enclosure against projected matrices at both correlation
signs, the saved witness and interval, invariances, and conservative failure.
The 23 integration tests check coverage, surviving explanations, work limits,
calibration, unit conversion, and damaged certificates. The verifier
reconstructs residual states separately but shares the bound arithmetic;
the preceding determinant tests supply independent exact dense algebra checks.

The `summary <data/spectral_reporting_v1_summary.json>`_ and compressed
`full archive <data/spectral_reporting_v1_diagnosis.json.gz>`_ preserve the
seven results, previous results, observations, interval proofs, source hashes,
and numerical diagnoses. The `running guide <running.rst>`_ gives commands.
ASV's production detector and the frozen studies remain unchanged. As before,
the model assumes stationary Gaussian noise of common variance; the inherited
t/F reporting cutoffs are numerical.

The next question is how much looseness remains in the spectral lower bounds
versus the tilted-density inequality itself. Use the new 51/64 witness to
separate those costs before changing the method again. A bound below 0.01
at one explanation still needs a complete replay to establish an additional
detection. The middle-change loss remains a separate calibration question.
