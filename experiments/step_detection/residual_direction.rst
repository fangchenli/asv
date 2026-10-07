Confidence from the pattern of residuals
========================================

The `indexed study <indexed_report.rst>`_ leaves a specific problem: a clean
8% slowdown can still be explained by correlation near 1. Most misses have
such a surviving explanation; they are not searches that ran out of time.
This note derives a confidence construction that removes the unknown plateau
levels and overall noise scale before asking about correlation.

The construction is valid and much cheaper in predictive terms on the saved
examples. It still retains the difficult positively correlated witness.
The next useful question is how conservatively we are converting directional
evidence into a confidence decision, rather than another change to level priors.

What we keep and what we remove
-------------------------------

Choose a candidate change location t. Subtract the arithmetic mean on each
side of it. The resulting residuals describe the pattern within the two
plateaus. Divide the whole residual vector by its length, using one common
divisor for both sides. This keeps the pattern and removes its amplitude.

At the correct location, changing either plateau level does not change these
residuals. Multiplying every noise deviation by the same positive amount does
not change their normalized pattern. Correlation still matters: persistent
noise produces different patterns from rapidly alternating noise.

This transformation is for correlation confidence only. Reporting still uses
the original readings to ask whether the later level exceeds 1.05 times the
earlier level. Removing the jump from that size test would remove its target.

The model and the residual coordinates
--------------------------------------

At a fixed candidate location, write the Gaussian model as::

    y = X_t beta + sigma * epsilon
    epsilon ~ N(0, R_rho),  (R_rho)_ij = rho**abs(i-j)

X_t has two indicator columns, one for each plateau. The vector beta contains
their levels, sigma is the common marginal noise standard deviation, and
-1 < rho < 1. Let n be the number of readings and d=n-2. We need n>=3 and
1<=t<n for the direction to exist; useful correlation information requires
at least d=2, hence n>=4.

Choose an n by d matrix U_t whose columns are an orthonormal basis for the
vectors with zero sum on each plateau. Thus U_t' U_t=I and U_t' X_t=0. Define::

    z = U_t' y
    C_t(rho) = U_t' R_rho U_t
    s = z / ||z||

Then z is exactly Gaussian with mean zero and covariance sigma**2*C_t(rho).
The mean parameters disappear by projection, not by treating their estimates
as if they were known. C_t is positive definite for every interior rho.
The coordinates z contain the same information as the centered residuals;
using d coordinates avoids their two linear constraints.

If z=0, the direction is undefined. A future detector should abstain rather
than manufacture directional evidence. This event has probability zero under
the nondegenerate continuous model but matters for exact or quantized data.

The exact density of the direction
----------------------------------

With respect to *uniform probability* on the unit sphere in d dimensions,
the density of s is::

    g_t,rho(s) = determinant(C_t(rho))**(-1/2)
                 * (s' inverse(C_t(rho)) s)**(-d/2)

To derive it, substitute z=r*s into the centered Gaussian density. Polar
coordinates contribute r**(d-1). Integrating the positive radius r leaves
the quadratic form to power -d/2. The powers of sigma cancel. Dividing by
the uniform sphere density cancels the remaining sphere-area constant.
At C=I, the result is g=1, fixing its normalization.

This is the angular central Gaussian distribution studied by
`Tyler (1987) <https://doi.org/10.1093/biomet/74.3.579>`_. The
`Riemann package documentation
<https://search.r-project.org/CRAN/refmans/Riemann/html/acg.html>`_ gives the
same uniform-measure convention. Multiplying C by a positive scalar leaves
g unchanged. Changing the orthonormal residual basis rotates both s and C
and also leaves g unchanged.

Choose a proper predictor for directions
----------------------------------------

The simplest predictor is uniform: q_t(s)=1. It predicts the residual
directions of independent noise and has no adjustable parameter. Its
weakness is that strongly correlated observations can be poorly predicted.

A second, illustrative predictor is a fixed mixture::

    q_t(s) = sum_j w_j * g_t,r_j(s)
    w_j >= 0,  sum_j w_j = 1

Every component is a normalized directional density, so their mixture is
proper. For the accompanying development calculation, specify five equally
weighted correlations -0.9, -0.5, 0, 0.5, and 0.9 before calculating any
directional results on the archived examples. They give a simple symmetric
control spanning weak and strong correlation. This is not a calibrated
choice or a claim that five components are sufficient.

Both predictors live on residual directions, not on complete histories.
They need no prior on beta or sigma and no external measurement unit.
Choosing the largest component after seeing s would generally not integrate
to one. Keeping the weighted sum is what makes the next argument work.

Coverage and its connection to reporting
----------------------------------------

For each candidate pair define::

    E(t,rho) = q_t(s_t) / g_t,rho(s_t)
    retain (t,rho) when E(t,rho) <= 1/delta

At the true location and correlation, the direction has density g, so::

    expectation(E(t_true,rho_true)) = integral(q_t_true(s)) = 1
    probability(E(t_true,rho_true) > 1/delta) <= delta

The second line follows from Markov's inequality. This is the density-ratio
principle also used in `Universal Inference
<https://pmc.ncbi.nlm.nih.gov/articles/PMC7382245/>`_, here applied to a
statistic with an exact nuisance-free distribution.

There is no independence requirement between confidence and reporting.
A false alert must either discard the true pair during confidence
construction, or reject it in the size/shape tests. With delta=0.01 and the
existing reporting allocation 0.04, the union bound is still 0.05. There is
no additional location penalty: an alert requires rejection of every pair,
and the proof only needs coverage at the fixed true location. For a constant
mean, any preselected location is a valid representation with equal levels.

These statements concern the same stationary Gaussian, common-variance,
at-most-one-change model as the preceding study. The construction removes
level and scale nuisance parameters; it does not remove model assumptions.

An identity that connects it to the existing solver
---------------------------------------------------

We need not build a dense residual covariance for every rho. Let::

    A_rho = (1-rho**2) * inverse(R_rho)
    Q_t(rho) = min_beta (y-X_t*beta)' A_rho (y-X_t*beta)
    S_t = Q_t(0)

A_rho is the existing tridiagonal precision numerator, and S_t is the ordinary
within-plateau sum of squares. The two-column Gram determinant is::

    D_t(rho) = determinant(X_t' A_rho X_t) = (1-rho)*B_t(rho)
    B_t(rho) = (1+rho) + (n-2)*(1-rho)
               + (t-1)*(n-t-1)*(1-rho)**3

The diagonal Gram entries are 1+(t-1)*(1-rho)**2 and
1+(n-t-1)*(1-rho)**2; the off-diagonal entries are -rho. This verifies D.
A block determinant and a Schur complement give::

    determinant(C_t) = determinant(R_rho)
                        * determinant(X_t' inverse(R_rho) X_t)
                        / determinant(X_t' X_t)
    s' inverse(C_t) s = Q_t(rho) / ((1-rho**2)*S_t)

Substituting determinant(R_rho)=(1-rho**2)**(n-1) yields::

    g_t,rho(s_t)**2
        = t*(n-t)*(1+rho)/B_t(rho) * (S_t/Q_t(rho))**(n-2)

This is the key implementation identity. The existing exact GLS polynomials
already supply Q=N/h. B is only cubic and is strictly positive on [-1,1]
for n>=3. Its derivative is -(n-3)-3*(t-1)*(n-t-1)*(1-rho)**2, so B is
nonincreasing. Powers need not be expanded into high-degree polynomials.

What happens near correlation 1
-------------------------------

Let V_t be the sum of squared successive differences *within* the plateaus,
omitting the difference across the candidate boundary. Then::

    Q_t(1) = V_t
    B_t(1) = 2
    limit(rho -> 1) g_t,rho(s_t)
        = sqrt(t*(n-t)) * (S_t/V_t)**((n-2)/2)

Q_t(1) follows by fitting the boundary jump in the difference-based limiting
precision. V_t>0 whenever S_t>0. The limit is finite and positive.

This changes a useful feature of the previous full-history rule. Its
stationary likelihood eventually vanishes as rho approaches 1 for these
nonconstant residuals, forcing confidence exclusion sufficiently close to
the endpoint. Here the common covariance shrinkage cancels when scale is
removed. A finite positive directional density remains. Near-1 exclusion
must come from the residual pattern, and is not automatic.

For intuition, larger V_t relative to S_t means more rapid variation within
the plateaus. Such a direction is less plausible under extreme persistence.
For the uniform predictor, the limiting evidence is exactly::

    E(t,1 limit) = (V_t/S_t)**((n-2)/2) / sqrt(t*(n-t))

For the mixture, multiply this by q_t(s_t). This separates the information
in the pattern from the choice of directional predictor.

What happens near correlation -1
--------------------------------

At rho=-1, R approaches v*v', where v alternates +1 and -1. For n>=3,
v is not constant on both plateaus, so U_t'v is nonzero. One residual
covariance direction remains of order one; the other d-1 shrink linearly
in 1+rho. For d>=2 and a direction not parallel to U_t'v, the density tends
to zero at rate sqrt(1+rho). A positive fixed predictor can therefore exclude
an endpoint neighborhood.

There is an exception: if the centered history is exactly proportional to
the centered alternating pattern, its density diverges at rate
(1+rho)**(-(d-1)/2). That direction supports extreme negative correlation.
It must not be certified away using the generic limit. When d=1 (n=3),
only two signed directions exist and each has probability one half for every
rho; the density relative to uniform is identically one. There is then no
correlation information in the direction.

For the generic case, the same conclusion follows from B_t(-1)>0 and
Q_t(-1)>0 in the scalar identity. In the alternating exception Q_t(-1)=0,
so substituting -1 into the displayed ratio would incorrectly create 0/0.

How a future interval certificate could work
--------------------------------------------

Let q_lower be a nonnegative lower bound on the fixed directional predictor
at this location. Confidence exclusion follows from::

    (delta*q_lower)**2 * B_t(rho) * N(rho)**d
        > t*(n-t)*(1+rho) * S_t**d * h(rho)**d

For an interval [left,right], use exact Bernstein bounds L_N for N from below
and H_h for h from above. A sufficient check is::

    L_N > 0, H_h > 0
    (delta*q_lower)**2 * B_t(right) * L_N**d
        > t*(n-t)*(1+right) * S_t**d * H_h**d

The same subdivision engine can combine this with the existing reporting
certificates. It must still cover every candidate location and all interior
correlations; testing only a saved witness or a grid cannot certify an alert.

There is an additional numerical advantage. For rational observations and
rational component correlations, each g_t,r_j**2 is rational. An integer
square root supplies a rigorous rational lower bound on each g_t,r_j, hence
on their weighted mixture. The confidence calculation can avoid floating
logarithms and exponentials entirely. This does not certify the numerical
quantiles used by the separate reporting tests.

What the archived witnesses show
--------------------------------

Keep delta=0.01, so exclusion needs log evidence above log(100)=4.605.
The five-component predictor was specified above before these calculations.
The first case is the earlier independent-noise 8% slowdown at location 25;
the second is the positive-correlation case identified in the latest study.
Both have 100 readings and noise standard deviation 0.05.

.. list-table:: Confidence evidence at existing witnesses
   :header-rows: 1

   * - Archived case
     - Candidate correlation
     - Indexed-jump log evidence
     - Direction, uniform predictor
     - Direction, fixed mixture
   * - Independent noise, seed 1400
     - 4095/4096
     - 9.545
     - 22.859
     - 21.259
   * - Correlation 0.7, seed 1500
     - 16383/16384
     - -45.458
     - -32.066
     - 1.603

The independent case's old witness is excluded by both directional predictors.
Exact rational bounds also exclude the entire interval from 4095/4096 to 1
at that location. Indexed-jump already excludes this witness, so this is
stronger confidence evidence, not a new recovered detection.

The correlated witness survives. The uniform predictor describes its pattern
poorly; the mixture raises log evidence to 1.603 but still falls short of
4.605. That is an evidence ratio of about 5, where this rule requires more
than 100. Rigorous upper and lower predictor bounds establish these decisions,
so this is not ambiguity from rounding. At the positive endpoint, mixture
log evidence tends to 1.604: moving even closer to 1 does not resolve it.
Both predictors retain the generating correlation at the generating location
on both examples.

These are checks of already inspected histories, not estimates of detection
power. No new observations or complete reporting searches were run. In
particular, eliminating one witness does not rule out other locations and
correlations.

The next mathematical lever: calibrate the direction itself
-----------------------------------------------------------

The density-ratio construction uses only the mean-one bound to set its
cutoff. We now know the entire distribution of the direction at a candidate
rho, without unknown levels or noise scale. That opens another route:
calibrate a directional statistic under that distribution directly.

For example, retain the ordering of the uniform-predictor statistic 1/g.
Write z=sigma*C**(1/2)*xi, where xi has independent standard normal
coordinates. Then::

    s' inverse(C) s = (xi'xi) / (xi'C xi)
    R_observed = 1 / (s' inverse(C) s)
    R_null = sum(lambda_j * xi_j**2) / sum(xi_j**2)

The lambda_j are the eigenvalues of C at the candidate location and
correlation. Larger 1/g corresponds to smaller R_observed. Its exact lower
tail is a weighted chi-square probability::

    p(t,rho) = P(sum((lambda_j-R_observed)*xi_j**2) <= 0)

Rejecting the candidate when this tail probability is at most delta gives
the same true-pair coverage argument, but uses a candidate-specific tail
instead of the common Markov cutoff. When C is a scalar multiple of I,
the statistic is constant and p=1; that candidate is retained. In the
nonconstant case the distribution is continuous. A common positive scaling
of C scales both lambda_j and R_observed, leaving the probability unchanged.

For the same uniform-predictor ordering, this exact tail test contains the
Markov rejection region: whenever E=1/g>1/delta, Markov's inequality gives
P(E_null>=E_observed)<=1/E_observed<delta. This is a rejection-region
comparison with the uniform directional rule, not a dominance claim over
the mixture or earlier full-history rules.

The next task is to derive reliable bounds for this weighted chi-square
tail and determine how to certify them across a correlation interval.
Ordinary floating quadrature at a few rho values would not suffice for the
existing continuous reporting guarantee. Keep the fixed mixture as a
control; do not adjust its components to this correlated example. A new
detector and a fresh frozen evaluation should follow that analysis.

The follow-up `directional tail certificate <directional_tail.rst>`_ now
provides a rational upper bound below 0.008194 across the saved correlated
witness's interval up to 1. It combines exponential tilting, a density bound,
and matrix enclosures; a complete reporting search remains the next step.

Files and verification
----------------------

`residual_direction.py <residual_direction.py>`_ is an algebra and diagnostic
helper, not a reporting detector. It evaluates the projected density through
an independent dense covariance calculation and through the exact rational
GLS identity, supplies the square-root and interval bounds, and reads two
existing examples with archive hash checks. The `saved diagnosis
<data/residual_direction_v1_diagnosis.json>`_ records both calculations,
predictor bounds, endpoint values, and source hashes. The `running instructions
<running.rst>`_ reproduce it.

All 32 analytical tests pass. They cover Gaussian radial integration,
directional density normalization, the determinant identity, level/scale/reversal invariance,
both endpoint behaviors including the alternating exception, n=3 degeneracy,
rational square-root enclosures, and a nontrivial interval certificate.
No frozen study sources, settings, observations, or decisions are changed.
