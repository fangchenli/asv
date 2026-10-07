Certifying the probability of a residual pattern
================================================

The saved positively correlated example now has a rigorous confidence
exclusion where the previous directional mixture failed. At candidate
location 25, the tail probability is below **0.008194 throughout
[16383/16384, 1)**. This is below the 0.01 confidence budget. The old witness
and every more persistent explanation at that location are excluded.

This is an interval certificate on an archived example. Other locations and
correlations still need checking before the history can produce an alert.
No fresh histories or detection-rate comparison were run.

Why calibrating the pattern helps
---------------------------------

The `residual-direction derivation <residual_direction.rst>`_ removes the
two plateau levels and overall noise scale. The remaining direction s has
a known distribution for each candidate location and correlation. The first
rule compared a predictive density with that distribution and used a common
mean-based cutoff. Its fixed mixture supplied evidence of about 5 on the
correlated witness, while the cutoff required more than 100.

Here we ask directly how often the candidate model produces a direction at
least this unfavorable to it. A numerical calculation puts that probability
near 0.00253 at the old witness. To use it in the continuous reporting search,
we need a rigorous upper bound over an interval, not just a small numerical
answer at one correlation.

The probability being bounded
-----------------------------

Keep the notation from the preceding derivation: U is an orthonormal residual
basis, C=U'R_rho U, and s=U'y/||U'y||. Define::

    r = 1 / (s' inverse(C) s)
    A = C/r - I
    W = sum(a_j * xi_j**2),  xi_j independent N(0,1)
    p(t,rho) = P(W <= 0)

The a_j are the eigenvalues of A. This is the exact lower tail of the
directional Rayleigh quotient. It uses the ordering of the *uniform*
directional predictor; it is not the tail distribution of the five-component
mixture statistic. All a_j>-1. Positive rescaling of C also rescales r,
leaving A and p unchanged.

At the true pair, this tail probability is a valid p-value. If the quotient
is nonconstant, its distribution is continuous. If C is proportional to I,
the quotient is constant and we set p=1. Rejecting only when a proven upper
bound is below delta=0.01 therefore preserves true-pair coverage. Adding the
existing 0.04 reporting budget keeps the per-history false-alert bound 0.05
under the declared stationary Gaussian, common-variance, one-change model.

A first bound: exponential tilting
----------------------------------

Choose a number tau with 0<tau<1/2. The upper restriction ensures
1+2*tau*a_j>0, since every a_j>-1. The Gaussian integral gives::

    M = expectation(exp(-tau*W))
      = product(1+2*tau*a_j)**(-1/2)
      = determinant(I+2*tau*A)**(-1/2)

On W<=0, the exponential is at least one, so p<=M. This is a Chernoff bound.
For the saved interval, the implementation proves M<=0.013289. That alone
is not small enough for the 0.01 budget.

Tilting also gives a useful identity. Define a new probability law by
weighting outcomes by exp(-tau*W)/M. Under it, the Gaussian coordinates have
variances 1/(1+2*tau*a_j), and hence::

    W under the tilted law = sum(b_j * eta_j**2)
    b_j = a_j / (1+2*tau*a_j),  eta_j independent N(0,1)
    p = M * expectation_tilted(exp(tau*W) * indicator(W<=0))

The ordinary Chernoff bound replaces that last expectation by one. We can
bound it more tightly by controlling the tilted density near zero.

Four covariance directions provide the missing factor
-----------------------------------------------------

Suppose at least four tilted coefficients b_j are at least b0>0. Each
Gaussian-square term contributes characteristic-function magnitude
(1+4*b_j**2*u**2)**(-1/4). Four such terms give::

    |characteristic_function_tilted(u)| <= (1+4*b0**2*u**2)**(-1)

Other independent terms have characteristic-function magnitude at most one.
Fourier inversion therefore bounds the whole tilted density f by::

    sup_w f(w) <= integral((1+4*b0**2*u**2)**(-1), u=-infinity..infinity)/(2*pi)
               = 1/(4*b0)

Consequently::

    expectation_tilted(exp(tau*W) * indicator(W<=0))
        <= integral(exp(tau*w)/(4*b0), w=-infinity..0)
        = 1/(4*tau*b0)
    p <= M * min(1, 1/(4*tau*b0))

This additional factor is enough for the saved interval. The construction
needs four suitable residual covariance directions, not four readings. If
there are fewer than four, or their coefficients cannot be bounded away
from zero, the code falls back to the ordinary Chernoff bound.

Why the tilt may be chosen after seeing the history
---------------------------------------------------

For a fixed observed threshold r, each admissible tau bounds the same tail
probability p. Taking the smallest of these upper bounds is valid. This is
different from maximizing an uncalibrated density ratio: we are bounding a
known probability function pointwise. Adaptive choices of tau, interval
width, and proof precision need no extra error budget.

The saved calculation uses the simple rational value tau=1/16 for every
interval. It is a development choice, not a tuned statistical prior. The
certificate remains valid for any admissible tilt that passes its checks.

Uniform bounds near correlation 1
---------------------------------

We can certify an interval without numerically tracking its eigenvalues.
Use first differences within each plateau as a nonorthonormal residual
basis B. There are d=n-2 columns, each with entries +1 and -1 at adjacent
readings, with the difference across the change omitted. Its Gram matrix::

    G = B'B

has two tridiagonal blocks, of sizes t-1 and n-t-1, with diagonal 2 and
off-diagonal -1. Its determinant is t*(n-t).

Normalize the covariance by 1-rho and put::

    K(rho) = B'R_rho B / (1-rho)
    K_ii = 2
    K_ij = -(1-rho)*rho**(abs(edge_i-edge_j)-1),  i != j

At rho=1 this has the continuous limit K=2I. If rho is in [l,u] within
[0,1], each off-diagonal magnitude is at most 1-l. A row-sum bound gives::

    K(rho) >= kappa*I,  kappa = 2-(n-3)*(1-l)

The matrix inequality means the difference is positive semidefinite. If
kappa<=0, this particular enclosure is unhelpful and returns no certificate.
It is intentionally a bound for a neighborhood of +1; it supplies no
certificate on negative correlations.

Choosing U=B*G**(-1/2) gives an orthonormal residual basis. In those
coordinates the normalized covariance is G**(-1/2)*K*G**(-1/2), so it is bounded below by
kappa*inverse(G). The observed quotient after the same normalization is::

    r_scaled(rho) = (1+rho)*S_t/Q_t(rho)

The preceding derivation already supplies Q=N/h and S_t. Exact Bernstein
bounds on N and h give, for a positive lower bound N_low::

    r_scaled(rho) <= (1+u)*S_t*h_high/N_low <= r_high

We round the right side *up* to a rational with denominator 2**32. This
keeps subsequent arithmetic smaller while preserving the upper bound.

A determinant lower bound from two short recurrences
----------------------------------------------------

The matrix used by the Chernoff bound is uniformly bounded below by::

    I+2*tau*A >= a*I + c*inverse(G)
    a = 1-2*tau > 0
    c = 2*tau*kappa/r_high > 0

Determinants preserve this ordering for positive-definite matrices, so::

    determinant(I+2*tau*A) >= D_low
    D_low = determinant(a*G+c*I) / determinant(G)

For each tridiagonal block, its size-m determinant is obtained from::

    H_0 = 1
    H_1 = 2*a+c
    H_m = (2*a+c)*H_(m-1) - a**2*H_(m-2)
    D_low = H_(t-1)*H_(n-t-1) / (t*(n-t))

Every quantity is rational. No computed eigenvalue or floating determinant
is part of this certificate.

Bounding the four required coefficients
---------------------------------------

For a plateau of length L, the corresponding Gram eigenvalues are
2-2*cos(j*pi/L), for j=1,...,L-1. The inequalities
2-2*cos(x)<=x**2 and pi<22/7 give rational upper bounds
(22*j/(7*L))**2. Let mu be the fourth smallest of these bounds across both
plateaus. At least four distinct eigenvalues of G are at most mu. Thus at
least four eigenvalues of A are at least::

    a0 = kappa/(r_high*mu) - 1

When a0>0, the increasing transformation a/(1+2*tau*a) gives::

    b0 = a0/(1+2*tau*a0)
    F = max(1, 4*tau*b0)
    p <= 1/(sqrt(D_low)*F)

If a0<=0 or d<4, use F=1. The final exclusion check needs no square root::

    D_low * F**2 > 1/delta**2

This strict rational inequality certifies every interior correlation in the
interval. The endpoint 1 itself is only a limit; no stationary model at
exactly 1 is added to the parameter space. Failed bounds remain inconclusive.

Results on the saved examples
-----------------------------

The independent example is the older 100-reading, 8% slowdown at location 25.
The correlated example has the same design with generating correlation 0.7
and is the remaining miss named in the indexed-study report.

.. list-table:: Proven upper bounds over complete intervals at location 25
   :header-rows: 1

   * - Example
     - Correlation interval
     - Probability upper bound, rounded upward
     - Below 0.01?
   * - Independent, seed 1400
     - [0.99, 1)
     - 0.005190
     - Yes
   * - Independent, seed 1400
     - [4095/4096, 1)
     - 0.000000472
     - Yes
   * - Correlation 0.7, seed 1500
     - [0.99, 1)
     - 1
     - No certificate
   * - Correlation 0.7, seed 1500
     - [0.999, 1)
     - 0.014857
     - No certificate
   * - Correlation 0.7, seed 1500
     - [16383/16384, 1)
     - 0.008194
     - Yes

For the last row, the ordinary Chernoff bound is 0.013289; the tilted-density
factor supplies the improvement needed for exclusion. Wider intervals lose
precision because the common matrix and residual bounds become looser.
A failed interval certificate does not mean that every point in it survives.

As a separate numerical check, characteristic-function inversion estimates
the correlated example's tail probabilities as 0.8012 at its generating
correlation 0.7, 0.04203 at 0.9, 0.003253 at 0.99, and 0.002531 at the old
witness. The inversion follows the quadratic-form formulation documented by
`qfratio <https://stat.ethz.ch/CRAN/web/packages/qfratio/vignettes/qfratio_distr.html>`_
and originating in `Imhof (1961) <https://doi.org/10.1093/biomet/48.3-4.419>`_.
These floating calculations and their quadrature error estimates are
diagnostics; the exclusion proof uses only the rational bounds above.
Very small numerical probabilities on the independent example suffer
cancellation and are not treated as accurate tail estimates.

How this fits a future detector
-------------------------------

The existing reporting search still needs to reject every location and
every correlation. This result supplies one continuous confidence route,
not a complete alert certificate. Next, integrate the route into a separate
reference and reconstruct its full certificates on the archived examples.
Retain unresolved regions as non-alerts and preserve the existing reporting
allocation and work limits.

There is a compatible fallback for the same p-value: the uniform directional
density g is also an upper bound on p by Markov's inequality. Either the
uniform-direction interval certificate or the new tail certificate may
therefore reject a region at the same delta without splitting the budget.
Taking the minimum of two bounds on the same p is safe. The numerical
bound implemented here does not automatically inherit all rejections of
the ideal exact-tail test, so that fallback is useful.

The five-component directional mixture and the previous full-history rules
use different statistics. Simply accepting rejection by any of those rules
at full delta would not have the same coverage proof; combining them needs
its own error-budget argument. Keep them as comparison methods for now.
Measured detection gains still require a new frozen evaluation on fresh
histories after implementation and analytical checks.

Files and verification
----------------------

`directional_tail.py <directional_tail.py>`_ implements the exact interval
bound and a separate floating tail diagnostic. The `saved results
<data/directional_tail_v1_diagnosis.json>`_ contain rational certificates,
source/archive hashes, and numerical checks for both archived histories.
See `running instructions <running.rst>`_ for reproduction commands.

The 28 new analytical tests check the difference covariance, Gram spectrum,
determinant recurrence, and tail bounds against closed-form F probabilities.
They also check outward rounding, the saved interval, scale/level/reversal
invariance, degenerate cases, and independent dense calculations. Together
with the 32 preceding directional tests, all 60 pass. The previous frozen
evaluation and its decisions remain intact.
