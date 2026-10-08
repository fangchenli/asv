Where the remaining eigenvalue approximation loses information
==============================================================

The saved explanation at correlation 0.796875 now has a rigorous probability
bound below **0.009568**, crossing the 0.01 cutoff. The bound remains below
**0.009631** on an interval around it. This improves the previous bound of
0.010687 without changing the tilt choices or probability inequality.

The key is to use the proposed split. This explanation places the change
after reading 1, leaving a plateau of 99 readings. Bounding the covariance
directions within that long plateau loses less information than a generic
bound that only knows two levels were removed.

This checkpoint implements and tests the mathematical bound. It does not
rerun the complete reporting search, so there is no additional detection
claim. The latest complete replay still has four alerts and three non-alerts
among seven saved cases; the fresh-study count remains 58/144.

Separating the losses before changing the method
------------------------------------------------

At the saved split 1 and rho=51/64, keep the existing tilt tau=1/8 and the
same density correction. Change only how the residual eigenvalues are
obtained. The resulting comparisons are:

.. list-table:: Counterfactual probability bounds at the saved witness
   :header-rows: 1

   * - Eigenvalue information
     - Probability bound or estimate
     - Directions used
   * - Previous generic Toeplitz comparison
     - 0.010687
     - 8
   * - Exact precision counts, full-history comparison
     - 0.010235
     - 12
   * - Exact precision counts, long-plateau comparison
     - 0.009568
     - 12
   * - Numerical eigenvalues of the projected residual covariance
     - 0.009546
     - 12
   * - Numerical estimate of the directional probability itself
     - 0.0043470
     - All residual coordinates

The new long-plateau result has a rational certificate. The table's other
counterfactual calculations use floating arithmetic to separate the costs;
the last two rows are diagnostic, not rigorous certificates. The `saved
diagnosis <data/sturm_v1_diagnosis.json>`_ includes exact point and interval
proofs as well as these comparisons at both existing tilts.

Improving the eigenvalue input is enough to cross the cutoff here. The new
bound is already close to the density bound using numerical projected
eigenvalues. The larger remaining gap to the numerical probability comes
from the density inequality itself. This tells us where further mathematical
effort would have the most effect at this particular explanation.

Why using the split helps
-------------------------

The previous `spectral refinement <spectral_refinement.rst>`_ bounds the
k-th largest eigenvalue of the residual covariance C=U'R U by the
(k+2)-th largest eigenvalue of the full covariance R. Removing the two
plateau levels costs two positions in this generic comparison.

There is another valid comparison. Consider residual contrasts supported
entirely inside one plateau, with coefficients summing to zero there.
Such contrasts are also valid residual directions for the full two-plateau
model. If that plateau has length L, its covariance is the ordinary AR(1)
matrix R_L, and only one level has to be removed. Therefore::

    lambda_k(C) >= lambda_k(one-plateau residual covariance)
                >= lambda_(k+1)(R_L).

Eigenvalues here are ordered from largest to smallest. The first inequality
follows by restricting the candidate subspaces in the eigenvalue min-max
principle. The second is codimension-one interlacing. We use the longer
plateau when it has at least k+1 readings.

At split 1, the long plateau has 99 readings. Its residual space is the
entire residual space of this proposed two-level model: the first level
absorbs the single first reading. This makes the comparison especially useful.
For balanced splits, the original full-history comparison may be stronger.
The implementation keeps the maximum of all valid lower bounds.

This argument restricts the directions used for an eigenvalue comparison.
It does not condition on the first reading or replace the reporting model.
The determinant and observed directional threshold still use the existing
full-history calculations.

Exact eigenvalue bounds from the precision matrix
-------------------------------------------------

For v=|rho|, write::

    R_L(v) = (1-v**2) * inverse(A_L(v)).

A_L is tridiagonal: its first and last diagonal entries are 1, its interior
diagonal entries are 1+v**2, and its adjacent off-diagonal entries are -v.
Thus a lower bound on the j-th largest covariance eigenvalue follows from
an upper bound on the j-th smallest precision eigenvalue.

The preceding Toeplitz comparison replaced the precision matrix with a
simpler one and then bounded a sine function. We can avoid those losses by
counting the eigenvalues of A_L directly with rational arithmetic.

For a proposed threshold h, form the leading-principal-minor determinants
of A_L-hI. They obey the recurrence::

    D_0 = 1
    D_1 = 1-h
    D_i = (diagonal_i-h)*D_(i-1) - v**2*D_(i-2).

The number of sign changes in D_0, ..., D_L is the number of precision
eigenvalues below h. This is the Sturm count. With nonzero minors, each
sign change corresponds to a negative pivot in an LDL decomposition;
congruence preserves the number of negative eigenvalues.

The implementation scales each rational minor by a positive common
denominator and updates its integer numerator. This preserves every sign
while avoiding repeated fraction simplification inside the recurrence. A
test compares the integer signs with the original rational recurrence.

Zero minors need care. When an interior minor is zero and v is nonzero,
its two neighboring minors have opposite signs by the recurrence. Skipping
the zero therefore counts one sign change. A zero final minor is also
skipped, giving the strict count below an exact eigenvalue. When v=0,
the matrix is the identity and the count is handled directly. Tests cover
exact eigenvalues and interior zeros, rather than perturbing them numerically.

Bisection starts from the known spectral enclosure::

    0 < eigenvalue_j(A_L) <= (1+v)**2.

Whenever at least j eigenvalues lie below the midpoint, it becomes the
new upper endpoint; otherwise it becomes the lower endpoint. After 32
iterations, the upper endpoint is a proved rational bound. The bracket
width is at most (1+v)**2 / 2**32. The calculation never rounds a numerical
eigensolver result and treats it as a proof.

Covering a correlation interval
-------------------------------

Let [v_low, v_high] enclose |rho| over the requested interval, and let
m=(v_low+v_high)/2. The maximum absolute row sum of A_L(v)-A_L(m) is at most::

    epsilon = (v_high**2-m**2) + 2*(v_high-m).

The first term bounds the diagonal change; the second covers the two
neighboring off-diagonal changes. This also bounds the operator norm, so::

    eigenvalue_j(A_L(v)) <= eigenvalue_j(A_L(m)) + epsilon.

If H_j is the certified eigenvalue upper bound at m, then throughout the
whole correlation interval::

    lambda_j(R_L(rho)) >= (1-v_high**2) / (H_j + epsilon).

Negative correlations have the same covariance eigenvalues as their
positive magnitudes, through an alternating diagonal sign change. The
same argument therefore covers both signs and intervals crossing zero.

The new `implementation <directional_sturm.py>`_ compares three residual
eigenvalue lower bounds: the previous Toeplitz bound, the full-history
Sturm bound with j=k+2, and the long-plateau Sturm bound with j=k+1.
Their maximum replaces the eigenvalue input to the existing density
formula. The tilts remain 1/16 and 1/8, and direction counts remain
4, 8, 12, and 16. Every bound still majorizes the same directional p-value.

The saved certificates
-----------------------

For split 1 and rho=51/64, the selected certificate uses twelve directions,
the long-plateau bound, and tau=1/8. Its point upper bound is approximately
0.009567494435, with the exact rational result in the archive. The interval::

    [52223/65536, 52225/65536]

has an upper bound below 0.009631. The interval has nonzero width, which
is necessary for eventual use in the continuous reporting search.

The second lost case, with its witness still at rho=7/8, does not cross
the cutoff. Its new point bound is about 0.04603 and its saved numerical
probability estimate is about 0.01682. As before, tightening an upper bound
on that same probability is unlikely to eliminate this second explanation.

Validation and next step
------------------------

The 28 tests cover exact Sturm roots and zeros, independent dense precision
eigenvalues, interval perturbation bounds, projected residual eigenvalues
at both correlation signs, split reversal, unit/level invariance, the saved
point and interval, conservative failure, and argument validation before
cache lookup. Eigenvalue bisections are cached because the same covariance
bound is reused at both tilts; the cache does not change the arithmetic.

The diagnosis checks its input archive and previous source hashes, then
recomputes the previous witness certificate before constructing the new
proofs. The `running guide <running.rst>`_ gives reproduction commands.
Frozen references and ASV's production detector are unchanged.

The subsequent `complete Sturm replay <sturm_reporting.rst>`_ integrates this
eigenvalue bound and verifies all seven saved cases. Its early negative
correlation case is unresolved at the default depth limit. Later experiments
below sharpen the Fourier density estimate at its 101/128 witness; further
subdivision then finds a different survivor at 99/128.

Keeping stronger directions in the density bound
-------------------------------------------------

At the 101/128 explanation, the four groups of positive tilted coefficients
have different lower bounds. The earlier calculation used the weakest bound
in all sixteen directions, discarding information from the first groups.

For group j, let b_j be a lower bound on the four coefficients at ranks
4j-3 through 4j. The four corresponding Gaussian-square factors bound the
characteristic function by ``(1 + 4*b_j**2*u**2)**-1``. Applying Hölder's
inequality to the m group factors and integrating in Fourier inversion gives::

    sup density <= c_m / (b_1 * ... * b_m)**(1/m)
    c_m = binomial(2*m-2, m-1) / 4**m.

When all b_j are equal, this is the existing m-group bound. When earlier
groups have stronger coefficients, their larger values raise the geometric
mean and lower the density bound. The implementation rounds that geometric
mean downward to a dyadic rational using integer nth roots. This keeps the
bound conservative with exact rational arithmetic.

At split 1 and rho=101/128, this changes the saved point bound from about
0.012103 to about 0.009224. On [51711/65536, 51713/65536], the bound is about
0.009284. Both are below the 0.01 cutoff. The complete replay with this
grouped bound is recorded in `the reporting results <sturm_reporting.rst>`_.
It moves the unresolved region toward |rho|=1 rather than changing the
overall alert count.

Integrating the unlike groups directly
--------------------------------------

Hölder is still a relaxation: it bounds the integral of a product using the
separate integrals of each group. We can instead integrate the product
itself. Write a_j=2*b_j for the lower bound of group j. The characteristic
function is bounded by::

    product_j (1 + a_j**2*u**2)**-1.

When the rates a_j are distinct, partial fractions give the exact Fourier
integral and hence this density upper bound::

    density <= (1/2) * sum_i a_i**(2*m-3)
                         / product_(j != i)(a_i**2 - a_j**2).

The terms can have different signs, but their exact rational sum is
positive. For equal or nearly equal group bounds, the implementation rounds
each b_j downward to a dyadic rational and separates collisions by one
dyadic unit. The resulting rates are distinct and no larger than the
original rates, so this remains conservative. Only four groups are used,
keeping the exact partial-fraction calculation small.

At rho=101/128, the point bound falls from the grouped-Hölder value about
0.009224 to **0.0087203**. On [51711/65536, 51713/65536], it is **0.0087751**.
Both use tau=1/8 and sixteen directions. These are exact rational
certificates rounded here for readability. The full replay excludes this
local explanation but the search finds a different survivor at split 2 and
rho=99/128 when its depth limit is raised to 20. The overall saved-case alert
count therefore remains unchanged.
