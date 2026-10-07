A stronger bound for the saved missed slowdown
==============================================

We can now rule out one explanation that stopped the previous detector.
The improvement comes from keeping the full noise covariance in a probability
bound. Its former approximation threw away too much information.

For the first of the two lost detections in the `fresh directional study
<directional_report.rst>`_, the old probability bound was 1. The new bound
is about **0.006072 at the saved explanation**, below the 0.01 rejection
cutoff. It remains below **0.007037 throughout a nonzero correlation
interval** around that explanation. These are rigorous upper bounds.

This is a mathematical development result on saved data. It removes a known
obstacle; establishing a new detection requires the complete search over
all possible explanations. The measured result remains 58/144 detections.

What explanation was stopping us?
---------------------------------

Imagine a history of 100 benchmark timings. The data were generated with an
8% slowdown after reading 25. The detector does not know where the change
happened or how strongly neighboring measurement errors are related.

It therefore considers alternatives. One says: the change happened after
the very first reading, and the measurement errors have correlation 0.875.
Under that alternative, a long rise in timings can partly be explained by
persistent noise. The first plateau has only one observation, so its level
is uncertain. As long as this explanation remains plausible, the detector
cannot confidently report a slowdown greater than 5%.

The confidence check asks whether the *pattern* left after removing the two
plateau levels is unusual under that proposed correlation. It returns a
directional p-value: smaller values provide stronger evidence against the
proposed explanation. The detector rejects it when a proved upper bound on
that p-value is below 0.01.

Both saved losses have this same proposed location and correlation. Their
actual generating correlation was -0.5. The proposed 0.875 is an alternative
the detector has to evaluate, not an estimate of the known simulation input.

Why the previous bound failed
-----------------------------

The previous `tail certificate <directional_tail.rst>`_ replaced a whole
covariance matrix by a single lower bound. This is useful very close to
correlation 1. For 100 readings at correlation 0.875, its coefficient is::

    kappa = 2 - (100 - 3) * (1 - 0.875) = -10.125.

A negative coefficient provides no useful positive lower bound on the
covariance. The code consequently returns probability bound 1. This says
the calculation is inconclusive; it does not say that the proposed model
fits the readings well.

The new calculation retains the complete covariance matrix. Its determinant
measures the combined contribution of all residual directions, avoiding
the loss caused by replacing them with one common coefficient.

The probability inequality
--------------------------

Fix a proposed split t and correlation rho. Let X contain the two plateau
indicator columns, and let U be an orthonormal basis perpendicular to them.
Multiplying the readings by U' removes both unknown plateau levels. There
are d = n - 2 remaining coordinates.

Write R for the AR(1) covariance matrix with entries rho**|i-j| and C = U' R U.
Normalize the residual vector z = U' y to s = z / ||z||. This removes the
unknown overall noise scale. Define::

    r = 1 / (s' C^-1 s)
    W = sum_j (lambda_j(C) / r - 1) * xi_j**2,

where the xi_j are independent standard normal variables. The directional
p-value is P(W <= 0), with value 1 for a constant directional density.
The earlier `directional derivation <directional_tail.rst>`_ establishes
this representation.

For any fixed 0 < tau < 1/2, the event W <= 0 implies exp(-tau W) >= 1.
Taking expectations gives::

    p <= E[exp(-tau W)]
       = det((1 - 2*tau) I + (2*tau/r) C)^(-1/2).

This is the Chernoff bound. The determinant appears because the normal
coordinates are independent: their exponential expectations multiply, and
the determinant is exactly that product over eigenvalues. We retain the
previous choice tau = 1/16; no search over tilts is needed for this result.

The old method bounded this determinant using the crude covariance
coefficient. We instead enclose the full determinant directly. This version
does not need the previous tilted-density refinement to exclude the first
saved witness.

How a large determinant becomes a small calculation
---------------------------------------------------

We could construct U and calculate a dense (n-2)-dimensional determinant.
The AR(1) model lets us avoid that and keep every operation rational.

For positive constants a and b, set M = a I + b R. Completing the normalized
columns of X with U gives an orthogonal basis. The block determinant identity
in that basis gives::

    det(U' M U) = det(M) * det(X' M^-1 X) / det(X' X)
               = det(M) * det(X' M^-1 X) / (t * (n-t)).

The second determinant is only two by two: there are two unknown plateau
levels. To calculate the first determinant and the inverse products,
use the tridiagonal AR(1) precision numerator A::

    R^-1 = A / (1-rho**2).

A has endpoint diagonal entries 1, interior diagonal entries 1+rho**2,
and neighboring off-diagonal entries -rho. With q = 1-rho**2, define::

    T = a A + b q I.

Since M = R T / q and det(R) = q**(n-1),::

    det(M) = det(T) / q
    M^-1   = (I - b q T^-1) / a.

T is tridiagonal and positive definite. An LDL decomposition calculates its
determinant as a product of positive pivots. Two tridiagonal solves supply
T^-1 X. Thus the full residual determinant needs a linear number of scalar
operations and a final two-by-two determinant. This also works for negative
correlations within the stationary range.

Proving the bound throughout an interval
----------------------------------------

A point calculation is insufficient for the complete detector: it must
cover every correlation, including values between sampled points.

The existing exact GLS polynomials express the minimized precision-weighted
residual error as Q = N/h. With S equal to the within-plateau Euclidean
residual sum of squares, the radius has the identity::

    r(rho) = (1-rho**2) * S / Q(rho).

Here Q uses A, rather than R^-1, which explains the factor 1-rho**2.
Bernstein polynomial coefficients enclose N and h over a whole interval.
When their lower bounds are positive, they give::

    r(rho) <= r_upper
           = upper(1-rho**2) * S * upper(h) / lower(N).

Use a = 1-2*tau and b = 2*tau/r_upper. Since C is positive definite and
r <= r_upper, replacing r by r_upper can only decrease the determinant.
It therefore gives a conservative probability upper bound.

The implementation carries out the tridiagonal calculation with intervals.
Every arithmetic operation rounds its lower endpoint down and its upper
endpoint up to a rational number with denominator 2**192. Dependency between
occurrences of rho can widen the enclosure, but cannot invalidate it.
An uncertain pivot, a denominator containing zero, or a nonpositive final
determinant lower bound produces an inconclusive result.

If L is the resulting positive determinant lower bound, the final decision
is the exact rational comparison::

    p**2 <= min(1, 1/L) < delta**2.

No numerical eigensolver, quadrature result, or square-root approximation
determines that decision. The displayed probability is rounded upward too.
The new routine accepts intervals strictly inside (-1, 1). The existing
endpoint certificates remain necessary for a complete search.

What the two saved cases show
-----------------------------

The new interval is [3583/4096, 3585/4096], centered at 7/8. Both rows use
the proposed split after reading 1. ``Early`` and ``middle`` identify the
actual simulated change locations, after readings 25 and 50 respectively.

.. list-table:: Saved-witness diagnosis
   :header-rows: 1

   * - Actual change
     - Old point bound
     - New point bound
     - New interval bound
     - Numerical p estimate
   * - Early
     - 1
     - 0.006072
     - 0.007037
     - 0.0003212
   * - Middle
     - 1
     - 0.082641
     - 0.095379
     - 0.01682

The displayed new bounds are rounded upward. The final column repeats
the previous floating quadrature estimates and is not a rigorous bound.
The `saved certificates <data/determinant_v1_diagnosis.json>`_ contain the
exact rational inequalities, input hashes, and source hashes.

The early case now has an interval certificate below 0.01. In the middle
case, even the estimated ideal directional p-value exceeds that cutoff.
Tightening an upper bound on the same p-value is therefore unlikely to
remove that particular explanation.

Would combining size and noise evidence solve the second case?
--------------------------------------------------------------

There is a useful independence result. At the true split and correlation
under the common-scale Gaussian model, the GLS estimate of the plateau
levels is independent of the residual vector. Whitening that residual
vector separates its length from its direction; those are independent too.

The directional p-value depends only on the whitened residual direction.
The one-sided size statistic depends on the fitted levels and the whitened
residual length. Consequently their p-values are independent at the true
pair. This remains true when the actual slowdown is below the threshold.

More explicitly, the GLS estimator is L y with
L = (X' R^-1 X)^-1 X' R^-1. Its covariance with z = U' y is proportional to
L R U = (X' R^-1 X)^-1 X' U = 0. Joint normality gives independence.
The whitened vector C^-1/2 z is isotropic normal, and its squared length
is the GLS residual sum of squares. These identities give the two
independences used above.

For independent valid p-values p_direction and p_size, define::

    v = p_direction * p_size
    p_combined = v * (1 - log(v)).

The correction appears because the product of two independent uniform
variables has cumulative distribution v*(1-log(v)). The raw product alone
would be too small to use as a calibrated p-value. Valid conservative
p-values give a conservative combined value by the same argument. The
formula is increasing in v, so certified upper bounds on the component
p-values would also give an upper bound on their combination.

One possible new rule would allocate 0.042 to this combined check and retain
0.008 for shape checks, preserving a total budget of 0.05. It would replace
the separate confidence/size rejection rule: appending it at full budget
would require a new error calculation. It also need not preserve every
rejection made by the old rule.

For the middle witness, however, the size p-value is about 0.9881. Its
estimated earlier and later levels are 10.641 and 10.386: that particular
explanation offers essentially no evidence for a greater-than-5% slowdown.
Combining it with the estimated directional p-value 0.01682 gives about
**0.08472**, above the proposed 0.042 cutoff. This natural combination does
not remove the second witness. That diagnosis helps us avoid implementing
a larger change on an unsupported expectation of fixing both losses.

Checks and next step
--------------------

The new `standalone implementation <directional_determinant.py>`_ has tests
against independently calculated exact dense determinants, at both signs
of correlation and several split locations. Tests cover interval arithmetic,
the saved witnesses, level/scale/reversal invariance, uninformative residual
directions, and conservative failure on wide intervals or low precision.
The `running guide <running.rst>`_ gives reproduction commands.

The subsequent `complete reporting replay <determinant_reporting.rst>`_ now
adds this certificate as another bound on the same directional p-value.
Taking the minimum of valid upper bounds requires no extra confidence
budget. Seven saved cases verify with no gained or lost alerts: the first
lost case finds another witness at correlation 13/16. The next step is a
sharper probability inequality, followed by complete-search verification.
The second lost case remains a modeling/calibration question.
