Choosing a shared noise floor
=============================

A shared floor can be defined as a lower bound on the noise scale in a
Laplace model. Profiling that model gives a finite score for perfect fits
and a clear interpretation for the floor. Its value must come from a
declared measurement or modeling assumption; the current relative weights
do not supply it automatically.

There is a second decision: how persistent the noise may be. Even with a
shared floor, the current correlation model can explain a lasting level
change as one noise shock. The example below proves that longer histories
can then make a real step less attractive to the selection score.

The experimental implementation in ``noise_model.py`` therefore requires
an explicit floor and complexity penalty. It implements the derived score
without changing ASV's production default. Its checks evaluate analytical
cases; no new benchmark campaign is needed for these conclusions.
The subsequent `persistence design <noise_persistence.rst>`_ adds an explicit
half-life cap and proves a bound on the discounted residual error.
The `exact independent reference <exact_reference.rst>`_ finds the global
optimum of this score for declared floor and beta inputs.

.. contents::
   :local:
   :depth: 1

Deriving the score from a noise bound
-------------------------------------

Fix a candidate partition with K segments and n retained observations.
Write S for its weighted absolute residual sum. For an independent-noise
model, S is the usual weighted fitting error D. For ASV's conditional
AR(1) model, S is the weighted innovation error minimized over rho.

With positive, dimensionless relative weights w[i], assume each innovation
has a Laplace scale b/w[i]. The negative log likelihood, after dropping
terms shared by all candidate partitions, is::

    L(b) = n * log(b) + S / b

Now specify a positive floor b_min, shared by every candidate, and constrain
b >= b_min. Since ``L'(b) = (n*b - S) / b**2``, its minimum occurs at::

    b_hat = max(S/n, b_min)

Substitute b_hat, divide by n, and subtract the common constant log(b_min).
With ``t = S / (n*b_min)``, the resulting loss is::

    g(t) = t                 if 0 <= t <= 1
           1 + log(t)        if t >= 1

The proposed selection score is::

    Q = beta * K + g(S / (n*b_min))

Below the floor, improving the fit reduces the score linearly. Above the
floor, the data determine the fitted scale and the reduction is logarithmic.
The two branches meet with equal value and slope at t=1. A perfect fit has
g(0)=0, so among exact fits a positive beta favors fewer segments.

This differs from adding a shared constant delta inside ``log(delta+S)``.
For delta = n*b_min, that alternative's excess loss is ``log(1+t)``.
It is smooth but is not the constrained maximum-likelihood calculation
above. Changing the score also changes the gain that competes with beta;
the existing coefficient ``4*log(n)/n`` remains a heuristic, not a newly
derived calibration.

The units and the role of weights
---------------------------------

Normalized weights are dimensionless, so b and b_min have the same units
as the observations. S is a sum across n observations. The corresponding
floor on the sum's scale is n*b_min, not b_min alone.

Converting observations from seconds to milliseconds multiplies S and
b_min by 1000, leaving t and Q unchanged. Adding a constant baseline leaves
residuals unchanged, so a fixed shared b_min preserves candidate rankings.
Multiplying every relative weight by q likewise requires multiplying b_min
by q to keep the same model.

For two candidates differing by one segment, the extra segment wins when
its reduction in g exceeds beta. If both residual sums lie below n*b_min,
that condition reduces to::

    S_merged - S_split > beta * n * b_min

Thus the floor determines how much absolute error reduction is enough to
justify another segment. It is a statistical assumption, not merely a way
to avoid evaluating log(0).

What ASV actually retains about uncertainty
-------------------------------------------

The local data path is:

* ``asv_runner.statistics.compute_stats`` estimates the sample median and
  its 99 percent confidence interval. The repository's
  ``test/test_statistics.py`` checks this interpretation.
* ``asv/_stats.py:get_weight`` converts the interval [a,b] to
  ``2 / abs(b-a)``. Missing, infinite, or zero-width intervals have no weight.
* ``asv/commands/publish.py`` supplies these weights to the graph.
* ``asv/graph.py:get_data`` averages repeated values and weights for each
  revision separately.
* ``asv/step_detect.py:detect_steps`` divides weights by their median and
  passes the normalized weights to the fitter.

When the raw weights are valid, their reciprocal gives a confidence-interval
half-width. But normalizing them removes their common absolute scale. If
every interval doubles in width, every raw weight halves and the normalized
weights remain the same. The fitter receives no separate copy of that scale.

Preserving raw widths would help, but would not finish the model. An interval
for an estimated median is different from a noise distribution for a single
run, and from correlated fluctuations between revision measurements. Averaging
weights is also not a general rule for computing the uncertainty of an average.
The independent sample sizes and relationships among repeated measurements
matter. In particular, zero-width intervals cannot establish zero physical
noise or a positive measurement-resolution bound.

The missing contract is therefore the meaning of b_min for the graph's final
observations, after reduction. We should carry absolute uncertainty and
measurement-resolution metadata through that reduction before treating them
as a calibrated floor. Simply setting b_min to a median interval width would
silently equate different kinds of uncertainty.

Why the observations alone cannot always supply a floor
-------------------------------------------------------

Suppose a data-only estimator T(y) is required to be both translation
invariant and homogeneous under positive unit scaling::

    T(y + c) = T(y)
    T(a*y) = a*T(y), a > 0

On a constant history, translation invariance gives ``T(c,...,c)=T(0,...,0)``.
Homogeneity gives ``T(0,...,0)=a*T(0,...,0)`` for every a > 0, so any finite
such estimator must be zero on constant histories. It cannot provide a
strictly positive floor everywhere.

A robust estimate from adjacent differences can be useful when its noise
assumptions hold, but it does not escape this result. Differences within
constant or quantized plateaus may vanish. Steps can contaminate them, and
correlation changes their distribution. A positive lower bound needs an
external scale or an explicit convention that sacrifices an invariance.

For now, the prototype requires b_min from the caller and supplies no
automatic fallback. This makes the uncalibrated choice visible rather than
hiding it in the selected partition's levels.

What an exact penalty path would guarantee
------------------------------------------

In the independent-error model, define ``h(D)=g(D/(n*b_min))``. This function
is increasing, concave, and differentiable, including at D=n*b_min. At an
optimum with error D_star, the tangent inequality gives::

    0 <= Q(P) - Q(P_star)
       <= beta * (K - K_star) + h'(D_star) * (D - D_star)

Therefore that optimum also minimizes a fixed-penalty Potts objective with::

    gamma_star = beta / h'(D_star)
               = beta * max(D_star, n*b_min)

This extends the earlier shared-additive-floor argument to the actual
constrained-likelihood score. Some segment counts may be skipped by the
penalty path without losing an optimal score. A practical path algorithm
still needs exact fixed-penalty fits, a sufficient penalty range, and
appropriate handling of ties. The identity does not prove that repeatedly
updating gamma from the latest fit converges globally.

The same tangent argument applies to S if candidate generation globally
minimizes ``S + gamma*(K-1)``. ASV's candidates instead minimize independent
error D approximately. Rescoring that candidate set with correlated S cannot
inherit the independent-error path guarantee.

A lasting step can look like a noise shock
------------------------------------------

Consider n=2m noiseless observations: m at 0, then m at 2a, with a > 0,
unit weights, and m >= 2. The one-segment median chosen by ASV is a.
Its residuals are m copies of -a followed by m copies of +a.

There are n-2 within-block transitions and one transition across the step.
For rho in [-1,1], the one-segment correlated cost is::

    S(rho) = a + (n-2)*a*(1-rho) + a*(1+rho)
           = a * (n - (n-3)*rho)

It is minimized at rho=1, with S=3a, independent of the length of either
block. The two-segment fit has S=0. In contrast, the one-segment independent
error is D=n*a and grows with the length of the history.

For sufficiently large n, both correlated candidates are in the floor's
linear regime. Taking ASV's beta=4*log(n)/n, the split wins exactly when::

    3*a / (n*b_min) > 4*log(n)/n
    3*a / b_min > 4*log(n)

The left side is fixed while the right side grows. Eventually the score
prefers one segment, even though more observations support the lasting change.
For a=1 and b_min=0.1, it prefers the split at n=32 and the merged fit at
n=4096. These are comparisons of the two specified candidates, not a claim
about the full detector's selected output. The tests calculate these scores
directly. They are a counterexample to assuming the floor alone fixes the
persistence ambiguity.

Excluding rho=1 with an open interval does not resolve it: rho can approach
one arbitrarily closely and the cost has the same infimum. A fixed cap
rho_max < 1 changes the one-segment cost per observation to::

    S(rho_max)/n = a*(1-rho_max) + 3*a*rho_max/n

Its limit is now positive, so this particular step eventually earns its
extra segment as beta tends to zero. But choosing that cap requires a
meaning for the maximum noise persistence. A per-observation cap also depends
on how measurements are ordered and spaced. ASV orders this series by
revision, which need not be the time order in which benchmarks were run.

Decisions and next implementation
---------------------------------

The shared-floor score is specified and implemented as an experimental pure
function. Its caller must provide a positive floor and a nonnegative beta.
Production keeps its current score while the following model choices are
resolved:

1. Define the floor's absolute scale for reduced graph observations. Preserve
   the uncertainty or resolution metadata needed to justify it; do not infer
   it from candidate levels.
2. Define how long noise may persist, or use independent noise as an explicit
   reference model. A correlation cap or prior needs a declared interpretation.
3. For that chosen model, establish the candidate-generation guarantee before
   comparing search algorithms. Exact independent fitting and a penalty path
   provide a tractable reference; jointly correlated fitting is a different
   optimization problem.

Use ``profile_loss(S, n, b_min)`` for the excess fit loss and
``selection_score(S, n, K, b_min, beta)`` for candidate comparison. Both live
in ``experiments/step_detection/noise_model.py``. ``test_noise_model.py`` checks
the constrained likelihood, scale invariance, numerical extremes, supporting
penalties, and the persistent-step counterexample.

All 48 analytical tests pass. Ruff, whitespace checks, and reStructuredText
parsing also pass. This step changes the experimental score and design
documents; it does not change the production detector or rerun the baseline.
