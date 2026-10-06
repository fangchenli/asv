Bounding noise persistence
==========================

A strict correlation cap gives a useful guarantee: correlated noise cannot
make a large unexplained residual history almost free. With equal weights
and abs(rho) <= r < 1, the minimized correlated error is between (1-r) times
the independent error and the independent error itself.

The experimental scorer now accepts this cap explicitly. A maximum shock
half-life, measured in retained observations, provides an interpretable way
to specify it. The cap is a modeling input; the derivation does not identify
a universal half-life for all benchmark histories.

This extends the `shared-floor model <shared_noise_floor.rst>`_. Production
continues to use its existing correlation domain and selection score while
the prototype makes noise scale and persistence explicit.

.. contents::
   :local:
   :depth: 1

What a half-life means here
---------------------------

In an AR(1) model, a disturbance's contribution after k observations is its
initial size multiplied by rho**k. Subsequent disturbances add their own
contributions. The magnitude of the initial contribution therefore halves
after H observations when::

    abs(rho)**H = 1/2
    abs(rho) = 2**(-1/H)

Specifying a maximum half-life H gives a cap r = 2**(-1/H). Negative rho
alternates the sign of the contribution; the same bound controls its magnitude.
H=0 denotes independent noise, with r=0.

.. list-table::
   :header-rows: 1

   * - Maximum half-life in observations
     - Correlation cap
   * - 0
     - 0
   * - 1
     - 0.5
   * - 2
     - about 0.7071
   * - 4
     - about 0.8409
   * - 8
     - about 0.9170

These are parameter conversions, not recommended settings. A cap of 0.99
is also a duration assumption, even if it is described merely as a numerical
margin below one. Choosing H makes that assumption visible.

The residual error bound
------------------------

For a fixed fitted history, define unit-weight errors::

    D = sum_{i=0}^{n-1} abs(e[i])
    S(rho) = abs(e[0]) + sum_{i=1}^{n-1} abs(e[i] - rho*e[i-1])

The reverse triangle inequality gives::

    S(rho) >= D - abs(rho) * sum_{i=0}^{n-2} abs(e[i])
           >= (1-r) * D

Also rho=0 is allowed and gives S(0)=D. Thus::

    (1-r)*D <= min_{abs(rho)<=r} S(rho) <= D

This is a bound on every residual sequence, not a statement based on selected
examples. A fixed r < 1 limits how much correlation can discount the error.
It does not say that correlated and independent fitting will choose the
same boundaries.

For positive unequal weights, write w_min and w_max for their extremes.
Applying the unweighted bound gives::

    (w_min/w_max)*(1-r)*D_weighted <= min S_weighted <= D_weighted

The bound weakens when the weight ratio grows. Uniform guarantees over
increasing history lengths therefore need a weight ratio bounded away from
zero, in addition to a strict correlation cap. Replacing the innovations
with differences of standardized residuals would define a different model;
the current model correlates raw residuals and weights their innovations.

What this proves for a lasting step
-----------------------------------

Consider the earlier noiseless example: equal blocks at 0 and 2a, n observations
in total, unit weights, and a > 0. Any one-level fit has independent error
at least n*a. Consequently, even if its level and rho are jointly optimized,
its correlated error is at least (1-r)*n*a.

Let g be the shared-floor profile loss and keep b_min > 0 fixed. Every
one-segment candidate has score at least::

    beta + g((1-r)*a/b_min)

The correct two-segment fit has zero residuals and score 2*beta. For fixed
r < 1, the g term is a positive constant. As beta=4*log(n)/n tends to zero,
the two-segment fit eventually beats every one-segment fit. Any fit with
more than two segments costs at least 3*beta, so cannot beat that exact
two-segment fit when beta > 0.

This establishes the eventual preference of this objective on the specified
noiseless history. A candidate-generation heuristic still has to supply the
correct fit; the proof does not give it that guarantee or establish general
statistical consistency on noisy data.

For the particular one-level median fit at a, the exact correlated error is::

    min S = a * (n - (n-3)*r), n >= 4

At r=1 it collapses to 3a. At any fixed r < 1 it grows with n. The analytical
test at n=4096, a=1, and b_min=0.1 confirms the old two-candidate comparison
prefers merging when r=1 and prefers splitting with an illustrative half-life
cap of four observations. Four is a test input, not an inferred default.

A cap and a soft prior make different promises
----------------------------------------------

A hard cap excludes longer persistence from the model. A prior can allow it
while charging for it, but simply adding a prior does not guarantee that a
lasting step becomes preferable.

To see why, keep the two-block example and consider positive rho near one.
Write delta=1-rho. Its one-level median fit has::

    S = a * (3 + (n-3)*delta)

In the floor's linear regime, the total fit loss is S/b_min. A prior whose
negative log penalty stays bounded as rho approaches one only adds a bounded
cost there. It cannot eventually outweigh the extra segment's total charge
``n*beta = 4*log(n)`` in this example.

Even a divergent prior needs a specified strength. For illustration, use
``penalty(rho) = -lambda*log(1-rho)`` on 0 <= rho < 1, with lambda > 0.
This corresponds to a density proportional to (1-rho)**lambda. In the linear
regime, minimizing S/b_min plus this penalty gives::

    delta_star = lambda*b_min / (a*(n-3))

For sufficiently large n this lies in (0,1), and the combined fit/prior cost
grows as ``lambda*log(n) + a bounded term``. For the exact two-segment fit,
rho=0 gives zero residual and prior cost. If lambda < 4, the one-segment
candidate still eventually wins their comparison. If lambda=4, the remaining
constants matter. A generic weak prior therefore need not cure the problem.

This calculation concerns a particular prior and the profiled shared-floor
objective, not all Bayesian models. A full stationary initial-state model
would require its own likelihood. The current conditional score gives the
first residual a Laplace term independent of rho.

Time ordering is part of the assumption
---------------------------------------

The implemented half-life is in retained observations: after filtering
missing values, two neighboring entries count as one transition. It is not
a half-life in commits or seconds. ASV's graph orders observations by revision;
that ordering can differ from the order and timing of benchmark execution.

With known measurement times and an appropriate positive-correlation model,
one might use a decay factor depending on elapsed time. But skipping several
AR(1) observations also combines several innovations. A sum of independent
Laplace innovations is generally not Laplace with the original scale.
Substituting a time-dependent rho into the current score therefore needs
a corresponding innovation model, not only timestamps.

The prototype deliberately specifies an observation-order model. A future
measurement-time model should preserve execution timestamps or batch labels
and define the innovation distribution across gaps.

Implementation and analytical checks
------------------------------------

``_fit_ar1`` now accepts an optional ``rho_max`` in [0,1]. It first finds the
closest-to-zero weighted-median minimizer in [-1,1], then projects that value
onto [-rho_max,rho_max]. Convexity makes this an exact constrained minimizer
in real arithmetic. Because the smaller interval contains zero, projection
also preserves the closest-to-zero tie rule. The default remains 1, so
production calls retain the preceding implementation's behavior.

``noise_model.correlation_cap(H)`` converts an explicit maximum half-life
to r. It rejects a finite H so large that floating-point rounding would
produce r=1, which would lose the strict bound. ``fit_correlated_score``
accepts a candidate's residuals and weights plus its segment count, floor,
beta, and correlation cap. It returns rho, the cap, the residual error, and
the shared-floor score. It holds the candidate's levels and boundaries fixed.

The tests compare the conditional fitter with rational-arithmetic minima
on 1500 small residual/weight/domain cases. Further tests check the contraction
bound with equal and unequal weights, the half-life conversion, zero-residual
ties, input validation, and the lasting-step calculation.

The combined detector, harness, and noise-model suite passes 184 tests, with
two skips for tests requiring a free-threaded Python build. Ruff, whitespace,
and reStructuredText parsing checks pass. No benchmark campaign was run.

What is ready and what comes next
---------------------------------

The prototype now has an explicit noise scale and an explicit persistence
limit, each with a mathematical interpretation. Independent noise is the
r=0 member of the same family. No data-derived cap or floor is claimed.

The next useful implementation is an exact independent-noise reference for
the shared-floor score. It can establish a global optimum for a declared
floor and beta, separating candidate-generation errors from the choice of
noise model. That provides a mathematical reference before evaluating the
correlated candidates or calibrating defaults on real histories.
