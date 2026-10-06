Where the correlation confidence set loses information
======================================================

The forward-and-reverse proposal does **not** fix the missed 8% example.
It improves prediction substantially, but the same near-1 correlation still
survives both confidence checks and both reporting tests. Averaging the two
evidence values gives a stronger combination with the same error bound;
that also leaves the example unresolved statistically.

This calculation identifies two costs worth addressing: weak predictive
likelihoods and information lost by conditioning on half the history.
The full-history stationary likelihood retains information that both
conditional directions discard. That is the next mathematical design to
develop. No new sensitivity experiment or detector change is made here.

Start with the missed example
-----------------------------

The `frozen study <ar1_report.rst>`_ includes a history with 100 readings,
independent Gaussian noise of standard deviation 0.05, and a change from
10 to 10.8 after reading 25. Its archived identifier is
``independent_constant-ar1-study-n100-d0.08-s0.005-early-seed1400``.

When given the true correlation, the reporting rule easily establishes a
slowdown above 5%. With correlation unknown, it retains location 25 and
rho=4095/4096, approximately 0.999756. This pair also fails both reporting
rejection tests. Removing that pair is a necessary first check for the
proposed confidence improvement; removing it alone would not prove an alert.

The following calculations reuse exactly those archived readings. They
diagnose an observed failure and are not an independent evaluation of a
new method.

What the likelihood ratio is paying for
---------------------------------------

For a proposed location and correlation, the confidence rule compares::

    E = q / U

Here q is the probability density assigned to the validation readings by
a predictor chosen from training data. U is the largest conditional density
the candidate model can assign after adjusting its unknown levels and noise
scale, or an upper bound on that density. A candidate is excluded when E
exceeds the reciprocal of its error budget.

A large U means the candidate explains the validation readings well. A
small q means our predictor explains them poorly. Either keeps E small.
Thus a candidate can survive because prediction is weak even when its
correlation is substantially less convincing than the true correlation.

To make that distinction visible, let U0 be the exact conditional maximum
at correlation zero and U1 its limit as correlation approaches 1. Then::

    log E1 = [log U0 - log U1] - [log U0 - log q]
                correlation evidence       prediction cost

This is an accounting identity for the saved example, not a new test.
The first term measures how much worse near-1 correlation explains the
readings than correlation zero. The second measures how much the chosen
predictor loses against the fitted zero-correlation model.

.. list-table:: The balance in natural-log density units
   :header-rows: 1

   * - Direction
     - Correlation evidence
     - Prediction cost
     - Resulting log E1
   * - Forward
     - 15.129
     - 29.292
     - -14.163
   * - Reverse
     - 11.580
     - 10.799
     - 0.780

Splitting the 0.01 confidence budget equally requires log E greater than
log(200), approximately 5.298, in either direction. Both fall short.
The endpoint 1 is used only to derive a limit. At the actual stationary
witness 4095/4096, the corresponding values are -14.169 and 0.774.

Why reversing helps, and why it is insufficient
-----------------------------------------------

In forward order, the training half contains the real step. The predictor
fits it as persistent noise and chooses correlation 0.95. In reverse order,
the first half lies entirely on the later plateau, and its fitted correlation
is 0.0454. The reverse predictor therefore describes the noise much better.

It still must predict an unknown future jump. The current density mixes
equally over 51 possibilities: no jump and each possible validation location.
Its jump magnitude also has a broad Gaussian prior. This spreads predictive
probability over many histories. These are legitimate choices for coverage,
but that does not make them efficient for detecting this particular change.

For comparison, inserting the known generating conditional density yields
limiting log E values of 14.786 forward and 10.380 reverse. Those exceed
5.298. This uses unavailable knowledge of the true levels, location, and
noise scale, so it is only a diagnostic. It shows that this conditional
likelihood comparison can reject the problematic correlations with a better
numerator. It does not show that a practical predictor can recover that gain.

How much does the enlarged candidate model cost?
------------------------------------------------

The current denominator fits three independent innovation means: one before
the proposed step, one for the transition reading, and one afterward. A
real two-level model links those means through its two levels, a and b.
Removing that extra freedom can only reduce U and increase E.

For validation indices i, write::

    z_i = y_i - rho*y_(i-1)
    c = (1-rho)*a
    d = b-a
    expected z_i = c + d*v_i

The coefficient v is zero before the step, 1 at the transition, and 1-rho
afterward. Fitting an intercept and this one column gives the exact
conditional profile over unrestricted levels.

There is also an explicit formula for the cost of enforcing the link. Let
L and R be the numbers of ordinary validation innovations before and after
the step, u and v their respective means, and z the transition innovation.
When both groups are nonempty::

    exact SSE - relaxed SSE
        = [(1-rho)*z + rho*u - v]**2
          / [(1-rho)**2 + rho**2/L + 1/R]

This follows by projecting the three unconstrained group means onto their
one linear constraint. It is nonnegative, as required for a smaller model.
If the change is already inside training, validation has only one plateau
intercept and this enlargement costs nothing.

That is exactly the forward situation in the saved example. The reverse
direction has two validation plateaus, but the enlargement has little cost:
at rho approaching 1, exact SSE is 0.188680 versus relaxed SSE 0.188150.
The resulting improvement in log E is only 0.0703, far below the roughly
4.59 needed by the relaxed check to reach its cutoff.

Requiring positive levels does not remove the saved witness either. At
rho=4095/4096, the exact forward fit can have a positive unobserved early
level and a validation level of 0.6012. The reverse fit has levels 10.0577
and 9.2799. Both attain the unrestricted conditional residual minimum with
positive levels. The forward fit's implausibly low level illustrates what
conditioning discards: its likelihood never scores the training observations.

Take the limit after profiling
------------------------------

For every rho below 1, the reparameterization by c and d is invertible.
As rho approaches 1, the original levels can grow while c stays nonzero.
Consequently, at an interior validation step the limiting exact fit retains
a common intercept for all non-transition differences. Its SSE is their
centered sum of squares.

Substituting rho=1 into the original a,b design first incorrectly deletes
that intercept. Profiling and taking the limit do not commute in that
parameterization. The analytical checks explicitly cover this distinction.

Combining two directions with a valid error bound
-------------------------------------------------

Let EF(t,rho) be forward evidence and ER(n-t,rho) be reverse evidence.
Stationary Gaussian AR(1) noise has the same distribution after reversal;
the mean step reverses direction and maps from t to n-t. The confidence
model allows either mean direction. Reporting still uses the original order.

At the true pair, each evidence value has expectation at most one. The
conditional expectation arguments use different training prefixes, but
each also gives an unconditional bound.

The originally proposed intersection keeps a pair when both EF <= 200
and ER <= 200. A union bound limits true-pair exclusion to 0.005+0.005=0.01.
The directions need not be independent.

There is a stronger combination at the same total budget::

    E_average = (EF + ER) / 2
    retain the pair when E_average <= 100

Its expectation is at most one by linearity, so Markov's inequality gives
the same 0.01 exclusion bound. Every pair excluded by the intersection is
also excluded by this average, because either value above 200 makes their
nonnegative sum exceed 200. The average can additionally reject two values
that are individually below 200 but sum to more than 200.

For the saved stationary witness, however, the log of the averaged evidence
is just 0.011 with relaxed denominators, or 0.081 with exact denominators.
Both are far below log(100), approximately 4.605. Thus even this stronger
combination, with the denominator tightened, cannot make this history alert.
Multiplying EF and ER would require a different validity argument; dependence
prevents simply multiplying their expectation bounds.

These calculations profile the nuisance parameters separately in each
direction. Requiring both directions to use the same levels and variance
could strengthen the check further, but still does not fix this example.
One common witness has original-order levels 10 and 10.8, rho=4095/4096,
and innovation variance 0.004562. Its forward and reverse log evidence
values are -13.812 and 1.251; the log of their average is 0.558, still below
4.605. These are candidate parameters demonstrating survival, not parameters
used to choose either predictive numerator. They imply much larger marginal
noise variance than the actual generator, which conditioning fails to exclude.

The information a full-history likelihood restores
--------------------------------------------------

The conditional calculation treats the training history as given. In
particular, it omits the stationary distribution of the first observation.
Reversing and intersecting two conditional confidence sets does not restore
the full stationary likelihood.

Write tau squared for innovation variance. The stationary first residual
has variance ``tau**2/(1-rho**2)``. For a complete history, let Q_t(rho) be
the minimum over plateau levels of::

    (1-rho**2)*(y_0-a)**2
      + sum(((y_i-mu_i) - rho*(y_(i-1)-mu_(i-1)))**2, i=1..n-1)

Here mu_i is a before t and b afterward. Profiling the unknown innovation
variance gives the full density maximum::

    U_full(t,rho)
        = sqrt(1-rho**2) * [2*pi*e*Q_t(rho)/n]**(-n/2)

The square-root factor matters. With positive limiting residual variation,
U_full tends to zero as rho approaches 1. The stationary model cannot keep
assigning substantial density to this fixed complete history while its
common-level uncertainty grows without bound.

For a fixed candidate step, the full profiled residual has limit::

    Q_t(1-) = sum((y_i-y_(i-1))**2 for i != t)

Unlike the conditional limit, these differences are not centered. To retain
a nonzero conditional drift, the levels would have to diverge as
1/(1-rho); the initial stationary term makes that increasingly costly.
For the saved example, this limit is 0.455346, strictly positive.

Suppose q_full is a proper, finite, positive joint density specified without
fitting it to the same complete history it scores. Then q_full/U_full has
expectation at most one at the true pair, by the same density-ratio argument.
Moreover it diverges as rho approaches 1 whenever the residual limit above
is positive. It therefore excludes a sufficiently small neighborhood of 1.
The conditional rule has no corresponding automatic endpoint exclusion.

This identifies information discarded by our construction; it is not an
impossibility result for unknown correlation. Nor does excluding a tiny
endpoint neighborhood guarantee useful detection: a surviving correlation
farther from 1 could still defeat the size test.

What to build next
------------------

Develop the full-history density construction before implementing another
reporting variant. Two parts need explicit treatment:

* Choose a proper joint predictor. One option is a proper density g for
  the training prefix multiplied by the existing conditional predictor:
  ``q_full = g(training) * q(validation | training)``. The density g must
  itself be proper; fitting a maximum likelihood density to that prefix
  does not establish this. A prior mixture or a sequential predictive
  construction could supply it. Units and level/scale assumptions need
  to be declared before evaluating fresh histories.
* Derive continuous confidence exclusion. The existing GLS polynomials
  already give ``Q=N/d``. For positive q_full, the strict exclusion inequality
  can be written as::

      (delta*q_full)**2 * (2*pi*e/n)**n * N**n
          - (1-rho**2)*d**n > 0

  This follows by squaring positive quantities in ``q_full/U_full > 1/delta``.
  Its degree grows with n, so the present low-degree quadratic machinery
  will need a new numerical design, possibly using certified bounds in log
  coordinates. The algebra alone is not a floating-point certificate.

Prediction efficiency remains essential. Adding a valid but diffuse density
for the prefix could spend more evidence than the restored likelihood gains.
Check the archived witness and deterministic examples first; freeze a new
evaluation only after those calculations justify a candidate.

Verification and sources
------------------------

`ar1_information.py <ar1_information.py>`_ reproduces the saved example from
its hashed archive. `The diagnostic output <data/ar1_information_v1.json>`_
records both directions, exact and relaxed fits, the saved witness, and the
full likelihood's endpoint behavior. No original frozen source or record is
changed.

The 28 analytical tests check the conditional reparameterization against
direct least squares, the relaxation formula, the order of limits, the full
likelihood against a stationary multivariate Gaussian density, reversal,
and the saved witness against independent GLS reporting calculations.

The general density-ratio confidence principle comes from
`Universal Inference <https://arxiv.org/abs/1912.11436>`_ by Wasserman, Ramdas,
and Balakrishnan. The AR(1) calculations, example diagnosis, and full-history
extension here are our derivations. They establish mathematical behavior
and implementation identities, not a measured sensitivity improvement.
