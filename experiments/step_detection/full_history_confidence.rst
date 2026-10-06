A full-history confidence reference
===================================

The new reference restores the stationary information identified in the
`information-loss analysis <ar1_information.rst>`_. It uses a proper joint
predictive density and checks the continuous correlation range with exact
rational interval bounds. It is implemented separately from the frozen
conditional reference and production detection.

The saved missed 8% slowdown still does not alert. Correlation 0.9999 is now
excluded at its true location, but 4095/4096, approximately 0.999756, remains
and fails both reporting tests. The construction restores useful evidence;
the first predictive density still spends too much of it.

This is a mathematical reference and a development diagnosis. No new
statistical evaluation or claim about average sensitivity is made here.

Choose a density before looking for a correlation
-------------------------------------------------

The confidence rule needs a proper density q for the entire history. It
cannot use a fitted maximum likelihood value in the numerator: such a value
is generally not a normalized density over possible histories.

We use a finite mixture of Bayesian predictive densities. First put half
the probability on no change and half on exactly one change. Distribute
the latter half equally over the n-1 possible locations. Within each model,
integrate out the plateau levels and a shared noise variance.

The predictive models use independent Gaussian noise. This does not impose
independence on the true data. Any proper density can serve as q in the
density-ratio confidence argument. Poor prediction under correlated noise
will, however, make the confidence set wider.

All the following calculations use readings divided by an externally
specified reference unit. The initial prior is fixed as follows:

* A plateau level, conditional on variance s squared, has a normal prior
  with mean zero and variance ``s**2/kappa``, where ``kappa=0.0001``.
  Different plateau levels are conditionally independent.
* The common variance has an inverse-gamma prior with shape ``alpha=0.5``
  and scale ``beta``. Its density is proportional to
  ``(s**2)**(-alpha-1) * exp(-beta/s**2)``.
* Mix equally over 25 choices ``beta = 0.5 * 2**(2*j)``, for integer j from
  -12 through 12.

Each prior is proper and every mixture weight is positive. The scale
mixture covers a broad range without selecting a variance from the scored
history. Its finite range is not a restriction on candidate noise variance:
each inverse-gamma component itself has support over all positive variances.

These are initial reference settings, chosen before the first calculation
of this density on the archived case and left unchanged after that check.
They have not been calibrated for detection quality. In particular, the
independent level priors give the jump a broad prior too. A more economical
model of the baseline and jump is a possible subsequent design.

Why the reference unit is an explicit input
-------------------------------------------

A proper prior on levels and noise introduces a scale. ``unit`` declares
that scale in the same measurement units as the observations. For example,
an external reference of one millisecond is ``unit=1`` for millisecond
readings and ``unit=0.001`` after converting those readings to seconds.
Scaling observations and unit together preserves the calculation.

The API requires this input rather than estimating it from the current
history. Choosing it from the same observations would make the numerator
data-dependent in a way the present coverage proof does not justify.
The archived diagnosis uses unit=1, an explicit development choice.

The integrated density in closed form
-------------------------------------

For each plateau j, let n_j be its length, ybar_j its sample mean, and
W_j its within-plateau sum of squared deviations. Define::

    D = sum(W_j + kappa*n_j/(kappa+n_j) * ybar_j**2)

For a particular partition and beta, integration gives::

    log q_partition,beta
        = -n/2 * log(2*pi)
          + 1/2 * sum(log(kappa/(kappa+n_j)))
          + log Gamma(alpha+n/2) - log Gamma(alpha)
          + alpha*log(beta)
          - (alpha+n/2)*log(beta+D/2)

The shrinkage term in D is the cost of levels far from their prior center.
The log factors involving plateau lengths account for integrating uncertain
levels rather than substituting their best fits. Mixing these normalized
densities over beta and partitions gives q. All possible locations contribute;
we do not take the largest component as the predictive density.

The component density is also a multivariate Student distribution. Tests
compare the formula against that independently evaluated density, and check
the mixture weights and invariance under reversing the history.

The confidence set and error bound
----------------------------------

For a candidate location t and correlation rho, the full stationary
likelihood, maximized over unrestricted plateau levels and innovation
variance, is::

    U(t,rho) = sqrt(1-rho**2) * [2*pi*e*Q_t(rho)/n]**(-n/2)

Q_t is the profiled residual cost using the stationary AR(1) precision
numerator. Allowing unrestricted levels makes U an upper bound for the
model with positive levels as well. At the true pair, U is at least the
true joint density p. Consequently::

    expectation(q/U) <= integral((q/p)*p) = integral(q) = 1

Retain a pair when ``q/U <= 1/delta``, with delta=0.01. Markov's inequality
gives a true-pair exclusion probability at most 0.01. The existing reporting
tests receive the remaining 0.04, with the same size/shape allocation and
cutoffs as before. Requiring every retained pair to fail a reporting test
therefore keeps the total per-history false-alert bound at 0.05 under the
declared stationary Gaussian, constant-variance, one-change model.

This argument uses the full density. It does not require splitting the
history or independence between confidence construction and reporting.

Certifying correlation intervals without a large polynomial
-----------------------------------------------------------

The existing exact GLS algebra supplies low-degree polynomials N and d
such that ``Q=N/d``, with d positive in the stationary domain. Cancel only
factors common to N and d; removing endpoint factors independently would
preserve their signs but change the residual cost.

Define a positive coefficient::

    C = (delta*q)**2 * (2*pi*e/n)**n

Then a pair is excluded precisely when::

    C*N(rho)**n > (1-rho**2)*d(rho)**n

Expanding these powers would create a polynomial whose degree grows with n.
Instead, on an interval [left,right], use exact Bernstein coefficients to
find a lower bound L for N and an upper bound H for d. Let R be the exact
maximum of ``1-rho**2`` on the interval. If::

    L > 0, H > 0, and C*L**n > R*H**n

the whole interval is excluded. This is a sufficient certificate; failure
to prove it calls for subdivision or another rejection route. The Gram
determinant identity establishes the denominator's positivity inside the
stationary range.

This avoids both a correlation grid and high-degree polynomial expansion.
The calculation still raises rational numbers to the nth power. Close to
1, R becomes small; for a positive limiting residual cost, that permits
certification of an entire endpoint neighborhood.

The search combines this confidence certificate with the existing exact
size and extra-step certificates. A midpoint can supply a surviving pair;
it cannot certify an interval. Only complete coverage of every location
and the open correlation range (-1,1) produces an alert. Exhausting 4096
cells or depth 16 returns an unresolved non-alert.

Predictive log densities, distribution cutoffs, and the initial log C are
computed in floating point. The implementation converts the exponential
of that computed log C downward to a rational number using 80-digit decimal
arithmetic. Interval certificates are exact for the resulting coefficient;
this does not certify rounding errors in the preceding density calculation.
Unrepresentable density/coefficient calculations disable confidence exclusion.

What happens to the saved missed slowdown
-----------------------------------------

The archived history has independent noise with standard deviation 0.05
and an 8% increase from 10 after reading 25 of 100. With unit=1, the new
predictive log density is 128.093. At the generating location:

.. list-table:: Confidence evidence for the same archived readings
   :header-rows: 1

   * - Assumed correlation
     - log(q/U)
     - Excluded at cutoff log(100)=4.605?
   * - 0
     - -26.783
     - No
   * - 0.9
     - -3.585
     - No
   * - 0.99
     - 1.864
     - No
   * - 4095/4096
     - 4.194
     - No
   * - 0.9999
     - 4.648
     - Yes

The reference also certifies the whole interval [0.99999,1) as excluded at
this location. Nonetheless, its complete reporting search returns the
same surviving pair, location 25 and rho=4095/4096, after 49 visited cells.
Independent GLS checks confirm that this pair fails both reporting tests.

At that pair, improving the predictive density by a factor greater than
``exp(4.605-4.194)``, about 1.51, would cross the confidence cutoff. This is
the cost of excluding one particular witness. Other locations and
correlations may still survive, so it is not a sensitivity forecast or a
sufficient condition for an alert.

The diagnostic also records two explicit execution fixtures with a 20%
step and variation ``0.1*sin(1.7*i)``. The 24-reading version retains a
near-1 witness. The 100-reading version produces a complete alert
certificate. The longer fixture was added to exercise successful certificate
reconstruction after the short fixture failed; these are not independent
evaluation histories or estimates of power.

Where this leaves the next design
---------------------------------

We now have a proper full-history numerator, a continuous confidence check,
and an independently verifiable reporting reference. The missing stationary
information is restored, but the predictive density remains the immediate
efficiency question.

Next, decompose its prediction cost into the level prior, jump/location
mixture, and noise-scale mixture. A natural model uses one common baseline
and a separate jump prior, instead of independently placing broad priors on
both levels. Derive its integrated density and the effect of its assumptions
before selecting another candidate. Any use of the saved history is
development work; measured gains require a subsequent frozen evaluation on
fresh observations. The present diagnosis does not justify replacing the
existing detector.

The follow-up `baseline/jump analysis <baseline_jump_prior.rst>`_ now supplies
that decomposition and identifies a removable location-mixture cost.
Candidate-specific proper predictors certify the saved missed slowdown with
either the original or revised level prior. This is still development evidence;
a fresh frozen comparison is the next step.

Code and verification
---------------------

The reference is `reporting_ar1_full.py <reporting_ar1_full.py>`_, with API
``evidence(values, unit=..., config=None, max_cells=4096, max_depth=16)``.
Results preserve the predictive prior, density, coefficient, residual ratios,
interval certificates, and any surviving witness.

`full_history_diagnosis.py <full_history_diagnosis.py>`_ reproduces the archived
case and the two explicit fixtures. `Saved diagnostics
<data/full_history_v1_diagnosis.json>`_ include their full reference outputs.
See `running instructions <running.rst>`_ for commands.

The 23 new analytical tests compare the integrated density to multivariate
Student densities, verify mixture weights and reversal, check GLS residual
ratios against direct stationary likelihoods, and reconstruct complete alert
certificates. They also cover a narrow surviving region that a grid misses,
the saved non-alert witness, unit conversion, numerical fallback, degenerate
data, and work limits. The mathematics extends the density-ratio argument
described in `Universal Inference <https://arxiv.org/abs/1912.11436>`_; the
specific predictive mixture and interval bounds here are our construction.
