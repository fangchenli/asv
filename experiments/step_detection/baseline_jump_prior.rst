Stop paying twice for the unknown change location
=================================================

The archived missed 8% slowdown now produces a complete alert certificate.
The decisive change is letting each candidate location use its own proper
predictive density. A better prior on jump size helps, but is not what
changes the decision in this example.

The confidence rule already checks every possible location. Its previous
predictor also spread probability across all those locations, paying a cost
of about 5.29 log-density units on the saved history. That second cost is
not required by the coverage proof. Removing it restores enough evidence to
exclude every pair that survives the size and shape tests.

This is a development result on a previously inspected history. The new
construction has not yet been evaluated for sensitivity on fresh histories.
The short 24-reading, 20% deterministic example still does not alert.

Account for the old prediction cost first
-----------------------------------------

The `full-history reference <full_history_confidence.rst>`_ assigns the
saved history log density 128.093. The best independent Gaussian fit at
its generating location has log likelihood 154.876. The gap is 26.783.
It can be separated exactly:

.. list-table:: Costs for the archived 100-reading, 8% slowdown
   :header-rows: 1

   * - Contribution
     - Log-density cost
   * - Integrating variance within the best noise-scale component
     - 2.502
   * - Averaging noise-scale components
     - 2.083
   * - Pulling levels toward the prior center
     - 3.931
   * - Integrating uncertain plateau levels
     - 12.979
   * - Averaging possible change locations and no change
     - 5.288
   * - Total
     - 26.783

These numbers are an accounting identity, not independently removable
penalties. Integration makes the predictor a proper density. Substituting
fitted levels or variance to avoid those costs would invalidate the current
argument unless replaced with another valid construction.

For the accounting, let B(S) be the variance-integrated log factor evaluated
at residual cost S. Let RSS be the least-squares residual, D the residual
plus the Gaussian-prior penalty, V the log determinant factor for integrating
the levels, and q_old the full mixture density. Then::

    log L_max - log q_old
        = [log L_max - max_beta B_beta(RSS)]
          + [max_beta B_beta(RSS) - B(RSS)]
          + [B(RSS) - B(D)]
          - V
          + [B(D) + V - log q_old]

The terms correspond to the rows above. In this example the generating
location dominates the mixture. Its weight is ``1/(2*(n-1))``, so the last
cost is effectively ``log(198)=5.288``. B(RSS) is used only for this
decomposition; without the level-integration factor it is not the predictive
density used by the detector.

A baseline and a jump are different quantities
----------------------------------------------

The original model places independent priors on the two levels a and b.
Instead write::

    m = (a+b)/2
    d = b-a
    a = m-d/2
    b = m+d/2

The midpoint m represents the overall level. The signed difference d
represents the jump. Using the midpoint rather than the earlier level keeps
the model symmetric under reversing the history: m stays the same and d
changes sign.

Changing coordinates alone must not change the density. If the old levels
each have prior variance ``s**2/kappa``, then::

    prior variance of m = s**2/(2*kappa)
    prior variance of d = 2*s**2/kappa

They are independent under that original Gaussian prior. Thus the equivalent
precisions are ``kappa_m=2*kappa`` and ``kappa_d=kappa/2``. The tests recover
the original density with exactly these choices.

The new prior keeps kappa_m unchanged and mixes over jump standard deviations
of 1, 2, 4, 8, 16, 32, or 64 times the noise standard deviation. All seven
receive equal weight. In comparison, the original jump standard deviation
is about 141 times the noise standard deviation. The new mixture assigns
more predictive probability to moderate changes while allowing either sign.

These scales were specified before calculating this variant on the saved
history and have not been adjusted to make it alert. They are development
choices, not calibrated defaults. The noise-scale mixture, location weights
in the global control, reference unit, and reporting budgets are unchanged.

Integrating the new prior
-------------------------

For a fixed location t, form a two-column design X. The first column is 1.
The second is -1/2 before t and +1/2 afterward. With
``Lambda = diag(kappa_m,kappa_d)``, define::

    G = X'X + Lambda
    beta_hat = inverse(G) * X'y
    D = ||y-X*beta_hat||**2 + beta_hat'*Lambda*beta_hat
    V = 1/2 * [log determinant(Lambda) - log determinant(G)]

The integrated log density is ``V + B_beta(D)``, with the same inverse-gamma
variance integration as the preceding reference. Averaging over its 25
noise scales and the seven jump scales gives a proper density r_t for each
fixed location. The nonnegative expression for D avoids subtracting two
large, nearly equal quadratic forms.

Independent multivariate Student-density calculations check this formula.
Under reversal, the corresponding location is n-t; the symmetric jump prior
ensures ``r_t(y) = r_(n-t)(reversed y)``.

Why a density for each candidate location is valid
--------------------------------------------------

For every fixed t, choose a proper density q_t over complete histories.
It may depend on t as a candidate parameter. It must not be fitted afterward
to the observed history in an unaccounted-for way. The density r_t above
integrates unknown levels and variance using the declared priors, so it meets
this requirement.

Let U(t,rho) be the stationary profile likelihood from the full-history
reference. The confidence set is now::

    C = {(t,rho): q_t(y) / U(t,rho) <= 1/delta}

At the fixed true pair (t_star,rho_star), U is at least the true density p.
The only numerator needed for true-pair coverage is q_t_star::

    expectation[q_t_star / U(t_star,rho_star)]
        <= integral[(q_t_star/p)*p]
        = integral[q_t_star]
        = 1

Markov's inequality therefore gives the same delta=0.01 coverage failure
bound. There is no union bound over the candidate locations. We need to
retain the true pair, and an alert requires rejecting every retained pair.
We are not selecting the largest density over t and using it as a common
numerator. That maximum would generally fail to integrate to one.

A constant mean also fits this argument: choose any fixed admissible t as
its representation with equal levels. No data-driven location selection is
needed for the proof. The reporting tests still receive 0.04, keeping the
same total 0.05 false-alert bound under the declared Gaussian AR(1) model.

Keep some probability on the original predictor
-----------------------------------------------

The implemented indexed density is::

    q_t = (q_old + r_t) / 2

Both terms integrate to one, so this is proper for each t. Retaining q_old
also guarantees ``q_t >= q_old/2`` on every history. In log-density terms,
the largest possible loss relative to q_old is log(2), approximately 0.693.
This is a density guarantee, not a bound on changes in detection rate.

For any two distinct normalized densities, one cannot dominate the other
everywhere: otherwise their integrals could not both be one. The retained
weight makes this tradeoff explicit instead of promising uniform improvement.

Only the confidence coefficient changes in the continuous search. Previously
it used a common ``C=(delta*q_old)**2*(2*pi*e/n)**n``. Now each location has
its own C_t built from q_t. The exact rational interval bounds and the
reporting size/shape certificates are unchanged. Each result records the
log density and coefficient separately for every visited location.

Separate the two changes on the saved history
---------------------------------------------

The diagnostic compares global and candidate-specific predictors, with
original and midpoint/jump priors. The indexed controls both retain half
the original global density. The global jump control retains half the
original density and half a global mixture using the new jump prior.

.. list-table:: Development comparison on the archived 8% slowdown
   :header-rows: 1

   * - Predictor
     - Log density at location 25
     - Full reporting result
   * - Original global mixture
     - 128.093
     - Non-alert
   * - Global mixture with new jump prior
     - 128.539
     - Non-alert
   * - Candidate-specific, original level prior
     - 132.694
     - Certified alert
   * - Candidate-specific, midpoint/jump prior
     - 133.444
     - Certified alert

Changing the jump prior in the global mixture raises its density by a factor
of 1.56. That exceeds the previously calculated 1.51 requirement for the
old witness. But another pair survives: location 25 and rho=16379/16384,
approximately 0.999695. Its size statistic is 1.811, below 1.873, and its
extra-step statistic is 6.333, below 16.924. Its confidence log evidence is
4.525, below log(100)=4.605. Removing one witness was not sufficient.

Both indexed variants cover every location and correlation with rejection
certificates after 121 visited cells. Each certificate uses 98 shape
intervals, ten size intervals, and two confidence intervals. The revised
jump prior raises density further, but candidate-specific prediction is
sufficient to change this example's decision even with the original priors.

The 100-reading deterministic 20% fixture also alerts with the new reference.
The 24-reading 20% fixture still retains a near-1 witness. A 24-reading 4%
fixture returns a non-alert at its actual location and correlation zero.
These are checks of mathematical behavior, not estimates of detection power
or false-alert rates.

What to evaluate next
---------------------

There is now a justified candidate for a new frozen comparison. Stop
adjusting priors on the saved example. Keep both indexed variants so that
the simpler original-prior model can win if the added jump mixture offers
little benefit.

Compare the conditional reference, the original full-history reference,
the two indexed references, and the known-correlation oracle on fresh paired
histories. Keep the reporting budget, external unit, priors, and search limits
fixed first. Include short histories, near-threshold changes, negative and
positive correlation, and early/recent locations. Count certified alerts,
surviving pairs, and unresolved searches separately. This comparison should
establish whether removing the location cost improves sensitivity broadly
and whether the jump prior contributes anything beyond that change.

Files and verification
----------------------

`reporting_ar1_jump.py <reporting_ar1_jump.py>`_ provides the indexed reference
and the two controls. `jump_prior_diagnosis.py <jump_prior_diagnosis.py>`_
reproduces the cost decomposition and saved examples. The `diagnostic output
<data/jump_prior_v1_diagnosis.json>`_ preserves all results and certificates.
See `running instructions <running.rst>`_ for commands.

The shared full-history search was factored into one engine so the variants
use identical confidence bounds and reporting tests. The original predictor's
API and results are preserved. Analytical tests check prior equivalence,
proper Student densities, mixture weights, reversal, the retained-density
bound, the cost identity, candidate-specific coefficient recording, units,
and the complete archived alert certificate.

The 26 new analytical tests accompany the existing reference checks. All
three saved alert certificates were independently reconstructed, including
complete interval coverage. All four saved non-alert witnesses passed exact
confidence membership checks and independent GLS reporting checks. The
original full-history reference reproduces its previously saved result.

The indexed coverage argument is an application of the density-ratio
principle described in `Universal Inference
<https://arxiv.org/abs/1912.11436>`_. The priors and their finite mixtures here
are our construction. As before, density and quantile calculations are
floating point; the subsequent interval certificates are exact for their
computed coefficients and cutoffs.
