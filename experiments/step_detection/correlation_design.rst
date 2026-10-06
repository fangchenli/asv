Exact fitting of residual correlation
=====================================

The correlation parameter can be fitted by a weighted median for each
candidate history. This replaces the inner golden-section search and fixes
its premature stopping on equal trial costs. The segment fitter, outer
penalty search, and existing noise floor retain their present objectives.

The implementation is ``_fit_ar1`` in ``asv/step_detect.py``. Both the Python
and C++ segment-fitting backends use this Python helper. The
`mathematical analysis <mathematical_analysis.rst>`_ gives the wider context;
this design specifies the decisions needed for the focused correction.

Objective and domain
--------------------

For a fixed candidate, e[i] is the measured value minus its fitted level.
The conditional correlation objective is::

    S(rho) = w[0] * abs(e[0])
             + sum_{i=1}^{n-1} w[i] * abs(e[i] - rho * e[i-1])

Use the closed interval [-1,1]. It includes the limiting values of the
documented stationary model and guarantees that a minimum is attained.
At the endpoints, the finite-history conditional score remains well-defined;
this does not assert that a stationary infinite-history process exists there.
There is no extra stationary-initial-distribution term in this objective.

This intentionally corrects the previous expanded bracket, approximately
[-4.236,4.236]. A fit outside [-1,1] no longer competes. Consequently, the
new conditional score may exceed a previous score obtained outside the
allowed domain, even though the new solver minimizes the stated objective.

Solution and ties
-----------------

For each nonzero previous residual, define::

    knot[i] = e[i] / e[i-1]
    mass[i] = w[i] * abs(e[i-1])

The variable part of S is the sum of ``mass[i] * abs(rho - knot[i])``.
Its minimizers form the weighted-median interval. Clip this interval to
[-1,1] and choose its member closest to zero. This chooses the least
persistence among equally good fits; it does not alter their score.
If there are no nonzero previous residuals, every rho ties and the result
is zero.

We can clip each knot before finding the median. For example, a knot q > 1
contributes ``mass * (q - rho)`` throughout the allowed domain. Replacing
q with 1 subtracts a constant from that term, preserving the minimizers.
The same argument applies below -1. The final score is always evaluated
from the original residuals, so those constants are retained.

Numerical implementation
------------------------

The helper accepts finite residuals and positive finite weights, as used
by the prepared fitting problem. It uses ordinary floating-point arithmetic;
the exactness claim concerns solving the conditional objective without an
iterative search tolerance.

* Compare absolute residuals before dividing. When the numerator is at least
  as large as the denominator, the clipped knot is simply -1 or 1 according
  to their signs. Otherwise the quotient has magnitude below one. A tiny
  denominator therefore never creates an overflowing ratio.
* Split each weight and absolute previous residual into a mantissa and
  power of two. Form their product after extracting the exponents, then
  normalize all masses by a shared power of two. This prevents overflow
  of large products and underflow when all products are uniformly small.
* Sort the clipped knots and find where their signed cumulative mass balance
  crosses zero. Binary search with freshly computed ``math.fsum`` balances
  reduces cancellation error compared with subtracting rounded prefix sums.
  If the balance is zero, the interval extends to the next knot.
* Zero previous residuals have no effect on the minimizing rho. Their
  constant contributions, the first residual, and transitions across segment
  boundaries all remain in the final S calculation.

Normalization cannot retain arbitrarily large relative dynamic ranges in
double precision. Masses that underflow after normalization are omitted;
rounding may affect ties. The final residual score can also overflow if the
input's weighted errors exceed the representable range. The helper avoids
unnecessary product overflow but does not redefine ASV's input domain.

Sorting and the signed-balance searches require O(n log n) time and O(n)
storage. This establishes a complexity bound, not a runtime advantage over
the previous small number of O(n) search evaluations.

Analytical acceptance cases
----------------------------

For unit weights and residuals [-1,0,1], S(rho) = 2 + abs(rho), so the answer
is zero and the score is 2. The previous search could stop at -1 or 1 with
score 3. The production integration test forces this candidate and checks
the actual automatic-selection score, including its existing floor.

Residuals [1,-1,-1] with unit weights give equal mass at -1 and 1. Every
rho in [-1,1] minimizes S, so the chosen value is zero. Residuals [4,1,1]
with weights [1,1,4] give equal mass at 0.25 and 1; the chosen value is 0.25.
Increasing the last weight to 5 moves the optimum to 1.

The tests also enumerate 375 small residual/weight combinations. For each,
an independent rational-arithmetic calculation evaluates the original S at
all its in-domain breakpoints, both endpoints, and zero. Those positions
contain a minimum and the closest-to-zero representative of any flat minimum.
Additional cases cover zero residuals, boundary optima, dangerous ratios,
overflowing or underflowing unscaled products, and rescaling invariance.

Noise model decision
--------------------

Retain the existing candidate-dependent floor for this correction. Replacing
it with a shared lower noise bound requires a definition of that bound from
measurement resolution or an explicit noise assumption. No such definition
has yet been established for all ASV benchmark types.

The proposed shared-floor model and its penalty-path guarantee remain the
next modeling question. Exact conditional correlation fitting does not make
the full segmentation problem jointly optimal, establish the outer search's
unimodality, or calibrate regression alerts.

The subsequent `shared-floor analysis <shared_noise_floor.rst>`_ derives that
score and shows that allowing rho arbitrarily close to one can hide a lasting
step. The closed correlation domain solves the specified conditional fit;
the persistence assumption remains a separate model decision.

Validation and historical results
----------------------------------

Use the analytical and existing regression tests before any new benchmark
campaign. The harness now times the production correlation helper directly,
records ``rho_fits`` and ``rho_fit`` time, and preserves optional reuse for
identical candidate fits in the grid and hybrid variants.

The step-detection and harness suite passed 103 tests, with two skips for
tests requiring a free-threaded Python build. The 375 rational-oracle cases
are included in that suite. All 36 graph and publishing integration tests
also passed. Ruff, whitespace, and documentation parsing checks passed.

The earlier baseline measured the old correlation search. Its results remain
historical evidence, not measurements of this implementation. Reproduce that
version from revision ``d526c2e``; the current branch uses the corrected fit.
