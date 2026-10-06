Mathematical analysis of step detection
=======================================

ASV combines three decisions: where to put boundaries, how many boundaries
to allow, and which changes to report as regressions. A better algorithm for
one decision need not improve the others. We can establish several useful
results before collecting more benchmark results:

* The fixed-penalty fitting problem admits exact dynamic programming and
  safe pruning with positive weights.
* The correlation parameter has an exact weighted-median solution. ASV's
  numerical search can stop at a nonoptimal value even on three observations.
* An exact penalty path suffices for a simpler, consistent selection score.
  The current correlation adjustment and candidate-dependent noise floor
  require additional reasoning.

The derivations below use the implementation in ``asv/step_detect.py`` at
``d526c2e``. Statements about proposed models are labeled as such. The
`ASV explanation <https://asv.readthedocs.io/en/latest/step_detection.html>`_
provides the statistical motivation; the formulas here follow the local code
where they differ. These are mathematical arguments and hand-worked examples,
independent of the previous experiment results.

.. contents::
   :local:
   :depth: 1

What the fitting objective asks for
-----------------------------------

Let y[i] be the measured value at observation i. A partition P divides those
observations into K consecutive segments. Within each segment, x[i] is a
single fitted level: the value used to represent that stretch of history.
The weights w[i] are positive and fixed during fitting.

Define the fitting error and objective as::

    D(P) = min over segment levels of sum_i w[i] * abs(y[i] - x[i])
    F_gamma(P) = D(P) + gamma * (K - 1)

The error rewards fidelity to the measurements. The charge per boundary
requires a reduction in error before a more complicated history is worthwhile.
Counting boundaries instead of their sizes is the Potts choice: a large
genuine jump pays the same boundary charge as a small jump.

For one segment, the best level is a weighted median. Moving the level to
the right increases error from observations on its left and decreases error
from those on its right. A minimum occurs where neither side has more than
half the total weight. If a whole interval of levels minimizes the error,
a deterministic choice within that interval gives the same D.

There is an exact detection threshold in a simple noiseless example. Suppose
one run has level a and total weight W_left, followed by level b with total
weight W_right. Fitting both runs with one level costs::

    min(W_left, W_right) * abs(b - a)

Fitting the two runs separately has zero error and costs one boundary charge.
The split wins exactly when::

    gamma < min(W_left, W_right) * abs(b - a)

This explains why duration and uncertainty matter as well as jump size.
Four unit-weight observations at 10 followed by four at 12 save
``4 * 2 = 8`` by splitting, so their critical penalty is 8. A short run
provides less evidence for its own level than a long run of the same height.
Noise and neighboring changes make the full problem less local.

Exact fitting and safe pruning
------------------------------

Write C(s,t) for the minimum weighted absolute error on observations with
indices s through t-1. The best cost for a prefix ending at t satisfies::

    B(0) = -gamma
    B(t) = min over s < t of [B(s) + C(s,t) + gamma]

Every partition has a last segment. Once its start s is chosen, everything
before it must be an optimal prefix; otherwise replacing that prefix would
improve the whole partition. This establishes the recurrence. The initial
minus gamma makes the final charge count boundaries rather than segments.

Weighted absolute error has an additional useful property::

    C(a,b) + C(b,c) <= C(a,c)

Two pieces can always reuse the best level of the combined piece. Allowing
them separate levels can only reduce their total error.

This gives a pruning proof. If at position t a previous candidate s satisfies::

    B(s) + C(s,t) >= B(t)

then for any later u::

    B(s) + C(s,u) + gamma
        >= B(s) + C(s,t) + C(t,u) + gamma
        >= B(t) + C(t,u) + gamma

Starting the last segment at s cannot beat starting it at t. We can discard
s while preserving an optimal objective. Equality can discard tied solutions;
preserving a particular boundary tie rule needs more care.

This is the zero-constant case of the pruning condition in
`PELT, Theorem 3.1 <https://arxiv.org/html/1101.1438v3>`_. The proof above
assumes unrestricted nonempty segments. With minimum lengths, t must become
a legal competitor before s is discarded; other constraints need their own
feasibility check. The expected linear-runtime result needs additional
assumptions and does not make interval costs free.

ASV's approximation instead restricts initial segment lengths, merges
segments, and shifts boundaries locally. Each accepted edit reduces or
preserves the objective, but that establishes improvement from its starting
partition, not a global optimum. A single boundary-adjustment pass also need
not reach a local optimum after neighboring boundaries move.

For an interval-cost implementation, maintain cumulative weight W and
weighted value sum T on each side of a weighted median m. Then::

    C = m * W_left - T_left + T_right - m * W_right

An order-statistic tree storing both sums can update a growing interval in
O(log n) time. Starting a new structure at each left endpoint yields an
O(n^2 log n) exact-cost construction with O(n) working storage when costs are
consumed immediately. Processing left endpoints in increasing order makes
B(s) final before its growing intervals update later prefix costs B(t).
This is a construction, not a speed claim for ASV's short windows.
The `L1 Potts paper <https://arxiv.org/abs/1207.4642>`_ gives
an O(n^2)-time, O(n)-space method for its discrete problem; adapting that
stronger bound to ASV's arbitrary weights requires a separate analysis.

What changing gamma can and cannot do
-------------------------------------

For a fixed partition, D and K are fixed. Its penalized cost is a straight
line in gamma, with slope K-1. The exact optimum is the lowest of these lines.
It is therefore piecewise linear and concave.

The number of segments in an exact optimum cannot increase as gamma grows.
To see this, compare optimizers (D1,K1) and (D2,K2) at gamma1 < gamma2::

    D1 + gamma1 * (K1 - 1) <= D2 + gamma1 * (K2 - 1)
    D2 + gamma2 * (K2 - 1) <= D1 + gamma2 * (K1 - 1)

Adding gives ``(gamma2 - gamma1) * (K2 - K1) <= 0``, so K2 <= K1.
Some segment counts may be skipped. The remaining boundaries need not be
a subset of those at the smaller penalty: fewer segments can require a
different arrangement of boundaries.

These results describe exact fitting, with consistent choices among ties.
They do not automatically transfer to ASV's merge heuristic. Moreover,
ASV does not minimize this lowest-line function to select gamma. It assigns
a different score to the fit produced at each gamma.

With a deterministic median and partition tie rule, that selection score is
constant wherever the selected fit is unchanged. Neither this fact nor
monotonicity of K establishes the unimodality required by golden-section
search. A logarithmic reparameterization does not add that guarantee.

The `CROPS method <https://arxiv.org/pdf/1412.3617>`_ uses intersections of
penalized-cost lines to recover the optimal tradeoffs across a penalty range.
It requires an exact fixed-penalty solver and suitable tie handling. It
cannot supply the same guarantee around an arbitrary approximate solver,
or guarantee that those candidates optimize a different correlated score.

When a penalty path is enough
-----------------------------

Consider the following simpler selection model, with beta > 0 and a positive
floor delta shared by every candidate::

    Q(P) = beta * K + log(delta + D(P))

Here candidate generation and scoring use the same fitting error D. We can
prove directly that a global optimum P_star is also an exact Potts optimum
for some penalty. Write its error as D_star and segment count as K_star.
Concavity of log gives, for any candidate P::

    log(delta + D) - log(delta + D_star)
        <= (D - D_star) / (delta + D_star)

Since P_star minimizes Q::

    0 <= Q(P) - Q(P_star)
       <= beta * (K - K_star) + (D - D_star) / (delta + D_star)

Multiplying by delta + D_star proves that P_star minimizes the fixed-penalty
objective at::

    gamma_star = beta * (delta + D_star)

This derivation is useful because it does not require every segment count
to occur on the path. Searching the exact tradeoffs can attain an optimal
Q even when some counts are skipped. A restricted penalty range must still
include a suitable gamma, and this identity alone does not give an iteration
with a global convergence guarantee.

ASV's current score is instead::

    Q_ASV(P) = beta * K + log(delta(P) + min_rho S(P,rho))
    beta = 4 * log(n) / n

Its candidates minimize independent error D, while S accounts for correlation
and delta depends on the fitted levels. Those substitutions break the proof's
premises. Even perfect exploration of the independent-error penalty path
would not establish a global optimum of Q_ASV. We should decide which scoring
model we want before investing in a guaranteed path-search implementation.

The correlation fit has an exact solution
-----------------------------------------

Fix a candidate history and write its residuals as e[i] = y[i] - x[i]. ASV
uses the following correlation-adjusted error::

    S(rho) = w[0] * abs(e[0])
             + sum_{i=1}^{n-1} w[i] * abs(e[i] - rho * e[i-1])

The parameter rho describes how much of one residual the previous residual
predicts. For every nonzero e[i-1], rewrite its contribution as::

    w[i] * abs(e[i] - rho * e[i-1])
        = w[i] * abs(e[i-1]) * abs(rho - e[i] / e[i-1])

Thus minimizing S is another weighted-median problem:

* Values: ``e[i] / e[i-1]``.
* Weights: ``w[i] * abs(e[i-1])``.
* Terms with zero e[i-1] are constant in rho and can be omitted when finding
  the minimizer. They remain in the final score.

For a closed allowed interval [L,U], clip a weighted median to that interval.
Convexity ensures this is a constrained minimizer. If all variable weights
vanish, every rho ties and zero is a natural choice when allowed. Sorting
gives an O(n log n) implementation; its practical cost must later be compared
with repeated O(n) objective evaluations. Very small denominators need careful
numerics, for example evaluating the piecewise-linear objective's subgradient
without forming unbounded ratios.

There is a concrete issue with the current search even though S is convex.
Take unit weights and measurements ``[9, 10, 11]`` fitted with one level, 10.
Then e is ``[-1, 0, 1]``, and::

    S(rho) = 1 + abs(rho) + 1 = 2 + abs(rho)
    S(-1) = S(1) = 3
    S(0) = 2

The inner ``golden_search`` first evaluates -1 and 1. Its default stopping
condition accepts equal function values as convergence, so it can return
an endpoint without investigating zero. Equality at two points of a convex
function does not imply a minimum between them has been found. This example
follows from the formula and stopping condition; no simulation is needed.

There is also a domain discrepancy. The call passes -1 and 1 with
``expand_bounds=True``. Substitution into ``golden_search`` gives an actual
bracket ``[-(2 + sqrt(5)), 2 + sqrt(5)]``, about [-4.236, 4.236]. The endpoints
-1 and 1 are initial evaluation points, not hard limits. The documented
stationary AR(1) model assumes abs(rho) < 1. Choosing a closed numerical
interval such as [-1,1], or a margin inside it, is a model decision; reproducing
the existing expanded bracket is a separate compatibility decision.

An exact conditional rho fit is therefore a focused first change to design.
It fixes this subproblem, but does not jointly optimize segment levels and
correlation: S includes terms crossing segment boundaries, so neighboring
levels are coupled. Independent weighted medians generally need not minimize
that joint objective.

Why the logarithm and floor need a model
----------------------------------------

Under independent Laplace errors with scales b/w[i], the negative log
likelihood, omitting terms fixed across fits, is::

    n * log(b) + D / b

For D > 0, minimizing over the unknown scale gives b = D/n. Substitution,
division by n, and removal of shared constants leave log(D). That is why
the logarithm appears: it comes from estimating the noise magnitude along
with the history. This calculation does not determine ASV's coefficient 4
or justify treating unknown discrete boundary locations as ordinary smooth
parameters in a BIC calculation.

ASV adds a floor to make perfect fits finite. Its raw floor is::

    K <= 2:  0.001 * min_j abs(level[j])
    K > 2:   0.1 * min_j abs(level[j+1] - level[j])

It then applies a lower bound of 1e-300. Because this floor changes between
candidates, creating a tiny fitted jump in a history with more than two
segments can reduce the floor as well as the residual error. If a K=3 fit
has smallest jump epsilon and zero residual error, its score is
``3 * beta + log(max(1e-300, 0.1 * epsilon))``. The floor can reward that
tiny jump more strongly as epsilon decreases, until the numerical bound
is reached. This describes the score; it does not assert that every such
candidate is visited by the fitter.

A proposed alternative is to specify one minimum noise scale b_min > 0
for the whole history and profile the same Laplace likelihood subject to
b >= b_min. The optimal scale and per-observation loss are then::

    b_hat = max(D/n, b_min)

    if D >= n * b_min:  loss = log(D/n) + 1
    if D <  n * b_min:  loss = log(b_min) + D / (n * b_min)

This produces a finite score for perfect fits from an explicit noise
assumption. It differs from adding a candidate-dependent floor inside log.
We would still need to justify b_min from measurement resolution or another
specified noise model, and decide whether independent errors are adequate.

Which transformations should preserve the answer
------------------------------------------------

The independent fixed-penalty problem has simple invariances. For a > 0,
replacing y by a*y+c multiplies D by a. Scaling gamma by a therefore preserves
optimal boundaries. Adding c alone preserves them without changing gamma.
Multiplying every weight by q > 0 also preserves boundaries if gamma is
multiplied by q. ASV normalizes measurement weights before fitting.

For a fixed candidate, adding c leaves residuals and S unchanged. But the
current floor for K <= 2 uses absolute levels, so the automatic score is
not translation invariant. For example, shifting levels 1 and 1.02 to 1001
and 1001.02 changes their floor from 0.001 to 1.001.

Positive unit rescaling multiplies both S and the raw floor by a. Away from
the hard 1e-300 bound, the score then changes by the shared constant log(a),
preserving candidate rankings. These are algebraic properties; finite
precision and numerical stopping tolerances can affect implementation behavior.

Whether translation invariance is desirable depends on the intended signal.
Absolute changes are preserved by adding a baseline; percentage changes are
not. For strictly positive timings, fitting log(y) is one possible relative
change model. If original uncertainty is approximately s, first-order error
propagation gives log-space uncertainty s/y, so weights proportional to 1/s
would become proportional to y/s before normalization. Zero or negative
measurements need a different model. Regression reporting has its own
percentage threshold and remains a separate decision.

Decisions before the next implementation
----------------------------------------

1. Specify the correlation domain and tie rule, then design the exact
   weighted-median rho fit with stable handling of zeros and tiny residuals.
   The three-observation example supplies an analytical acceptance case.
2. Choose whether to retain the current heuristic score or adopt a noise
   model with a shared floor. Establish its units, invariances, and treatment
   of exact fits before choosing a search algorithm.
3. For an exact fitting replacement, write the pruning rules for every
   supported segment constraint and specify the weighted interval-cost
   structure. Separate the correctness proof from runtime assumptions.
4. Apply penalty-path guarantees only after checking that candidate generation
   and selection satisfy their premises. Use experiments afterward to measure
   speed, false alerts, and detection sensitivity, which these proofs alone
   cannot determine.

The earlier `baseline report <baseline_report.rst>`_ remains a record of the
existing heuristic. The next work should resolve these mathematical choices
before extending that experimental comparison.
