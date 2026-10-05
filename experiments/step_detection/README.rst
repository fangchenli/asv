How ASV step detection works
============================

ASV estimates an underlying performance history as a sequence of constant
plateaus, then decides which upward changes should be reported as regressions.
The useful separation is between the signal model, the algorithm that fits it,
the rule that chooses its complexity, and the reporting policy.

This experiment directory records the implementation at source revision
``d33754e129c025beb5c2ca440c3c280433b264f7``. The companion
`improvement plan <improvement_plan.rst>`_ defines experiments and acceptance
criteria. It includes executable reproductions of two initial findings.

The path from measurements to reports
-------------------------------------

The publishing pipeline connects the components as follows::

    Stored benchmark values and confidence intervals
        |
        | commands/publish.py and _stats.get_weight
        v
    Graphs grouped by benchmark and environment parameters
        |
        | Graph.get_data: sort revisions and average repeated entries
        | Graph.detect_steps: dispatch each parameter series
        v
    step_detect.detect_steps
        |
        | filter missing observations and normalize weights
        v
    solve_potts_autogamma
        |
        | try penalties; share an interval-cost cache
        +--> solve_potts_approx
        |        |
        |        +--> solve_potts with a maximum segment length
        |        |        |
        |        |        +--> C++ RangeMedian or Python L1Dist
        |        |
        |        +--> merge_pieces and adjust boundaries
        |
        | score candidates using complexity and correlated residuals
        v
    Segments in observation positions
        |
        | graph._compute_graph_steps maps positions to revisions
        v
    regressions plugin filters history and applies reporting thresholds
        |
        | step_detect.detect_regressions
        v
    Regression records and Atom feed

`Publishing <../../asv/commands/publish.py>`_ reads each benchmark result and
its statistics. `Graph <../../asv/graph.py>`_ sorts entries by revision and
averages values and weights separately when several entries share a revision.
Parameterized benchmarks are split into individual series for detection.
Publishing can distribute these series across worker processes.

The detector operates on ordered observations, not elapsed time. After missing
values are removed, two neighboring observations may represent distant commits
or measurement dates. Neither the fit nor its correlation model uses those
gaps as distances.

What the weights mean
---------------------

`get_weight <../../asv/_stats.py>`_ returns ``2 / abs(ci_99_b - ci_99_a)``
for an available, finite, nonzero confidence interval width. A wider interval
therefore gives an observation less influence. Missing statistics, infinite
interval endpoints, and zero widths yield a missing weight.

``detect_steps`` removes ``None`` and NaN observations. It also removes
observations whose provided weight is nonpositive. For retained observations,
missing weights become 1 after normalization; provided weights are divided by
the median of the nonmissing, non-NaN input weights. A missing or zero median
falls back to 1. The implementation collects that median from the original
weight list, rather than only the retained observations.

These are relative uncertainty weights, not inverse variances. The loss is
absolute deviation, so inverse uncertainty scale is the relevant analogy to a
heteroscedastic Laplace model. The procedure estimates overall residual noise
from the history instead of taking individual confidence intervals as a
complete noise model.

The Potts objective
-------------------

Let ``y[i]`` be the measured value, ``w[i]`` its positive relative weight,
and ``x[i]`` the fitted underlying value. For a fixed penalty, ASV fits::

    loss(x) = sum_i w[i] * abs(y[i] - x[i]) + gamma * J(x)

    J(x) = number of positions where x[i] != x[i - 1]

The first term rewards agreement with observations. The second charges a
fixed amount for each change, encouraging neighboring fitted values to agree.
This is the one-dimensional Potts formulation for a signal with few jumps.
It is also described as an L0 penalty on successive differences. Unlike a
total-variation penalty, the jump term counts changes rather than summing their
magnitudes. See the `Potts reconstruction paper
<https://arxiv.org/abs/1304.4373>`_ for the connection to sparse recovery.

The documentation sometimes charges ``gamma * K`` for K segments. Since
``K = J + 1`` for a segmentation with distinct adjacent levels, this differs
by a constant and selects the same optimum at fixed gamma.

For a fixed interval, the level minimizing weighted absolute error is a
weighted median. This makes the fitted level less sensitive to extreme
measurements than a mean under squared error. Robustness is limited: an
extreme point can still justify a separate short segment when the reduction
in loss pays for the extra boundaries.

For example, with unit weights::

    Measurements    10 10 10 10 | 12 12 12 12
    One plateau     11 11 11 11   11 11 11 11    loss = 8
    Two plateaus    10 10 10 10 | 12 12 12 12    loss = gamma

The two-plateau solution wins for gamma below 8; there is a tie at 8. More
generally, two noise-free adjacent plateaus of lengths a and b and separation
delta save ``min(a, b) * abs(delta)`` in absolute error when split. This
explains how duration and magnitude work together as evidence for a change.

Exact optimization for a fixed penalty
--------------------------------------

Use half-open intervals and define::

    C(s, t) = min_mu sum_{i=s}^{t-1} w[i] * abs(y[i] - mu)

    F(0) = -gamma
    F(t) = min_{0 <= s < t} [F(s) + C(s, t) + gamma]

The recurrence tries every possible start of the final segment. Each choice
combines the best fit of the prefix with one new segment. Saving the winning
start at each endpoint permits backtracking to recover the full partition.
The initial value avoids charging for the first segment.

``solve_potts`` in `step_detect.py <../../asv/step_detect.py>`_ implements this
dynamic program, subject to its interval constraints. Its returned ``right``
positions are exclusive, while the internal ``mu(left, right)`` and
``dist(left, right)`` methods take an inclusive right endpoint. Keeping these
conventions separate is essential when comparing solvers.

Where C++ helps
---------------

``get_mu_dist`` selects `RangeMedian <../../asv/_rangemedian.cpp>`_ when the
extension is importable, otherwise ``L1Dist`` provides the Python fallback.
The extension performs both weighted-median calculations and the inner
dynamic-programming loop. Python still controls approximation, penalty
selection, and regression reporting.

For an uncached interval, C++ copies its value-weight pairs, sorts them,
locates the cumulative-weight midpoint, and sums weighted absolute deviations.
It caches the median and cost together in a ``std::map`` keyed by endpoints.
The fallback memoizes medians and distances in dictionaries and periodically
clears a large cache during its dynamic-programming loop. The C++ cache has
no analogous size limit in this revision. Automatic penalty selection reuses
one interval-cost object across candidate penalties.

The dynamic program has quadratically many candidate intervals in the
unrestricted case. That does not make the present implementation quadratic:
an uncached interval of length m requires sorting, costing O(m log m).
Summing that work over all intervals gives an O(n^3 log n) upper bound for
straightforward unrestricted evaluation, with potentially O(n^2) cached
entries. The ``solve_potts`` docstring's O(n^2 log n) claim therefore should
not be read as a bound established by this interval-cost implementation.

The approximation used in normal detection
------------------------------------------

``solve_potts_approx`` initially restricts segments to at most 20 observations
when ``min_size < 10``. Otherwise it uses ``min_size + 50``. The resulting
partition can have artificial boundaries solely because of this restriction.

``merge_pieces`` repeatedly considers neighboring segments and removes a
boundary when the merged interval costs no more than the separate intervals
plus gamma. Its scan uses an early-exit heuristic, so this is not a global
search over all possible merges. It then makes one pass over surviving
boundaries, trying offsets within the initial maximum segment length, and
recomputes the segment levels and costs.

The final segments can be much longer than 20 observations. The fixed window
limits the initial dynamic program, not the output plateau length. Merging
and boundary adjustment also evaluate long intervals, so the docstring's
linear-time description is not a demonstrated worst-case bound for the
whole procedure. The initial short-window pass does explain its practical
speed advantage.

Choosing the penalty and estimating noise
-----------------------------------------

``solve_potts_autogamma`` scales candidate penalties by the absolute-error
cost of fitting the entire series with one plateau, using 1 if that cost is
zero. It searches on a logarithmic penalty scale. Each candidate produces
a segmentation through the approximate solver.

For a candidate fit, let ``e[i] = y[i] - fitted_level[i]``. The score uses::

    S(rho) = w[0] * abs(e[0])
             + sum_{i=1}^{n-1} w[i] * abs(e[i] - rho * e[i - 1])

    score = (4 * log(n) / n) * K + log(sigma_0 + S(rho))

Here K is the segment count and the implementation approximately minimizes
S over rho in [-1, 1]. Residual correlation continues across fitted segment
boundaries. Positive correlation lets a sequence of similar residuals count
as less surprising than independent errors would suggest.

This correlation model only ranks candidate segmentations. Their medians and
boundaries are fitted with the independent weighted absolute-error objective;
ASV does not jointly optimize the correlated model over levels and boundaries.

The noise floor ``sigma_0`` prevents a zero residual sum from producing
negative infinity after the logarithm. The current implementation sets it to::

    More than two segments:  0.1 * smallest adjacent fitted-level difference
    One or two segments:     0.001 * smallest absolute fitted level
    Lower bound in all cases: 1e-300

Both the outer penalty search and the inner correlation search use
``golden_search``. The outer score is constant while a candidate partition
stays unchanged and can have multiple competing minima. Golden-section
search therefore has no general global-optimum guarantee for that objective.

The `statistical explanation <../../docs/source/step_detection.rst>`_ motivates
the score using Laplace noise and a Bayesian information criterion. Its
prefactor, approximation path, and noise floor are practical choices. The
score is not a posterior probability or a calibrated significance test.

From fitted steps to reported regressions
-----------------------------------------

``detect_steps`` maps filtered positions back to original observation
positions and returns tuples::

    (left, right, fitted_median, minimum_observed_value, error)

Bounds are half-open. The error is the segment's weighted absolute-error sum
divided by its retained observation count, not by its total weight. It is a
residual scale estimate, not a confidence interval for the fitted median.
Missing observations need not belong to a segment, and the final included
position is always a retained observation.

``graph._compute_graph_steps`` translates bounds into revision coordinates.
The `regressions plugin <../../asv/plugins/regressions.py>`_ selects the
configured history and threshold, then calls ``detect_regressions``. The
function's default relative threshold is zero; the plugin supplies 0.05 when
no configured threshold matches. These are different defaults at different
layers.

``detect_regressions`` scans backward through the steps. It compares earlier
levels with later best levels, accounting for residual error and the relative
threshold. It conditionally accepts short intervals and can retract an
apparent regression when the earlier level explains a brief dip. Its default
``min_size=2`` is a reporting rule, not a minimum plateau size in fitting.
After graph mapping, segment widths are revision spans, which can differ
from numbers of measured observations.

The reported value after a regression is the best subsequent fitted value,
not necessarily the next plateau. Consequently, a temporary slowdown that
fully recovers need not remain a reported regression. The plugin maps
revision bounds to commit hashes, distinguishes a single responsible commit
from a range, and writes the regression records and feed.

Reading and validation guide
----------------------------

Read ``detect_steps`` and ``solve_potts_autogamma`` first to understand the
public path, then ``solve_potts_approx``, ``merge_pieces``, and ``solve_potts``.
The C++ implementation explains the computational costs; the graph and
regressions modules explain the user-visible behavior.

`Existing tests <../../test/test_step_detect.py>`_ cover clear steps, noise,
weights, missing values, correlated signals, and reporting rules. The
`existing benchmark <../../benchmarks/step_detect.py>`_ measures a repeated
two-level signal. The improvement plan broadens this into separate checks
of optimizer correctness, statistical detection quality, and runtime.

The `online explanation
<https://asv.readthedocs.io/en/latest/step_detection.html>`_ is useful context,
but differs from this implementation in details: it describes replacing zero
weights rather than dropping their observations, and describes a different
segment-count cutoff for the noise-floor formulas. Experiments should pin
the source revision and state which behavior they evaluate.
