ASV step detection implementation reference
===========================================

The `conceptual guide <README.rst>`_ explains why ASV uses fitted levels,
absolute error, change penalties, and a noise score. This reference connects
those choices to the source code, with exact formulas and executable examples.
Its example positions start at zero to match Python; the guide numbers
versions from one.

Suppose a benchmark takes about 10 milliseconds for several commits, then
about 12 milliseconds for later commits. Individual measurements fluctuate:
one run takes 9.9 milliseconds, another 10.1. ASV tries to work out whether the
underlying performance changed, where that change happened, and whether the
slowdown should appear in its regression report.

It does this by fitting flat stretches to the measurements. A flat stretch is
called a plateau or segment. The boundary between two different plateaus is
a step, also called a change point. The main idea is to explain the data with
a small number of plateaus, allowing ordinary measurement noise within each.

Start with the worked example below, then read the four questions that the
implementation answers. The later sections explain the mathematics and map
it to the code. No prior knowledge of the Potts model or Bayesian statistics
is needed for the example.

This experiment directory records the implementation at source revision
``d33754e129c025beb5c2ca440c3c280433b264f7``. The companion
`improvement plan <improvement_plan.rst>`_ defines experiments and acceptance
criteria. It includes executable reproductions of two initial findings.

.. contents:: On this page
   :local:
   :depth: 2

A small example from beginning to end
-------------------------------------

Imagine these timings, in milliseconds, for eight consecutive commits::

    Observation position   0     1     2     3  |  4     5     6     7
    Measured timing       10   10.1   9.9   10  | 12   12.1  11.9   12
    Fitted timing         10    10    10   10  | 12    12    12    12

Positions start at zero because the implementation uses Python list indices.
For now, every measurement has equal weight. The vertical line is the change
we would like the detector to find; it is not supplied as an input.

First, consider fitting just one plateau. The median of all eight values is
11, so one possible fit assigns 11 to every observation. Add up the absolute
differences from that fitted level::

    abs(10 - 11) + abs(10.1 - 11) + ... + abs(12 - 11) = 8

Now fit two plateaus, splitting before position 4. The first four values have
median 10 and the last four have median 12. The total absolute error is::

    First plateau:   0 + 0.1 + 0.1 + 0 = 0.2
    Second plateau:  0 + 0.1 + 0.1 + 0 = 0.2
    Total:                              0.4

The split improves the fit by 7.6. But better fit alone cannot decide the
answer: assigning each observation its own plateau would reduce the error
all the way to zero, including the tiny 0.1 millisecond fluctuations.

ASV therefore charges for each added boundary. Call that charge gamma. At
gamma = 1, the one-plateau fit costs 8, while the two-plateau fit costs
``0.4 + 1 = 1.4``. The two-plateau fit wins. At gamma = 10, it costs 10.4,
so the one-plateau fit wins. Gamma determines how much evidence a new
boundary must earn.

For gamma = 1, the exact solver returns::

    right endpoints:  [4, 8]
    fitted levels:    [10, 12]
    interval errors:  [0.2, 0.2]

An endpoint is exclusive: 4 means the first segment contains positions
0, 1, 2, and 3. The next segment starts at 4 and ends just before 8.

In normal use, ASV chooses gamma automatically. For this example, the
automatic path also returns these two plateaus. It converts each interval
error to an error per retained observation: ``0.2 / 4 = 0.05``. That number
describes residual scatter around the fitted plateau.

Finally, the reporting stage decides whether to call the step a regression.
At a 5 percent reporting threshold, the earlier level of 10 gives a threshold
of 0.5. The increase of 2 exceeds both that threshold and the residual scatter
of 0.05, so this example produces a regression between positions 3 and 4.
This stage receives fitted steps, not the original noisy measurements.

Try the whole example from the repository root:

.. code-block:: python

    import runpy

    # Load just this source file using its Python fallback.
    sd = runpy.run_path("asv/step_detect.py")
    y = [10, 10.1, 9.9, 10, 12, 12.1, 11.9, 12]

    fixed = sd["solve_potts"](y, [1] * len(y), gamma=1)
    print("Fixed penalty:", fixed)

    steps = sd["detect_steps"](y)
    print("Automatic steps:", steps)

    report = sd["detect_regressions"](steps, threshold=0.05)
    print("Regression report:", report)

The report is ``(12.0, 10.0, [(3, 4, 10.0, 12.0)])``: latest level 12,
best level 10, and a reported change from 10 to 12 between positions 3 and 4.
Printed errors may look like ``0.1999999999999993`` because of floating-point
arithmetic. In an installed ASV environment, ``from asv import step_detect``
is the normal entry point; the standalone example avoids installation and
deliberately uses the Python implementation.

Four questions answered by different parts of the code
------------------------------------------------------

It helps to separate four decisions that can otherwise look like one large
algorithm:

1. **What shape are we looking for?** Flat stretches with occasional jumps.
   This is the signal model. It expresses a belief about performance history.
2. **For a given price per jump, where should the boundaries go?**
   ``solve_potts`` answers this optimization question; normal detection uses
   its approximate wrapper, ``solve_potts_approx``, for speed.
3. **What price per jump should we use?** ``solve_potts_autogamma`` tries
   candidate prices and compares the resulting fits. This is model selection.
4. **Which fitted jumps matter to the user?** ``detect_regressions`` applies
   reporting rules, including the relative slowdown threshold and treatment
   of later recovery.

These stages have different notions of success. An exact optimizer can find
the best answer for a poorly chosen gamma. Automatic selection can choose a
plausible step that is too small to report. A faster C++ implementation can
compute the same answer sooner without making it statistically more reliable.

Words and symbols used below
----------------------------

* **Observation:** one benchmark value in the history, usually a summary of
  repeated timing measurements for a commit. It is not each individual repeat.
* **Segmentation or partition:** a division of that history into contiguous
  segments. These two words mean the same thing here.
* **Level:** the single fitted value assigned to a segment.
* **Residual:** measured value minus fitted level. At position 1 in the example,
  it is ``10.1 - 10 = 0.1``.
* **Cost, loss, or objective:** a numerical rule for ranking fits; smaller
  is better according to that rule. It is not automatically a probability.
* **n:** number of retained observations. **K:** number of segments.
  **J:** number of changes; normally ``J = K - 1``.
* **y[i]:** measured value. **x[i]:** fitted value. **w[i]:** influence of
  observation i in the fitting error. **mu:** a candidate segment level.
* **gamma:** charge per boundary in the fitting problem. **rho:** strength
  of the relationship between neighboring residuals in the noise score.

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

For example, a graph might have measurements at revision numbers
``[100, 103, 110, 111]``. The fit sees a list with positions ``[0, 1, 2, 3]``.
A split before position 2 is later mapped to a change somewhere between
revisions 103 and 110. Without intermediate measurements, the detector cannot
identify a unique responsible commit inside that range.

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
absolute deviation, so weights are proportional to inverse uncertainty scale,
not its square. The corresponding statistical model allows observations to
have different noise scales, a property called heteroscedasticity. The
procedure estimates overall residual noise from the history instead of taking
individual confidence intervals as a complete noise model.

For a numerical example, confidence interval widths of 0.2 and 2 give raw
weights of 10 and 1. The same absolute fitting error therefore costs ten times
as much at the more precise observation. Dividing all weights by their median
preserves that ratio while making the overall weight scale easier to handle.
The weights express confidence in measurements; they do not say which commits
are important or how large a regression would be acceptable.

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

The name Potts refers to the preference for neighboring states to agree.
For this application, the states are the fitted performance levels. You do
not need to simulate a physical system to solve it: ASV directly minimizes
the fitting cost. The essential feature is that a jump from 10 to 11 and a
jump from 10 to 100 pay the same boundary charge. The observations determine
whether either jump saves enough fitting error to be worthwhile.

The phrase L0 counts how many differences are nonzero. L1, used for the
measurement error, means summing absolute differences. These describe two
different parts of the same objective: L1 measures how far the fit is from
the data, and L0 counts how often the fitted level changes.

For a fixed interval, the level minimizing weighted absolute error is a
weighted median. This makes the fitted level less sensitive to extreme
measurements than a mean under squared error. Robustness is limited: an
extreme point can still justify a separate short segment when the reduction
in loss pays for the extra boundaries.

To see why the median appears, take one segment with values ``[10, 10, 100]``.
At level 10, its absolute error is ``0 + 0 + 90 = 90``. At the mean, 40, its
absolute error is ``30 + 30 + 60 = 120``. Moving the plateau toward the outlier
makes the absolute-error fit worse. A squared-error fit would instead choose
the mean because squaring gives the large error much more influence.

A weighted median uses influence rather than just counts. Sort observations
by value, add their weights, and find the point where the accumulated weight
reaches half the total. For values ``[10, 11, 20]`` with weights ``[1, 4, 1]``,
the total weight is 6 and the midpoint is 3. The cumulative weights are
``[1, 5, 6]``, so the median is 11. Its cost is
``1 * 1 + 4 * 0 + 1 * 9 = 10``. With a tie, an interval of levels can be
equally good; ASV can return the midpoint of two neighboring values.

For example, with unit weights::

    Measurements    10 10 10 10 | 12 12 12 12
    One plateau     11 11 11 11   11 11 11 11    loss = 8
    Two plateaus    10 10 10 10 | 12 12 12 12    loss = gamma

The two-plateau solution wins for gamma below 8; there is a tie at 8. More
generally, two noise-free adjacent plateaus of lengths a and b and separation
delta save ``min(a, b) * abs(delta)`` in absolute error when split. This
explains how duration and magnitude work together as evidence for a change.

That duration argument is useful for recent regressions. A small slowdown
seen for only one new commit has little accumulated evidence. If subsequent
commits stay slow, the same change can become easier to detect even though
its size has not grown. Automatic selection also changes with history length,
so this is intuition about the fixed-penalty problem, not a detection guarantee.

Exact optimization for a fixed penalty
--------------------------------------

There are many ways to split a history: each gap can either contain a boundary
or not. With n observations there are ``2 ** (n - 1)`` such choices if all
segment lengths are allowed. Trying every complete partition quickly becomes
impractical.

Dynamic programming avoids repeating the same work. Suppose we already know
the cheapest fit for every shorter prefix. To fit the first t observations,
we only need to try each possible start of the final segment. Everything
before that start has already been solved.

Use half-open intervals and define::

    C(s, t) = min_mu sum_{i=s}^{t-1} w[i] * abs(y[i] - mu)

    F(0) = -gamma
    F(t) = min_{0 <= s < t} [F(s) + C(s, t) + gamma]

The recurrence tries every possible start of the final segment. Each choice
combines the best fit of the prefix with one new segment. Saving the winning
start at each endpoint permits backtracking to recover the full partition.
The initial value avoids charging for the first segment.

Here ``C(s, t)`` means the cost of fitting just one plateau to positions s
through t minus one. ``F(t)`` means the cheapest cost of fitting the first t
observations with any allowed number of plateaus. The minimum asks us to try
all eligible s values and keep the cheapest choice.

For a complete numerical example, use ``y = [1, 1, 3, 3]``, unit weights, and
gamma = 1. The shorter-prefix answers are::

    F(0) = -1  initialization; cancels the first boundary charge
    F(1) =  0  [1]
    F(2) =  0  [1, 1]
    F(3) =  1  [1, 1] | [3]

To compute F(4), try each start of its last segment::

    Last start s    Prefix F(s)    Last cost C(s,4)    + gamma    Total
         0               -1                 4             1        4
         1                0                 2             1        3
         2                0                 0             1        1
         3                1                 0             1        2

The winning start is 2. The last plateau is ``[3, 3]``, and the best prefix
ending there is ``[1, 1]``. Remembering those winning starts recovers the
partition ``[1, 1] | [3, 3]`` with total cost 1. This recovery step is called
backtracking. The algorithm is exact for the stated cost and allowed segments;
that does not imply it knows the true history that generated noisy data.

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

Think of ``RangeMedian`` as an object answering two repeated questions:
what is the best flat level for this interval, and how much error remains?
For the first four observations in the opening example, ``mu(0, 3)`` is 10
and ``dist(0, 3)`` is approximately 0.2. A cache saves those answers so later
questions about the same interval do not need another sort. The answers do
not depend on gamma, which is why different penalty trials can share them.

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

In this notation, O(n^2) describes work that grows roughly with the square
of history length, ignoring constant factors and smaller terms. Doubling n
then suggests about four times as much work. C++ can make each operation
cheaper, but it does not remove the growth in the number or size of intervals.
That is why ASV combines a compiled inner loop with an approximation.

The approximation used in normal detection
------------------------------------------

``solve_potts_approx`` initially restricts segments to at most 20 observations
when ``min_size < 10``. Otherwise it uses ``min_size + 50``. The resulting
partition can have artificial boundaries solely because of this restriction.

For example, a perfectly flat history of 100 observations cannot initially
be represented by one interval when intervals are limited to 20. The first
pass must divide it into at least five pieces. The merging stage can then
join those pieces: there is no extra fitting error, and each removed boundary
saves gamma. This is how the method can discover long plateaus despite its
small initial window.

``merge_pieces`` repeatedly considers neighboring segments and removes a
boundary when the merged interval costs no more than the separate intervals
plus gamma. Its scan uses an early-exit heuristic, so this is not a global
search over all possible merges. It then makes one pass over surviving
boundaries, trying offsets within the initial maximum segment length, and
recomputes the segment levels and costs.

The merge decision can be written in ordinary arithmetic. Suppose neighboring
segments have fitting errors 2 and 3, and gamma is 1. Keeping the boundary
costs ``2 + 3 + 1 = 6`` locally. If fitting their union costs 5.5, merging
saves 0.5. If it costs 7, the boundary is kept. After merging, moving a surviving
boundary a few positions left or right can further reduce fitting error.

These operations improve the current partition locally. They cannot freely
replace it with every other possible partition. A better result might require
several coordinated changes that no individual merge or boundary shift finds.
That is the reason to compare this approximation with the exact solver at the
same gamma, as the improvement plan proposes.

The final segments can be much longer than 20 observations. The fixed window
limits the initial dynamic program, not the output plateau length. Merging
and boundary adjustment also evaluate long intervals, so the docstring's
linear-time description is not a demonstrated worst-case bound for the
whole procedure. The initial short-window pass does explain its practical
speed advantage.

Choosing the penalty and estimating noise
-----------------------------------------

The fixed-penalty solver answers a conditional question: if a new boundary
costs gamma, what is the cheapest partition? It does not answer what gamma
should be. ASV wraps that solver in a second decision process::

    Try a gamma
        -> get one fitted history
        -> inspect its residuals and number of segments
        -> assign that history a model-selection score
    Compare the histories found at different gammas
        -> keep the lowest-scoring history encountered

There are therefore two costs to keep separate. The inner fitting cost
contains ``gamma * number_of_jumps``. The outer selection score below uses a
different complexity penalty and a noise estimate. Comparing inner costs at
different gammas would change the measuring rule each time; the outer score
provides a common rule for comparing their fitted histories.

``solve_potts_autogamma`` scales candidate penalties by the absolute-error
cost of fitting the entire series with one plateau, using 1 if that cost is
zero. It searches on a logarithmic penalty scale. Each candidate produces
a segmentation through the approximate solver.

Searching on a logarithmic scale means varying gamma by ratios instead of
fixed additions, like trying 0.1, 1, and 10 rather than 1, 2, and 3. Those are
illustrative values, not the actual sequence tried by the implementation.
This helps explore penalties over different orders of magnitude. Scaling
them by the one-plateau error also adapts the search to the size of the data.

For a candidate fit, let ``e[i] = y[i] - fitted_level[i]``. The score uses::

    S(rho) = w[0] * abs(e[0])
             + sum_{i=1}^{n-1} w[i] * abs(e[i] - rho * e[i - 1])

    score = (4 * log(n) / n) * K + log(sigma_0 + S(rho))

Here K is the segment count. The correlation search first evaluates rho at
-1 and 1, but ``expand_bounds=True`` expands its bracket to about
[-4.236, 4.236]; these initial values are not hard limits. Residual correlation
continues across fitted segment boundaries. Positive correlation lets a
sequence of similar residuals count as less surprising than independent
errors would suggest.

The `mathematical analysis <mathematical_analysis.rst>`_ derives an exact
weighted-median solution for rho and shows why this stopping condition
can return a nonoptimal value even for this convex subproblem. The current
branch now implements that correction; see the
`correlation-fit design <correlation_design.rst>`_. This reference's search
description follows its original pinned revision.

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

Why the score contains a logarithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Bayesian discussion in the original documentation is motivation for the
selection rule; it is not necessary to implement the fixed-penalty solver.
A shorter route to understanding the logarithm is to start with independent
Laplace noise. Laplace noise is a probability model in which large absolute
errors become exponentially less likely. Allowing a different uncertainty
at each observation gives a term proportional to ``exp(-w[i] * abs(e[i]) / b)``,
where b is the unknown overall noise scale.

With fixed weights, the negative log likelihood, ignoring terms that do not
change between candidate fits, is::

    n * log(b) + D / b
    where D = sum_i w[i] * abs(e[i])

For positive D, the best noise-scale estimate is ``b = D / n``. Substituting
it gives ``n * log(D / n) + n``. Add a complexity charge of
``c * K * log(n)``, divide by n, and drop constants shared by all candidates::

    log(D) + c * K * log(n) / n

This explains the shape of ASV's score. The implementation uses c = 4,
replaces D with a residual sum adjusted for correlation, and adds a noise
floor. This simplified derivation does not establish that 4 is the uniquely
correct coefficient. It also does not account rigorously for searching over
unknown boundary locations, which is part of the difficulty discussed in the
original documentation.

What correlated noise means
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Independent noise means the last measurement's error does not help predict
the next error. In practice, background machine load might make several
neighboring measurements all a little slow. ASV uses rho to ask how much of
today's residual can be predicted from the preceding residual. AR(1) is the
name for this model with one previous value.

For residuals ``[1, 1, 1, 1]`` and unit weights, rho = 0 gives a sum of 4.
At rho = 0.8, the first residual still costs 1, but each later residual costs
``abs(1 - 0.8 * 1) = 0.2``. The total becomes 1.6. Thus a persistent residual
pattern can be explained partly by correlated noise rather than extra steps.
This is an illustration of the residual score, not a full fitted dataset.

That flexibility is useful but ambiguous: a real performance change also
creates a persistent pattern if the fitted history misses it. The complexity
penalty, candidate fits, and noise score must therefore be evaluated together.
Changing the correlation score can select a different candidate, but it cannot
create a boundary that none of the candidate fits contains.

Why a noise floor is needed
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If every observation gets its own segment, all residuals vanish. Without a
floor, the score contains ``log(0)``, which tends to negative infinity. No
finite charge for extra segments can compete with that. The floor makes the
reward for a perfect fit finite.

The current floor also depends on the fitted levels. With two plateaus at
1 and 1.02, it is ``0.001 * 1 = 0.001``. Shift both plateaus upward by 1000,
and it becomes ``0.001 * 1001 = 1.001`` even though their separation remains
0.02. Now eliminating that small residual error matters much less inside the
logarithm. This explains why the baseline-shift example in the improvement
plan is worth testing rather than assuming automatic selection is invariant
to adding a constant.

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

Two examples show why finding steps and reporting regressions are separate.
Assume each plateau has several observations and negligible residual error::

    Fitted levels         What the history says
    10 -> 12 -> 10        A slowdown occurred, then fully recovered
    10 -> 14 -> 12        A slowdown occurred, then partly recovered

At a 5 percent threshold, the first history produces no remaining regression
with these assumptions. The second reports a lasting change from 10 to 12,
located at the original rise to 14. The peak slowdown and the best subsequent
performance answer different questions.

Changing the reporting threshold does not refit the plateaus. It changes
which fitted changes are worth reporting. Gamma influences segmentation;
the relative threshold influences reporting. They are not interchangeable
sensitivity controls.

What each function gives to the next
------------------------------------

Use this compact map when navigating the source:

* ``get_weight(stats)`` gives publishing a relative weight from timing
  uncertainty. ``Graph.get_data()`` gives detection ordered revision, value,
  and weight triples.
* ``detect_steps(y, w)`` filters inputs and asks ``solve_potts_autogamma``
  for a fitted history. It eventually returns position bounds, fitted levels,
  observed minima, and residual error per observation.
* ``solve_potts_autogamma(y, w)`` asks ``solve_potts_approx`` for fits at
  different gammas. It returns exclusive endpoints, levels, interval error
  sums, and the selected gamma.
* ``solve_potts_approx(y, w, gamma)`` asks ``solve_potts`` for a restricted
  partition, then passes it to ``merge_pieces`` for local improvements.
* ``solve_potts(y, w, gamma)`` asks the interval-cost object for costs and
  medians. With the extension available, it also delegates the dynamic
  programming loop to that object. Its output is three parallel lists:
  exclusive endpoints, levels, and interval error sums.
* ``graph._compute_graph_steps`` replaces observation bounds with revision
  bounds. The plugin then calls ``detect_regressions`` with fitted segments
  and a reporting threshold to obtain the user-facing regression records.

The median and interval-cost calculations are reused throughout, but each
outer layer answers a broader question than the layer it calls.

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
