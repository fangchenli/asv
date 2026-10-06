An exact reference for the shared-floor score
=============================================

``exact_reference.fit_independent`` finds the globally best partition for
the independent-noise shared-floor score, given a positive noise floor and
a nonnegative complexity penalty. It returns the chosen fit and the best
fit at every feasible segment count. This gives us an inspectable reference
before changing candidate searches or comparing noise models.

The implementation uses dynamic programming over segment counts. It does
not choose gamma or assume that ASV's outer search visits enough candidate
fits. It reuses either the Python or C++ weighted interval-cost backend;
the partition search is new and is checked independently by enumerating
every partition of tiny inputs.

Why one fit per segment count is enough
---------------------------------------

For a partition P with K segments, let D(P) be its minimum weighted absolute
error. The declared score is::

    Q(P) = beta*K + g(D(P)/(n*b_min))

Here g(t)=t below one and g(t)=1+log(t) above one, as derived in the
`shared-floor analysis <shared_noise_floor.rst>`_. For fixed K, beta*K
is constant and g increases with D. Therefore the partition with the least
error at that K also has the least score at that K.

Define D_K as that minimum error. Finding every feasible D_K and selecting
the smallest ``beta*K + g(D_K/(n*b_min))`` gives a global optimum over the
declared family of partitions. The collection of (K,D_K) values is the
segment-count frontier returned by the implementation.

For example, take unit weights and measurements [10,10,12,12], with
b_min=0.5 and beta=0.6:

.. list-table::
   :header-rows: 1

   * - Segment count K
     - Minimum error D_K
     - Score
   * - 1
     - 4
     - 0.6 + 1 + log(2), about 2.293
   * - 2
     - 0
     - 1.2
   * - 3
     - 0
     - 1.8
   * - 4
     - 0
     - 2.4

One level costs four units of error. Two levels remove that error. Additional
segments cannot improve the fit but still pay their complexity charge, so
the two-segment solution wins. The floor and beta here illustrate the score;
they are not recommended defaults.

The dynamic program
-------------------

Let C(s,t) be the best weighted absolute error on observations s through
t-1. A weighted median supplies the optimal level for that interval.
Let E(k,t) be the minimum error for the first t observations in exactly
k segments. Initialize and update it as follows::

    E(0,0) = 0
    E(0,t) = infinity for t > 0

    E(k,t) = min over legal s of [E(k-1,s) + C(s,t)]

Every k-segment partition has some last segment [s,t). Replacing its prefix
with an optimal (k-1)-segment prefix cannot increase the error. Conversely,
each legal transition joins a valid prefix to a valid final segment. These
two facts prove the recurrence by induction. D_K is E(K,n).

For minimum length m and maximum length M, legal transitions require both::

    m <= t-s <= M
    (k-1)*m <= s <= (k-1)*M

Unreachable prefix states are skipped. This supports minimum and maximum
segment lengths without inserting invalid intermediate fits. A set of
constraints may have no feasible partition; for example, five observations
cannot be partitioned into segments of exactly three observations each.
The implementation raises an error for that case.

Each winning transition stores its last start s. Backtracking reconstructs
the endpoints, levels, and interval errors at every feasible K. Equal-error
transitions retain the earliest last start at each state; prefix ties follow
the same rule recursively. Equal final scores prefer fewer segments. Interval
level ties use the selected backend's weighted-median convention.

K counts partition intervals, including adjacent intervals with equal levels
if length constraints force them. With unrestricted maximum length and a
positive beta, such redundant boundaries cannot improve the selected score.

Inputs and output
-----------------

From Python at the repository root::

    from experiments.step_detection.exact_reference import fit_independent

    result = fit_independent(
        y=[10, 10, 12, 12],
        w=[1, 1, 1, 1],
        noise_floor=0.5,
        beta=0.6,
        backend="python",
    )

    result["selected"]["right"]       # [2, 4]
    result["selected"]["values"]      # [10, 12]
    result["selected"]["score"]       # 1.2
    result["frontier"]                # One fit for each feasible K

Each fit includes ``segments``, exclusive endpoints ``right``, segment
``values``, interval ``costs``, total ``error_sum``, and ``score``. The result
also records the floor, beta, backend, input size, and length constraints.
Use ``backend="native"`` to select the locally built C++ interval backend.

Inputs must be nonempty finite observations with equally many positive finite
weights. The solver accepts prepared sequences; it does not remove missing
values, normalize weights, map revision coordinates, or apply regression
reporting. The caller's floor must use the same weight scale as its input.
``harness.prepare`` remains available for obtaining ASV's prepared observations
and their original index mapping.

Once the frontier is available, another floor or beta can rescore its D_K
values without rerunning the partition search. Changing the weights or segment
constraints requires a new frontier. Correlated rescoring is possible, but
does not make these independent-error candidates globally optimal for a
correlated model.

Cost and numerical scope
------------------------

There are O(n^2) prefix/count states, each considering up to n final starts.
The partition dynamic program therefore takes O(n^3) time. Storing parents,
interval costs, and all reconstructed fits takes O(n^2) space; the forward
cost calculation uses two rows.

Both current interval backends sort values for uncached intervals. Across
all intervals, preprocessing has an O(n^3 log n) upper bound. This reference
prioritizes an independently checkable result over a fast implementation.
It rejects inputs longer than 200 observations to make that scope explicit.

The optimization is exact in the combinatorial sense. Its costs and scores
use floating-point arithmetic, which can affect very close comparisons or
ties. Nonfinite interval costs, accumulated costs, or scores are rejected. If a
rescaling is needed, scale the floor consistently with the data or weights.
The test oracle uses rational arithmetic to check the selected objectives
on small integer-valued cases.

Verification
------------

The oracle enumerates all partitions and tries every observed value as each
segment's level. It therefore depends on neither ASV's weighted-median
routine nor this dynamic program. Tests check every length-four sequence
from {-1,0,2}, two weight patterns, four segment-length constraint choices,
and two floor/penalty pairs: 1296 cases on each of the two backends.

For each case, the tests verify every feasible segment count's error, returned
levels and interval costs, coverage, length constraints, and the globally
selected score. Other checks cover known steps, ties, single observations,
infeasible constraints, unit scaling, baseline shifts, weight scaling, and
invalid inputs. All 48 tests pass, including the 2592 exhaustive comparisons.
Each exhaustive case also checks its optimum's supporting penalty against
every feasible segment count and the sufficient penalty range derived below.
Ruff, whitespace, and reStructuredText parsing checks also pass.

This is analytical correctness validation. No new benchmark campaign or
production-default change is part of this implementation.

How this supports a faster solver
---------------------------------

The shared-floor tangent argument gives a supporting penalty for an optimum
with error D_star::

    gamma_star = beta * max(D_star, n*b_min)

Let D_one be the unrestricted one-level weighted error, even if the chosen
length constraints disallow that single interval. For any feasible partition,
each segment could use that same level; optimizing its levels only reduces
error. Hence D_star <= D_one, giving the sufficient penalty range::

    beta*n*b_min <= gamma_star <= beta*max(D_one, n*b_min)

This range does not depend on a guessed logarithmic search bracket. An exact
penalty-path implementation can use it, with consistent constraint and tie
handling, and compare its selected score against this reference. The full
frontier can also reveal segment counts that the path skips. Their absence
alone is not an error: the proof guarantees an optimum on the exact path,
not that every count must occur there.
