Step detection improvement plan
===============================

The goal is to detect lasting performance changes more reliably while keeping
analysis practical for long benchmark histories. We should evaluate three
layers separately: finding the best fit at a fixed penalty, selecting the
penalty and noise model, and reporting useful regressions.

The `conceptual guide <README.rst>`_ explains why those layers exist. The
`implementation reference <implementation_details.rst>`_ maps them to the
code at revision ``d33754e129c025beb5c2ca440c3c280433b264f7``.

Start with a reproducible evaluation harness, then compare penalty searches.
Those experiments will show whether a larger investment in the solver or
noise model is justified.

The initial harness, bounded search comparison, and scaling sample are now
implemented. The `baseline report <baseline_report.rst>`_ records the findings;
`running the experiments <running.rst>`_ explains how to reproduce them.
The broader corpus, internal cache profiling, and replacement solver experiments
below remain follow-up work.

.. contents:: On this page
   :local:
   :depth: 1

What we already know
--------------------

Two small examples identify behavior worth measuring:

* At gamma = 1, the approximate solver can select three segments with cost
  16.52, while the exact solver selects one with cost 16.05. The approximation
  reduces residual error to 14.52 but pays 2 for its extra boundaries.
* Two constant runs at 1 and 1.02 are separated by automatic selection, while
  the same runs at 1001 and 1001.02 are merged. The current noise floor depends
  on absolute fitted levels.

The first example establishes an optimization gap. Its generated history
actually contains a change halfway through, so the exact solver's lower
objective does not establish better recovery of the true history. We need
both objective comparisons and detection-quality measurements.

The second example raises a modeling choice: adding a baseline preserves the
absolute change but reduces its relative size. Decide which behavior is
intended before changing the floor.

Source inspection also suggests investigating penalty search and interval
costs. The outer score need not have the single minimum shape assumed by
golden-section search, and the C++ backend repeatedly sorts uncached intervals.
Their practical impact remains to be measured.

Reproduce the initial observations
----------------------------------

Run this block with Python from the repository root. ``runpy`` deliberately
loads the standalone module without a package context, selecting its Python
fallback and avoiding installation dependencies. This reproduces algorithmic
behavior; it is not a C++ performance measurement. Expected results below
refer to the pinned source revision.

.. code-block:: python

    import math
    import random
    import runpy

    sd = runpy.run_path("asv/step_detect.py")
    assert sd["_rangemedian"] is None

    rng = random.Random(7)
    y = [
        round(rng.gauss(1 if i < 35 else 1.3, 0.3), 2)
        for i in range(70)
    ]
    w = [1] * len(y)
    gamma = 1.0

    exact = sd["solve_potts"](y, w, gamma)
    approx = sd["solve_potts_approx"](y, w, gamma)

    def objective(solution):
        right, values, distances = solution
        return sum(distances) + gamma * (len(right) - 1)

    print("exact:", exact[0], objective(exact))
    print("approximate:", approx[0], objective(approx))
    assert exact[0] == [70]
    assert approx[0] == [16, 19, 70]
    assert math.isclose(objective(exact), 16.05)
    assert math.isclose(objective(approx), 16.52)

    boundaries = []
    for offset in [0, 1000]:
        y = [offset + 1] * 30 + [offset + 1.02] * 30
        result = sd["solve_potts_autogamma"](y, [1] * len(y))
        boundaries.append(result[0])
        print("offset:", offset, "boundaries:", result[0])
    assert boundaries == [[30, 60], [60]]

Phase 1 measures the baseline
-----------------------------

Build a harness that runs identical inputs through detector variants and saves
comparable outputs. Record the source revision, interpreter, backend, seed,
input arrays, parameters, fitted levels and boundaries, selected gamma,
objective, runtime, and peak memory. Save generated inputs as well as seeds.

Use three kinds of evidence:

* **Calculation correctness:** compare interval costs with direct calculations;
  enumerate every partition of tiny histories to independently verify the
  exact solver. Include weights, ties, restricted positions, and segment-size
  constraints. Compare objectives when different partitions tie.
* **Optimization quality:** compare exact and approximate solvers at the same
  gamma on tractable histories. Record objective gaps and boundary differences.
* **Detection quality:** use histories with known changes to measure missed
  changes, false boundaries, and location error. Evaluate final regression
  alerts separately from fitted boundaries.

Include flat histories, isolated steps, several changes, recent changes,
short dips, recovered slowdowns, unequal plateau lengths, and gradual drift.
Vary noise magnitude, outliers, unequal uncertainty, correlation, quantization,
missing observations, and revision gaps. Run multiple seeds, including held-out
seeds for evaluating settings selected during development.

Predeclare a boundary-matching tolerance, match each true and detected boundary
at most once, and report actual localization error alongside match rates.
Measure both false boundaries per history and the probability of any false
boundary. Report uncertainty intervals for aggregate detection rates.

Run Python and C++ backends. Start performance measurements at 100, 1000, and
10000 observations, using the exact oracle only where tractable. Include long
flat histories as well as frequent changes. Record interval lengths, candidate
counts, cache hits, and time spent in fitting, merging, correlation scoring,
and penalty search to identify the main costs.

The deliverable is a reproducible baseline report. Resolve disagreements with
the tiny-history oracle before using the exact solver to judge replacements.

Phase 2 compares penalty searches
---------------------------------

Keep the candidate fitter, correlation score, and noise floor fixed. Compare
the existing golden-section search with a logarithmic grid, then refine near
penalties that produce different partitions. Use the same effective range and
an explicit solver-call budget; test wider ranges separately.

This isolates the question of whether search overlooks useful candidate fits.
Save the sampled penalties, partitions, and scores so their relationship is
visible. Deduplicate repeated partitions before correlation scoring, and retain
the current search's best candidate in a hybrid comparison.

A dense grid provides a useful diagnostic for short histories, though narrow
regions can still fall between samples. Measure score improvement, detection
quality, solver calls, and wall time. Favor a search method only if its gains
on held-out histories justify its computational cost.

Once an exact fixed-penalty solver is practical, consider
`CROPS <https://arxiv.org/abs/1412.3617>`_ to explore optimal segmentations across
a continuous penalty range. Its guarantee concerns the underlying penalized
objective. ASV's correlation-adjusted selection score remains a separate
criterion evaluated on those candidate fits.

Phase 3 improves fixed penalty fitting
--------------------------------------

There are two costs to attack: the number of intervals considered and the
work needed to evaluate each interval. Measure changes to each separately
before combining them.

For interval costs, compare sorting each uncached interval with maintaining
weighted order statistics as an interval grows. Cumulative weights locate
the median; weighted value sums can help calculate absolute-error cost.
Benchmark at the interval sizes actually requested, since sorting small
windows may be faster despite worse asymptotic scaling. Measure memory and
cache reuse across penalties as well as cold-run time.

For candidate search, evaluate `PELT <https://arxiv.org/abs/1101.1438>`_. It
prunes starts of segments that can be proven unable to win, preserving an
exact solution under its conditions. Weighted absolute-error costs with
positive weights satisfy the relevant splitting inequality::

    C(a, b) + C(b, c) <= C(a, c)

Here C is the minimum error of fitting one level to an interval. Allowing
separate levels for two pieces cannot make the fit worse. Check the full
pruning rule, constraints, and tie handling against the oracle. PELT's
favorable linear scaling is conditional, and interval evaluation still
contributes to total runtime.

The `L1 Potts paper <https://arxiv.org/abs/1207.4642>`_ describes an exact
quadratic-time, linear-space algorithm for its discrete problem. Review whether
its data structure and guarantees extend to ASV's arbitrary weights.

Accept solver changes when they match the verified objective, respect segment
constraints, and meet measured runtime and memory budgets across both sparse
and frequent change regimes. Keep experimental solvers selectable while their
practical limits are being established.

Phase 4 evaluates the noise assumptions
---------------------------------------

Choose the intended invariances first. Converting seconds to milliseconds
should preserve the interpretation of a history. Adding a constant baseline
preserves absolute changes but alters percentage changes. The fixed-penalty
fit, automatic selection, and reporting threshold have different roles in
these transformations and should be tested separately.

Compare the current noise floor with a floor shared across candidate fits,
using a robust noise estimate or measured uncertainties. Score the same saved
candidates first, then rerun the full search to measure its interaction with
the changed score. Include perfectly fitted data, quantized readings, near-zero
levels, and the switch in the current formula between two and three segments.

For strictly positive timing data, explore fitting logarithms if relative
changes are the target. Changes from 10 to 12 and from 100 to 120 then have the
same size, log(1.2). Transform uncertainty weights consistently and establish
behavior for zero-valued and non-timing benchmarks before generalizing it.

Compare independent-noise scoring with the current correlation adjustment.
A jointly fitted correlated model would couple residuals across boundaries
and needs a separate algorithmic design. Consider run order or batch metadata
where available: measurement-time correlation may differ from commit order.

Select settings on training simulations and evaluate sensitivity and false
alerts on held-out histories. The useful result is a measured tradeoff between
missed changes and false detections, rather than merely a lower selection score.

Phase 5 checks reporting and integration
----------------------------------------

Replay promising fits through graph coordinate mapping and regression
postprocessing. Cover missing revisions, brief fast dips, recovered slowdowns,
configured history filters, and reporting thresholds. Compare the final
regression records as well as the fitted steps.

For example, ``10 -> 12 -> 10`` fully recovers, while ``10 -> 14 -> 12``
retains a slowdown. A boundary between sparsely measured revisions should
remain a range when the evidence cannot identify one responsible commit.

Use the baseline report to set numerical quality and resource budgets before
choosing a default. A speed improvement should preserve the intended fitting
and reporting behavior; a statistical change should show its sensitivity
versus false-alert tradeoff.

Land changes in reviewable units: evaluation harness, focused correctness
fixes, penalty search, interval costs, then larger solver or model changes.
Run the relevant step-detection, graph, and publishing tests for implementation
changes, plus the required repository checks before integration. Preserve
comparison results and the baseline implementation during the experiments.

First deliverable
-----------------

Implement the phase 1 harness and produce its baseline report, then run the
bounded grid-search comparison from phase 2. Together these should establish
whether penalty selection, approximation error, or interval evaluation is the
most useful next target.
