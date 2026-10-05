Step detection improvement plan
===============================

Improve ASV step detection by measuring three things separately: whether the
optimizer minimizes its stated objective, whether the selected segmentation
captures real changes, and whether the reporting policy produces useful
regression alerts. A faster optimizer can still select a poor noise model;
a lower fitting objective does not by itself establish better alerts.

The `implementation walkthrough <README.rst>`_ describes the baseline at
``d33754e129c025beb5c2ca440c3c280433b264f7``. This plan proposes experiments;
it does not change the production detector or its defaults.

Initial evidence and open questions
-----------------------------------

Two small reproductions establish behaviors worth measuring:

* At a fixed penalty, the approximate solver can have a larger objective
  than the exact solver. A 70-observation example below yields three segments
  and cost 16.52, versus one segment and cost 16.05 for the exact solver.
* Adding a constant baseline can change automatic segmentation. Two constant
  runs at 1 and 1.02 are separated, while the same runs at 1001 and 1001.02
  are merged. The automatic noise floor depends on absolute fitted levels.

These examples establish possibility, not frequency or practical severity.
The first is expected from an approximation. The second requires a modeling
decision: baseline dependence may be appropriate for relative measurement
accuracy, but should be explicit and calibrated.

Source inspection raises further questions to evaluate:

* Does golden-section search miss useful penalty regions on its discontinuous,
  potentially non-unimodal score? No failing search example is established here.
* How much time and memory are spent sorting long merged intervals and caching
  costs across penalties?
* Does scoring correlation only after independent-noise fitting miss changes
  or interpret machine drift as changes?
* How do missing observations, revision gaps, and short recent plateaus affect
  the final regression records?

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

Phase 1 establishes an evaluation baseline
------------------------------------------

Build a small harness in this directory before changing the detector. Record
the source revision, interpreter, backend, random seed, input arrays, all
parameters, fitted boundaries and levels, objective, selected gamma, runtime,
and peak memory. Save generated inputs as well as seeds so examples survive
changes in generators. Keep large generated results out of version control.

Use four evaluation layers:

1. Verify interval medians and costs against direct calculations on small
   weighted examples, including ties and repeated values.
2. Enumerate all partitions of very short series as an independent oracle for
   ``solve_potts``. Compare objectives and constraint satisfaction, allowing
   equivalent partitions when optima tie. Verify restricted positions and
   minimum and maximum segment sizes separately.
3. Compare approximate and exact solutions at identical fixed penalties on
   short and medium histories. Report absolute and normalized objective gaps,
   extra or missing boundaries, runtime, and memory.
4. Evaluate automatic selection and final regression reporting separately
   against known changes and desired reporting outcomes.

Run both the Python fallback and compiled extension. Agreement on objectives
is more important than identical tie-breaking. Instrument interval requests,
cache hits, interval lengths, candidate counts, and time spent in fitting,
merging, correlation fitting, and penalty search. Measure cold runs as well
as reuse across penalties inside one automatic fit.

The synthetic suite should cross several dimensions without requiring a
full Cartesian product:

* Histories with no changes, one change, several changes, short dips, temporary
  slowdowns followed by recovery, and small changes near the newest revision.
* Gaussian and Laplace noise, occasional extreme values, unequal uncertainty,
  AR(1) correlation of both signs, and shifts in machine background load.
* Constant data, repeated or quantized values, zero and near-zero levels,
  unequal plateau lengths, and positive timing series over several scales.
* Missing results and weights, nonpositive weights, sparse revision sampling,
  and unequal gaps in observation order and measurement time.
* Slowly drifting signals, where a step model is deliberately misspecified.

Use the exact oracle only where tractable. For speed and memory scaling,
start at 100, 1000, and 10000 observations, then increase the size of methods
that remain practical. Include long nearly constant histories as well as
frequent changes; the current repeated two-level benchmark covers only a
narrow shape of input.

For statistical quality, report false boundaries per history, probability
of any false boundary, missed changes, and localization error with a stated
matching tolerance. Separately report false and missed regression alerts,
including the expected treatment of recovered slowdowns. Sweep step size and
available observations after a change to measure sensitivity near the latest
commit. Use multiple seeds and uncertainty intervals for aggregate rates.

Acceptance for this phase is reproducibility of the two examples, agreement
between the exact solver and the independent oracle on valid small cases,
and a saved baseline report. Any oracle disagreement must be understood before
using that solver to judge replacements.

Phase 2 evaluates penalty search
--------------------------------

Keep the current segment cost, approximate solver, correlation score, and
noise floor fixed. Compare the existing search against a logarithmic grid
with refinement around changes in the returned partition. Start with the same
effective search range and an explicit budget of solver calls; evaluate wider
ranges as a separate experiment. Deduplicate identical partitions before
repeating correlation scoring.

On short histories, use a dense penalty grid as a diagnostic baseline. Record
each candidate partition and its score so plots show the plateaus and jumps
in the search objective. A finite dense grid is still an approximation and
must not be labeled a global oracle. The current search should remain one
candidate in a hybrid experiment, making its best observed score available.

Measure score improvement, changes in detection quality, solver-call count,
and wall time. Lower scores alone do not justify a default change. Prefer a
candidate only if it improves the predeclared detection metrics at an
acceptable runtime cost on held-out seeds and histories.

If an exact fixed-penalty solver is practical, evaluate
`CROPS <https://arxiv.org/abs/1412.3617>`_, which finds optimal segmentations
across a continuous penalty range. Its guarantees depend on exact solutions
to the underlying penalized problem. It neither guarantees every possible
segment count appears on that path nor makes ASV's distinct correlation score
globally optimal over all segmentations.

Phase 3 evaluates faster exact fitting
--------------------------------------

Profile interval evaluation before choosing a data structure. Candidate
approaches include maintaining weighted order statistics while an interval
grows, with cumulative weights and weighted value sums to obtain both its
median and absolute-error cost. Compare that against the current sort-and-cache
method at the actual interval sizes requested. Short windows may favor the
simpler implementation despite worse asymptotic behavior.

Evaluate `PELT <https://arxiv.org/abs/1101.1438>`_ as an exact boundary-pruning
algorithm. For positive weights and independent segment levels, splitting
an interval cannot increase its minimized absolute-error cost::

    C(a, b) + C(b, c) <= C(a, c)

This supplies the zero-constant cost inequality relevant to PELT pruning.
Check the complete pruning rule, indexing, ties, and segment constraints
before implementation. PELT's favorable linear scaling requires conditions
on the problem; long segments and expensive weighted-median queries can
still make it costly. Report actual candidate counts and interval-query
costs rather than inferring total linear runtime from the algorithm's name.

The `L1 Potts paper <https://arxiv.org/abs/1207.4642>`_ also describes an exact
quadratic-time, linear-space algorithm for its discrete problem. Review its
data structure and whether its guarantees extend to ASV's arbitrary weights
before treating it as an implementation option.

Acceptance requires agreement with the verified objective oracle, valid
constraints, and measured runtime and memory across both sparse and frequent
change regimes. Evaluate compiled and fallback behavior. Keep any new solver
selectable in experiments until its practical limits and tie behavior are
understood.

Phase 4 calibrates the noise model
----------------------------------

First decide which invariances the model should have. The fixed-penalty
absolute-error fit is invariant to adding a constant to all measurements.
Multiplying values by a positive factor preserves boundaries when gamma is
scaled by the same factor; scaling all weights also requires scaling gamma.
Automatic selection and relative regression thresholds introduce additional
scale choices. Test these layers separately instead of imposing one invariant
on the entire pipeline.

Compare the current noise floor with a floor shared across candidate
partitions, derived from a robust residual-scale estimate or measured
uncertainties. Candidate-derived floors change both the noise model and the
relative scores of partitions. Pay particular attention to the switch between
two and three segments and to exactly fitted or quantized data. A principled
prior on noise scale is another option, but its assumptions and consequences
must be explicit.

For strictly positive timing data, test fitting logarithms when relative
changes are the desired target. Transform uncertainty weights consistently;
do not reuse original-scale weights without justification. Establish behavior
for zero values and non-timing benchmarks before proposing a general default.

For correlated noise, first evaluate the existing scoring-only treatment
against independent noise. A jointly fitted correlated objective couples
residuals across boundaries and cannot simply replace the independent segment
cost in the existing dynamic program. Treat that as a separate algorithmic
project. Consider recorded run order or batch information when available,
because commit order need not represent the temporal correlation mechanism.

Choose model settings on training simulations and evaluate on held-out seeds
and representative histories. Report false alerts and sensitivity together.
Do not describe the selected information score as a p-value or posterior
probability of regression.

Phase 5 validates reporting and prepares integration
----------------------------------------------------

Replay candidate fits through graph coordinate mapping and regression
postprocessing. Exercise missing revisions, recovered slowdowns, brief fast
dips, configured history filters, and relative reporting thresholds. Record
changes to both fitted steps and user-visible regression records. A deliberate
change in reporting semantics needs its own rationale.

Before proposing a default change, agree on numerical quality and resource
budgets using the baseline report. Preserve raw comparison results so a
different tradeoff can be evaluated later. A speed-only change should preserve
the intended objective and reporting semantics; a statistical change should
show the sensitivity versus false-alert tradeoff explicitly.

Land improvements in reviewable units: evaluation harness, documented behavior
and focused correctness fixes, penalty search, interval-cost improvements,
then larger solver or model changes. Run the repository's relevant step,
graph, and publishing tests for implementation changes, plus broader required
checks before integration. Keep the existing implementation available for
controlled comparisons during the experiment.

Concrete first deliverable
--------------------------

The first implementation should produce the evaluation harness and baseline
report from phase 1, followed by the bounded grid-search comparison in phase
2. This tests a small, isolated change before investing in a new C++ solver.
The report should identify which bottleneck dominates and whether a better
optimizer or a better model-selection rule is likely to improve actual alerts.
