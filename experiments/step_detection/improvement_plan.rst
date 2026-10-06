Step detection improvement plan
===============================

ASV tries to tell lasting slowdowns apart from ordinary variation in timing
measurements. We want to find out whether it can do that more reliably and
whether the calculations can run faster.

Start with `the beginner guide <README.rst>`_ if benchmarks, fitted values,
or penalties are unfamiliar. The first section below explains the improvement
plan in everyday terms. The later sections are working notes for implementing
the experiments and contain more technical detail.

These are proposed experiments. The existing ASV detection behavior has not
been changed. The `implementation reference <implementation_details.rst>`_
describes the source code being studied at revision
``d33754e129c025beb5c2ca440c3c280433b264f7``.

.. contents:: On this page
   :local:
   :depth: 1

The plan in everyday terms
--------------------------

First we need examples where we know the answer. Imagine making up a history
of measurements ourselves. We decide that the usual timing is 10 milliseconds
for the first 50 versions and 12 milliseconds for the next 50. Then we add
small variations, such as 0.1 above or below the usual time, to make the
readings resemble real measurements.

Now we can run ASV on those readings and ask whether it finds the change we
put there. We should also make a history with a usual timing of 10 throughout,
again with small variations. On that second history, there is no lasting
change to find. A useful detector needs to handle both kinds of example.

We will measure three basic outcomes:

* Did it miss a slowdown that we deliberately put in the history?
* Did it invent a lasting change where we only added ordinary variation?
* When it found a change, how close was its location to the one we put there?

We will also measure how long the analysis takes and how much computer
memory it uses. Finding the same answers faster can still be valuable when
there are many benchmarks and many program versions.

Experiment 1 checks the current behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the current ASV code on many examples and save the results. Use easy
examples, hard examples, histories with no change, and histories with several
changes. Repeat them with different small variations so one lucky set of
readings does not decide the result.

This gives us something to compare changes against. Without it, we could
make the detector better on one example while making many other cases worse
without noticing.

Experiment 2 tries different change penalties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The beginner guide shows that the penalty for adding a step changes the
chosen explanation. ASV already tries several penalties automatically.
We want to check whether it sometimes overlooks a useful choice.

Keep the fitting calculation the same and try penalties more systematically.
Then check whether the resulting descriptions find more of our known changes
without inventing more false ones. Also count the extra calculation time:
trying a thousand penalties instead of ten may be too expensive for the
benefit it provides. Those numbers illustrate the tradeoff, not ASV's actual
number of trials.

Experiment 3 checks the shortcut used to find segments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a particular penalty, ASV normally uses a shortcut to find a good division
into segments. The shortcut can miss a division with a better score. We have
already found a small example where that happens.

On short histories, compare the shortcut with a slower calculation that
finds the smallest possible score. This tells us how often the shortcut
makes a difference. Next, try ways to make the more thorough calculation
faster, such as reusing partial calculations instead of repeating them.

Remember that the smallest score only means best under the chosen rule.
We must still check whether its chosen steps match the changes we put into
our examples.

Experiment 4 checks how ordinary variation is treated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some measurements vary independently. Others may all be a little slow during
a period when the computer is busy. ASV tries to account for that distinction.
We want to test whether its rules sometimes explain away a real slowdown,
or mistake background variation for a software change.

There is also an existing example where adding the same large number to
every timing changes the detected steps. This keeps the absolute increase
the same but makes the percentage increase smaller. We need to decide what
behavior is useful for benchmark reporting and test the rule against that
decision.

Experiment 5 checks the report users see
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, run the promising changes through the whole reporting process.
Check a lasting slowdown, a slowdown that fully recovers, and a slowdown
that only partly recovers. Also check histories with missing measurements,
where ASV can identify a range of possible versions rather than one exact
version responsible for the change.

A new method is worth adopting if the evidence shows a useful improvement,
such as fewer missed slowdowns at a similar false-alarm rate, or the same
answers with less computation. We should decide which tradeoff matters
before selecting a winner.

The rest of this document turns those five experiments into a technical
work plan. The beginner guide and this overview provide the main ideas;
the implementation reference supplies the mathematical and programming details.

What we are trying to improve
-----------------------------

Think of ASV as trying to explain a noisy picture of performance history.
There are three different ways the result can disappoint us:

* It can do a poor job solving the mathematical problem we gave it. Another
  partition may have a lower cost at the same gamma. This is an optimization
  problem, and changing the solver may help.
* It can solve that problem well, but the chosen noise assumptions or penalty
  may make it prefer the wrong explanation of the data. This is a modeling
  or model-selection problem. Faster C++ alone will not fix it.
* It can find a real step, but report something unhelpful, such as a slowdown
  that later fully recovered. This is a reporting-policy problem.

For example, imagine we generated 100 observations with a true change from
10 to 12 after observation 50, then added random noise. If ASV finds a boundary
at 52, that is a localization error of two positions. If it finds no boundary,
that is a missed change. If it also finds a boundary at 20, that extra boundary
is a false detection. These judgments use the history we deliberately
generated, rather than just the detector's own cost function.

The proposed order of work is:

1. Build a repeatable way to measure current behavior.
2. Test whether searching more carefully for gamma helps.
3. Test whether we can fit candidates exactly at an acceptable cost.
4. Revisit how noise is modeled and scored.
5. Check the effect on the reports users actually see.

This order makes the causes of improvements easier to identify. Changing
the optimizer, noise floor, and reporting threshold together could produce
different answers without telling us which change helped or hurt.

Terms used in the experiments
-----------------------------

An **evaluation harness** is a small program that runs the same datasets
through several detector variants and saves comparable results. It is not a
new detection algorithm. A **baseline** is the current implementation used
as the reference for those comparisons.

An **oracle** is a deliberately simple, trusted way to compute the answer
for small inputs, even if it is too slow for real histories. Here we can
try every possible partition of a tiny list and keep the cheapest one.
The oracle tells us the optimum of a cost function, not the true history
behind noisy measurements.

**Synthetic data** are observations we generate ourselves so the underlying
changes are known. A **seed** makes random generation repeatable. **Held-out
data** are examples we reserve for evaluation rather than use while choosing
settings. Otherwise we risk tuning the detector to the examples we keep seeing.

An **acceptance criterion** states what evidence would justify continuing
with an idea. It prevents a visually appealing plot or a single favorable
example from becoming the whole argument for a change.

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

In the first example, the exact solution pays 16.05 in residual error and
zero boundary charges. The approximate solution pays 14.52 in residual error
and two boundary charges of 1, for 16.52 in total. It fits the observations
more closely, but the improvement is not enough to pay for its extra steps.
That is precisely the tradeoff the fixed-penalty objective is supposed to
make.

The generated data actually change their mean halfway through that example.
The exact solver's preference for one segment therefore does not prove it
found the true generating history. It proves only that the approximation
missed the cheapest solution at this particular gamma. This distinction is
why the plan measures both optimization error and statistical detection error.

In the second example, the absolute change remains 0.02 after adding 1000,
but its relative size shrinks from 2 percent to about 0.002 percent. A detector
of absolute changes might want the same answer; a detector designed around
relative measurement accuracy might not. We need to decide which behavior
is intended before labeling the difference a bug.

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

This phase asks: what does the current implementation do well, where does it
fail, and which part consumes the time? Its output is a comparison report,
not a new default algorithm.

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

A small first batch can use 100 observations with a known step at position
50, repeat that dataset with different noise samples, then repeat the entire
batch with no step. The changed histories measure sensitivity; the flat
histories measure false alarms. A method that always reports a step would
look excellent on the first batch and fail immediately on the second.

For boundary matching, specify a tolerance before examining results. For
example, a detected boundary within two positions of a true boundary might
count as a match, with each true and detected boundary matched at most once.
Report the actual position error too, so that a tolerant matching rule does
not hide imprecise localization. The choice of two is illustrative; suitable
tolerances depend on the histories and intended use.

Keep the output readable: one row per method and scenario can summarize
missed changes, false changes, typical location error, time, and memory.
Save the detailed per-run records behind that summary so surprising results
can be reproduced and inspected.

Phase 2 evaluates penalty search
--------------------------------

This phase asks whether ASV is overlooking good candidate histories because
it tries too few penalties or tries them in the wrong places. The candidate
fitting algorithm stays the same, so any differences can be traced to search.

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

Why might search matter? As gamma rises, a fit can stay exactly the same for
a while, then suddenly lose a boundary. Its outer score can therefore look
like a staircase, with several low regions, rather than one smooth bowl.
Golden-section search efficiently narrows a region when the objective falls
toward one minimum and then rises. That shape is not guaranteed here.

A grid experiment simply tries a list of penalties spaced by a fixed ratio.
If several adjacent trials produce the same partition, they give us the same
history to score. Refinement tries more penalties near a transition between
different partitions. This is easy to inspect, but a coarse grid can miss a
narrow useful region, and a dense grid can be slow.

CROPS means Changepoints for a Range of Penalties. Its purpose in this plan
is to explore changes in the optimal partition as the boundary charge varies,
rather than spending many trials rediscovering the same partition. It depends
on having an exact solver for each penalty; using the current approximation
does not automatically inherit that guarantee.

Phase 3 evaluates faster exact fitting
--------------------------------------

This phase asks whether we can remove approximation error without making
publishing impractically slow. There are two separate costs to attack:
how many candidate intervals we examine, and how much work each interval
requires. Improving only one may leave the other as the main bottleneck.

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

The intuition behind incremental interval statistics is simple. Suppose we
already know the sorted values and weighted sums for positions 0 through 99.
To extend the interval through position 100, we would like to insert one new
observation and update those statistics, rather than sort all 101 values
again. A weighted order-statistics data structure maintains enough information
to find the halfway point in cumulative weight. Efficient insertion, median
lookup, and error calculation must all be considered together.

PELT stands for Pruned Exact Linear Time. Pruning means proving that some
candidate starts of future segments can never beat another candidate, then
stopping consideration of those starts. This differs from limiting every
segment to 20 observations: the aim is to discard only choices proven
unnecessary for the optimum. The proof conditions and interval-cost work
still matter, which is why the name alone is not a runtime guarantee for ASV.

Treat these as two experiments first: faster interval calculations with the
same search, and fewer candidate intervals with the same cost calculations.
Then combine them if both are useful. This identifies which idea provides
the gain and whether their memory costs interact badly.

Phase 4 calibrates the noise model
----------------------------------

This phase asks whether the score favors the right kinds of histories.
Calibration here means choosing settings by measured detection behavior
across many examples, rather than because one setting looks convincing on
one graph. It does not turn the current score into a formal significance test.

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

An invariance is a change to the input that we believe should leave the
answer unchanged. Converting timings from seconds to milliseconds changes
their numerical scale but not which commits became slow. Adding a constant
baseline is a different operation: it preserves absolute changes while
altering percentage changes. Keeping those examples separate prevents a
useful unit-conversion property from being confused with a modeling choice.

To isolate the noise-floor question, score the same saved candidate
segmentations with the current floor and with each proposed replacement.
Then compare their selected histories against the known generating histories.
Only after that should we rerun the full search, where a changed score may
also alter which penalties are tried.

For logarithmic timing data, a change from 10 to 12 and a change from 100 to
120 have the same difference: ``log(1.2)``. That makes relative changes easier
to express. It also changes the measurement-error model, so it is a separate
statistical experiment rather than a harmless preprocessing step.

Phase 5 validates reporting and prepares integration
----------------------------------------------------

This phase asks whether an improvement in the fitting experiments makes ASV
more useful in practice. The detector is one part of a publishing pipeline;
users see commit ranges and regression records, not just lists of boundaries.

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

For example, a fit of ``10 -> 12 -> 10`` contains a real upward step, but the
last plateau shows full recovery. A fit of ``10 -> 14 -> 12`` retains a lasting
slowdown. Reporting tests must encode that difference. Likewise, a boundary
between two observations may correspond to a range of commits if measurements
are sparse. A numerical improvement in boundary placement should not create
false certainty about the exact responsible commit.

Concrete first deliverable
--------------------------

The first implementation should produce the evaluation harness and baseline
report from phase 1, followed by the bounded grid-search comparison in phase
2. This tests a small, isolated change before investing in a new C++ solver.
The report should identify which bottleneck dominates and whether a better
optimizer or a better model-selection rule is likely to improve actual alerts.
