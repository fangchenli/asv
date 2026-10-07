Step detection improvement plan
===============================

The goal is to detect lasting performance changes more reliably while keeping
analysis practical for long benchmark histories. We should evaluate three
layers separately: finding the best fit at a fixed penalty, selecting the
penalty and noise model, and reporting useful regressions.

The `conceptual guide <README.rst>`_ explains why those layers exist. The
`implementation reference <implementation_details.rst>`_ maps them to the
code at revision ``d33754e129c025beb5c2ca440c3c280433b264f7``.

Continue from the `mathematical analysis <mathematical_analysis.rst>`_ and
`shared noise-floor design <shared_noise_floor.rst>`_. Establish the objective,
guarantees, and noise assumptions before extending the experiments. Exact
conditional correlation fitting is implemented. The experimental score now
requires explicit noise-scale and persistence inputs; production defaults
have not been chosen. An `exact independent-noise reference <exact_reference.rst>`_
now establishes the global optimum for a declared floor and beta.

The `before/after comparison <version_comparison.rst>`_ now measures the
correlation correction on 180 synthetic histories, six public NumPy histories,
and four scaling cases. All fitted steps and alerts stay the same. Native
length-100 runs take about 11% less time at the median; longer runs show no
consistent gain and use more memory. Exact checks of seven representative
boundary-error cases find no fixed-penalty optimization gap.

The `near-threshold study <threshold_report.rst>`_ now evaluates that harder
statistical comparison. Settings selected on 360 development histories were
frozen before 600 held-out histories. Bounded persistence reduces missed
above-threshold histories from 147/240 to 71/240, while below-threshold alert
histories rise from 0/240 to 10/240. The common-pool control leaves the current
detector's decisions unchanged. This is a scoring tradeoff, with no demonstrated
need for a faster solver on these small histories.

The `component study <ablation_report.rst>`_ now separates those changes on
1200 development and 2400 fresh histories. Lowering the penalty with the
current floor creates one segment per observation in 2370 of 2400 histories.
The shared-floor formulation prevents that failure. With settings calibrated
to a 5% development false-alert budget, shared-floor/cap-1 has 298 misses
and 24 below-threshold alerts; shared-floor/cap-0.5 has 306 and 23, respectively,
each out of 960 eligible histories per error type. The cap provides little
advantage after calibration. Actual held-out rates remain essential because
some settings exceed their development budgets.

The `robustness study <robustness_report.rst>`_ now evaluates 7680 fresh paired
histories without changing those settings. At the inherited 5% budget, the
shared floor reduces misses from 1696 to 1347 and below-threshold alert
histories from 208 to 163 versus the calibrated current floor, each error
type measured over 3840 eligible histories. Its observed false-alert rate
exceeds 5% with positive outliers and falling noise amplitude. Missing data
also changes the sensitivity/false-alert tradeoff.

The `reporting uncertainty derivation <reporting_uncertainty.rst>`_ now defines
the claim as a positive lower bound on ``later - 1.05*earlier``. A mathematical
reference uses simultaneous median bounds to keep every plausible one-change
explanation, including uncertainty about the boundary. It gives a conservative
per-history false-alert bound for independent observations with positive
medians and at most one true change. Analytical tests verify its binomial
ranks and confidence-set construction. It is not a production default.

The `joint calibration study <reporting_report.rst>`_ now compares those
rules on 2400 fresh independent histories, with the shared-floor detector
held fixed. Fair-sign calibration was frozen before evaluation. Joint
calibration detects 126 of 960 above-threshold histories versus 73 for the
original reference, adding one false alert among 1440 null histories.
Exactly 5% is included in the null. Both rules miss every recent change.
Standalone and gated evidence agree throughout, identifying the confidence
construction itself as the sensitivity limit in this study.

The `direct threshold study <direct_report.rst>`_ now tests every fixed-location
null using size and shape checks, then requires all locations to be rejected.
Its finite-sample error argument assumes independent Gaussian noise with a
common variance. With analytical cutoffs frozen before 2400 fresh histories,
the direct gate detects 530 of 960 positive histories versus 119 for the
joint-sign gate, with 10 versus zero false alerts among 1440 null histories.
It detects 225 of 480 recent changes, where the sign references detect none.
Nine false alerts occur under Laplace noise, outside its Gaussian guarantee.

The `variance and correlation study <direct_stress_report.rst>`_ now keeps
that rule fixed across 2400 versions of 400 fresh paired Gaussian histories.
Direct-gate false alerts rise from 3/240 in the independent constant-variance
control to 22/240 with rising variance, 22/240 with correlation, and 33/240
with both. True detections in the last condition fall from 124/160 to 97/160.
Both size and shape checks reject valid below-threshold explanations.

The `known-covariance reference <known_covariance.rst>`_ is now derived and
implemented for Gaussian noise with known covariance up to an unknown
overall scale. It transforms both the plateau design and observations,
uses GLS for the size contrast, and compares nested weighted residuals
for the shape scan. Orthogonal Gaussian projections give the same t/F
cutoffs and per-history error bound. Deterministic tests recover the old
test at identity covariance and check independent GLS calculations.

The `covariance oracle study <covariance_report.rst>`_ now evaluates 2400
fresh paired histories with their generating covariance supplied. Gated
false alerts fall from 85 to eight, while true detections fall from 585 to
499. Independent falling variance improves both counts; recent detections
there rise from 23/80 to 46/80. Rising correlated noise loses half its
detections, from 96/160 to 48/160, while false alerts fall from 35/240 to
1/240. All 400 identity-control decisions agree, and no numerical abstentions
occur. Supplying a variance transition tied to the mean boundary is additional
oracle information, so this is not a practical covariance-estimation result.

The `unknown-correlation reference <unknown_correlation.rst>`_ now handles
constant marginal variance with unknown stationary AR(1) correlation and an
unknown mean step. A conditional predictive-likelihood construction produces
a joint confidence set for correlation and location with failure budget 0.01.
Reporting receives the remaining 0.04. Exact polynomial interval certificates
cover the continuous range -1<rho<1; work limits produce unresolved non-alerts.
Analytical checks verify conditional densities, GLS equivalence, and the
certificate construction.

The `unknown-correlation study <ar1_report.rst>`_ now evaluates 360 fresh
versions of 120 paired histories with settings frozen in commit ``601bed9``.
The rule detects 16 of 144 true slowdowns versus 103 for the oracle at the
same reporting budget, with zero false alerts among 216 nulls for both.
All searches resolve: there are 16 certified alerts and 344 explicit surviving
explanations. The confidence sets are too broad. A follow-up mathematical
diagnosis explains 72 of the 87 lost oracle detections through correlations
arbitrarily close to 1, where uncertainty in the baseline defeats the
percentage-change test. Increasing the search budget cannot fix these cases.

The `information-loss analysis <ar1_information.rst>`_ now derives that
forward/reverse intersection and a stronger averaged-evidence combination.
Both retain the archived missed 8% history's near-1 witness, even with exact
conditional nuisance fits. The reverse predictor improves substantially, but
its prediction cost still absorbs most of the evidence against persistence.
The candidate-model relaxation has little effect in this example. Conditioning
also drops stationary information: the full-history likelihood contains a
square-root determinant factor that rules out correlations sufficiently close
to 1 when the limiting residual cost is positive.

The `full-history reference <full_history_confidence.rst>`_ now uses a proper
mixture over no change, one change, plateau levels, and noise variance. Exact
Bernstein bounds on low-degree GLS polynomials certify continuous confidence
exclusion without expanding their nth powers. Its reference unit must be
specified externally. The saved case excludes correlation 0.9999 at the true
location but still retains 4095/4096, which fails both reporting tests.
A short deterministic 20% step also remains inconclusive; a longer fixture
supplies a complete alert certificate. These are development diagnostics.

The `baseline/jump analysis <baseline_jump_prior.rst>`_ now decomposes that
cost and derives a midpoint/jump prior. More importantly, it proves that
each candidate location can use its own proper predictive density without
paying a mixture penalty over locations. True-pair coverage uses only that
pair's indexed density. The new prior alone raises the saved example's
density by 1.56 times but merely moves the surviving witness. Indexed
prediction produces a complete alert certificate with either level prior.
The short 20% fixture remains inconclusive.

The `fresh indexed study <indexed_report.rst>`_ now holds those settings fixed
on 360 histories. Indexed-jump detects 28 of 144 true slowdowns, compared with
20 for conditional, 25 for global full-history, 27 for indexed-original, and
114 for the matching oracle. All five have zero false alerts among 216 nulls.
Indexed-jump gains ten detections and loses two versus conditional. The new
jump prior adds only one detection beyond indexing with the original prior.
Positive correlation remains difficult: all four unknown-correlation methods
detect only one of 48 true changes. All archived inputs, certificates,
witnesses, and summaries were verified.

The `residual-direction analysis <residual_direction.rst>`_ now derives that
next construction. Projecting out both plateau levels and normalizing the
residual vector removes their nuisance parameters exactly. The angular
Gaussian density supplies a valid confidence rule and reduces to the existing
GLS cost times a cubic factor. Both correlation endpoints are analyzed,
including the exceptional alternating residual direction near -1. At +1,
the density has a finite positive limit, so endpoint exclusion is no longer
automatic. Rational square-root bounds permit exact confidence arithmetic.

On the saved positively correlated miss, a fixed five-component directional
mixture increases log evidence at the old witness from -45.458 to 1.603,
still below the 4.605 exclusion cutoff. Its limit at +1 also remains below
the cutoff. This establishes a remaining limitation without generating new
histories or running a new detection study.

The `directional tail analysis <directional_tail.rst>`_ now derives an exact
interval certificate. Exponential tilting bounds the weighted chi-square
tail; a density bound using four covariance directions tightens it enough
to reject the saved correlated witness. Matrix bounds and tridiagonal
determinant recurrences prove a probability below 0.008194 throughout
[16383/16384,1) at location 25. All 60 directional and tail tests pass.
This removes the endpoint obstruction at that location, not every remaining
explanation of the history.

The `complete directional reference <directional_reporting.rst>`_ now
integrates both confidence certificates with reporting alpha=0.04 and the
existing work limits. All ten archived and deterministic runs pass full
certificate verification. The new tail route turns the saved correlated
8% miss into a certified alert; the uniform-only control retains the old
witness. The 4% fixture remains a non-alert.

The `fresh directional comparison <directional_report.rst>`_ now freezes
that reference and evaluates 360 new histories. Combined confidence detects
58/144 true slowdowns, versus 47 for uniform-only and 27 for indexed-jump;
all three have zero false alerts among 216 null histories. The tail bound
adds 11 detections, all in 100-reading positively correlated histories.
Relative to indexed-jump, combined confidence gains 33 and loses two.
All 2160 unknown-correlation outcomes are verified and archived.

The matching oracle detects 108/144. Of combined confidence's 86 misses,
84 retain explicit explanations and two hit a work limit. Its two losses
to indexed-jump both retain location 1 and correlation 0.875. Numerical
tail estimates at those witnesses are 0.0003212 and 0.01682, on opposite
sides of the 0.01 confidence cutoff. These estimates are diagnostic and
do not replace rigorous certificates.

Next, derive interval bounds useful at moderate correlations, then examine
joint calibration of size and noise-pattern evidence. The first loss
suggests bound tightening can remove a witness; the second suggests that
the same tail test would still retain its witness after exact evaluation.
Use these saved cases for development and retain all frozen comparisons.
Any revised statistical rule requires a new freeze and fresh histories.
Keep predictive mixtures as separate controls unless a new coverage
argument explicitly combines them with directional evidence.

Unknown variance ratios and noise-transition locations remain subsequent
mathematical tasks; use the completed studies as diagnostics.
Laplace behavior near exactly 5%, especially on short plateaus, remains
a separate distributional question.
Preserve observations and usable counts; fitted step tuples alone are
insufficient. Multiple changes, recovery semantics, weights, and repeated
publication remain separate tasks. Continue investigating scale estimation
under varying noise and gaps. Production defaults remain unchanged.

The initial `baseline report <baseline_report.rst>`_ retains its historical
results; `running the experiments <running.rst>`_ explains reproduction.

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

Phase 0 establishes the mathematical design
-------------------------------------------

The `analysis <mathematical_analysis.rst>`_ provides derivations for fixed-penalty
dynamic programming, weighted-cost pruning, penalty-path structure, conditional
correlation fitting, noise-scale profiling, and transformation invariances.
Use those results to resolve three design questions:

* **Correlation:** the weighted-median correction is implemented with an
  explicit domain and numerical handling. For residuals [-1, 0, 1], the
  exact score is ``2 + abs(rho)``; the old equal-value stopping rule could
  return rho near -1 instead of the optimum 0. The corrected fit also
  replaces the expanded bracket with [-1,1].
* **Selection score:** the constrained-Laplace shared-floor score is derived
  and implemented in an experimental helper with an explicit floor and beta.
  Its independent-error version admits an exact penalty-path guarantee.
  The floor still needs an absolute measurement-scale interpretation, and
  correlation near one can hide a lasting step even with this new score.
  The `persistence design <noise_persistence.rst>`_ adds an explicit half-life
  cap and proves a lower bound on how much error correlation can remove.
* **Exact fitting:** specify pruning and tie rules under each supported segment
  constraint, together with a weighted interval-cost data structure. Establish
  correctness before evaluating computational savings.

The `correlation-fit design <correlation_design.rst>`_ is now implemented:
a weighted median constrained to [-1,1], with ties resolved toward zero and
analytical acceptance tests. It replaces the premature-stopping search.
The `shared-floor analysis <shared_noise_floor.rst>`_ establishes the new
score, its invariances, and a persistence counterexample. Production retains
its existing floor until absolute-scale metadata and the intended persistence
model have been specified. Those decisions come before wider search experiments.

The bounded-persistence prototype implements a chosen model without assigning
automatic defaults. The exact independent-noise reference now computes the
minimum error at every feasible segment count and selects the global shared-floor
optimum. It is verified by exhaustive partition enumeration on both backends.
Use it to check a faster exact penalty path next; the derivation supplies
a sufficient penalty range. Uncertainty metadata and persistence calibration
remain separate tasks.

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

After resolving phase 0, keep the chosen candidate fitter, correlation score,
and noise floor fixed within a comparison. Compare
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

Existing deliverables
---------------------

The initial harness, baseline, bounded search comparison, mathematical
analysis, correlation correction, experimental noise model, and exact
reference are available. The controlled version comparison adds fresh timing,
real-history replay, and exact diagnoses of remaining boundary errors.
The near-threshold study adds a frozen development/evaluation split and
measured alert tradeoffs for independent and bounded-persistence scoring.
The component study adds controlled factor comparisons and frozen calibration
at three false-alert budgets. The frozen robustness study measures how those
settings behave under less regular noise and missing readings. Statistical
meaning and calibration of near-threshold alerts now precede a new default.
