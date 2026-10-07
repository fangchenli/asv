Why ASV detects steps this way
==============================

A benchmark runs a repeatable task to measure performance. ASV runs benchmarks
across versions of a program so we can see when performance changes. For a
timing benchmark, an increase means the task became slower.

The challenge is that timings fluctuate even when the code performs the same
way. ASV looks for lasting changes beneath that variation. Its central idea
is to describe the history with a few constant levels, choosing a new level
only when the measurements give enough reason to do so.

The `mathematical analysis <mathematical_analysis.rst>`_ derives what we can
establish before further experiments and identifies an exact solution for
the correlation fit. The existing `harness instructions <running.rst>`_ and
`baseline findings <baseline_report.rst>`_ record the initial experiments.
The `correlation-fit design <correlation_design.rst>`_ describes the first
correction derived from that analysis and its analytical tests.
The `shared noise-floor analysis <shared_noise_floor.rst>`_ derives a new
score and shows why noise persistence needs a separate modeling decision.
The `persistence bound <noise_persistence.rst>`_ gives that decision an
explicit half-life parameter and proves what a strict bound guarantees.
The `exact independent reference <exact_reference.rst>`_ now finds the global
shared-floor optimum and exposes the best fit at every feasible segment count.
The `before/after comparison <version_comparison.rst>`_ measures the correction:
identical decisions across 190 histories, faster short-history detection,
and no consistent speedup on long histories.
The `near-threshold study <threshold_report.rst>`_ tests statistical changes:
bounded persistence finds more small slowdowns, with an explicit increase in
below-threshold alerts. Its settings were frozen before evaluating fresh cases.
The `component comparison <ablation_report.rst>`_ isolates why: the shared
floor prevents severe overfitting when the segment penalty is reduced. The
correlation cap offers little advantage after calibrating false-alert rates.
The `robustness study <robustness_report.rst>`_ keeps those settings fixed
under outliers, gaps, heavy tails, and changing noise. The overall improvement
persists, but false-alert rates vary by condition; uncertainty about change
size and location is the next modeling question.
The `reporting uncertainty derivation <reporting_uncertainty.rst>`_ now explains
how to ask whether a slowdown convincingly exceeds 5%. An independent-noise
reference accounts for every possible single-change boundary and provides
a conservative error bound. It also shows why three post-change readings
and correlated noise require more care than the current residual allowance.
The `joint calibration study <reporting_report.rst>`_ then calibrates the
interval checks together and compares reporting rules on fresh independent
histories, with the shared-floor detector held fixed and exactly 5% changes
included among the null cases.
The `direct threshold test <direct_report.rst>`_ then checks every
below-threshold explanation using size and shape tests. Its derivation
accounts for the unknown boundary under a stronger assumption: independent
Gaussian noise with one common variance. A fresh evaluation compares it
with both references and treats Laplace noise as a stress check.
The `variance and correlation study <direct_stress_report.rst>`_ then holds
that test fixed and shows where its assumptions matter. Rising noise and
correlation increase false alerts through both the size and extra-step
checks, motivating a covariance model shared by the two tests.
The `known-covariance derivation <known_covariance.rst>`_ supplies that
reference: generalized least squares corrects both tests when relative
variances and correlations are known. Analytical checks recover the
independent-noise test and verify the projection identities behind the
error bound. Estimating covariance remains a separate problem.
The `covariance oracle study <covariance_report.rst>`_ measures that reference
on fresh paired histories: gated false alerts fall from 85 to eight, while
true detections fall from 585 to 499. Quiet recent plateaus benefit, while
noisy correlated plateaus lose sensitivity. The next question is uncertainty
in the covariance itself.
The `unknown-correlation reference <unknown_correlation.rst>`_ now constructs
a joint confidence set for AR(1) correlation and change location. It checks
the continuous stationary correlation range using exact polynomial interval
certificates, with unresolved regions producing a non-alert. Its error budget
includes uncertainty in the confidence set.
The `unknown-correlation study <ar1_report.rst>`_ now measures its cost:
16 of 144 true slowdowns detected, versus 103 for the oracle at the same
reporting budget, with zero false alerts for both. Every search resolves;
broad confidence regions cause the loss. A mathematical diagnosis explains
why retained correlations near 1 undermine percentage-change evidence.
The `information-loss analysis <ar1_information.rst>`_ checks that next
candidate on a saved missed slowdown. Combining both directions still keeps
the problematic correlation, even after tightening the candidate likelihood.
The derivation separates prediction cost from discarded stationary
information and motivates a full-history confidence construction.
The `full-history reference <full_history_confidence.rst>`_ now integrates a
proper predictive mixture and certifies confidence exclusion over correlation
intervals. It removes sufficiently persistent explanations but still retains
the saved 8% slowdown's witness, motivating a closer look at predictive cost.
The `baseline/jump analysis <baseline_jump_prior.rst>`_ now separates those
costs and finds an unnecessary location-mixture penalty. A proper density
for each candidate location avoids that penalty with the same coverage
argument. Both indexed variants certify the archived missed 8% slowdown;
changing the jump prior alone does not.
The `fresh indexed comparison <indexed_report.rst>`_ now measures the gain:
28 of 144 true slowdowns detected, versus 20 for the conditional reference,
25 for global full-history prediction, and 27 for indexed-original. All four
have zero false alerts among 216 null histories; the matching oracle detects
114. Of indexed-jump's 116 misses, 115 retain an explicit noise explanation.
The `residual-direction derivation <residual_direction.rst>`_ now removes
plateau levels and overall scale before constructing correlation confidence.
It establishes coverage, a connection to the existing GLS polynomials, and
both endpoint limits. A fixed directional mixture still retains the saved
positively correlated witness.
The `directional tail bound <directional_tail.rst>`_ now supplies that
calibration near correlation 1. Exact rational inequalities put the saved
witness's tail probability below 0.82% throughout its interval up to 1,
crossing the 1% confidence cutoff. The `complete directional reference
<directional_reporting.rst>`_ now combines both bounds in the full search.
It certifies the saved correlated 8% slowdown, with every interval checked;
the uniform-only control retains the old witness. The `fresh directional
comparison <directional_report.rst>`_ now detects 58/144 true slowdowns,
versus 27/144 for indexed-jump on the same inputs, with zero false alerts
among 216 null histories for both. The tail bound adds 11 detections beyond
uniform directional confidence. Most remaining misses retain an explicit
explanation. The `full-determinant derivation <directional_determinant.rst>`_
now certifies a correlation interval around one of the two lost detections'
witnesses: its probability bound falls from 1 to below 0.007037. A product-based
combination of size and directional evidence retains the second saved witness.
The `complete determinant search <determinant_reporting.rst>`_ now replays
seven saved cases with every result verified. Decisions are unchanged:
removing the first lost case's old witness reveals another at correlation
0.8125. The `spectral refinement <spectral_refinement.rst>`_ now tightens that
bound from 0.032029 to 0.006590 using eight residual directions and a second
tilt. A verified seven-case replay retains the same decisions and finds
another surviving correlation at 0.796875. The remaining approximation gap
is the next mathematical question. The `split-aware eigenvalue derivation
<sturm_refinement.rst>`_ now separates that gap and uses exact precision
eigenvalue counts within the long plateau. It certifies the 0.796875
explanation's neighborhood below 0.009631. The `complete Sturm replay
<sturm_reporting.rst>`_ verifies all seven saved outcomes with no changed
decisions and finds another survivor at 0.7890625. Numerical projected
eigenvalues still leave its density bound above 0.01, directing the next
work toward a sharper probability inequality.

The isolated `Lean pilot <lean/README.rst>`_ now proves the scalar interval
rules and checks the final exact arithmetic of one saved certificate. Its
`determinant trace <lean/determinant_trace.rst>`_ connects those rules to all
1,723 nodes of the saved determinant calculation, throughout the correlation
interval. The matrix identity, spectral bound, and probability argument
remain separate proof obligations; detection results are unchanged.

.. contents:: On this page
   :local:
   :depth: 1

What we are trying to fit
-------------------------

Consider timings in milliseconds for eight successive program versions::

    Version          1      2      3      4      5      6      7      8
    Measured time   10     10.1    9.9   10     12     12.1   11.9   12

The first four readings cluster around 10; the last four cluster around 12.
A useful description would assign a typical time of 10 to the first group
and 12 to the second::

    Fitted time     10     10     10     10     12     12     12     12

**Fitted** means the value chosen to describe the measurements. Here the fit
expresses the idea that performance changed once, while the smaller deviations
around 10 and 12 are measurement variation. Choosing the levels and where they
change is the fitting problem.

Each constant stretch is a segment or plateau. The change between two levels
is a step. ASV assumes this shape is a useful description of benchmark history:
most commits leave a benchmark roughly unchanged, while some introduce jumps.

Why measure absolute error
--------------------------

To compare descriptions, we need a measure of how well each agrees with the
readings. ASV adds the absolute distance between each measurement and its
fitted value. Both above and below count as disagreement, so they cannot
cancel one another.

Suppose we describe all eight versions with one level, 11. The distances are::

    Measurement     10   10.1   9.9   10   12   12.1   11.9   12
    Distance to 11   1    0.9   1.1    1    1    1.1    0.9    1

    Total error = 1 + 0.9 + 1.1 + 1 + 1 + 1.1 + 0.9 + 1 = 8

That 8 is the accumulated mismatch of the proposed description. For the
alternative with levels 10 and 12, the distances are much smaller::

    Distance to fit  0    0.1   0.1    0    0    0.1    0.1    0

    Total error = 0 + 0.1 + 0.1 + 0 + 0 + 0.1 + 0.1 + 0 = 0.4

So the split explains the observations better: it removes 7.6 of fitting
error. This is why the absolute-error sum appears in the algorithm. It
quantifies how much the evidence favors one description over another.

Absolute error also limits the influence of extreme timings compared with
squared error. An unusually slow measurement contributes in proportion to
its distance, rather than the square of that distance.

For a fixed segment, the level that minimizes absolute error is a **median**.
This explains why ASV uses medians: they follow from the error measure. For
example, the readings 10, 10, and 100 have total absolute error 90 around their
median of 10, versus 120 around their mean of 40. The isolated high reading
pulls the mean farther from the two ordinary readings.

In the eight-reading example, the middle two sorted values are 10.1 and 11.9.
Their midpoint, 11, is ASV's one-segment level. Any value between those middle
two gives the same absolute-error total. The two separate groups instead have
medians 10 and 12.

Why charge for changes
----------------------

If agreement with the readings were the only goal, copying every reading
would give zero error. The fit would claim a performance change at every
version, including each small fluctuation.

ASV adds a penalty for each change in fitted level. With a penalty of 1,
three possible descriptions have these scores::

    Description                 Fitting error    Changes    Total score
    One level, 11                     8             0          8
    Level 10, then level 12            0.4           1          1.4
    Copy each reading                 0             7          7

The two-level description wins: one substantial change explains most of the
variation, and fitting the remaining small fluctuations would cost too much.

The penalty is called **gamma**. Increasing it favors fewer changes. In this
example, gamma = 10 makes the two-level score 10.4, so the one-level score of
8 wins. An added change must reduce fitting error by enough to pay its penalty.

This also explains why persistence matters. A small slowdown repeated over
many versions can accumulate enough error under the old level to justify
a step. The same slowdown seen once provides less evidence.

The connection to the Potts model
---------------------------------

The fitting problem we have just built is a one-dimensional Potts problem::

    score = sum of weighted absolute errors + gamma * number of changes

The Potts idea is that neighboring fitted values are encouraged to agree.
A disagreement costs a fixed amount, regardless of the jump's size. That
encourages a small number of constant stretches while allowing their levels
to be determined by the data.

In mathematical terminology, this combines an L1 data error with an L0
penalty on successive differences: L1 sums the sizes of errors, and L0 counts
how many changes occur. The distinction matters because the penalty controls
how often the level changes, rather than directly penalizing how far it jumps.

The weights account for unequal measurement uncertainty. A precise reading
has more influence than an uncertain one. ASV derives relative weights from
timing confidence intervals; with weights, the best level for a segment is
a weighted median. The example above gives every reading equal influence.

How the fitting algorithm finds the steps
-----------------------------------------

ASV must choose both the boundaries and the levels. Once boundaries are
chosen, medians give the levels. The remaining problem is to find a good
division of the ordered measurements into segments.

For a fixed gamma, the exact solver uses dynamic programming. It considers
each possible start of the final segment and combines that segment's cost
with the best solution for everything before it. Reusing the answers for
shorter prefixes avoids trying every complete division independently.

Normal detection uses an approximation for speed. It first solves a restricted
problem with short segments, then merges neighboring segments when doing so
reduces the score, and adjusts the surviving boundaries. The initial segment
limit is normally 20 observations; merging allows the final plateaus to be
much longer. These local improvements can miss the global best division.

The C++ module accelerates interval medians, interval errors, and the dynamic
programming loop. It caches interval answers so different trials can reuse
them. Python controls the approximation, penalty selection, and reporting.

Why choosing gamma needs another score
--------------------------------------

The fitting solver answers: given this price per change, what description
has the lowest score? We still need to choose the price.

ASV tries several gamma values and compares their fitted histories using a
separate score. This creates two layers:

* The inner fit balances absolute error against the specified price per jump.
* The outer selection decides which resulting history has a useful level of
  detail, accounting for its remaining noise.

A fixed numerical gamma would be awkward across benchmarks with different
units and noise levels. The outer score lets ASV estimate noise from each
history while selecting how many segments to keep.

Why the outer score uses a logarithm
------------------------------------

The logarithm comes from estimating an unknown noise scale. To see the
connection, start with independent Laplace noise, a statistical model whose
likelihood leads to absolute error. Write D for the total weighted absolute
error, n for the number of observations, and b for the unknown noise scale.
Its negative log likelihood, dropping terms shared by the candidate fits, is::

    n * log(b) + D / b

The best scale for a given fit is b = D / n. Substituting it gives::

    n * log(D / n) + n

After dividing by n and dropping terms common to every candidate, the part
that depends on fit quality is log(D). Add a charge for model complexity,
and we obtain the shape of ASV's selection score::

    complexity charge + log(remaining error)

The practical consequence is that error reductions are judged relatively.
Reducing error from 100 to 50 gives the same improvement in the logarithmic
term as reducing it from 10 to 5. That is useful when different benchmarks
have different timing scales. The complexity charge still determines whether
that improvement is worth adding another segment.

ASV uses a complexity charge of ``4 * log(n) / n`` per segment and makes two
further adjustments to the remaining-error term.

First, it accounts for **correlation** between neighboring residuals, where a
residual is a measurement minus its fitted level. If several readings are
slow because the machine was busy, their errors can resemble one another.
The score allows some of that pattern to be explained as related noise,
reducing the pressure to turn every persistent fluctuation into another step.
This adjustment ranks the candidate fits; the inner solver still uses the
ordinary weighted absolute-error objective.

Second, it adds a **noise floor** before taking the logarithm. Otherwise a
fit that exactly copies every reading would have zero residual error and an
unbounded advantage from log(0). The floor makes that reward finite. ASV's
current floor depends on the fitted levels, so its choice affects which
histories are selected.

The coefficient, floor, and approximate search are practical choices around
this statistical motivation. Their effect on detection quality is a central
part of the improvement plan.

Why reporting regressions is a separate step
--------------------------------------------

A fitted history describes both speedups and slowdowns. ASV then decides which
upward changes should appear in the regression report, considering their
size, residual variation, and later recovery.

Our example changes from 10 to 12 milliseconds, a 20 percent slowdown. At a
5 percent reporting threshold, it is large enough to report and well above
the variation within each plateau.

A history of ``10 -> 12 -> 10`` tells a different story: a slowdown followed
by full recovery. It need not remain in the regression report. A history of
``10 -> 14 -> 12`` still has a lasting slowdown, so ASV can report the remaining
change from 10 to 12 at the original upward boundary.

Gamma controls the fitted history. The reporting threshold controls which
changes in that history matter enough to display.

How the pieces connect
----------------------

The complete path is::

    Benchmark results and uncertainty estimates
        -> order results by program version and prepare weights
        -> try a gamma
        -> fit segments using medians, errors, and change penalties
        -> compare fits from different gammas using the noise score
        -> map boundaries back to versions of the program
        -> apply regression reporting rules

The measurements establish the evidence; the Potts objective defines the
tradeoff; the solver searches for a fit; automatic selection chooses its
complexity; and reporting turns fitted changes into actionable records.

For function names, the exact recurrence, and executable examples, see the
`implementation reference <implementation_details.rst>`_. The
`improvement plan <improvement_plan.rst>`_ describes how to evaluate changes
to the solver, penalty search, and noise model separately.

This explanation describes source revision
``d33754e129c025beb5c2ca440c3c280433b264f7``.
