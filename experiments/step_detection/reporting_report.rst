Does joint calibration make threshold evidence more useful?
===========================================================

Joint calibration improves the conservative reference, but its sensitivity
is still too low for a production reporter. It detects 126 of 960 true
above-threshold slowdowns, compared with 73 for the original reference.
It adds one false alert among 1440 null histories. Both rules miss every
recent-change case in this evaluation.

The result supports joint calibration as a better mathematical reference.
It also shows that demanding a confidence set for the whole step function
is too costly for the short histories we want to analyze.

What changed
-------------

The original reference gives a small separate error allowance to each
interval. The new method calibrates all the interval checks together,
accounting for their overlap. Both methods keep every plausible one-change
explanation and ask whether all of them imply a slowdown exceeding 5%.

Calibration used fair-coin sign sequences, not benchmark outcomes. The
`protocol <reporting_protocol.rst>`_, code, and resulting cutoffs were committed
as ``2302698`` before generating the evaluation histories. No settings were
changed after seeing the results.

For each of the two lengths, 8191 calibration sequences supplied the maximum
imbalance across intervals. The cutoff is their 7832nd sorted maximum. Its
exact binomial calibration-failure bound is about 0.004922 per length. With
at least 99% confidence over both calibrations, the frozen rules have at
most 5% per-history false-alert probability under the independent, at-most-
one-change model. This includes uncertainty about the boundary. It does not
certify a particular revision or cover correlated observations.

The first finite interval bound now appears at 11 readings for length-40
histories, versus 16 with the original reference; at length 100, the change
is from 18 to 12. These are properties of the calibrated constraints, not
minimum sample sizes that ensure an alert.

What was compared
------------------

The evaluation contains 2400 fresh histories: 1440 at or below 5%, and 960
above 5%. The null group includes 480 histories with exactly a 5% slowdown.
Half use independent Gaussian noise and half independent Laplace noise of
the same variance. Each complete case has a distinct random stream. All
methods receive identical observations within a case.

The shared-floor detector retains the inherited ``shared-r1-c2`` settings.
Its fitted segments and existing reports are saved once per history. The
two evidence gates keep an existing alert only when their confidence rule
also establishes an above-threshold initial-to-final change. A gate can
remove an alert but cannot add one. Standalone evidence results show what
the rule could establish without that detector requirement.

The certified claim is about change size across the history. Existing
reported locations are retained as detector output, not certified locations.
The model allows at most one true change even if the fitted detector happens
to produce several segments. This does not implement ASV's full semantics
for multiple lasting changes and recoveries.

Overall results
----------------

The miss columns each have 960 true 6%/8% histories. The false-alert columns
each have 1440 true 0%/4%/5% histories.

.. list-table:: Reporting comparison
   :header-rows: 1

   * - Method
     - True slowdowns detected
     - Missed slowdowns
     - False alerts
   * - Production default
     - 589 (61.4%)
     - 371 (38.6%)
     - 146 (10.14%)
   * - Shared floor, existing reporter
     - 806 (84.0%)
     - 154 (16.0%)
     - 253 (17.57%)
   * - Shared floor, original evidence gate
     - 73 (7.6%)
     - 887 (92.4%)
     - 0
   * - Shared floor, joint evidence gate
     - 126 (13.1%)
     - 834 (86.9%)
     - 1 (0.07%)

Joint calibration gains 53 positive histories over the original reference
and loses none. Relative to the existing shared-floor reporter, however,
it removes 252 false alerts while also losing 680 true detections. This is
a large sensitivity cost, not a general improvement over current reporting.

Standalone evidence and gated evidence agree on every history for both
references. In this evaluation, the fixed detector does not limit the new
rules: every history they certify already has an existing alert.
Neither confidence-set construction returns ``incompatible_model`` here;
the missed detections have plausible explanations that prevent certification.

Where the gains occur
---------------------

Each row below has 240 positive histories, pooled across Gaussian and Laplace
noise. Counts are detections, not misses.

.. list-table:: Sensitivity by history length and change location
   :header-rows: 1

   * - History
     - Existing shared-floor reporter
     - Original evidence gate
     - Joint evidence gate
   * - 40 readings, middle change
     - 217
     - 0
     - 4
   * - 40 readings, recent change
     - 161
     - 0
     - 0
   * - 100 readings, middle change
     - 233
     - 73
     - 122
   * - 100 readings, recent change
     - 195
     - 0
     - 0

In Gaussian histories the reference-to-joint detection count rises from
31 to 58 out of 480 positive histories. In Laplace histories it rises from
42 to 68 out of 480. The joint rule detects 41 of the 480 true 6% changes
and 85 of the 480 true 8% changes; the original reference detects 26 and 47,
respectively. It remains weak even on changes comfortably above 5%.

For a concrete gain, take
``gaussian-n100-d0.06-s0.005-middle-seed601``. The fixed detector splits at
observation 50 and fits levels 10.0070 and 10.6005. The original reference
keeps splits 33 through 67 and gives a lower bound on D of -0.02855. Joint
calibration narrows the plausible splits to 39 through 61 and raises the
lower bound to 0.02719. The latter clears zero and supports the alert without
pretending that the boundary is known exactly.

The one joint false alert occurs on
``gaussian-n100-d0.05-s0.05-middle-seed608``: a true 5% middle change with
noise standard deviation 0.5. Its lower excess is only 0.00326. No flat or
4% history is alerted by either reference. Exactly-5% alerts number 143
for production, 219 for the existing shared-floor reporter, zero for the
original reference, and one for joint calibration, each out of 480 histories.

Why the remaining conservatism matters
--------------------------------------

The rule constructs a confidence set for the whole step function, then asks
whether every member has a sufficiently large change. That is stronger than
what a size-only alert ultimately needs. A wrong boundary can leave a short
plateau whose level is poorly bounded, and that explanation can prevent an
alert even when the detector's preferred fit looks convincing.

The recent-change design leaves only four later readings at length 40 and
ten at length 100. Both are shorter than the first finite joint interval
bound. If the true boundary survives the constraints on its earlier segment,
its later level remains unbounded above and bounded below only by zero.
That surviving explanation prevents an above-threshold conclusion. Better
calibration alone does not fix this lack of information in the chosen
confidence construction.

This differs from the existing reporter's task. It can alert when a fitted
increase is just above 5%, including histories whose true increase is exactly
5%. The new rule asks for evidence that the underlying increase exceeds 5%.
Counting exactly 5% as null therefore makes the existing reporter's false-alert
rates here larger than in earlier studies that omitted those histories.
Those rates should not be compared directly across the different designs.

Next decision
--------------

Keep joint calibration as a checked improvement to the mathematical reference.
It is not yet a useful replacement for the existing reporter. The next
statistical design should test the composite null ``D <= 0`` more directly,
while still allowing unknown boundaries, instead of first demanding a
confidence set for the entire function. If stronger noise assumptions are
needed to obtain useful power on short histories, state and test them
explicitly.

Do not raise this cutoff's error allowance or adjust its statistic using
these evaluation histories. A new procedure needs a new derivation, frozen
calibration, and fresh evaluation. Keep production defaults unchanged and
retain both references to check the next method.

Evidence and reproduction
--------------------------

The `running instructions <running.rst>`_ give the freeze and evaluation
commands. ``data/reporting_v1_frozen.json`` preserves every calibration
maximum, the selected cutoffs, inherited detector, design, and source hashes.
``data/reporting_v1_results.json`` records results and archive checksums.
The ``reporting_v1_*_inputs_*.jsonl.gz`` and
``reporting_v1_*_records_*.jsonl.gz`` archives preserve all numerical inputs
and per-history outputs. Infinite bounds are stored as explicit strings.

The reporting reference, joint calibration, and pipeline have 44 passing
tests. Together with the existing robustness, component, threshold, and
exact-reference checks, 132 tests pass. Exact checks enumerate small sign
sequences, verify calibration ranks against direct binomial sums, and recover
the original reference when its cutoff is supplied. Production code remains
unchanged.
