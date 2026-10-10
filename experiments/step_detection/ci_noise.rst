Investigating large CI timing spikes
====================================

A CI job can be 50% slower even when the code has not changed. The most
useful first step is to run the old and new code in that same job. Their
ratio can remove a shared machine slowdown. Independent jobs then tell us
whether the remaining difference repeats.

This note derives that approach before running experiments. The accompanying
`exact reference <ci_noise.py>`_ checks the algebra and probability budgets.
It does not yet collect timings or change ASV's reporting behavior. We have
not inspected the user's actual CI data, so both suite-wide and isolated
spikes remain possible explanations.

What one slow result can tell us
--------------------------------

Suppose yesterday's result was 10 ms and today's is 15 ms. Two explanations
fit exactly:

* The code became 50% slower on an unchanged machine.
* The code stayed unchanged and today's environment made it 50% slower.

With only those two measurements, the distinction cannot be recovered by a
better step detector. A longer history supplies evidence about persistence,
but a sustained environment change can still resemble a sustained code
regression. A fresh measurement of the old code supplies different evidence.

For example, a real 5% regression under a shared 50% machine slowdown gives::

                         Old code       New code       New / old
    Normal job           10 ms          10.5 ms         1.05
    Slow job             15 ms          15.75 ms        1.05

Comparing 15.75 ms with the historical 10 ms suggests a 57.5% regression.
The within-job comparison still identifies 5%.

ASV's absolute-error Potts fit reduces the influence of outliers within a
plateau, but the fit also chooses the plateaus. In a simple unweighted history
that is otherwise exactly 10 ms, fitting an isolated interior 15 ms reading
as its own plateau removes 5 ms of absolute error and adds two boundaries.
With a fixed penalty ``gamma`` per boundary, that is worthwhile when
``5 > 2*gamma``. At the newest reading it costs only one boundary. Thus an
absolute-error fit can still interpret a large spike as a step; this example
concerns the fitting objective, before ASV's later regression reporting.

A model that separates the causes
---------------------------------

For positive timing measurements, write

.. math::

   Y_{jvt} = \theta_v S_j \exp(g_j(t) + e_{jvt}).

Here ``j`` identifies a CI job, ``v`` is the old or new version, and ``t``
identifies a measurement slot. ``theta`` is the runtime we want to compare.
``S`` is a slowdown shared by both versions throughout the job. ``g`` is
drift during the job, and ``e`` is the remaining measurement disturbance.

Taking logs makes the contributions additive:

.. math::

   z_{jvt} = \mu + \delta\,1\{v=B\} + u_j + g_j(t) + e_{jvt}.

``A`` means baseline and ``B`` means candidate. The code effect is
``delta = log(theta_B / theta_A)``. Subtracting baseline log time from
candidate log time removes ``u_j``. This is why logarithms are useful here:
we can remove a multiplicative disturbance by subtraction, while retaining
the percentage comparison through ``exp(delta) - 1``.

The cancellation is exact for any shared positive ``S_j``; its size need not
be small or Gaussian. It does not remove a disturbance affecting only the
candidate, an additive delay, or different responses of the two versions to
the machine. Those effects require additional measurements or assumptions.

Why order matters
------------------

Running all baseline measurements before all candidate measurements lets
machine drift imitate a regression. Consider four equally spaced slots in
the order ``A B B A``. Give each slot one summarized log timing and compute

.. math::

   D_j = (z_1 + z_2 - z_0 - z_3)/2.

This is the difference between the candidate's average log timing and the
baseline's average log timing. A constant job effect cancels. Linear drift
``g(t) = a + bt`` also cancels because ``1 + 2 = 0 + 3``. Equivalently,
``exp(D_j)`` is the ratio of the two versions' geometric mean timings.

Actual benchmark slots need not be equally spaced: builds, setup, and other
benchmarks can take different amounts of time. In clock time, the required
condition is ``t1 + t2 = t0 + t3``. Record timestamps so we can check this;
the order alone does not establish it. Curved drift also leaves a remainder.

Randomizing between ``ABBA`` and ``BAAB`` adds a separate protection. Suppose
the four slot disturbances are fixed independently of which version we
assign to each slot, and the version effect is the same ``delta`` in every
slot. Conditional on these disturbances, the two orders give

.. math::

   D_j = \delta + H_j \quad\hbox{or}\quad D_j = \delta - H_j,

with equal probability. This follows by reversing the signs of every slot
disturbance when we reverse the version labels. No Gaussian assumption or
independence between slots is needed. At the boundary ``delta = tau``, at
most one order produces ``D_j > tau``. Below that boundary, the same bound
holds. A coin flip controls the order; it does not guarantee a small error
in any particular job.

The assignment-independent disturbance assumption matters. A candidate that
heats the CPU differently, changes cache state for later slots, or changes
their start times can violate it. A fixed warmup/setup procedure and recorded
timestamps help diagnose these interactions. With several blocks in one job,
a fixed weighted mean of block contrasts preserves symmetry if all block
orientations are independently randomized under the same assumption. The
job still contributes one contrast to the analysis below.

What independent jobs buy us
----------------------------

Choose a slowdown threshold before collecting confirmation data. For a 5%
threshold, ``tau = log(1.05)``. Count jobs whose contrast strictly exceeds
``tau``. A value exactly at the threshold stays in the job count and does
not count as an exceedance.

Under the clean model above, when the true code effect is at or below the
threshold, each job has exceedance probability at most one half. Alternatively,
we can assume this sign condition directly; then the conclusion concerns
the clean distribution of paired contrasts, and a connection to the code
effect still needs a measurement model.

For ``n`` independent jobs and ``k`` observed exceedances, a conservative
one-sided tail bound is

.. math::

   P(K \ge k) \le \sum_{i=k}^n {n\choose i}(1/2)^n.

This is a sign-test reference. It uses whether a job exceeds the threshold,
so one enormous spike contributes just one exceedance. It gives up information
about the magnitude of each excess. That makes it a useful conservative
baseline, rather than a claim of optimal sensitivity.

Two out of three exceedances have tail probability ``1/2`` at the null
boundary. Even three out of three have probability ``1/8``. Five out of
five reach ``1/32``, below a one-sided 5% budget. These are the strongest
conclusions this particular sign rule can obtain at those job counts;
stronger noise assumptions can support different tests.

More measurements inside one job can improve its contrast, but cannot be
counted as more independent jobs. If 100 measurements share one disturbance,
their agreement mostly tells us that the disturbance persisted.

Allowing some jobs to violate the model
---------------------------------------

Now suppose each independent job has contamination probability at most
``epsilon``. Conditional on being clean, its exceedance probability is at
most one half. A contaminated job can produce any result, including always
exceeding the threshold. Consequently,

.. math::

   p_0 = (1-\epsilon)/2 + \epsilon = (1+\epsilon)/2,
   \qquad
   P(K\ge k) \le \sum_{i=k}^n {n\choose i}p_0^i(1-p_0)^{n-i}.

The second inequality follows because independent Bernoulli trials with
success probabilities at most ``p0`` are dominated by ``n`` trials at
``p0``. This derivation is our explicit contamination extension to the
sign reference. It allows arbitrarily large corrupted timings, but still
requires independence across jobs and a declared contamination bound.

For illustration, with ``epsilon = 0.10``, the null exceedance probability
can be as high as 0.55. At a fixed one-sided 5% budget:

.. list-table:: Minimum exceedances needed
   :header-rows: 1

   * - Independent jobs
     - Clean model
     - At most 10% contamination probability
   * - 3
     - Unattainable
     - Unattainable
   * - 5
     - 5
     - Unattainable
   * - 6
     - 6
     - 6
   * - 10
     - 9
     - 9
   * - 20
     - 15
     - 16

Five contaminated-model exceedances have bound ``0.55**5 = 0.0503284375``;
six have bound ``0.55**6 = 0.027680640625``. The 10% figure is an assumption
for this calculation, not an estimate of the user's CI and not a selected
production default. It is a per-job probability bound, not a promise that
every small batch contains at most 10% bad jobs. Systematic bias in every
candidate measurement is outside this model.

This budget is for one benchmark, one fixed pair of commits, and one fixed
confirmation batch. To support suite-wide alerts, allocate the error budget
across benchmarks or use a separately justified multiple-testing procedure.
Repeatedly inspecting results and stopping when they pass this fixed-sample
test also needs a different rule. An initial noisy run can trigger a fixed
fresh confirmation batch whose decision uses only that batch; keep the
trigger run for diagnosis. Failed or missing jobs need a predeclared handling
rule so failures do not silently select favorable measurements.

Where this connects to ASV
--------------------------

The following observations refer to repository revision ``f3ec685``:

* `continuous.py <../../asv/commands/continuous.py>`_ runs both commits and
  passes their results to the comparison command.
* `run.py <../../asv/commands/run.py>`_ alternates commit traversal direction
  across interleaved rounds. With two commits and two rounds this gives an
  ``ABBA`` or ``BAAB`` sequence. The initial orientation is deterministic;
  actual measurement times are not forced to be equally spaced.
* `results.py <../../asv/results.py>`_ appends raw sample lists and recomputes
  pooled statistics. This sample representation does not retain per-sample
  job, round, and assignment identities needed by this proposed analysis.
* `compare.py <../../asv/commands/compare.py>`_ combines a point-estimate
  ratio threshold with ``_stats.is_different``. With sufficient raw samples,
  `the statistics helper <../../asv/_stats.py>`_ uses Mann--Whitney U;
  otherwise it checks confidence-interval overlap. This is not a paired
  job-level test that the slowdown itself exceeds the practical threshold.

ASV therefore already has useful collection machinery. A possible addition
would preserve pairing metadata and analyze independent job contrasts. A
first prototype should save a separate manifest containing commit IDs,
benchmark/parameter identity, environment and runner identity, job and block
IDs, assigned order, timestamps, warmup/settings, raw samples, and failures.
It can reuse existing runs without immediately changing ASV's result format.

If many benchmarks spike together, their pattern is evidence worth inspecting
for a shared job effect. Automatically dividing by the suite's median change
would assume most benchmarks are unaffected; a broad real regression could
then disappear. Rerunning the old commit provides a more direct control.
An isolated benchmark can also benefit from pairing, but persistent spikes
in its ratios point toward benchmark-specific disturbances or real effects.

How this changes the next step
------------------------------

Our recent AR(1) certificates address uncertainty about correlation in a
Gaussian historical model. They do not establish the contamination bound
above. The earlier `robustness study <robustness_report.rst>`_ also found
that positive outliers weakened a different formulation's calibration.

For this CI question, prioritize measurement design and a job-level reference
before refining the correlation endpoint calculation:

1. Inspect a representative slow job and its raw samples: identify whether
   spikes are within a sample set, across a whole benchmark, or suite-wide.
   Check whether baseline and candidate were both measured in the same job.
2. Build a separate collector/manifest for balanced, randomized comparisons.
   Decide the target claim, confirmation size, missing-job handling, and
   suite-wide error budget before using results to make alerts.
   Start with A/A comparisons: run the same commit under both labels. Any
   measured difference then comes from the collection setup or environment.
   This diagnoses baseline instability and order effects; a small A/A batch
   does not establish a universal contamination rate or cover interactions
   caused by genuinely different code.
3. Freeze a study comparing historical absolute timings, ASV's existing
   comparison, and paired job evidence. Include shared 50% slowdowns,
   isolated spikes, drift, unequal slot durations, version-dependent machine
   effects, and dependent jobs. Include no change, exactly-at-threshold
   changes, and larger regressions. Assumption violations should be labeled
   as stress cases, not included in a claimed guarantee.
4. Evaluate false alerts, missed regressions, inconclusive outcomes, and CI
   cost on fresh data, then decide whether ASV integration is useful.

The current deliverable establishes cancellation identities and exact
fixed-sample budgets. Detection power and actual CI contamination remain
unmeasured. The implemented tail calculation supports an alert or insufficient
evidence; a three-way decision also needs a separately budgeted upper bound
to establish that a slowdown is below the practical threshold.

Reproduce and read further
--------------------------

Print the exact-reference table and run the deterministic checks::

    .venv/bin/python -m experiments.step_detection.ci_noise
    .venv/bin/python -m pytest -p no:rerunfailures \
        experiments/step_detection/test_ci_noise.py -q

The checks cover shared-scale cancellation, balanced drift, random-order
symmetry, failure under unequal spacing and dependent repeats, and exhaustive
enumeration of clean/contaminated outcomes. No CI workload or Monte Carlo
experiment is run by these commands.

The `ASV tuning guide <https://asv.readthedocs.io/en/stable/tuning.html>`_
describes sampling controls for short- and long-time variation. NIST's
`randomized block design discussion
<https://www.itl.nist.gov/div898/handbook/pri/section3/pri332.htm>`_ explains
the general idea of grouping measurements to remove nuisance effects.
Its `sign-test reference
<https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/signtest.htm>`_
describes the binomial sign calculation. Here we explicitly retain ties as
non-exceedances for a conservative one-sided bound and use the randomized
assignment argument above to connect signs to the code effect.
