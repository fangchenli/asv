Running the step detection experiments
======================================

The harness compares fitting and penalty-selection methods while retaining
ASV's production scoring rule. Start with the `baseline findings
<baseline_report.rst>`_ for the results of the first run.

Setup and a short run
---------------------

From the repository root, create an environment and build the extension::

    uv venv .venv --python python3
    uv pip install --python .venv/bin/python -e '.[test,dev]' ruff==0.13.3

Check the reporting mathematics
-------------------------------

The `reporting derivation <reporting_uncertainty.rst>`_ has a standalone,
standard-library reference for independent observations with positive
medians and at most one true change. Its checks enumerate sign patterns and
level constraints; they do not run a new statistical experiment::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_reference.py -q

Inspect the known-boundary and unknown-boundary calculations directly::

    .venv/bin/python - <<'PY'
    from experiments.step_detection.reporting_reference import (
        known_boundary_lower,
        single_change_evidence,
    )

    for count in (20, 40):
        values = [10] * count + [12] * count
        print("readings on each side:", count)
        print("known split, lower excess:", known_boundary_lower(values, count))
        print("unknown split:", single_change_evidence(values))
    PY

Both fixed-split lower bounds are 1.5. The unknown-boundary reference returns
``insufficient`` for twenty readings per side and ``evidence`` for forty.
``feasible_splits`` uses retained-observation indices and represents uncertainty
about a single boundary. An empty confidence set returns
``incompatible_model`` and never becomes an alert. Inputs are limited to
200 retained observations; remove missing readings explicitly and do not
treat weighted or correlated readings as independent replicates.

Evaluate the jointly calibrated reporting rule
----------------------------------------------

The `reporting protocol <reporting_protocol.rst>`_ calibrates interval
constraints jointly using fair signs, then freezes them before evaluating
new Gaussian and Laplace histories. Reproduce the held-out evaluation::

    .venv/bin/python -m experiments.step_detection.reporting_study run \
        --frozen experiments/step_detection/data/reporting_v1_frozen.json \
        --output experiments/step_detection/runs/reporting_replay

The runner verifies source, protocol, and native-extension hashes. On a
different build, regenerate the calibration to a new file after verifying
the inherited component-study provenance; do not overwrite the recorded
freeze or call changed-code results an exact replay. The original freeze
command was::

    .venv/bin/python -m experiments.step_detection.reporting_study freeze \
        --inherited experiments/step_detection/data/ablation_v1_frozen.json \
        --output experiments/step_detection/data/reporting_v1_frozen.json

That command refuses to overwrite an existing file. It uses sign sequences
only, without generating evaluation histories. The frozen JSON preserves
all calibration maxima, selected ranks, cutoffs, and failure bounds. The
run saves numerical inputs, fits, existing reports, both confidence-set
results, summaries, and condition breakdowns. Infinite lower bounds are
serialized as the explicit string ``-infinity``.

Run the analytical and pipeline checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_reference.py \
        experiments/step_detection/test_reporting_signs.py \
        experiments/step_detection/test_reporting_study.py -q

Evaluate the direct Gaussian threshold test
-------------------------------------------

The `direct-test protocol <direct_protocol.rst>`_ derives fixed critical
values from Gaussian t and F distributions, without a calibration simulation.
It reuses the previous detector and reporting references on new histories::

    .venv/bin/python -m experiments.step_detection.direct_study run \
        --frozen experiments/step_detection/data/direct_v1_frozen.json \
        --output experiments/step_detection/runs/direct_replay

The original freeze command was::

    .venv/bin/python -m experiments.step_detection.direct_study freeze \
        --inherited experiments/step_detection/data/reporting_v1_frozen.json \
        --output experiments/step_detection/data/direct_v1_frozen.json

Existing freeze files and output directories are never overwritten. The
runner checks the inherited hashes, new source/protocol hashes, native
extension, and SciPy version. Reproduction requires the matching environment;
a changed build needs separately recorded provenance and calibration.

Each result includes the existing fit and all previous evidence decisions.
In ``direct``, array entry i refers to retained-observation split i+1:
``size_t`` tests the relative increase, ``lack_of_fit_f`` tests the largest
improvement from an extra boundary, and ``best_extra_split`` identifies that
boundary. ``null_splits`` contains locations not rejected by either route.
It is not a confidence interval for the change location. Infinite statistics
in deterministic cases are serialized as explicit strings.

Run the new analytical and pipeline checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_direct.py \
        experiments/step_detection/test_direct_study.py -q

Run the original harness
-------------------------

Run a small comparison::

    .venv/bin/python experiments/step_detection/harness.py \
        --output experiments/step_detection/runs/quick \
        --scenarios flat step weak --sizes 100 --seeds 0 1

Use a new output directory for each run. The harness refuses to overwrite
existing results. Generated run directories are ignored by Git.

For a Python-only run, use ``--backends python``. That backend loads the local
source without importing the ASV package and requires only the Python standard
library. The native backend loads the locally built C++ extension. Both use
isolated detector modules, so search instrumentation does not modify an ASV
module imported elsewhere in the process.

Reproduce the initial experiments
---------------------------------

The recorded baseline uses the old correlation search at revision ``d526c2e``.
To reproduce it, check out that revision in a separate worktree and follow
the setup above there. Running these commands on the current branch instead
measures the corrected weighted-median correlation fit.

The baseline uses 15 scenario families, 100 observations, four seeds, both
backends, and all five methods::

    .venv/bin/python experiments/step_detection/harness.py \
        --output experiments/step_detection/runs/baseline

The scaling sample uses flat and single-step histories at 1000 and 10000
observations::

    .venv/bin/python experiments/step_detection/harness.py \
        --output experiments/step_detection/runs/scaling \
        --scenarios flat step --sizes 1000 10000 --seeds 0 \
        --methods approximate current --timeout 120

``--max-exact-size`` defaults to 200; the unrestricted exact method is skipped
above that size. ``--timeout`` is a per-worker limit in seconds. Failed workers
or timeouts stop the run; already completed records remain in ``records.jsonl``.
A final report is written only after all requested measurements succeed.

Compare the original and corrected detector
-------------------------------------------

The version comparison replays identical inputs through revisions ``d526c2e``
and ``123405a``. Both use the same locally built extension. Run from the
repository root with Git history available::

    .venv/bin/python -m experiments.step_detection.compare_versions \
        --output experiments/step_detection/runs/version_comparison

Defaults are 15 synthetic families with seeds 0 through 11 at length 100,
six saved public NumPy histories, and flat/step histories at lengths 1000
and 10000. Each version/backend/case gets three fresh processes. Version
order alternates; timings include input preparation, step detection, graph
coordinate mapping, and the 5 percent regression-reporting rule. Imports
and diagnostic tracing are excluded. Accuracy counts each history once.

``--seeds``, ``--backends``, and ``--repeats`` narrow the run;
``--no-real`` and ``--no-scaling`` omit those corpora. ``--before`` and
``--after`` select other Git revisions. A run saves source snapshots,
input arrays, source and extension hashes, every timing sample, paired
summaries, and untimed candidate traces for changed fits.

``data/numpy_snapshot.json`` is the frozen real-data input, so ordinary
runs require no downloads. ``fetch_numpy.py`` records the six named
operations, environment, source URLs, hashes, and commit/date mappings.
The published graphs contain no uncertainty weights; these replays use
unit weights and have no labeled change-point truth. Refreshing the
snapshot changes the experiment and requires explicitly moving the old
file before running::

    .venv/bin/python -m experiments.step_detection.fetch_numpy

Validate the comparison with::

    .venv/bin/python -m pytest experiments/step_detection/test_compare_versions.py -q

Inspect the five weak-change boundary errors and representative dip/outlier
cases against the exact reference, scoring its candidates with the actual
production closure from each revision::

    .venv/bin/python -m experiments.step_detection.inspect_comparison \
        experiments/step_detection/runs/version_comparison

``--cases`` selects other saved cases of at most 200 retained observations;
``--backend`` chooses the interval backend. This is an untimed diagnostic.
It saves ``inspection-native.json`` by default. Only one independent-error
minimizer is retained at each segment count, so rescoring that frontier does
not establish a global optimum for the correlated score.

See `the version comparison <version_comparison.rst>`_ for measured results
and their interpretation.

Run the near-threshold study
-----------------------------

The `protocol <threshold_protocol.rst>`_ specifies the simulation and setting
selection before evaluation. The `results <threshold_report.rst>`_ report
the observed sensitivity/false-alert tradeoff. Run development first::

    .venv/bin/python -m experiments.step_detection.threshold_study development \
        --output experiments/step_detection/runs/threshold_development

This evaluates 360 histories and writes ``frozen.json`` with one independent
and one bounded-persistence configuration. Evaluate the 600 fresh histories
using that file::

    .venv/bin/python -m experiments.step_detection.threshold_study heldout \
        --frozen experiments/step_detection/runs/threshold_development/frozen.json \
        --output experiments/step_detection/runs/threshold_heldout

Every output directory must be new. Each phase saves all inputs, per-history
results, aggregate results, and breakdowns. Held-out execution checks the
frozen protocol, source, design, and native-extension hashes before generating
cases. In another environment, run both phases to create a freeze matching
that build. The archived original configuration deliberately pins its original
extension binary. Use ``--backend python`` in both phases to avoid the extension.

Validate the runner and its exact reference::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_threshold_study.py \
        experiments/step_detection/test_exact_reference.py -q

Compressed input and result snapshots live in ``data/threshold_v1_*.jsonl.gz``.
For example, inspect the saved held-out results without rerunning fitting::

    import gzip
    import json

    with gzip.open(
        'experiments/step_detection/data/threshold_v1_heldout_records.jsonl.gz', 'rt'
    ) as stream:
        rows = [json.loads(line) for line in stream]

Run the component and calibration study
---------------------------------------

The `component protocol <ablation_protocol.rst>`_ separates the segment
penalty, correlation cap, and floor formulation. It calibrates each scoring
family to the same development false-alert budgets, then measures the actual
rates on fresh histories. Run both phases in the same environment::

    .venv/bin/python -m experiments.step_detection.ablation_study development \
        --output experiments/step_detection/runs/ablation_development

    .venv/bin/python -m experiments.step_detection.ablation_study heldout \
        --frozen experiments/step_detection/runs/ablation_development/frozen.json \
        --output experiments/step_detection/runs/ablation_heldout

Development uses 1200 histories and evaluation uses 2400, with random seeds
disjoint from the earlier study. The held-out run evaluates the eight fixed
factorial settings and the frozen calibrated choices. It rejects changed
protocol, code, extension, or calibration-grid hashes. ``--backend python``
works when supplied to both phases. Each output directory must be new.

Validate the scoring controls and calibration with::

    .venv/bin/python -m pytest experiments/step_detection/test_ablation_study.py -q

See `the component results <ablation_report.rst>`_ for interpretation. Archived
inputs and results use ``data/ablation_v1_*.jsonl.gz``, readable in the same
way as the preceding threshold-study snapshots.

Run the frozen robustness study
--------------------------------

The `robustness protocol <robustness_protocol.rst>`_ reuses the component
study's full-correlation-cap settings at all three calibration budgets.
There is no new configuration search. Freeze settings from a completed
component-development run in the same environment::

    .venv/bin/python -m experiments.step_detection.robustness_study freeze \
        --calibration experiments/step_detection/runs/ablation_development/frozen.json \
        --output experiments/step_detection/runs/robustness_frozen.json

Then run the paired control and stress cases::

    .venv/bin/python -m experiments.step_detection.robustness_study run \
        --frozen experiments/step_detection/runs/robustness_frozen.json \
        --output experiments/step_detection/runs/robustness

The original study used ``data/robustness_v1_frozen.json``; it pins its original
sources and extension build. As with previous studies, a different build
requires a matching calibration/freeze rather than changing the archived
hashes. Supply ``--backend python`` consistently when using that backend.

The study saves all 7680 inputs and per-history results, condition breakdowns,
and paired gains/losses against the Gaussian control. Missing readings are
filtered through production, and fits are mapped back to original revision
coordinates before reporting. Validate these mechanisms with::

    .venv/bin/python -m pytest experiments/step_detection/test_robustness_study.py -q

See `the robustness results <robustness_report.rst>`_ for findings. Archived
``data/robustness_v1_*.jsonl.gz`` files preserve the input and result snapshots.

What the methods compare
------------------------

* ``exact`` and ``approximate`` use the same fixed gamma,
  ``3 * single_segment_error * log(n) / n``. This is the approximate solver's
  default rule, evaluated on the prepared observations.
* ``current`` runs the actual ``solve_potts_autogamma`` implementation, recording
  the penalties, fits, and scores it evaluates.
* ``grid`` replaces only the outer search, using exactly as many candidate
  fits as ``current`` used for that input. It covers the expanded bracket from
  the production golden-section search, then refines between sampled penalties
  that yield different partitions, favoring lower-scoring regions.
* ``hybrid`` runs current search first, then adds a grid/refinement budget equal
  to its candidate count. It keeps the best score encountered across both stages.

Grid and hybrid retain the production objective closure, including its noise
floor, complexity penalty, and residual correlation calculation. Repeated fits
reuse the fitted correlation parameter. The final residual sum still comes
from production code. The tests verify each recorded grid score against a
fresh production evaluation forced to use the same gamma.

Grid and hybrid require ``current`` in the requested methods to establish the
per-input budget. In the initial baseline, current used ten candidate fits,
so grid used ten and hybrid used twenty. Their runtime differences include
both candidate selection and reuse of correlation fits.

Correctness and detection metrics
---------------------------------

Each backend first passes 200 small-input checks against exhaustive partition
enumeration. The oracle evaluates possible segment levels directly from the
observed values, independently of ASV's median and dynamic-programming code.
It checks weights, ties, interval constraints, and subranges. The harness stops
before benchmarking if these checks disagree.

Boundary matching maximizes the number of one-to-one matches within a fixed
tolerance, then minimizes their total location error. The default tolerance is
two input positions, including positions with missing values. For missing-data
cases, truth boundaries are projected to the next retained observation.
Regression alerts are scored after mapping segments into revision coordinates
and applying ASV's existing 5 percent reporting rule.

Synthetic truth includes both physical step positions and explicit expected
regression positions. For example, a brief dip has two true boundaries but
no expected regression alert. Drift has no discrete boundary truth and is
excluded from detection accuracy counts. Even seeds form the development
split and odd seeds form the held-out split. No parameters were tuned during
the initial comparison.

Timing and artifacts
--------------------

Each method/input/backend combination runs in a fresh subprocess. Timings
cover the fit and search, including instrumentation; they exclude interpreter
startup, input preparation, and the separate fixed-gamma scale calculation.
Peak resident memory covers the entire worker, including interpreter overhead
and native allocations. It is unavailable on platforms without ``resource``.

Phase times cover dynamic programming, merging, candidate fitting, and the
correlation fit. Candidate fitting contains the dynamic-programming and
merge phases, so those times overlap. Internal C++ interval queries and cache
hits are not instrumented in this first version.

The corrected helper is recorded as ``rho_fits`` and ``rho_fit`` time.
Historical artifacts instead have ``rho_objective_evaluations`` and
``rho_search`` time from the numerical search.

Each run writes:

* ``manifest.json``: command, settings, interpreter, platform, Git state, and
  SHA256 hashes of detector and harness sources.
* ``inputs.json``: complete generated values, weights, revisions, seeds, and
  expected changes and alerts.
* ``oracle.json``: correctness-check counts and any disagreements.
* ``records.jsonl``: one result per measured method, including every sampled
  gamma, candidate fit, score, timing, memory, and detection metrics.
* ``summary.json``: aggregate metrics and paired method comparisons, including
  descriptive Wilson intervals for the rate of histories with false detections.
* ``report.rst``: a readable summary of the results.

Keep the full artifacts when investigating a changed boundary. They distinguish
an optimizer finding a lower score from a method finding a more accurate change.

Checks
------

The `exact independent reference <exact_reference.rst>`_ has a Python API
and exhaustive correctness checks::

    .venv/bin/python -m pytest experiments/step_detection/test_exact_reference.py -q

It returns the best fit at each feasible segment count and the global
shared-floor optimum. Inputs are already prepared values and weights; the
floor and beta are explicit. Histories are limited to 200 observations.

The shared-floor prototype has analytical tests that do not run the
benchmark harness::

    .venv/bin/python -m pytest experiments/step_detection/test_noise_model.py -q

See `the shared-floor design <shared_noise_floor.rst>`_ for its required
noise-scale and complexity-penalty inputs. Production uses its existing score.

The same test file covers the `bounded-persistence scorer <noise_persistence.rst>`_.
For a fixed candidate, it can be used from Python at the repository root::

    from experiments.step_detection.noise_model import (
        correlation_cap, fit_correlated_score,
    )

    result = fit_correlated_score(
        residuals=[-1, -1, -1, -1, 1, 1, 1, 1],
        weights=[1] * 8,
        segments=1,
        noise_floor=0.1,
        beta=0.2,
        rho_max=correlation_cap(4),
    )

These values illustrate the API; the half-life, floor, and beta are caller
inputs rather than calibrated defaults. The scorer imports ASV's production
correlation helper, so use the environment from the setup section.

Run the experiment tests and the production step-detection tests::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_harness.py test/test_step_detect.py -q

The experiment tests live outside the default production test directory and
must be named explicitly. The native cases are skipped if the extension is
not built; a benchmark run requesting native measurements fails in that case.
