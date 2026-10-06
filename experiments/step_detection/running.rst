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

The shared-floor prototype has analytical tests that do not run the
benchmark harness::

    .venv/bin/python -m pytest experiments/step_detection/test_noise_model.py -q

See `the shared-floor design <shared_noise_floor.rst>`_ for its required
noise-scale and complexity-penalty inputs. Production uses its existing score.

Run the experiment tests and the production step-detection tests::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_harness.py test/test_step_detect.py -q

The experiment tests live outside the default production test directory and
must be named explicitly. The native cases are skipped if the extension is
not built; a benchmark run requesting native measurements fails in that case.
