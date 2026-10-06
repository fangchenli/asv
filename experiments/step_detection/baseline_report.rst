Initial step detection findings
===============================

The first implementation provides a reproducible harness and an experimental
penalty search. This sample does not justify replacing the production search.
The correctness checks also found a single-segment subrange bug, now fixed.

For the next work, the `mathematical analysis <mathematical_analysis.rst>`_
takes priority over extending these experiments. It identifies an exact
correlation-fit solution and separates search guarantees from modeling choices.

These measurements precede the weighted-median correlation correction.
Revision ``d526c2e`` preserves the measured implementation and harness.

Scope
-----

The baseline contains 60 generated histories: 15 scenario families, 100
observations, and seeds 0 through 3. Each history was evaluated with five
methods on Python and C++, giving 600 runs. Four drifting histories have no
discrete change labels; accuracy counts therefore cover 56 histories per method.

Even seeds were designated development data and odd seeds held out. The
methods were specified before running the comparison; no settings were tuned
on these results. This is a small synthetic sample, not a real-project replay.

Correctness
-----------

* Each backend passed 200 exhaustive-partition comparisons, including positive
  weights, ties, segment-size limits, and restricted position ranges.
* Before the fix, 47 checks per backend failed because the one-segment shortcut
  returned the whole input instead of the requested subrange. The shortcut now
  uses the requested bounds. Dedicated production tests cover both backends.
* Python and C++ produced identical final partitions and fitted levels for all
  300 baseline method/input pairs. Scaling runs also agreed on final boundaries.
* The documented approximation-gap and baseline-offset examples are covered
  by executable tests.

The experiment and production step-detection tests finished with 76 passes
and two skips for tests requiring a free-threaded Python build. Graph and
publishing integration checks also passed after installing Mercurial and
putting the virtual environment's commands on ``PATH``.

Penalty search
--------------

Current search used ten candidate fits per history. Grid used the same budget;
hybrid retained current search and added ten grid/refinement trials. Both
alternatives evaluated the actual production score. Repeated candidate fits
reused correlation-search results.

Across the 60 histories on each backend:

* Grid had no strictly lower selection scores, one higher score, and one changed
  boundary. The other selected fits were unchanged.
* Hybrid had the same selected fits and scores as current search in all cases.
* Current, grid, and hybrid all reported the 36 expected regression alerts with
  zero additional alerts under the chosen reporting policy and tolerance.

The changed grid fit is useful evidence that selection score and truth recovery
are different measures. On ``weak-n100-seed3``, the true boundary is at position
50. Current chose 54 with score 2.654456; grid chose 52 with score 2.656737.
Grid localized the change more closely while receiving a slightly worse score.
The step is below the reporting threshold, so neither produced an alert.

With the fixed two-position matching tolerance, the aggregate C++ results are:

.. list-table::
   :header-rows: 1

   * - Method
     - Missed boundaries
     - Extra boundaries
     - Missed alerts
     - Extra alerts
     - Median fit time in ms
   * - exact
     - 10
     - 16
     - 0
     - 2
     - 3.856
   * - approximate
     - 9
     - 14
     - 0
     - 2
     - 0.607
   * - current
     - 10
     - 10
     - 0
     - 0
     - 3.204
   * - grid
     - 9
     - 9
     - 0
     - 0
     - 3.878
   * - hybrid
     - 10
     - 10
     - 0
     - 0
     - 5.386

A detected boundary outside the tolerance counts as both a missed true
boundary and an extra detected boundary. Counts include physical changes
such as brief dips that intentionally have no expected regression alert.

Exact and approximate use the same fixed penalty rule; the other methods
select a penalty automatically. Their accuracy differences include that
model-selection difference, not just optimizer quality.

Fixed penalty fitting
---------------------

At the same fixed gamma, the approximate objective was higher than the exact
objective in two of 60 histories per backend: correlated noise at seed 0
(gap 0.874985) and drift at seed 2 (gap 0.130472). A third partition differed
with the same objective within numerical tolerance. No approximation scored
below the verified exact objective.

The approximate solver was substantially faster in this small-input comparison.
The next exact-solver experiment should preserve that practical advantage
while measuring how often objective gaps change useful detections.

Scaling sample
--------------

These runs used seed 0 for flat and single-step histories at 1000 and 10000
observations, on both backends, with fixed-penalty approximation and current
automatic selection. The following table shows automatic selection:

.. list-table::
   :header-rows: 1

   * - History
     - Backend
     - Seconds
     - Process peak MiB
   * - flat, n=1000
     - python
     - 0.141
     - 34.3
   * - flat, n=1000
     - native
     - 0.037
     - 30.3
   * - flat, n=10000
     - python
     - 2.895
     - 115.3
   * - flat, n=10000
     - native
     - 0.591
     - 46.6
   * - step, n=1000
     - python
     - 0.192
     - 34.1
   * - step, n=1000
     - native
     - 0.033
     - 30.2
   * - step, n=10000
     - python
     - 2.632
     - 115.5
   * - step, n=10000
     - native
     - 1.273
     - 45.3

At 10000 observations, C++ merging and boundary adjustment took about
0.23 seconds for the flat history and 0.53 seconds for the step history.
The corresponding dynamic-programming phases took about 0.21 and 0.36 seconds.
This supports profiling long interval calculations and boundary adjustment
before assuming the dynamic-programming loop is the sole bottleneck.

Timings are instrumented, single-machine observations. They exclude process
startup and input preparation. Peak memory is the whole fresh worker process,
including Python and native allocations. Internal C++ cache-hit statistics
and repeated timing distributions remain future work.

Next experiment
---------------

Keep the current search as the default. Add harder and real benchmark histories
to the evaluation corpus, and profile interval queries during merging. Test an
incremental weighted-median cost implementation or exact pruning against this
baseline before adopting either. The grid comparison can remain available as
a diagnostic when a particular history exposes a search failure.

Reproduction
------------

See `running the experiments <running.rst>`_ for setup, commands, method budgets,
and artifact definitions. The baseline and scaling commands there reproduce
the inputs and comparisons; timing values depend on the machine.

Interpreter: ``3.14.5 (main, May 10 2026, 19:20:57) [Clang 22.1.3 ]``.

Platform: ``macOS-26.6.2-arm64-arm-64bit-Mach-O``.

Parent Git revision: ``f29d047c247e9c2485cd5c83c70bbbaba49887be``.

The measured working tree included the new harness and the subrange correction.
The run manifests recorded these exact source hashes:

* Detector: ``046af1c93b263016c02f1c1f34e61895083166d61d107ab78b850130c4904982``.
* Harness: ``25fb37356dfc874a7b4cf6e9ca749ad1bf1922f413f546ecf10bf5aad73db22a``.

Full input arrays, candidate traces, paired comparisons, and run manifests are
written under ``runs/`` when the commands are executed. Those generated artifacts
are excluded from Git; this report preserves the first-run findings.
