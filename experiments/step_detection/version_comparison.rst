Original versus corrected correlation fitting
=============================================

The correction preserves every fitted step and reported alert in this sample.
It reduces typical runtime on short synthetic histories, but does not give
a consistent speedup on longer histories. We have a mathematically correct
correlation fit; we have not yet demonstrated better change detection.

Measured results
----------------

We measured 190 histories with two versions, two backends, and three repeats:
2280 fresh worker processes. Within each backend, all 190 before/after fits
and full alert records are identical. The two backends also agree on boundary
and alert locations; some accumulated errors differ at floating-point precision.

On the 180 length-100 synthetic histories, the detection counts are:

.. list-table:: Counts per version and backend
   :header-rows: 1

   * - Measure
     - Original
     - Corrected
   * - Expected lasting-slowdown alerts found
     - 108 / 108
     - 108 / 108
   * - Extra alerts
     - 0
     - 0
   * - True boundaries matched within two observations
     - 151 / 180
     - 151 / 180
   * - Missed / extra boundaries
     - 29 / 27
     - 29 / 27
   * - Total position error among matched boundaries
     - 6 observations
     - 6 observations

The 12 drift histories are excluded from these accuracy counts. Of the 29
missed boundaries, 24 belong to deliberately brief, one-observation dips.
Of the 27 extra boundaries, 24 isolate an injected outlier. Neither produces
a lasting-slowdown alert. The remaining errors involve five weak-change
histories: three misplaced boundaries and two missed changes. Thus boundary
recovery and useful regression reporting measure different behavior.

The six NumPy histories have 784--802 observations, measured during 2015--2021
(the sorting history starts in 2016). Both versions produce the same fits,
with one to nine segments, and no regression alerts on any of them. This
sample checks real-data compatibility; it supplies no labeled accuracy result
and no positive real-world alert case.

Runtime depends on history length and backend:

.. list-table:: Runtime in milliseconds
   :header-rows: 1

   * - Corpus / backend
     - Median original
     - Median corrected
     - Median paired change
   * - 180 short synthetic / native
     - 3.159
     - 2.822
     - 11.2% less time
   * - 180 short synthetic / Python
     - 12.765
     - 12.454
     - 3.3% less time
   * - 6 NumPy histories / native
     - 29.654
     - 31.296
     - 2.3% more time
   * - 6 NumPy histories / Python
     - 126.426
     - 123.638
     - 3.0% less time

The last column first divides corrected by original time for each history,
then takes the median ratio. It is not the ratio of the two corpus medians.
The middle half of short-history native ratios is 0.840--0.943; for Python
it is 0.924--1.013. These are exploratory timings from one workstation with
three samples per combination, so differences of a few percent need more
repeats before being treated as stable improvements.

For the native backend, the scaling cases show where the short-history gain
stops carrying over:

.. list-table:: Native scaling, median of three processes
   :header-rows: 1

   * - History
     - Original time
     - Corrected time
     - Original / corrected peak RSS
   * - Flat, 1000 observations
     - 36.907 ms
     - 37.206 ms
     - 27.55 / 27.81 MiB
   * - Step, 1000 observations
     - 32.992 ms
     - 35.121 ms
     - 27.30 / 27.52 MiB
   * - Flat, 10000 observations
     - 584.267 ms
     - 577.411 ms
     - 44.12 / 47.31 MiB
   * - Step, 10000 observations
     - 516.926 ms
     - 559.928 ms
     - 43.30 / 46.20 MiB

The long step takes 8.3% more time and about 2.9 MiB more process memory.
With the Python backend, the same case goes from 2.393 to 2.508 seconds
(4.8% more time), and 114.53 to 116.48 MiB. There is no measured memory saving.

What the exact reference tells us
---------------------------------

There are no changed decisions to inspect. Instead, we examined all five
weak histories with boundary errors, plus one dip and one outlier example.
For both versions, all seven selected fits have zero measured objective gap
against the exact solver at the selected penalty. They also attain the
smallest independent absolute error at their selected segment count.

That narrows the explanation: these examples do not fail because the
approximate solver returned a worse fixed-penalty solution.

Consider ``weak-n100-seed7``. Its true level rises from 10 to 10.2 halfway
through. The corrected detector actually considers a two-segment fit, but
assigns it score 2.598534; the single-segment score is 2.591155. Lower wins,
so it keeps one level. Seed 9 behaves similarly: 2.567175 for two segments
versus 2.514249 for one. Making the inner correlation minimization exact
does not reverse either preference.

We also rescored the exact reference's best independent-error fit for every
segment count, using each revision's actual production scoring function.
None improves on the selected score in these seven examples. This is a
test of a useful candidate set, not a global solution of the correlated model.

Seed 3 shows why that distinction matters. Boundaries at observations 52
and 54 tie on independent error. The reference retains 52, within two
observations of the true change at 50. Production chooses 54, whose corrected
correlated score is about 0.00221 lower. An exact independent-error solution
can therefore have better location accuracy while scoring worse under the
correlated model. Tie handling and the scoring assumptions matter too.

What to do next
----------------

Keep the correlation correctness fix, but describe its demonstrated benefit
as correctness plus a short-history runtime improvement. A broader accuracy
claim would go beyond these results.

The next statistical comparison should put lasting slowdowns near the
5 percent reporting threshold, varying noise level, history length, and noise
persistence. The current weak changes are only 2 percent, while many positive
alert cases are much larger than the threshold. Compare the existing score
with the explicit shared-floor/persistence model on that harder corpus,
selecting settings on development histories and freezing them before testing
fresh histories. Measure missed and extra final alerts together.

For performance, profile the corrected correlation helper on long histories
before redesigning it. Its sorted weighted-median calculation and careful
summation allocate additional arrays; any faster implementation must retain
the existing analytical and extreme-value tests. Exact reference results
remain the check on solver changes, rather than evidence that a lower score
necessarily identifies a truer performance change.

What this comparison measures
-----------------------------

The original detector is revision ``d526c2e``; the corrected detector is
``123405a``. The production difference is the conditional correlation fit:
the original numerical search becomes an exact weighted-median solution
on [-1,1], with more careful floating-point summation. Both revisions include
the earlier subrange fix. The C++ interval backend, approximate segmentation,
outer penalty search, noise floor, and reporting threshold are unchanged.

The shared-floor model, persistence cap below one, and independent reference
remain experimental. These measurements do not evaluate them as new defaults.

The comparison covers:

* The original 60 synthetic histories: 15 families, 100 observations, seeds 0--3.
* Another 120 histories from the same families, seeds 4--11, without tuning
  either implementation. These extend the random samples, not the noise models.
* Six full histories from the public `NumPy benchmark site
  <https://pv.github.io/numpy-bench/>`_: dot-product application, array creation,
  range creation, sorting, matrix multiplication, and shuffling. All use the
  same i7 machine and Python 3.7/Cython 0.29.21 environment. The source URLs,
  values, revision numbers, and commit/date mappings are saved in
  ``data/numpy_snapshot.json``. These are historical measurements without
  labeled change points; changed decisions are not counted as accuracy gains.
* Four scaling cases: flat and single-step histories at lengths 1000 and 10000.

Each case runs through both the Python and native interval backends. For
each version/backend/case, three fresh processes measure the public detection
and reporting pipeline. Imports are excluded. Version order alternates and
the per-case median is used. Peak RSS includes the entire worker process,
including interpreter memory. It measures process footprint, not just the
correlation helper's allocations.

Boundary matches allow two observation positions of error, using one-to-one
matching. The location-error sum counts matched boundaries only; missed and
extra boundaries are reported separately. Changes inside missing-data gaps
map to the next usable observation. Alerts use the existing 5 percent rule.
Continuous drift has no discrete truth and is excluded from accuracy counts.

Untimed diagnostics preserve the actual outer-search candidates and chosen
penalties. For the seven inspected histories, the exact solver checks the
chosen penalty, and the independent reference checks the smallest absolute
error possible with the chosen number of segments. The latter uses only its
segment-count frontier: it does not claim an optimum for the correlated score.

Reproduction
------------

See `running the comparison <running.rst>`_ for the command and environment.
``data/comparison_v1.json`` records the paired results, all individual time/RSS
samples, and source provenance; ``data/comparison_v1_inspection.json`` stores
the diagnostic candidates and scores. Complete
worker records and source snapshots are saved locally under
``runs/version_comparison_v1``. Neither detector was tuned during this run.
The machine used CPython 3.14.5 on macOS 26.6.2, arm64. Measurements were
collected on October 5, 2026, local time.

The comparison and reference validation passed 63 tests. The scoring replay
is checked against both original and corrected production closures, including
the known equal-endpoint example where the old correlation search stops early.
