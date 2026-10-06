Can we detect smaller slowdowns more reliably?
==============================================

Yes, the experimental scoring rules detect substantially more of the
near-threshold slowdowns in this simulation. They also produce more alerts
below the requested threshold. The bounded-persistence model offers the more
conservative tradeoff of the two experimental models; these results do not
yet justify a new production default.

The baseline here is our branch's current detector, including the exact
correlation correction. The `previous comparison <version_comparison.rst>`_
found that correction preserved all decisions on its easier sample. This
experiment evaluates changes to the scoring model itself.

Results on fresh histories
--------------------------

We chose settings using 360 development histories, then committed the code,
protocol, and frozen settings as ``c31f7f6`` before running 600 held-out
histories. Of those 600, 240 have slowdowns above 5 percent, 240 are flat or
have a 4 percent slowdown, and 120 sit exactly at 5 percent. The last group
is reported separately and does not count toward either error rate.

.. list-table:: Held-out alert results
   :header-rows: 1

   * - Method
     - Missed above-threshold histories
     - Alerted below-threshold histories
     - Above-threshold alerts located within two observations
   * - Current detector
     - 147 / 240 (61.3%)
     - 0 / 240 (0%)
     - 91 / 240
   * - Independent shared floor
     - 34 / 240 (14.2%)
     - 41 / 240 (17.1%)
     - 170 / 240
   * - Bounded persistence
     - 71 / 240 (29.6%)
     - 10 / 240 (4.2%)
     - 153 / 240

The bounded model detects 169 of the 240 above-threshold histories, compared
with 93 for the current detector. That is 76 additional detected histories
and about 52 percent fewer misses, at the cost of ten additional
below-threshold histories receiving an alert. No previously detected
above-threshold history loses its alert in this sample.

The location column is stricter than detecting an affected history: an
alert has to point within two observations of the true change. Sixteen
bounded-model detections fail that location test; two current-detector
detections do. Detecting a slowdown and locating the responsible revision
remain separate quality measures.

All ten bounded-model below-threshold alerts occur on real 4 percent
slowdowns; none occur on flat histories. Here a false alert means crossing
the user's 5 percent reporting threshold incorrectly, even if a smaller
slowdown really happened. The independent model alerts on 29 of the
4 percent histories and 12 flat histories. Of its 41 below-threshold alerts,
37 occur with correlated noise.

At exactly 5 percent, the current, independent, and bounded methods alert
on 13, 59, and 40 of 120 histories respectively. These frequencies describe
behavior at the decision boundary; they are not counted as successes or errors.

Where the gains and costs occur
-------------------------------

The bounded model improves sensitivity across both tested noise-correlation
settings. The difference is especially large at the middle noise level:

.. list-table:: Current detector versus bounded persistence
   :header-rows: 1

   * - Noise standard deviation, relative to the original level
     - Current misses / positive histories
     - Bounded misses / positive histories
     - Bounded false-alert histories / negative histories
   * - 0.5%
     - 1 / 80
     - 0 / 80
     - 0 / 80
   * - 2%
     - 68 / 80
     - 11 / 80
     - 4 / 80
   * - 5%
     - 78 / 80
     - 60 / 80
     - 6 / 80

At 6 percent slowdown, misses fall from 78 to 40 out of 120 histories.
At 8 percent, they fall from 69 to 31. Recent changes, with only four or
ten observations after the change, remain difficult: misses fall from
76 to 37 out of 120. Full breakdowns by length, noise, correlation, location,
and change size are saved with the results.

The independent model fits many more noise fluctuations as steps. Across
all 600 histories it produces 395 extra fitted boundaries, compared with
88 for the bounded model and 6 for the current detector. These are fitted
boundaries, not final alerts. This helps explain why its higher sensitivity
comes with substantially more false alerts.

What the mathematical controls tell us
--------------------------------------

We give each scoring rule the same candidate pool: every candidate visited
by production, plus the exact minimum-independent-error fit at every segment
count. The independent shared-floor optimum is represented in that pool.
The correlated models have no global-optimality guarantee from this reference.

Rescoring the common pool with the unchanged current scoring rule produces
identical boundaries and alerts to production on all 600 held-out histories.
All 147 current-detector misses have a single fitted level. In this sample,
the main limitation is already present at model selection, before alert
postprocessing. A more thorough search of these independent-optimal candidates
does not resolve it.

One history makes the issue concrete: ``n100-d0.06-s0.02-r0-middle-seed102``.
The true level rises from 10 to 10.6 at observation 50, with independent noise
whose standard deviation is 0.2. Compare two candidate descriptions:

.. list-table:: Scores for that history; lower is preferred within each column
   :header-rows: 1

   * - Description
     - Current fitted residual correlation
     - Current score
     - Bounded-model score
   * - One level, about 10.255
     - 0.736
     - 3.14999
     - 2.01406
   * - Two levels, about 9.957 and 10.670, split at 50
     - -0.019
     - 3.18351
     - 1.86514

When one level describes this history, residuals tend to stay negative before
the change and positive afterward. The current model can interpret that
pattern as persistent noise, even though the generated measurement noise
is independent. Its score favors one level. The bounded experimental model
limits that explanation and chooses the two-level description.

The experiment changes several ingredients together: the shared floor,
the correlation bound, and the segment penalty. This example illustrates
their combined effect. It does not establish how much of the overall gain
comes from each ingredient. Scores from different columns also have different
normalizations and should not be compared numerically across models.

How settings were chosen
------------------------

The full `protocol <threshold_protocol.rst>`_ was written before either phase.
It covers lengths 40 and 100, middle and recent changes, five change sizes,
three noise levels, and Gaussian noise with correlation zero or 0.7. There
are 120 configurations per seed. Development uses seeds 0--2; evaluation
uses seeds 100--104. Random draws are paired across configurations, so the
600 histories are not 600 independent samples. Flat histories are repeated
under the two location labels. These are a small factorial simulation,
not estimates of error rates across real benchmark projects.

We tested nine independent and eighteen bounded settings. Each family was
selected by minimizing false-alert rate plus missed-alert rate on development
data, with predetermined tie rules. That gives the two error rates equal
weight; it does not encode a production user's cost of an unnecessary alert.

Both selected models use a segment penalty of ``2*log(n)/n`` and a shared
floor equal to half the median absolute adjacent difference. The bounded
model additionally caps absolute residual correlation at 0.5, corresponding
to a maximum half-life of one observation. The current segment penalty is
``4*log(n)/n``. The synthetic noise parameters are never passed to the scorers.

The same production alert-reporting function, including the 5 percent rule
and its noise-dependent threshold, is used for all methods. Settings were
frozen before evaluating fresh cases and were not adjusted afterward.

Recommendation
--------------

Continue with bounded persistence as an experimental candidate. The results
show a useful sensitivity/false-alert tradeoff, rather than an improvement
on both measures at once. Production remains unchanged.

The next comparison should change one ingredient at a time: correlation cap,
segment penalty, then floor. Compare sensitivity at a matched false-alert
rate, rather than letting a more permissive detector win by accepting more
false alerts. Freeze those choices before testing additional seeds and
histories with outliers, missing observations, and changing noise levels.

Reproduction and evidence
-------------------------

See `running the experiments <running.rst>`_. The frozen configuration is
``data/threshold_v1_frozen.json``. ``data/threshold_v1_results.json`` contains
the development grid, selected settings, held-out summaries, breakdowns,
source hashes, and the scored example above. Four ``threshold_v1_*.jsonl.gz``
files preserve every generated input and every per-history method result.
They can be read with Python's ``gzip.open`` and ``json.loads``.

The held-out run records freeze commit ``c31f7f6`` and verifies its protocol,
Python sources, and native-extension hashes before execution. Both phases
used CPython 3.14.5 on macOS arm64 with the native interval backend. The
experiment and exact-reference checks passed 65 tests, including agreement
between Python and native experiment results on representative cases.
