From a confidence bound to a complete reporting decision
========================================================

The new confidence bound now works inside a complete reporting reference.
It certifies the saved 8% slowdown with positively correlated noise that
the previous indexed reference missed. Turning off the new tail bound
retains the old explanation near correlation 1 and produces no alert.

This completes the implementation checkpoint. These are development
examples we already studied; the next step is a frozen comparison on fresh
histories to measure how often the improvement helps.

What the search must establish
------------------------------

Suppose a timing history seems to rise from one level to another. Before
reporting a slowdown above 5%, we check whether a smaller change plus noise
could explain the measurements. We do not know where the change happened
or how strongly neighboring measurement errors resemble each other.

Each candidate explanation therefore has two coordinates:

* A location t: readings before it share one underlying level; readings
  after it share another.
* A correlation rho: zero means independent errors, positive values mean
  neighboring errors tend to move together, and negative values mean
  they tend to alternate.

The search considers every location between readings and every stationary
correlation, -1 < rho < 1. Three kinds of evidence can rule out an explanation:

* **Size:** even allowing for uncertainty, its later level exceeds its
  earlier level by more than 5%.
* **Shape:** two constant levels with that noise model do not adequately
  explain the history; an additional boundary improves the fit too much.
* **Noise pattern:** after removing the two levels and overall noise
  amplitude, the remaining pattern is too unusual for that correlation.

An alert requires ruling out every candidate. One surviving candidate is
enough to withhold it. A search that runs out of its work allowance also
withholds the alert.

Where the new mathematics enters
--------------------------------

The `residual-direction derivation <residual_direction.rst>`_ supplies the
noise-pattern test without choosing priors for the plateau levels or noise
amplitude. Its directional density g describes how plausible a residual
direction is under a candidate model. The calibrated probability p asks how
often that model generates a direction with density this low or lower.

The implementation has two upper bounds on that same probability::

    p <= min(1, g)
    p <= tilted_tail_upper_bound

The first follows because integrating a density no larger than g over part
of a sphere of normalized area one gives probability at most g. It is cheap
and often sufficient. The `tilted tail bound <directional_tail.rst>`_ gives
a tighter answer near strong positive correlation, where the first bound
can leave an explanation alive.

Either bound below the 1% confidence cutoff rules out the noise pattern.
They share that budget because both bound the same p-value: taking their
minimum still gives an upper bound on p. This argument does not combine the
older five-component directional mixture or indexed predictive tests;
those use different statistics and remain separate comparisons.

For a true candidate under the stationary Gaussian, common-variance,
at-most-one-change model, confidence exclusion has probability at most 1%.
The inherited size and shape checks use a further 4% reporting budget.
Their union gives the existing nominal 5% per-history false-alert bound.
Searching many candidates does not multiply this budget: a false alert
requires rejecting the true candidate as well as all the others.

How the implementation connects the pieces
------------------------------------------

``reporting_directional.evidence`` is the new entry point. It reuses the
exact GLS fits from ``reporting_ar1.Models`` for the size and shape checks.
Those fits also provide the residual polynomials needed by the confidence
bounds. There is no second floating-point fit inside the confidence search.

For each candidate location, the search starts with the whole correlation
interval. It attempts a uniform-direction confidence certificate, then a
tail certificate for nonnegative intervals, then size or shape certificates.
If none covers the interval, it splits the interval in half and continues.
A midpoint that survives every check can end the search with a non-alert.

An interval certificate proves an inequality everywhere in that interval.
It does not merely sample correlations on a grid. The stored intervals for
an alert cover [-1,1] without gaps at every location; endpoint polynomial
limits establish coverage arbitrarily close to the stationary endpoints.
The endpoints themselves are not fitted stationary models.

The tail calculation retains its fixed tilt 1/16 and 32-bit outward rounding
for the radius upper bound. The existing limits remain 4096 visited cells
and subdivision depth 16. Unsupported tail intervals retain their candidates
unless another check rejects them.

``directional_diagnosis.verify`` reconstructs each certificate from the
measurements. It independently builds the confidence polynomials from
centered observations, checks every reported inequality and interval, and
requires complete coverage for alerts. For a surviving explanation, it also
checks size and shape with a separate dense GLS calculation.

The confidence and polynomial inequalities use exact rational arithmetic
for the supplied observations. The inherited t/F critical values are still
computed numerically; the exact reporting certificates are conditional on
those cutoffs. The confidence allocation uses the exact binary value of
the supplied calibration rather than a rounded display value.

What the complete checks found
------------------------------

Both variants use a 5% reporting threshold and the same calibration and
work limits. Each row below was run with and without the tail bound.

.. list-table:: Complete decisions on saved histories and deterministic fixtures
   :header-rows: 1
   :widths: 46 27 27

   * - History
     - Uniform bound only
     - Uniform + tail bound
   * - Saved independent noise, 8% slowdown, 100 readings
     - Alert, 115 cells
     - Alert, 115 cells
   * - Saved correlated noise, 8% slowdown, 100 readings
     - No alert, 53 cells
     - Alert, 125 cells
   * - Deterministic 20% slowdown, 24 readings
     - Alert, 39 cells
     - Alert, 39 cells
   * - Deterministic 20% slowdown, 100 readings
     - Alert, 113 cells
     - Alert, 113 cells
   * - Deterministic 4% slowdown, 24 readings
     - No alert, 12 cells
     - No alert, 12 cells

The saved correlated history has generating correlation 0.7, an early
change at location 25, and seed 1500. The uniform-only search stops at
location 25 and rho=16383/16384. The combined search certifies 112 intervals:
98 through shape, 12 through size, and two through the new tail bound.
Together they cover every location and correlation. This closes the gap
between rejecting the one old witness and proving a complete alert.

The deterministic histories have baseline 10, a change halfway through,
and variation ``0.1*sin(1.7*i)``. They check search behavior; they do not
measure error rates under the Gaussian model. The 4% example retains the
midpoint location with independent noise as a possible explanation.

All ten results passed certificate reconstruction. The compressed
`diagnosis archive <data/directional_reporting_v1_diagnosis.json.gz>`_ stores
the observations, source hashes, decisions, witnesses, and full certificates.
Its SHA256 is::

    138ebb469320dc0c682b9d2f20e6325de4e02c76d8224d25d909e73cf1ff96b2

Using and checking the reference
--------------------------------

From the repository root::

    from experiments.step_detection.reporting_directional import evidence

    result = evidence(values)
    control = evidence(values, use_tail=False)
    print(result['status'], result['has_alert'])

Inputs are 8 to 200 finite observations. ``certified_alert`` is the only
alert status. ``surviving_explanation`` records a candidate retained by the
implemented bounds; a more accurate probability calculation might still
exclude it. ``unresolved`` records an exhausted work limit, and
``insufficient_variation`` records a degenerate residual fit.

Reproduction and test commands are in the `running guide <running.rst>`_.
The existing reference modules and ASV's production detector are unchanged.

Next: measure the gain on fresh histories
-----------------------------------------

Freeze this implementation, its calibration, work limits, and evaluation
design before generating new histories. Compare uniform-only, combined
directional, indexed-jump, conditional, and known-covariance references on
the same inputs. Keep the oracle's reporting budget at 4% for comparison.

Use the previous study's range of history lengths, change locations, noise
amplitudes, and negative/zero/positive correlations, with new seeds. Count
exactly 5% changes among the null cases. Record paired gains and losses,
false alerts, surviving explanations, and unresolved searches, including
breakdowns by correlation and change location.

The previous fresh indexed study detected 28/144 above-threshold histories
with zero false alerts among 216 null histories. Those counts still describe
that version. The present checks establish a working mathematical mechanism;
the fresh comparison will establish its practical gain.
