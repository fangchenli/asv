Frozen evaluation with unknown AR(1) correlation
================================================

The `reference <unknown_correlation.rst>`_ combines a joint confidence set
for correlation/location with continuous polynomial rejection certificates.
This first evaluation asks how much sensitivity remains without oracle
correlation, and whether misses arise from surviving explanations or
unfinished certification. Freeze code, settings, and work budgets before
generating observations.

Design
------

Use stationary Gaussian AR(1) noise with constant marginal variance and
rho=-0.5, 0, or 0.7. Keep all three versions of each base history paired
through the same independent Gaussian initial draw and innovations. They
share the mean, noise amplitude, and change location. Different base histories
use distinct SHA-256-derived random streams with the prefix ``ar1-study``.

The Cartesian product is:

* Lengths 40 and 100.
* Increases 0%, 4%, 5%, 6%, and 8% from baseline 10.
* Marginal noise standard deviations 0.05 and 0.2.
* Early, middle, and recent changes, at n/4, n/2, and n-max(4,n/10).
* Fresh seed labels 1400 and 1401.

This gives 120 base histories and 360 correlated versions. Each correlation
condition has 48 positive and 72 null histories, including 24 exactly at 5%.
Exactly 5% remains null because the target exceeds 5%. Aggregate versions
are dependent. Counts in the smaller subdivisions are descriptive and cannot
establish precise false-alert rates. The mathematical guarantee is separate
from the empirical comparison.

The first half is the reference's training prefix. Early changes fall inside
training, middle changes occur at validation entry, and recent changes leave
only four or ten later readings. All cases satisfy the unknown-correlation
reference's declared model. No covariance transition reveals the mean boundary.

Comparisons and fixed settings
------------------------------

Reuse the complete frozen covariance-study pipeline, preserving production,
shared-floor reporting, sign references, the old independent-noise direct
test, and the known-correlation oracle. Add:

* The oracle with reporting alpha=0.04 and the same 80/20 size/shape allocation
  as the unknown-correlation rule. Its statistics are unchanged; only its
  critical values differ from the alpha=0.05 oracle. This isolates the
  allocation cost from the cost of covariance uncertainty.
* The unknown-correlation rule, both standalone and gating shared-floor
  reporting. Use total alpha=0.05, confidence alpha=0.01, reporting alpha=0.04,
  threshold 0.05, max_cells=4096, and max_depth=16. Do not raise work limits
  or change the predictor after inspecting evaluation data.

Do not include a plug-in estimator in this first comparison. The oracle at
the matching reporting budget supplies the immediate reference needed to
interpret uncertainty losses. The old sign references only have their
independence guarantee in the rho=0 condition.

Freeze the inherited configuration and all new source/protocol hashes,
NumPy/SciPy/Python versions, six covariance-matrix hashes, analytical
calibrations, and work budgets. Commit before generating evaluation histories.
The six oracle matrices depend only on length and rho. Save them once and
link each input to its matrix; only the oracle receives these matrices.

Diagnostics
-----------

Save full numerical inputs, all reporting outputs, confidence quadratics,
certificates, witnesses, work counts, and elapsed unknown-correlation time.
Separate certified alerts, surviving explanations, unresolved searches, and
insufficient-variation outcomes. Count every non-alert outcome as a miss
on positive histories; do not silently omit unfinished searches.

For each history, solve the saved quadratic inequalities to describe the
confidence region at every main split. Compute the width of the union over
splits, the width at the generating split, and whether the union permits
correlations arbitrarily close to one. Decimal root calculations at precision
80 are descriptive summaries only; reporting still uses exact interval
certificates. A singleton confidence interval has width zero but is nonempty.

Check whether the true (split,rho) pair is retained using exact rational
evaluation of the stored quadratic and the frozen floating-point bound.
At that pair, record confidence exclusion, size rejection, shape rejection,
both, or neither. A certified alert must reject this pair by some route.
For flat histories, the nominal split is a valid designated witness.

Report detection/false alerts and paired gains/losses against both oracle
budgets and the old direct test, by condition and by change location. Include
confidence coverage failures, widths, near-one admissibility, status counts,
and work/runtime summaries. Widths over locations describe the projected
joint set; they are not new rejection rules.

Validation and interpretation
-----------------------------

Before freezing, test the paired AR recurrence, shape registry, design counts,
reuse of the entire inherited pipeline, recalibration of oracle statistics,
quadratic-root summaries, and true-pair diagnostics using test-only seeds.
During evaluation, require agreement with the old direct test in the identity
control, gate subsets, valid witness rejection for every alert, and complete
numeric artifacts. Four or more workers may evaluate independent histories;
preserve input order and record worker/thread settings.

This reference uses exact rational arithmetic and is slower than the known-
covariance calculation. The smaller initial design limits total computation
while retaining the main statistical contrasts. It is not a runtime contest.
Report uncertainty and unresolved work explicitly. Keep completed histories
as diagnostics; any change motivated by this evaluation needs another freeze
and fresh data before a new quality claim.
