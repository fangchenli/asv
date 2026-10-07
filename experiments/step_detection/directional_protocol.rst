Frozen evaluation of directional confidence
===========================================

Evaluate the `complete directional reference <directional_reporting.rst>`_
on fresh histories. Freeze and commit the implementation, protocol,
calibration, and work limits before generating evaluation observations.
The saved correlated example is development evidence only.

Design and denominators
------------------------

Use stationary Gaussian AR(1) errors with one common marginal variance and
at most one change in the mean. Retain correlations -0.5, 0, and 0.7;
lengths 40 and 100; changes 0%, 4%, 5%, 6%, and 8% from baseline 10;
marginal noise standard deviations 0.05 and 0.2; and early, middle, and
recent locations. Locations are n/4, n/2, and n-max(4,n/10).

Use new seed labels 1600 and 1601 and the SHA-256 random-stream prefix
``directional-ar1-study``. Three correlation versions of a base history
share their Gaussian initial draw and innovations. Different base histories
have different streams. Before freezing, tests use seed label 12345 only.

There are 120 base histories and 360 correlated versions: 144 slowdowns
above 5% and 216 null histories, including 72 exactly at 5%. Each correlation
condition contains 48 positives and 72 nulls. The correlation versions are
dependent. Rates and small subgroup differences are descriptive; the study
does not precisely measure small differences in false-alert probability.

Comparisons
------------

Reuse the complete frozen indexed-study pipeline. Add uniform directional
confidence alone and uniform confidence combined with the tilted tail bound.
Both run the same complete size/shape search. Main comparisons are combined
versus uniform-only, each directional variant versus indexed-jump and
conditional confidence, and each versus the known-correlation oracle with
reporting alpha=0.04. Retain global full-history and indexed-original results
as secondary controls, together with the inherited reporting baselines.

Record standalone evidence and evidence gated by the same shared-floor
reporter. The primary reporting tables use gated decisions, with any
standalone differences stated explicitly. Only the oracle receives the
generating covariance. Other methods receive the observations and frozen
calibration; generating parameters are used afterward for diagnostics.

Fixed choices
--------------

Keep the 5% slowdown threshold, total alpha=0.05, confidence alpha=0.01,
and reporting alpha=0.04 with its existing 80/20 size/shape allocation.
Both directional bounds concern the same p-value and share the confidence
budget. Do not combine their rejections with predictive-mixture rejections.
All unknown-correlation methods keep max_cells=4096 and max_depth=16.
Unresolved searches are non-alerts and count as misses on positive histories.

The new tail bound keeps tilt=1/16 and 32-bit outward rounding of its radius
bound. The inherited predictors keep external unit=1 and all their previous
priors. Numerical t/F cutoffs are unchanged; rational reporting certificates
remain conditional on those cutoffs. No parameter or work limit may be
chosen from the evaluation outcomes.

Freeze source and protocol hashes, the inherited configuration chain,
analytical calibrations, six covariance matrices, and Python/NumPy/SciPy
versions. Use six worker processes with numerical library thread counts
fixed to one. Preserve deterministic result order and per-method elapsed
times; these timings describe the exact reference implementation under
parallel load, not production detector performance.

Validation and outputs
-----------------------

Before freezing, check stream separation and paired AR(1) innovations,
design counts and exact-threshold labeling, inherited pipeline reuse,
generating-pair diagnostics even after early termination, gate subsets,
comparison arithmetic, and rejection of altered frozen settings. Retain
the completed mathematical and interval-certificate tests.

Save all inputs, covariances, decisions, calibrated bounds, witnesses,
interval certificates, work counts, and generating-pair rejection routes.
For each directional method, compute point bounds at the generating pair
independently of whether the search reached that location. An alert must
reject that pair. A flat history's designated location represents equal
levels on both sides and is still a valid pair.

After evaluation, regenerate every input exactly and reproduce summaries.
Reconstruct all directional interval certificates, including partial ones,
and check surviving witnesses against separate dense GLS fits. Use the
existing verifier for inherited certificates and witnesses. Preserve
unresolved cases and all six unknown-correlation statuses. Archive numeric
inputs and full results in compressed shards with hashes and round-trip
checks. Verification caches must match complete inputs, records, and sources.

Report detection and false-alert counts by correlation, length, change,
noise amplitude, and location. Report paired gains and losses as well as
net differences. Separate surviving explanations from incomplete searches
when diagnosing misses. No changes to methods or settings are permitted
in response to these observations; subsequent changes require a new freeze
and fresh evaluation data.
