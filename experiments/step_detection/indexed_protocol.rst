Frozen evaluation of candidate-specific AR(1) prediction
========================================================

This study evaluates the `candidate-specific construction
<baseline_jump_prior.rst>`_ on fresh histories. Freeze its sources, priors,
external unit, analytical reporting cutoffs, and work limits before generating
observations. Compare the two indexed variants separately; the saved example
does not establish that the jump prior improves detection beyond indexing.

Design and denominators
-----------------------

Use the same declared Gaussian, constant marginal variance, one-change AR(1)
model as the previous unknown-correlation study. Correlations are -0.5, 0,
and 0.7. The Cartesian product contains lengths 40 and 100; changes of 0%,
4%, 5%, 6%, and 8% from baseline 10; marginal noise standard deviations 0.05
and 0.2; early, middle, and recent locations; and seed labels 1500 and 1501.

Locations are n/4, n/2, and n-max(4,n/10). Use the new SHA-256 random-stream
prefix ``indexed-ar1-study``. The three correlation versions of each base
history share their Gaussian initial draw and innovations. Distinct base
histories use distinct streams, and no previous observations are reused.
Tests use separate seed labels, never evaluation histories before the freeze.

There are 120 base histories and 360 correlated versions: 144 positive
histories exceeding 5% and 216 null histories, including 72 exactly at 5%.
Each correlation condition has 48 positives and 72 nulls. Paired versions
are dependent. Subgroup counts and rates are descriptive; this design is
not sized to distinguish small changes in false-alert probability precisely.

Comparators
-----------

Reuse the full preceding AR(1) study pipeline, preserving production,
shared-floor reporting, sign references, the independent-noise direct test,
both known-correlation oracle budgets, and the conditional AR(1) reference.
Add three methods:

* The original full-history predictor with a global mixture over locations.
* Candidate-specific prediction with the original independent-level prior.
* Candidate-specific prediction with the midpoint/jump prior.

Each new method is measured alone and as a gate on the same shared-floor
reporter. The oracle at reporting alpha=0.04 is the main sensitivity
reference; it receives the true correlation, unlike the practical candidates.
Do not add the global jump-prior control to this evaluation: the two indexed
variants isolate whether the new prior helps after candidate-specific
prediction. The development diagnostic retains that earlier control.

Fixed parameters
----------------

Use total alpha=0.05, confidence alpha=0.01, reporting alpha=0.04 with the
existing 80/20 size/shape split, and a 5% slowdown threshold. Keep max_cells
at 4096 and max_depth at 16 for every unknown-correlation reference. Unresolved
searches are non-alerts and count as misses on positive histories.

Use external unit=1 for every full-history predictor. Keep the original
mean precision 0.0001, inverse-gamma shape 0.5, and the 25 noise-scale
components with exponents -12 through 12. The midpoint precision is 0.0002.
The new jump-scale mixture is uniform over 1, 2, 4, 8, 16, 32, and 64 noise
standard deviations. Each indexed density retains half the original global
density. Do not choose priors, the unit, or budgets from evaluation results.

Freeze the inherited configuration chain, new source/protocol hashes,
analytical calibrations, Python/NumPy/SciPy versions, and the six oracle
covariance-matrix hashes. Commit the freeze before generating histories.
Only the oracle receives the covariance matrices. Six worker processes may
evaluate independent histories, with numerical library thread counts fixed
to one and recorded in the manifest. Preserve deterministic output order.

Outputs and validation
----------------------

Save every numeric input, the covariance registry, inherited outputs, all
new densities, coefficients, certificates, surviving witnesses, cell counts,
and elapsed reference times. Record the density and confidence coefficient
by location for indexed methods. Partial certificates never count as alerts.

For every reference, check whether the generating (location,correlation)
pair is retained. When an indexed search stops before that location, compute
its fixed density separately for diagnostics; this does not change its
decision. Record the generating pair's confidence/size/shape rejection route.
Every certified alert must reject that designated pair. For a flat history,
the nominal location remains a valid representation with equal levels.

Report detections and false alerts by correlation, length, change magnitude,
noise amplitude, and location. Report paired gains/losses for full-history
versus conditional, each indexed variant versus full-history and conditional,
indexed-jump versus indexed-original, and all three new variants versus the
matching oracle. Separate certified alerts, surviving explanations, unresolved
searches, insufficient variation, and numerical confidence fallback.

Do not estimate full correlation-region widths from a grid. The saved interval
certificates and exact true-pair checks are the confidence diagnostics for
these nonquadratic regions. The inherited conditional method retains its
previous quadratic-region summaries.

Before freezing, test the new paired recurrence and design counts, inherited
pipeline reuse, indexed true-pair diagnostics after early termination, gate
subsets, summary arithmetic, and freeze rejection after a settings change.
After evaluation, regenerate inputs, reproduce summaries, check archive
hashes and covariance matrices, reconstruct alert certificates, and check
surviving witnesses against independent GLS. Keep unresolved cases visible.
No method, setting, or work budget may be changed in response to these results;
such a change requires a new freeze and new evaluation data.
