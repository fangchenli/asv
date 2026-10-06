Component comparison at fixed development false-alert budgets
=============================================================

This protocol precedes development and evaluation. It separates three
changes bundled in the previous threshold study: segment penalty, residual
correlation cap, and noise-floor formulation. Production remains unchanged.

Use the previous threshold study's 120 configurations per seed: lengths 40
and 100, changes 0/4/5/6/8 percent, noise standard deviations 0.5/2/5 percent,
correlations zero/0.7, and middle/recent change locations. Development uses
new seeds 200--209 (1200 histories); evaluation uses seeds 300--319 (2400).
Exactly 5 percent remains unclassified. Random draws are paired across
configurations; duplicated flat-location cases remain in the factorial
design. Counts are descriptive, not independent binomial trials.

Use the same pool for every score: production's visited candidates plus
the exact minimum-independent-error fit at every segment count. Preserve
production itself as a control. The correlated objectives do not have a
global-optimality guarantee from this pool.

The fixed factorial comparison contains eight settings:

* Segment penalty c*log(n)/n with c equal to 2 or 4.
* Absolute residual correlation capped at 0.5 or 1.
* Current floor formulation or the shared-floor formulation.

The current formulation is beta*K + log(sigma_0(fit) + S), with sigma_0
copied exactly from production. The shared formulation is beta*K + g(S/(n*b)),
with the existing profiled Laplace loss g and b equal to half the median
absolute adjacent difference (at least 1e-12). Thus the floor factor includes
the profiling formula, not just replacement of a scalar constant. S is the
conditional innovation loss minimized on the declared correlation interval.

Compare settings differing in one factor while holding the other two fixed.
Do not select a factorial winner using the held-out results.

For calibration, consider four score families: each floor formulation crossed
with each cap. Search penalty coefficients 1, 2, 3, 4, 6, 8, 12, 16 within
each family. At each development false-alert budget of 0, 1, and 5 percent,
choose the coefficient giving the fewest missed positive histories subject
to meeting that budget. Break ties by fewer false alerts, then larger penalty.
Mark a family/budget infeasible if no coefficient qualifies. These budgets
are experimental comparisons, not proposed product policy.

Save the selected configurations and all protocol/source/extension hashes
and commit the freeze before generating held-out inputs. Evaluation runs
only the eight fixed factorial settings, the frozen calibrated choices,
and production. Keep the reporting function and 5 percent threshold fixed.

Report misses and below-threshold alert rates separately, localization within
two observations, flat versus 4 percent false alerts, and breakdowns by
noise, correlation, length, and change location. Compare paired histories
for wins and losses. At matched development budgets, always show the actual
held-out false-alert rates; calibration cannot guarantee identical future
rates. Do not interpolate or adjust penalties using held-out outcomes.

A useful conclusion identifies which component changes sensitivity at a
comparable false-alert rate, and where that comparison fails. Preserve inputs
and per-history outputs, including disappointing or infeasible results.
