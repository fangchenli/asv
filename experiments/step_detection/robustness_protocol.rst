Frozen shared-floor robustness study
====================================

This protocol precedes evaluation. Reuse the current-floor/cap-1 and
shared-floor/cap-1 configurations selected at the 0%, 1%, and 5% development
false-alert budgets in the component study. Include unchanged production.
Do not adjust penalties, caps, scale estimation, or reporting thresholds.
The nominal budgets are inherited calibration targets, not guaranteed rates.

Use new seeds 400--409. Each condition covers lengths 40/100, changes
0/4/6/8 percent, marginal noise scales 0.5/2/5 percent of level 10,
correlations zero/0.7, and middle/recent changes. Exactly 5% is omitted
because this study measures errors away from the decision boundary. There
are 960 histories per condition, 480 positive and 480 negative, for 7680
histories in total. Random draws are paired across conditions where possible,
and across levels, locations, lengths, and scales. Counts are descriptive;
the histories are not independent trials.

Conditions
----------

* Gaussian control: stationary Gaussian AR(1) noise, approximated with a
  256-observation warmup from zero before retaining the measured history.
* Heavy tails: replace Gaussian innovations with Student t noise with three
  degrees of freedom, scaled to variance one. Use the same AR coefficient,
  innovation-variance normalization, and warmup. At nonzero correlation the
  marginal distribution is an AR-filtered heavy-tailed distribution.
* Positive outliers: Gaussian control plus an independent 3% chance per
  reading of a positive spike equal to ten times the declared marginal
  noise standard deviation. Spikes do not change the true underlying level.
* Random missing: Gaussian control with each reading independently missing
  with probability 30%, always retaining the first and last readings.
* Gap at change: Gaussian control with a contiguous 10%-of-history gap
  centered at the true or nominal change position. Retain endpoints. The
  nominal position is still used for flat histories.
* Variance up: multiply the second half of the Gaussian control's noise
  deviations by three, independently of the true change position.
* Variance down: multiply the first half of the Gaussian control's noise
  deviations by three. The underlying level is unaffected in both cases.
* Outliers plus missing: combine the paired positive-outlier and random-missing
  mechanisms. Missingness and outlier masks use independent random streams.

The true scale/correlation parameters generate inputs only; scorers do not
receive them. For variable-variance cases, the reported base scale is before
the factor of three. Student t innovations use a Gaussian draw divided by
the square root of a sum of three squared independent Gaussian draws.

Filtering, fitting, and scoring
-------------------------------

Use ASV's actual input filtering and weight normalization, with unit weights.
Build the same common pool of production candidates and exact independent
minima at every segment count on retained observations. Compute correlation,
the adjacent-difference scale, and the segment penalty using that retained
sequence, just as the original methods do. No gap-aware correlation correction
or variance-specific weight is added.

Map fits back to the original revision indices before calling the unchanged
production regression reporter at 5%. Boundaries inside a missing gap map
to the next retained observation for evaluation. Report two-observation
localization tolerance in original observation coordinates. Do not claim
exact revision recovery inside gaps. Endpoint retention makes before/after
levels observable, although recent changes may have very few usable readings.

Freeze the protocol, runner, inherited settings, and dependency hashes in a
commit before generating evaluation inputs. Record all generated inputs and
per-history results. Compare conditions separately, with false-alert histories,
missed positive histories, localization, segment counts, and scale estimates.
Separate flat-history false alerts from below-threshold 4% alerts. Report paired
wins and losses against the corresponding Gaussian control. Preserve failures;
do not use these histories to select a replacement configuration.

The decision is whether the existing calibration survives these changes in
measurement noise. Any new estimator or reporting rule motivated by failures
must be evaluated later on new cases, rather than fitted to this test set.
