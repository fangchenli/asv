Near-threshold detection experiment
===================================

This protocol is written before running development or held-out cases.
The question is whether a shared noise floor and a bound on noise persistence
improve final regression alerts near ASV's 5 percent reporting threshold.
Production defaults remain unchanged.

Histories have 40 or 100 observations, starting at level 10. A lasting change
occurs halfway through or with only the last 10 percent of observations
remaining (at least four observations). Its size is 0, 4, 5, 6, or 8 percent.
Noise has marginal standard deviation 0.5, 2, or 5 percent of the original
level, with Gaussian AR(1) correlation 0 or 0.7. The initial noise is drawn
from the stationary distribution; innovation variance preserves the declared
marginal variance. Each history has unit weights and consecutive revisions.

There are 120 configurations per seed. Development uses seeds 0--2, giving
360 histories. Held-out evaluation uses seeds 100--104, giving 600 histories.
Identical random draws are shared across change sizes and locations within
each size/noise/correlation/seed combination for paired comparisons. The two
location variants coincide for flat histories; both remain in the balanced
factorial design and are not independent evidence.

Zero and 4 percent changes should produce no 5 percent alert. Six and
8 percent changes should produce an alert. Exactly 5 percent is reported
separately as an ambiguous threshold case, excluded from tuning and error
rates. We report alert frequency there without calling either outcome wrong.

Candidate fits and methods
--------------------------

For each history, collect the actual production search's candidates and
the exact independent-error minimizer at every segment count. Score the
union, removing duplicate boundaries and levels. The methods are:

* Production: the current detector, including its existing penalty search.
* Current score / common pool: the unchanged scoring formula evaluated over
  the shared candidate pool. This isolates candidate search from scoring.
* Independent shared floor: the experimental score with correlation zero.
* Bounded persistence: the experimental score with a strict correlation cap.

The independent score's optimum is represented in the pool because, at a
fixed count, it increases with independent absolute error. The correlated
scores have no corresponding global guarantee. The reference retains one
independent minimizer per count; other tied fits can score differently under
correlation. Fitted levels are independent weighted medians in all methods.

For all experimental settings, estimate a scale from the median absolute
adjacent difference. Use floor factors 0.25, 0.5, and 1 times that scale, with
a minimum floor of 1e-12 in these synthetic units. This is a declared heuristic,
not a known measurement uncertainty. Persistent noise and genuine changes
can affect the estimate. The floor is shared by every candidate for a history.

Segment penalties are c*log(n)/n for c = 2, 4, or 8. Independent scoring uses
zero correlation. Bounded scoring uses maximum half-lives 1 or 4 observations,
giving caps 2**(-1/H). These make nine independent and eighteen bounded
configurations. Neither the true noise scale nor true correlation is supplied
to the scorers.

Selection and evaluation
------------------------

Choose one independent and one bounded configuration using development data.
Minimize the sum of false-alert history rate among below-threshold histories
and missed-alert history rate among above-threshold histories. A history is
positive if it produces any final alert. Break ties by lower false-alert rate,
then higher segment penalty, then higher floor, then shorter half-life.
Save the settings and source/protocol hashes before generating held-out cases.
Held-out execution reads this frozen file and rejects changed code or protocol.

Report both rates, their separate flat and 4 percent components, and missed
alerts broken down by change size, noise, length, correlation, and change
location. Boundary matches use a two-observation tolerance and one-to-one
matching; these location metrics are separate from history-level alert rates.
Also report fitted segment counts and alert rates for exactly 5 percent.
Every method uses the same production regression-reporting function at 5
percent, including its noise-dependent threshold and short-segment handling.

No settings are adjusted using held-out results. A lower selection score
alone is not an accuracy improvement. A useful recommendation must consider
misses and false alerts together. These simulations can identify limitations
and tradeoffs; they cannot establish accuracy on unlabeled real histories.
