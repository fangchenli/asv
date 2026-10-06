When is a slowdown convincingly above 5%?
=========================================

An estimated 6% slowdown and convincing evidence of a slowdown above 5%
are different things. The estimate might move below 5% if we repeated the
measurements. We need to measure that uncertainty before making the stronger
claim.

The proposed reporting question is: **do the plausible explanations of the
history all imply a slowdown greater than 5%?** A detector supplies a useful
description of the data. A reporting rule decides whether the evidence is
strong enough for an alert.

This note derives a conservative reference with a provable error bound for
independent observations and at most one true change. The accompanying
`reference implementation <reporting_reference.py>`_ checks every possible
split. It does not replace ASV's reporter. Its purpose is to establish what
we can guarantee before trying more powerful methods or running a new study.

.. contents:: On this page
   :local:
   :depth: 1

Start with the three-reading example
------------------------------------

In the `robustness study <robustness_report.rst>`_, one history really changed
from 10 to 10.4, a 4% slowdown. After missing readings were removed, its fitted
levels were 9.904 from 28 readings and 10.543 from just three readings. The
estimated increase was 6.45%, so the existing reporter alerted.

At a 5% threshold, the later level must exceed ``1.05 * earlier level``.
The fitted amount above that requirement is only::

    10.542673 - 1.05 * 9.903976 = 0.143499

That 0.1435 is the quantity whose uncertainty matters. The earlier level is
also estimated, so treating it as an exact baseline would omit uncertainty.

ASV currently compares the jump with the larger of its percentage threshold
and the segments' mean absolute residuals. A residual measures variation
between readings and their fitted level. It does not measure how precisely
we know the level. Three readings and thirty readings can have similar
residual variation but very different precision.

This explains why changing the segmentation penalty alone cannot settle the
reporting question. We need both a plausible fit and evidence about its size.

Define the claim before choosing a test
---------------------------------------

Let ``a`` be the true typical timing before a change and ``b`` afterward.
Because ASV fits absolute error, we take these typical timings to be population
medians. Both are positive. Let ``d = 0.05`` be the reporting threshold and
define::

    D = b - (1 + d) * a

The claim is ``D > 0``. The null hypothesis, meaning the explanations under
which an alert would be wrong, is ``D <= 0``. It includes flat histories,
improvements, 4% slowdowns, and exactly 5% slowdowns. Testing only for a
nonzero change would answer the wrong question.

Choose an error allowance ``alpha`` separately from ``d``. For example,
``alpha = 0.05`` asks for at most a 5% probability of an erroneous alert under
the declared null model. The two numbers happen to match here; one describes
the size of a slowdown and the other describes reliability of the decision.

For this reference, the guarantee concerns one complete history with at most
one change. It does not cover repeated checks as new results arrive, or the
chance of any false alert across thousands of benchmarks. Those would need
their own error allocation. It also does not establish which exact revision
caused a change.

If the boundary were already known
----------------------------------

First suppose someone chose the boundary before seeing these observations,
and each side really has a constant median. We can place an upper confidence
bound ``U_a`` on the earlier median and a lower bound ``L_b`` on the later one.
Then::

    lower bound on D = L_b - (1 + d) * U_a

Using an upper bound for the earlier timing is deliberate: a larger baseline
makes a percentage slowdown harder to establish. Alert only when this lower
bound is positive. Giving each bound failure probability at most ``alpha/2``
limits the probability that either is wrong to ``alpha``. This last step does
not require independence between the bounds; constructing each bound below
does require independent readings within its segment.

There is an exact way to bound a median without choosing Gaussian or Laplace
noise. Sort the ``m`` observations as ``x[1] <= ... <= x[m]``. For independent
continuous observations sharing a median, the number below that median has
a Binomial(m, 1/2) distribution. Choose the largest rank ``k`` satisfying::

    P(Binomial(m, 1/2) < k) <= q

Then ``x[k]`` is a lower bound and ``x[m-k+1]`` an upper bound, each with
failure probability at most ``q``. The
`NIST quantile discussion <https://nvlpubs.nist.gov/nistpubs/TechnicalNotes/NIST.TN.2119.pdf>`_
describes order-statistic bounds obtained from binomial probabilities.
Here we use discrete ranks directly, without interpolation.

Why this works: a lower bound exceeds the true median only when fewer than
``k`` readings fall below it. The binomial tail controls exactly that event.
The upper-bound argument is symmetric. Independent observations may have
different spreads while sharing the same median. Atoms or ties at a common
median make these bounds conservative rather than invalid.

For just three readings, even the minimum is above the true median with
probability ``(1/2)^3 = 1/8``. Therefore this order-statistic construction
cannot supply a finite data-derived lower bound with 95% coverage. Under the
positive-timing assumption, its fallback is zero, which cannot establish
the desired slowdown. This is a limitation of this distribution-free
construction, not a ban on inference from three readings under stronger
distributional assumptions.

At ``q = alpha/2 = 0.025``, the first usable rank appears at six observations:
``2**(-5) > 0.025`` but ``2**(-6) <= 0.025``. Six is only the point where
a finite endpoint becomes possible; it does not guarantee a narrow interval
or sufficient evidence for a 5% alert.

What noise and correlation do to precision
------------------------------------------

A large-sample calculation helps explain the scale of the problem. It is
not a substitute for the exact short-segment bounds above.

For independent observations with density ``f`` at their median ``a``, the
empirical fraction below ``a`` has variance ``1/(4*m)``. Moving the proposed
median by a small amount ``h`` changes that fraction by approximately
``f(a)*h``. Dividing the fluctuation by this slope gives::

    Var(estimated median) approximately 1 / (4*m*f(a)**2)

The density matters because a densely populated center pins down a median
more tightly. For Gaussian noise of standard deviation ``sigma``, this becomes
``pi*sigma**2/(2*m)``. For Laplace noise of scale ``b``, it becomes ``b**2/m``.
Consequently, dividing ASV's residual average by ``sqrt(m)`` is not a
distribution-independent standard error. The Gaussian and Laplace factors
already differ, even before outliers or correlation enter.

For independent Gaussian plateaus at a known boundary::

    Var(estimated D) approximately
        (pi/2) * (sigma_after**2 / m_after
                  + (1+d)**2 * sigma_before**2 / m_before)

In the three-reading example, the generator's known ``sigma = 0.2`` gives
an approximate standard error of 0.153. A one-sided Gaussian lower calculation
would be ``0.1435 - 1.645*0.153 = -0.1082``. Even this optimistic calculation
does not clear zero. It uses the known generating scale and a large-sample
formula at a sample size of three, so it is an illustration, not a calibrated
confidence result for that history.

Correlation changes the calculation itself. Define
``s_i = 1/2 - indicator(reading_i <= its true median)``. Linearizing the median
as above gives, for a stationary plateau::

    estimated median - true median approximately sum(s_i) / (m*f)

    Var(estimated median) approximately
        [m/4 + 2*sum((m-h)*Cov(s_i, s_(i+h)), h=1..m-1)] / (m*f)**2

Repeated noise excursions make neighboring signs agree, so extra readings
add less information. The covariance of these signs is what a median needs;
the familiar AR(1) effective-sample-size formula for a mean is not generally
the right formula here.

For a declared Gaussian AR(1) noise model with correlation ``rho**h``, the
Gaussian quadrant identity gives ``Cov(s_i,s_j) = asin(rho**h)/(2*pi)``.
One way to obtain that identity is to integrate the quadrant probability's
derivative ``1/(2*pi*sqrt(1-rho**2))`` from correlation zero, where the quadrant
probability is 1/4. The long-plateau variance inflation for the median is::

    1 + (4/pi) * sum(asin(rho**h), h=1..infinity)

At ``rho = 0.7`` it is about 4.109, compared with 5.667 for the corresponding
mean. Both calculations assume stationary, weakly dependent noise and long
segments. Neither justifies replacing three readings by a fractional count
inside an exact binomial confidence calculation.

There can also be correlation across the boundary. The variance of the
contrast includes::

    Var(estimated b) + (1+d)**2 * Var(estimated a)
        - 2*(1+d)*Cov(estimated a, estimated b)

With original observation positions ``t_i``, a Gaussian AR(1) model must use
``rho**abs(t_i-t_j)`` across gaps. Compressing retained observations silently
changes that model. ASV's fitted residual correlation and its experimental
correlation cap are fitting choices, not confidence bounds on the true
correlation. They cannot simply be inserted here to obtain a guarantee.

Why a selected boundary needs more work
---------------------------------------

The detector examines many explanations and chooses one that fits well.
Random fluctuations help some boundaries look more convincing than others.
Treating that selected boundary as if it had been specified beforehand
ignores the search that made it look convincing.

`Jewell, Fearnhead and Witten <https://arxiv.org/abs/1910.04291>`_ develop
inference after changepoint detection, including methods for L0 segmentation.
Their work addresses testing for a mean change after selection. Applying
that machinery to ASV's weighted medians, correlation fitting, and relative
5% threshold would require another derivation.

There is a second issue: an incorrectly placed split can put two real levels
inside one proposed segment. A confidence formula assuming one common median
then describes the wrong model. Correcting only for the number of tested
boundaries does not repair that assumption.

A reference that includes the unknown boundary
-----------------------------------------------

We can handle both issues by keeping all plausible one-change explanations.
The following construction is deliberately simple enough to audit.

Assume independent observations, positive plateau medians, unit weights,
and at most one true change. A missingness mask must be fixed or independent
of the measurements. Count retained readings; a revision span is not a
sample size. The one-change assumption describes the underlying history,
not merely the number of segments selected by a detector.

For ``n`` retained readings there are ``M = n*(n+1)/2`` nonempty contiguous
intervals. Give each endpoint of each interval failure allowance
``q = alpha/(2*M)``. Construct its median interval ``C[I]`` using the ranks
above. If there is no usable rank, use ``[0, infinity)`` from positivity.

Some intervals cross the true change and do not have a common median. We
make no coverage claim for those intervals. For every interval wholly inside
a true plateau, however, its interval covers that plateau's median with
failure probability at most ``2*q``. Adding those probabilities over at most
``M`` intervals bounds the chance of any such failure by ``alpha``. Overlapping
intervals need not be independent for this addition to hold.

Now define which explanations are plausible:

1. For a proposed plateau ``S``, intersect ``C[I]`` over every subinterval
   ``I`` inside ``S``. Call the intersection ``J[S]``. A constant level for
   ``S`` must satisfy all these constraints.
2. Keep every split ``t`` for which both ``J[before t]`` and ``J[after t]``
   contain positive levels. Also keep a no-change explanation if ``J[whole
   history]`` contains a positive level.
3. For every surviving split, compute
   ``lower(J[after]) - (1+d)*upper(J[before])``. For a no-change explanation,
   ``D = -d*a <= 0``, so it prevents an alert.
4. Alert only if the smallest contrast across all surviving explanations
   is positive. If none survive, report model incompatibility, not an alert.

Why this provides the promised bound: with probability at least ``1-alpha``,
all true constant subintervals cover their true levels. On that event the
true one-change explanation survives, wherever its boundary lies. Taking
the smallest contrast over a set containing the truth cannot exceed the
true contrast. Therefore a history with ``D <= 0`` is falsely alerted with
probability at most ``alpha``. This reasoning includes boundary uncertainty;
it never conditions on the detector's preferred split.

The reference returns all feasible splits in retained-observation coordinates.
They represent uncertainty, not a list of detected changes. Mapping them
back to revisions would preserve gaps where the exact change is unobservable.
The evidence concerns the final versus initial plateau, not ASV's full
multiple-change and recovery semantics.

This is related in spirit to
`SMUCE <https://arxiv.org/abs/1301.7212>`_, which combines tests across intervals
with piecewise constant estimation and confidence sets. Our binomial union
bound is a separate, simpler construction; it does not implement SMUCE or
inherit its power results.

The price of that guarantee
---------------------------

For ``n = 40`` and ``alpha = 0.05``, there are 820 intervals and each tail
receives about 0.00003049. A finite rank first appears at length 16. For
``n = 100``, it first appears at length 18. These lengths apply to individual
interval bounds, not a guarantee that plateaus of those lengths suffice.

The reference makes the cost visible even on noiseless examples:

* Twenty readings at 10 followed by twenty at 12 give evidence at a known
  split, but not under this unknown-boundary construction. Some surviving
  explanations have too short a plateau to bound its level.
* Forty readings at 10 followed by forty at 12 do pass. All feasible
  explanations have contrast at least ``12 - 1.05*10 = 1.5``. Multiple
  boundaries remain possible; evidence for size does not identify one commit.

This is a correctness reference, not a competitive default. A less wasteful
allocation across interval lengths or a properly calibrated multiscale
statistic could improve it. Such changes need their own coverage argument.

Even with a known boundary, strong claims near 5% can require substantial
data. For equal independent Gaussian plateaus with sigma 0.2, a baseline of
10, and a true 6% increase, the contrast above 5% is only 0.1. The large-sample
median calculation suggests about 82 readings on each side for 80% power
at a one-sided 5% error allowance. That approximation assumes known noise
and no boundary search. It explains why a stricter interpretation of alerts
will sacrifice some sensitivity; it cannot promise fewer false alerts
without any additional misses.

What to implement and evaluate next
-----------------------------------

The mathematical reference and its analytical tests are now implemented.
No new performance experiment or calibration has been run. Production
fitting and reporting remain unchanged.

The next implementation should expose estimated change size separately from
evidence for exceeding the threshold. Start with the independent, one-change
case and compare the existing reporter with this reference while holding the
shared-floor detector fixed. Preserve the actual observations and usable
counts alongside fitted steps; the current step tuple alone lacks what
inference needs. Keep the conservative reference available to check more
powerful alternatives.

Before extending a guarantee to correlated data, specify either a covariance
model with uncertainty accounted for, or explicit dependence bounds. With
unrestricted dependence, every reading could share the same noise excursion;
the number of readings alone cannot establish increasing precision. Likewise,
ASV's relative weights are not automatically independent replicate counts.
For fixed weights and independent identically distributed observations, a
linearized weighted median has variance
``sum(w**2)/(4*f**2*sum(w)**2)``; data-dependent precision weights need more
care. The reference accepts only unweighted observations.

For a practical calibrated reporter, include the entire selection and
reporting pipeline in calibration. The null family must include slowdowns
at or below 5%, varying change locations, noise levels, and dependence.
Simulating only a flat history or holding the selected boundary fixed would
miss the failure we are investigating. Resampling fitted residuals and
rerunning selection is a possible approximation, not an automatic proof
of validity for an unknown boundary and composite null.

Freeze those choices before a new evaluation, include exactly 5% cases, and
report both false alerts and missed 6%/8% changes. The previous stress set
remains evidence motivating the design, not new validation data. Extending
the target to multiple changes, recoveries, repeated publication, and many
benchmarks comes after the one-change claim is working.
