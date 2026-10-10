Diagnosis of unresolved fresh searches
=======================================

The `fresh evaluation <sturm_fresh_report.rst>`_ left nine positive candidate
searches unresolved: seven stopped at the depth limit and two at the cell
limit. This follow-up examines the saved terminal interval for each search.
It does not change settings, decisions, or the original fresh-study result.

What the check shows
--------------------

At the exact midpoint of each of the nine saved intervals, the point
certificate proves that the tail probability is below the 0.01 cutoff. This
means the unresolved outcomes are not evidence that those particular midpoint
values have large tail probabilities. The interval certificate must cover
every correlation value in the interval, though, and its upper bound can be
too loose to prove that the whole interval is below the cutoff.

I then split each saved interval once at its midpoint and certified the two
halves separately. Both halves certified in five cases. In four cases, at
least one half still did not certify. Thus more subdivision can resolve some
of the current cases, while four cases also expose looseness that persists
after a simple extra split. Those four need a tighter interval probability
bound, a more informative partition, or both.

Where the four failures occur
-----------------------------

The four remaining failures are all positive-correlation intervals near
``rho = 1``. Both tested tilts fail while enclosing the determinant, before
the Fourier-density correction is evaluated. The determinant is positive at
each interval midpoint, so the issue is the interval enclosure rather than a
nonpositive point determinant.

Increasing the dyadic arithmetic precision from 192 to 512 bits did not make
these determinant enclosures positive. On one representative 100-reading
case, subdividing a failing half into as many as 1,024 equal pieces also left
every piece's determinant lower bound nonpositive. This makes simple
precision increases or brute-force subdivision poor next steps. The likely
source is dependency growth in the interval tridiagonal LDL recurrence near
the highly correlated endpoint.

This is a local diagnosis of the nine saved witnesses, not a complete replay
with a larger work budget. A successful midpoint or half-interval check does
not itself establish that the full search would finish within its cell and
depth limits. In one case, the saved witness interval certifies, while other
parts of the search remain unresolved.

The numerical tail values included in the archive are diagnostic estimates;
they are not used as proofs. All claims above use exact rational certificate
outputs.

Reproduction and data
---------------------

Run ``python -m experiments.step_detection.sturm_unresolved_diagnosis
--output PATH`` from the repository root. The script reads the frozen v2
inputs and records, refuses to overwrite an existing output, and records the
input archive hashes. The committed `diagnosis data
<data/sturm_fresh_v2_unresolved_diagnosis.json>`_ preserves all nine
certificates and the numerical diagnostics.

The next mathematical step is to control coefficient growth in a determinant
enclosure near ``rho = 1``. A larger global search budget is unlikely to help
until this determinant bottleneck is addressed.

Continuant prototype
--------------------

An experimental `continuant implementation
<determinant_continuant.py>`_ writes the projected determinant as an exact
rational polynomial, then bounds it around the interval midpoint. It matches
the existing point determinant on small rational examples, and the tests
check its lower bound against independent dense determinants at several
points in an interval.

On one saved 40-reading hard case, a shared-radius calculation with the
continuant and paired Bernstein bounds certifies the full witness interval
without subdivision. The largest p-value upper bound is about 0.00878. The
proof also applies the rank-group density correction on the full interval.
The `saved full-interval result
<data/sturm_fresh_v2_continuant_full_interval_n40.json>`_ preserves the exact
bound, and the `partition script <continuant_partition_diagnosis.py>`_
reproduces it with
``--divisions 1``. The earlier eight- and 16-piece calculations also
certified, but are unnecessary for this witness. This is a post-hoc proof of
the saved witness only; it does not revise the frozen study record.

The exact polynomial and Bernstein setup takes tens of seconds for this one
case. It avoids interval subdivision, though the setup cost remains high.
At that checkpoint, the other three hard cases had not been checked with
this method. It was not used in the frozen fresh-study results or production
detector.

I attempted the same full-interval calculation on a saved 100-reading hard
case. Exact polynomial construction had not completed after 120 seconds, so I
stopped the run before it could produce a certificate. This is a clear
scaling limit of the rational-polynomial implementation. That motivated the
outward-rounded coefficient recurrence described below.

Shared-denominator arithmetic
-----------------------------

The polynomial recurrence now stores each polynomial as integer coefficients
over one shared denominator, reducing common factors after arithmetic. The
Bernstein conversion consumes that representation directly instead of first
turning every coefficient back into a separate ``Fraction``. Exact tests show
that the direct conversion matches the Fraction route, and the full n=40
certificate reproduces the same exact probability bound. This reduces some
rational-normalization work, but the n=100 full-interval attempt with this
backend still did not finish within 120 seconds. Shared denominators help
with normalization; they do not prevent the coefficients themselves from
growing too large.

Outward-rounded coefficient intervals
--------------------------------------

The next version rounds each coefficient outward to a fixed dyadic scale
after every polynomial operation. Each coefficient is now an interval, so
rounding cannot invalidate the enclosure. The recurrence still uses the
power basis, which avoids the costly Bernstein-basis arithmetic, and converts
the final coefficient intervals to Bernstein bounds with outward rounding.
Small exact tests verify both coefficient containment and the transformed
Bernstein bounds.

This resolves one n=100 witness over its entire saved interval. At 128 bits,
the determinant lower bound is about 2808.39 and the rank-group probability
upper bound is ``763783093383/140737488355328``, about 0.00543. The `saved
certificate <data/sturm_fresh_v2_interval_poly_n100_middle.json>`_ records
the exact result, and the `polynomial diagnostic
<interval_polynomial_diagnosis.py>`_ reproduces its numerical bounds. The same
interval succeeds at 192 and 256 bits with the same reported lower bound.
This is a post-hoc proof of one saved witness and does not alter the frozen
study record or production detector.

I also tested a genuinely stubborn 100-reading witness with an 8% early
change. The full-interval bound is positive, but its best tail-probability
upper bound is about 0.02346, above the 0.01 cutoff. Splitting it into two
equal pieces barely changes the worst piece (about 0.02345). Four pieces
improve the first three local bounds, but the last quarter still has an upper
bound about 0.02345. The reproducible `four-piece diagnostic
<data/sturm_fresh_v2_interval_poly_n100_hard_early_d4.json>`_ records these
results. These runs retained the radius bound from the original interval
while subdividing the determinant calculation. They therefore do not show
that the probability inequality itself has reached a limit.

Recomputing the radius after subdivision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The radius compares the observed residual direction with the candidate noise
covariance. It varies with correlation. In the early-change witness, the
original interval's radius upper bound is about 0.00041794. Over just the
last quarter it can be bounded by about 0.00037296. Reusing the larger value
is safe, but prevents subdivision from recovering much of the lost precision.

The diagnostic now recomputes this bound on each piece. Both the polynomial
determinant and the density correction receive that piece's radius. For
correlation rho, the existing radius formula is
``(1 - rho**2) * ss * den(rho) / num(rho)``. Here ``ss`` measures the
observed residual size, and ``num`` and ``den`` are the model's existing
residual-ratio polynomials. Bounding those factors on a smaller interval
gives a valid local radius upper bound. No new probability inequality is
needed: a smaller valid upper bound strengthens both the determinant and
the positive weights used by the density correction.

The overall report must also cover every piece. It chooses the smallest
probability upper bound across the tested tilts within each piece, then
reports the largest of those piecewise bounds. The previous top-level
``p_upper`` incorrectly reported the best piece; ``max_cell_p_upper`` and
the all-pieces status were already conservative. Schema version 2 makes both
probability fields describe the full interval and identifies the worst piece.
This is a maximum, rather than a sum, because the pieces are alternative
values of a model parameter; the certificate must hold for whichever value
is the true one.

For the early-change witness, the `two-piece result
<data/sturm_fresh_v2_interval_poly_n100_early_local_d2.json>`_ improves the
worst bound to about 0.0147004, still above the cutoff. With eight pieces,
all pieces certify. The `complete eight-piece record
<data/sturm_fresh_v2_interval_poly_n100_early_local_d8.json>`_ has worst upper
bound ``175736591167/17592186044416``, approximately 0.00998947. That is only
about 0.00001053 below the cutoff, but the exact rational comparison is
strict. The best tilt is 1/16 in every piece. This resolves the saved witness;
it does not rerun or change the original search's outcome.

The same settings (eight pieces, 128-bit polynomial coefficients, and local
192-bit radius bounds) give these results for all three hard n=100 witnesses:

.. list-table:: Full-interval upper bounds after local-radius subdivision
   :header-rows: 1

   * - Change location
     - Worst probability upper bound
     - All pieces below 0.01?
   * - Early
     - 0.00998947
     - Yes
   * - Middle
     - 0.01161550
     - No
   * - Recent
     - 0.00550398
     - Yes

These displayed decimals are rounded upward. The `middle-case record
<data/sturm_fresh_v2_interval_poly_n100_middle_local_d8.json>`_ preserves the
failure, and the `recent-case record
<data/sturm_fresh_v2_interval_poly_n100_recent_local_d8.json>`_ preserves the
second new certificate. Combined with the earlier n=40 certificate, three of
the four stubborn witnesses were resolved at this checkpoint. No full search was
rerun, so the fresh study's alert and unresolved counts are unchanged.

Cancelling the shared stationary factor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The middle case exposes another dependency loss. The factor ``1 - rho**2``
appears in both the covariance representation and the radius. Bounding those
quantities separately combines a small factor from one correlation with a
large factor from another. The resulting ratio is more pessimistic than any
single correlation requires.

Write ``q = 1 - rho**2``. The AR(1) covariance matrix is
``R = q * inverse(A)``, where A is tridiagonal with endpoint diagonals 1,
interior diagonals ``1 + rho**2``, and adjacent entries ``-rho``. U removes
the two fitted plateau levels, so the residual covariance is ``C = U' R U``.
The observed radius is ``r = q * k``, where
``k = ss * den(rho) / num(rho)``. Therefore, at every interior correlation::

    C / r = U' inverse(A) U / k.

The stationary factor cancels exactly. This changes the order of the
calculation, without changing the statistic or the probability inequality.
The new opt-in route bounds k directly and uses the covariance ``R/q``.
If K is an upper bound for k, the determinant is bounded using::

    det((1 - 2*tau) I + (2*tau/K) U' inverse(A) U).

Here tau is the existing exponential tilt. The continuant recurrence now
uses a constant shift ``2*tau/K`` in ``(1 - 2*tau) A + (2*tau/K) I``.
Previously the shift was ``(2*tau/radius_upper) * q``, while radius_upper
already included a bound on q.
The projection determinant identity and outward-rounded Bernstein bound
otherwise stay the same.

The density correction must use the same normalization. A precision
eigenvalue upper bound H gives a normalized covariance eigenvalue lower
bound ``1/H``. Dividing this by K bounds the required eigenvalue-to-radius
ratio directly. The old route used a covariance bound containing the minimum
q and a radius bound containing the maximum q. Removing that artificial
ratio strengthens both parts of the certificate.

Exact small-matrix checks cover the normalized determinant, and dense
covariance checks cover the normalized radius and eigenvalue bounds. At a
fixed correlation, the normalized and original calculations agree exactly
when supplied matching radii. The rank-group loop is also capped at the
longer plateau's residual dimension; requesting a higher rank previously
raised an exception when all available directions remained positive.

The option is ``--stationary-normalized``. Schema version 3 records this
choice explicitly. When enabled, all ``radius_upper`` values mean bounds on
``r/q`` rather than r, and both the determinant and density correction use
those units. The default continues to use the original covariance and radius.

With 128-bit polynomial coefficients, this normalization certifies all four
hard witnesses over their full saved intervals, using one piece each:

.. list-table:: Normalized full-interval certificates
   :header-rows: 1

   * - Readings and change location
     - Probability upper bound
     - Selected tilt
   * - 40, early
     - 0.00743787
     - 1/8
   * - 100, early
     - 0.00876185
     - 1/16
   * - 100, middle
     - 0.00987724
     - 1/16
   * - 100, recent
     - 0.00482682
     - 1/16

Displayed decimals are rounded upward. The exact certificates are saved for
the `40-reading early case
<data/sturm_fresh_v2_interval_poly_n40_early_normalized.json>`_ and the
100-reading `early
<data/sturm_fresh_v2_interval_poly_n100_early_normalized.json>`_, `middle
<data/sturm_fresh_v2_interval_poly_n100_middle_normalized.json>`_, and `recent
<data/sturm_fresh_v2_interval_poly_n100_recent_normalized.json>`_ cases. The
middle case previously failed with eight local-radius pieces at 0.01161550;
cancelling the shared factor brings its single-piece bound below 0.01.
The exact middle-case bound is ``2780194237603/281474976710656``.

All four certificates concern saved intervals. A complete search must also
cover every other correlation interval at every proposed step location.

Full-search fallback
--------------------

The normalized certificate is now available to the complete experimental
search through ``reporting_sturm.evidence(values, use_normalized=True)``.
``directional_normalized.py`` owns the radius calculation and combines the
polynomial determinant with the existing Sturm density correction. The
standalone interval diagnostic uses the same radius calculation.

The search first tries its existing inexpensive routes. If none can certify
a cell and subdivision has reached depth 17, it can try the normalized
certificate. The cell must lie inside ``0 <= rho < 1``. A successful proof
closes that cell and the search continues; a failed proof leaves the history
unresolved. The default remains ``use_normalized=False``.

The replay fixes the added work at eight calls per history, using the same
two tilts, 128-bit polynomial coefficients, and 192-bit radius bounds. It
retains the 4,096-cell search limit. Reaching a limit with uncovered
intervals produces no alert. The confidence cutoff comes directly from the
existing calibration, without rounding it upward to a decimal fraction.

Each successful fallback is saved as a ``direction_normalized`` certificate.
The verifier rebuilds it from independently centered data, checks exact
equality of the proof, and checks interval coverage. Alert results still
require complete coverage of ``(-1, 1)`` at every candidate step location.
Checks cover changed units, altered probabilities, invalid endpoints,
disabled bounds, missing coverage, and inconsistent work counters.

The nine-case replay is post-hoc: these histories were selected because the
frozen study left them unresolved. Its outcomes describe search completion
on those cases; they are not a new estimate of detection or false-alert
rates. The `experiment guide <running.rst>`_ contains the reproduction
command and output format.

Next mathematical route: differences at the endpoint
----------------------------------------------------

The normalized polynomial certificate requires ``right < 1``. A search can
therefore clear its interior cells and still stop at the final cell ending
at 1. The following derivation suggests an endpoint certificate; it is not
implemented or used in the replay above.

For a proposed split s, take adjacent differences within each plateau,
omitting the difference across the proposed step. Write D for this
``(n-2) x n`` difference matrix. Its rows span exactly the residual space:
constants on either plateau disappear under D.

Let ``q = 1-rho**2`` and ``T = D R_rho D' / q``. Unlike the original
covariance, this matrix has a simple finite limit as rho approaches 1:

* Every diagonal entry is ``2/(1+rho)``.
* Differences at edge positions i and j have covariance
  ``-(1-rho)/(1+rho) * rho**(abs(i-j)-1)`` for distinct edges.
* Thus T approaches the identity matrix as rho approaches 1.

For every ``rho in [l, 1)`` with ``l >= 0``, each diagonal is at least 1
and each off-diagonal magnitude is at most ``(1-l)/(1+l)``. There are at
most ``n-3`` off-diagonal entries per row, so a row-sum eigenvalue bound gives
``T >= beta I`` in the positive-semidefinite order, where
``beta = 1 - (n-3)*(1-l)/(1+l)``. Use this only when beta is positive.

This also gives a cheap determinant bound. Set ``G = D D'``. It has two
tridiagonal blocks, each with diagonal 2 and adjacent entries -1, and
``det(G) = s*(n-s)``. If U is an orthonormal residual basis and ``B = D U``,
then ``G = B B'`` and ``T = B C B'`` for ``C = U' R_rho U / q``. Therefore,
for positive a and b,

``det(a I + b C) = det(a G + b T) / det(G)
                  >= det(a G + b beta I) / det(G)``.

The inequality follows from positive-semidefinite ordering of the two
positive-definite matrices. The determinant on the right uses two scalar
tridiagonal recurrences, with no factor tending to zero at rho=1. The same
matrix ordering gives covariance eigenvalue lower bounds from beta and the
eigenvalues of G, so the density correction can use matching units.

Exact small-matrix checks confirm the difference identity, the row-sum
bound, and the two-block determinant recurrence, including its limit at 1.
The next implementation should certify the normalized radius and matching
density correction over the endpoint cell, then test complete coverage.
If it works, testing wider cells before reaching
the depth limit could also reduce repeated neighboring calls. These are
specific follow-ups to the replay; they do not change its recorded outcomes.

Taylor-model follow-up
----------------------

I also tested a degree-limited Taylor model. It keeps a few exact coefficients
around the interval midpoint and encloses all discarded terms with a
conservative remainder bound. The model passes independent exact determinant
checks on small matrices. On the 40-reading witness, degrees 2 and 4 are too
loose; degree 6 produces a determinant-ratio lower bound around 1605, but the
full-interval probability bound is still about 0.0116, above the 0.01 cutoff.
Degree 8 tightens the determinant only slightly. The paired Bernstein ratio
does prove the full interval; the Taylor remainder bound is simply too loose.

The Taylor model also does not solve the scaling problem: degree 2 took about
23 seconds on one 100-reading witness and still returned no positive lower
bound. Since its remainder gets very large in the continuant recurrence, more
Taylor terms cost more without an evident route to a useful certificate. The
small exact tests are retained, but this is not currently a candidate for the
search implementation. Those negative results motivated the outward-rounded
recurrence and normalization described above.

I also tried computing the continuants directly in the Bernstein basis to
avoid the power-to-Bernstein conversion. A Fraction-based recurrence ran for
more than a minute without finishing. I then implemented the Bernstein
recurrence with shared-denominator integers. It finished the 40-reading case
in about 43 seconds and produced the same determinant lower bound as the
power-basis path, which takes about six seconds end to end on this witness.
The direct recurrence therefore costs more without tightening the proof. It
was removed; the representation change alone does not control coefficient
growth.
