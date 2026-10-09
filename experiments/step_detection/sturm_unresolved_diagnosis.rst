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
bound, and `continuant_partition_diagnosis.py`_ reproduces it with
``--divisions 1``. The earlier eight- and 16-piece calculations also
certified, but are unnecessary for this witness. This is a post-hoc proof of
the saved witness only; it does not revise the frozen study record.

The exact polynomial and Bernstein setup takes tens of seconds for this one
case. It avoids interval subdivision, though the setup cost remains high.
Its behavior on the other three hard cases is not yet known. It has not been
used in the frozen fresh-study results or production detector.

I attempted the same full-interval calculation on a saved 100-reading hard
case. Exact polynomial construction had not completed after 120 seconds, so I
stopped the run before it could produce a certificate. This is a clear
scaling limit of the current rational-polynomial implementation. The next
algorithm work needs to reduce or avoid global exact coefficient growth
before applying this proof to the other large cases.

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
search implementation. The useful next direction is to test the full-interval
Bernstein bound on the other three hard cases and reduce coefficient growth.

I also tried computing the continuants directly in the Bernstein basis to
avoid the power-to-Bernstein conversion. The small exact coefficients matched
the power-basis result after degree elevation, but the rational recurrence had
not completed on the 40-reading case after more than a minute and was
stopped. This basis change alone does not control the large exact fractions,
so it is not kept as implementation code.
