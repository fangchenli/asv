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

The next mathematical step is a different enclosure for the determinant near
``rho = 1``. A centered or polynomial/continuant formulation should reduce the
dependency in the LDL recurrence. It must be checked against exact point
values and existing interval certificates before using it in a new frozen
evaluation. A larger global search budget is unlikely to help until this
determinant bottleneck is addressed.

Continuant prototype
--------------------

An experimental `continuant implementation
<determinant_continuant.py>`_ writes the projected determinant as an exact
rational polynomial, then bounds it around the interval midpoint. It matches
the existing point determinant on small rational examples, and the tests
check its lower bound against independent dense determinants at several
points in an interval.

On one saved 40-reading hard case, it produces a positive determinant lower
bound of about 449 on a half-interval where interval LDL returns no positive
bound. Applying the existing density correction gives a p-value upper bound
of about 0.036, still above the 0.01 cutoff. The polynomial calculation took
about 28 seconds for this one interval in the local environment. The approach
therefore shows that a different enclosure can recover positivity, but it is
not yet tight or fast enough to replace the current search. It has not been
used in the fresh-study results or production detector.

The next step is to tighten and optimize the local polynomial bound, then
recheck all four stubborn cases and the other five cases before any fresh
evaluation. This prototype's exact tests pass, but the large-case result is
only an exploratory measurement.
