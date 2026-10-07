The complete search with split-aware eigenvalue bounds
======================================================

The Sturm bound is integrated and all seven saved cases verify. There are
**no gained or lost alerts**: four previous alerts remain alerts, both lost
detections remain non-alerts, and the below-threshold example remains a
non-alert. No search reaches a work limit.

The first lost detection moves to another explanation, at correlation
**0.7890625**. The new diagnosis indicates that tighter eigenvalue estimates
alone will not remove this explanation with the current probability formula.
The next mathematical change should address that formula.

What is now implemented
------------------------

The `Sturm derivation <sturm_refinement.rst>`_ gives rigorous residual
eigenvalue lower bounds using exact tridiagonal eigenvalue counts and the
longer proposed plateau. The `reporting reference <reporting_sturm.py>`_
uses those bounds in the existing two-tilt density certificate.

The search retains the uniform-direction and endpoint-capable tail routes,
plus the existing size and shape checks. It applies the Sturm certificate
to stationary interior intervals and to candidate surviving points. An
excluded point prompts further subdivision unless a whole-interval proof
is available. A point exclusion never substitutes for interval coverage.

The confidence allocation remains 0.01, with 0.04 for reporting. The
tilts remain 1/16 and 1/8; direction counts remain 4, 8, 12, and 16.
Work limits remain 4096 visited cells and depth 16. These calculations
refine upper bounds on the same p-value and require no extra confidence
allocation. Reporting's t/F cutoffs remain numerical, as in the previous
references, under the stationary Gaussian common-variance assumptions.

The first lost detection
------------------------

This history has a generated 8% slowdown after reading 25. The detector
must also consider explanations placing the change elsewhere. Its previous
survivor proposed a change after reading 1 and correlation 51/64.

The Sturm bound rejects that explanation, both at the point and throughout
its saved neighborhood. The complete search then finds another survivor
at split 1 and rho=101/128, or 0.7890625. It passes the size and shape checks
and has a probability upper bound of 0.012103, above the 0.01 cutoff.
The search visits 12 cells, versus 11 in the preceding spectral replay.

This explains why several successful local improvements have not yet
changed the history's final decision. Each improved bound removes a
particular explanation. An alert requires removing every explanation,
including the next one found after subdivision.

Where the new survivor loses accuracy
-------------------------------------

At this new witness, keep the tilt tau=1/8 and the existing density formula.
The diagnosis compares the following eigenvalue inputs:

.. list-table:: New witness at rho=101/128
   :header-rows: 1

   * - Calculation
     - Bound or numerical estimate
   * - Current rigorous Sturm certificate
     - 0.012103
   * - Same density formula using numerical projected eigenvalues
     - 0.012076
   * - Numerical directional tail probability
     - 0.0055051

The first row is rounded upward from a rational certificate. The last two
rows are numerical diagnostics. The tail estimate has a quadrature error
estimate of about 1.4e-9, which is not a rigorous error certificate.

The full numerical decomposition checks both current tilts and all four
direction counts. Its smallest bound using the projected eigenvalues is
still about 0.012076. Replacing the eigenvalue lower bounds with more
accurate ones therefore appears insufficient with those fixed choices.
The much larger gap is between the density formula and the directional
probability itself.

This gives us a reason to stop concentrating on eigenvalue accuracy at
this witness. The current formula controls the tilted density by replacing
a selected group of coefficients with their weakest lower bound and then
using a global density maximum. A sharper probability bound should retain
more of the coefficient information or bound the negative tail more directly.
Any such change still needs a uniform interval proof.

The other six cases
--------------------

The second lost detection retains split 1 and rho=7/8. Its Sturm probability
bound is about 0.04603, while the numerical directional probability estimate
remains about 0.01682. Even ideal evaluation of the same directional test
would therefore appear to retain this explanation at the current cutoff.

The archived independent and positively correlated 8% examples remain
certified alerts. Deterministic 20% changes with 24 and 100 readings also
remain certified alerts, and the deterministic 4% example keeps its middle
split with zero correlation. All four alert partitions and their routes
are identical to the preceding replay.

None of the seven runs needs a whole-interval certificate through the new
Sturm route: inherited routes cover the four alerts, and the early lost
case reaches its new surviving point first. Separate saved proofs and tests
cover the valid Sturm interval around the removed 51/64 explanation.

Checks and reproducibility
---------------------------

The `verifier <sturm_verification.py>`_ reconstructs residual states
independently of the reporting fits, checks new point and interval proofs,
and delegates inherited routes to the previous verifier. Every alert must
cover every split's full correlation domain without gaps or overlaps.
Surviving explanations are independently checked with dense GLS.

The new 24 integration tests cover complete coverage, survivors, unit
conversion, calibration, degenerate data, work limits, and damaged proofs
or metadata, including changed eigenvalue precision. Disabling the Sturm
route recovers the original directional search exactly after removing the
additional metadata. The preceding mathematical checkpoint provides 28
tests of the Sturm arithmetic and eigenvalue bounds. The detector and
verifier share the bound arithmetic; those mathematical tests provide
separate matrix and exact-root checks.

The `summary <data/sturm_reporting_v1_summary.json>`_ contains all seven
decisions and the new-witness decomposition. The compressed
`full archive <data/sturm_reporting_v1_diagnosis.json.gz>`_ retains inputs,
previous and current results, certificates, source hashes, and numerical
diagnoses. The archive builder checks the inputs and previous results
against the preceding archive, re-verifies every result and saved witness
interval, and checks compression round trips. Commands are in the
`running guide <running.rst>`_.

The replay finished before a server interruption. Its saved result files
and source hashes were checked on resumption; the search was not repeated.
This is development evidence on saved cases. The fresh-study detection count
remains 58/144, and ASV's production detector is unchanged.

Next, derive a tighter tilted-density or quadratic-tail bound using the
101/128 witness to distinguish the benefit from additional eigenvalue work.
Check its interval certificate before another complete replay. The second
lost history remains a separate calibration question.
