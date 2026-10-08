Lean pilot: checking the certificate arithmetic
===============================================

This project uses Lean to check interval arithmetic, the complete arithmetic
trace of one saved determinant calculation, its connection to the actual
residual covariance matrix, and the certificate's final probability cutoff.
It does not change ASV or improve the detection rate. Read the
`trace walkthrough <determinant_trace.rst>`_ for the arithmetic and the
`matrix connection <matrix_connection.rst>`_ for the determinant identity.

Why this helps
--------------

Our experimental detector asks whether noise could plausibly explain an
observed history. A certificate can rule out one proposed explanation by
putting its directional tail probability below a cutoff. For the saved
example, that cutoff is 1%.

The calculation involves values that we cannot represent exactly, so the
Python reference carries a lower and an upper bound instead. Every operation
must keep the true value between those bounds. For example, rounding a lower
bound upward could make the final probability bound too small and wrongly
exclude an explanation. Lean proves that our mathematical rounding rule
always rounds in the safe direction, for every input covered by the theorem.

Tests check selected examples. A Lean theorem checks the argument for all
values satisfying its stated assumptions. This is particularly useful for
sign changes, interval division, and long chains of inequalities.

What is proved
--------------

``StepDetection/Interval.lean`` defines intervals over real numbers and proves
that addition, subtraction, multiplication, reciprocal, division, and
squaring enclose their true results. Multiplication uses all four endpoint
products. Squaring uses zero as its lower bound when the interval crosses
zero. Reciprocal and division require the denominator interval to be entirely
positive or entirely negative.

It also proves that rounding the lower endpoint down and the upper endpoint
up preserves enclosure. Taking the rounding scale to be ``2^bits`` gives
the dyadic rounding used by ``directional_determinant.Arithmetic``. These
results compose: apply an operation's enclosure theorem, then apply
``outward_sound`` or ``dyadic_sound`` to its result.

These definitions describe the mathematical rules used in the Python
reference. The Python implementation and the equivalence between its integer
floor-division code and these real-number definitions have not been formally
verified. For the saved trace, Lean instead checks directly that each recorded
pair of endpoints encloses the exact operation's result. This also checks
that the concrete rounding was safe, without trusting Python's rounding code.

``StepDetection/Trace.lean`` proves that enlarging an interval preserves its
enclosure. The generated ``StepDetection/DeterminantTrace/`` modules apply
this fact and the operation rules to 1,723 nodes, covering the saved
determinant calculation for every correlation in its input interval.

The `matrix connection <matrix_connection.rst>`_ then proves that the recorded
pivots, solves, and final expression compute the intended residual determinant.
It includes the AR(1) covariance identity and the positivity needed for
matrix inverses. ``SavedResidual.lean`` combines this with the enclosure.

``StepDetection/Certificate.lean`` proves two connecting facts:

* Positive lower bounds on the determinant and density correction give an
  upper bound on their reciprocal product.
* If a squared probability bound lies below the square of a nonnegative
  upper bound, and that upper bound is below the cutoff, the probability
  is below the cutoff too.

``StepDetection/SavedCertificate.lean`` applies those facts to exact fractions
from ``../data/sturm_v1_diagnosis.json``. The selected case is
``negative-directional-ar1-study-n100-d0.08-s0.02-early-seed1601``: the interval
around correlation ``51/64``, using tilt ``1/8`` and 12 directions.

The final calculation
---------------------

The statistical derivation supplies a bound of the following form::

    p² ≤ 1 / (D × c²)

Here ``p`` is the directional tail probability for the proposed noise
explanation. ``D`` is a determinant from the covariance calculation. ``c``
is the density correction that tightens the bound using several residual
directions. Both enter the denominator: larger valid lower bounds on them
make our probability upper bound smaller.

The archived certificate contains lower bounds on ``D`` and ``c``. Lean
checks that they are positive, that the saved squared upper bound is exactly
their reciprocal product, and that the displayed bound rounds that result
upward. Its final exact comparison is::

    1355334578761 / 140737488355328 < 1 / 100

The left side is approximately 0.009630231395. All proof calculations use
the original fractions, including the much larger determinant and correction
fractions; the decimal is only for reading.

The theorem ``saved_model_cutoff`` says: if the true determinant and correction
are at least the saved lower bounds, and the displayed statistical inequality
holds, then ``p < 1/100``. Those three conditions are explicit theorem
hypotheses. Lean has checked the implication, not established those conditions
for the benchmark history.

The new ``DeterminantTrace.saved_enclosure`` theorem proves the archived
enclosure for the full recorded arithmetic expression throughout the
correlation interval. ``cutoff_with_trace`` connects it to the 1% comparison,
leaving the correction bound and the probability inequality as hypotheses.
``residual_determinant_eq_trace`` now identifies the expression with the
intended matrix determinant for any orthonormal basis perpendicular to the
two plateau indicators. ``cutoff_with_residual_matrix`` uses that identity
in the final probability comparison.

How the pieces connect
-----------------------

The research pipeline currently has these proof boundaries:

1. The saved benchmark history defines residuals and covariance matrices.
   This construction remains in Python and the mathematical writeups.
2. Determinant calculations and Sturm eigenvalue counts produce the lower
   bounds. Lean checks the complete determinant arithmetic trace, including
   every reciprocal's nonzero-denominator condition, and proves its matrix
   identity. Lean proves the generic residual-radius inequality and Bernstein
   coefficient enclosure rule, then checks the archived radius polynomials
   and their bounds for the saved example. The projected-covariance identity
   behind this radius formula, generic history derivation, and spectral
   correction remain outside Lean.
3. The Gaussian argument relates those quantities to a probability bound.
   It remains outside Lean.
4. The final scalar arithmetic converts the bound into a comparison with
   1%. This is checked here for one archived certificate.

Excluding this explanation alone does not establish a new slowdown alert.
The complete search must exclude every relevant below-threshold explanation.
The later `Sturm replay <../sturm_reporting.rst>`_ still retains a different
explanation at correlation ``101/128``. Its reported decisions are unchanged.

Reproduce the check
--------------------

The project pins Lean ``v4.33.1``, Mathlib commit
``0df444a360eaa60ab8c11dca51a86af692955474``, and all transitive dependencies
in ``lake-manifest.json``. Install `Elan <https://lean-lang.org/install/>`_
if needed. From this directory, fetch the matching public Mathlib build
cache once::

    lake exe cache get

Use the repository's Python environment, which provides NumPy and SciPy
for importing the frozen calculation. From this directory run::

    ../../../.venv/bin/python verify.py

If Elan is installed but absent from your PATH::

    ../../../.venv/bin/python verify.py --lake "$HOME/.elan/bin/lake"

The verifier checks all three generated exports against the archive, runs
``lake build``, and audits the foundational and final theorems' transitive
axiom dependencies. It also checks two deliberately broken proofs: a zero
probability upper bound and a zero enclosure for an intermediate determinant
pivot. A third check doubles the AR(1) covariance entries and requires Lean
to reject the resulting matrix identity. Dependencies and build output stay
under the ignored ``.lake/`` directory. No statistical experiments are run.

To see the individual checks::

    ../../../.venv/bin/python export_certificate.py --check
    ../../../.venv/bin/python export_trace.py --check
    ../../../.venv/bin/python export_matrix.py --check
    lake build
    lake env lean Audit.lean

``export_certificate.py`` reads the named case and records the archive's
SHA-256 hash in the generated source. Regenerate with
``python3 export_certificate.py`` after an intentional input change. The
exporter is ordinary Python; Lean verifies the emitted constants and
theorems. The source comparison ties those constants back to the saved
JSON without claiming a verified JSON parser.

The axiom audit permits only Lean's usual ``propext``, ``Classical.choice``,
and ``Quot.sound``. There are no admitted proofs or project-specific axioms.
Arithmetic uses proof-producing tactics, not ``native_decide``. See Lean's
`proof validation documentation
<https://lean-lang.org/doc/reference/latest/ValidatingProofs/>`_ for the trust
boundary of kernel-checked proofs.

What to formalize next
-----------------------

The next connection is the radius certificate: verify that the archived
radius upper bound follows from the benchmark residuals, and prove that
using it in the determinant is conservative. The spectral correction and
Gaussian probability inequality are further proof obligations.

For the statistical research, the next task is still a sharper tail bound
for the surviving ``101/128`` explanation. Lean can check the new algebra as
it develops; it does not supply that inequality or demonstrate improved
detection power on its own.
