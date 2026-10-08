Connecting the saved calculation to its matrix determinant
===========================================================

The determinant trace checks that each arithmetic operation encloses its
answer. This extension supplies the missing mathematical connection: the
answer is the determinant of the covariance matrix used by the residual
test, after its saved shift and scaling.

The main theorem, ``residual_determinant_eq_trace``, covers the same saved
100-observation example, proposed split 1, and full correlation interval.
It states that the recorded expression equals::

    det(a I + b U' R U)

Here ``R[i,j] = rho^|i-j|`` is the AR(1) covariance matrix. ``U`` contains
98 orthonormal directions perpendicular to the two plateau indicator
columns. Multiplying by U' removes the two unknown plateau levels.
``a = 3/4`` and ``b`` is the exact coefficient from the archived radius
bound and tilt ``1/8``.

The theorem applies to any such U; it does not depend on a numerical basis
chosen by Python. ``saved_residual_enclosure`` combines the equality with
the trace proof to enclose this matrix determinant in the archived interval.
The determinant enclosure itself is no longer a hypothesis.

The chain of equalities
-----------------------

The implementation avoids constructing a dense 98-by-98 determinant.
Lean now checks the reasons that its smaller calculation gives the same
answer.

1. **Covariance and precision.** The tridiagonal matrix A has endpoint
   diagonal entries 1, interior entries ``1 + rho²``, and neighboring
   entries ``-rho``. Lean proves ``A R = (1-rho²) I`` by checking the
   endpoint and interior row formulas. It also proves ``det(A) = 1-rho²``.
   This identifies the inverse with the intended AR(1) covariance.
2. **The shifted tridiagonal matrix.** Let ``q = 1-rho²`` and
   ``T = a A + b q I``. Lean verifies that the saved diagonal and
   off-diagonal expressions define this T.
3. **Pivots.** The recorded pivot and multiplier equations imply
   ``T = L diag(pivots) L'``. The lower triangular matrix L has ones on
   its diagonal, so its determinant is 1. Consequently the product of
   the saved pivots equals ``det(T)``. Their interval bounds prove that
   every pivot is positive.
4. **Solves.** The recorded forward and backward equations imply
   ``T Y = X``, where X contains the two plateau indicators. Positivity
   makes T invertible, so the recorded Y really equals ``T^-1 X``.
5. **The two-by-two matrix.** Write ``M = a I + b R``. Matrix algebra
   gives ``det(M) = det(T)/q`` and
   ``M^-1 X = (X - b q Y)/a``. Lean checks that the saved sums over the
   two plateaus assemble exactly ``X' M^-1 X``.
6. **The residual determinant.** A block-matrix identity gives::

       det(U' M U) = det(M) * det(X' M^-1 X) / det(X' X).

   Lean proves this identity, including the normalization for columns of
   X that are not unit vectors. For this saved split, the plateau lengths
   are 1 and 99, so ``X' X = diag(1, 99)`` and its determinant is 99.
   This explains the denominator in the final scalar calculation.

The positive-definiteness proofs establish the invertibility used along
this chain. They are derived from positive pivots, the stationary correlation
interval, and the positive shift coefficients.

Where the proofs live
----------------------

``Tridiagonal.lean`` proves the factorization, determinant, and solve rules
for arbitrary finite tridiagonal systems. ``MatrixIdentity.lean`` proves
the block determinant identity. ``PrecisionIdentity.lean`` proves the
shift and inverse formulas, and ``MatrixReindex.lean`` handles changing
matrix index types without changing the result.

``export_matrix.py`` reads the frozen calculation's pivot and solve arrays
while collecting the existing trace. The generated ``SavedMatrix/`` modules
check their equations, products, and plateau sums. The exporter does not
modify the calculation or its saved evidence.

``AR1.lean`` identifies the covariance entries, ``SavedPositive.lean`` proves
positivity, and ``SavedResidual.lean`` assembles the final equality and
enclosure. The correlation input, the two coefficients, and the recorded
fractions are the same ones checked by the previous pilot.

What this establishes for the probability bound
-----------------------------------------------

The radius used to scale the covariance has the form
``(1-rho²) * ss * denominator / numerator``. Here ``ss`` is the squared
length of the observed residual direction before normalization, and
``numerator / denominator`` is its residual quadratic form under the AR(1)
precision numerator. For a nonnegative correlation interval, Lean proves
that lower-bounding ``numerator`` and upper-bounding ``denominator`` gives
the radius bound used by the certificate. It uses ``1-rho² ≤ 1+right`` on
``0 ≤ rho ≤ right ≤ 1``. This is the inequality in
``residualRadius_le_of_ratio_bounds``. For the saved interval, Lean uses the
tighter ``1-rho² ≤ 1-left²`` in
``residualRadius_le_of_interval_ratio_bounds``.

For the archived split, the numerator is cubic and the denominator is linear.
The Lean module ``SavedRadius.lean`` checks their exact power-to-Bernstein
identities, proves the coefficient bounds, and combines them with the radius
inequality to recover the saved radius upper bound over the whole interval.
``export_radius.py`` reconstructs those coefficients from the frozen history
and checks their provenance hashes. This closes the radius step for this
formula and interval bound for this example; it does not yet prove the
general link from the projected-covariance radius to this precision-quadratic
formula, or derive the residual polynomials in Lean from arbitrary histories.

``cutoff_with_residual_matrix`` now takes a probability inequality expressed
using the actual residual matrix determinant. It derives the determinant
lower bound and combines it with the saved scalar arithmetic to conclude
``p < 1/100``.

The probability inequality and the density-correction lower bound remain
explicit hypotheses. The next useful connection is the radius certificate:
prove that the saved radius upper bound is valid for the benchmark residuals
and that using it in the determinant is conservative. The spectral correction
and Gaussian probability argument are further proof obligations.

This extension covers one archived certificate throughout its correlation
interval. It does not prove the full Python detector correct or improve
its measured detection rate. The separate statistical research still needs
a sharper probability bound for the surviving correlation ``101/128``.

Reproduce and challenge the proof
----------------------------------

From the repository root, use the same pinned Lean setup as the
`pilot guide <README.rst>`_::

    .venv/bin/python experiments/step_detection/lean/verify.py \
        --lake "$HOME/.elan/bin/lake"

The verifier checks all three generated exports, builds the proofs, and
audits the final theorems' transitive axioms. In addition to the previous
two corruption checks, it doubles every covariance entry in a temporary
copy of the AR(1) proof. Lean must reject that incorrect matrix identity.

The larger finite solve proof has an explicit computation budget above
Lean's default. This permits more proof-checking work; it does not admit
unproved statements or change the permitted axioms.
