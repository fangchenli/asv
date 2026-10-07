From arithmetic rules to a saved determinant calculation
=========================================================

Lean now checks the full arithmetic chain behind one saved determinant
enclosure. Previously, our final certificate theorem assumed that enclosure
was valid. We now have a proof that the recorded expression lies inside it
for every correlation in the saved interval.

The calculation is the one from the first `Sturm diagnosis
<../sturm_refinement.rst>`_: 100 observations, proposed split 1, and correlation
between ``52223/65536`` and ``52225/65536``. It uses tilt ``1/8`` and the saved
radius bound. No new benchmark histories or detector decisions are involved.

What Lean follows
------------------

Imagine a worksheet where each row contains a calculation, its inputs, and
an interval that should contain the answer. Lean checks each row using the
proofs for the earlier rows.

The first matrix pivot is a small example. The calculation forms::

    rho²
    1 - rho²
    b × (1 - rho²)
    a + b × (1 - rho²)

``rho`` can be any real number in the input interval. ``a`` and ``b`` are
fixed positive coefficients from the saved calculation. At each row, Lean
proves that the recorded lower endpoint is no greater than the exact result,
and that the upper endpoint is no smaller. The last row is the first pivot.
Before taking its reciprocal, Lean also proves that its interval excludes
zero. The same process checks the rest of the calculation.

Every generated node has three declarations:

* ``b8`` is an interval of exact rational endpoints, viewed as real numbers.
* ``v8 rho`` is the exact mathematical expression for the first pivot.
* ``h8`` proves that ``b8`` contains ``v8 rho`` whenever ``rho`` is in the
  input interval.

The value expressions use exact real arithmetic. Rounding appears only in
the recorded bounds. Lean verifies that those bounds are wide enough; it
does not assume that Python rounded correctly.

The complete trace
-------------------

The exporter runs the unchanged ``determinant_interval`` function with a
recording arithmetic object. It records 1,927 calls and shares repeated
expressions, leaving 1,723 distinct nodes:

* One correlation input and five constants.
* 203 additions, 502 subtractions, and 909 multiplications.
* 102 reciprocals and one square.

Division becomes a reciprocal followed by multiplication, preserving the
original rounding after each operation. The calculation covers the
tridiagonal factorization, both linear solves, the two-by-two matrix built
from those solves, and the final determinant expression.

The generated proofs are split into 27 files so each build step is small.
``DeterminantTrace/Result.lean`` connects the last node to the archived bounds.
The final interval is approximately::

    [1602.9753070094475, 1627.4754116404897]

The actual theorem uses the archived fractions. Its statement covers every
``rho`` in the input interval, rather than checking a grid of sample values.
``cutoff_with_trace`` then supplies this lower bound to the existing
certificate theorem. Given the correction bound and the probability
inequality, it concludes that the directional probability is below 1%.

What still connects this to statistics
---------------------------------------

Lean proves that the recorded arithmetic expression stays inside the saved
interval. To call that expression the intended matrix determinant, we still
need a matrix proof of the recurrence and the projection identity. The
existing derivation explains these identities; they are not yet Lean theorems.

The saved radius bound, the spectral density-correction bound, and the
Gaussian probability inequality also retain their existing proof boundaries.
This work checks the determinant arithmetic, not the entire statistical
argument or the complete search over possible explanations.

The exporter and its connection to the Python routine are ordinary Python.
It checks frozen source hashes, checks the history archive hash, and compares
both an ordinary run and a recorded run with the saved result. Regeneration
checks every emitted file. These checks provide provenance; a formal proof
of the Python-to-Lean translation remains outside the scope.

Checks and reproduction
-------------------------

Follow the dependency setup in the `Lean guide <README.rst>`_. From the
repository root::

    .venv/bin/python experiments/step_detection/lean/verify.py \
        --lake "$HOME/.elan/bin/lake"

This regenerates the expected source in memory, checks it against the saved
files, builds the proofs, and audits their axioms. It also changes the first
pivot's enclosure to ``[0, 0]`` in a temporary proof file. Lean must reject
that false intermediate claim. The previous corrupted-probability-bound
check remains in place.

To intentionally regenerate the trace sources::

    .venv/bin/python experiments/step_detection/lean/export_trace.py

The manifest alongside the generated modules records the inputs, node
counts, archive hash, and final enclosure. All frozen research sources and
archived evidence remain unchanged.
