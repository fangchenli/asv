The complete search with the determinant bound
==============================================

The full-covariance determinant bound is now part of a separate complete
reporting reference. Replaying seven saved histories produces **no additional
alerts and no lost alerts**. All seven results pass certificate verification.

The first lost detection illustrates why removing one explanation was not
enough. The new bound rejects its old explanation at correlation 0.875,
but the search finds another at 0.8125. Both propose a change after the
first reading. The second lost detection keeps its original explanation.

These are development cases, including four previous alerts, the two lost
detections, and a below-threshold example. They establish integration behavior;
they do not update the fresh-study result of 58/144 detected slowdowns.
ASV's production detector is unchanged.

What changed in the search
--------------------------

The `mathematical checkpoint <directional_determinant.rst>`_ supplies an
upper bound on the same directional p-value used by the existing confidence
checks. The `new reporting reference <reporting_determinant.py>`_ takes the
minimum of the uniform, tilted-tail, and full-determinant bounds. All three
bound the same probability, so this needs no additional confidence budget.

An alert still requires ruling out every below-threshold explanation:
every split after readings 1 through n-1, and every correlation in (-1, 1).
A candidate can be rejected by its residual pattern, slowdown size, or
extra-step evidence. The total nominal false-alert budget remains 0.05,
with 0.01 allocated to confidence and 0.04 to reporting under the inherited
common-variance Gaussian AR(1) assumptions.
As before, the size and shape cutoffs are computed numerically; the
confidence certificates use rational inequalities.

The search keeps the previous order of operations. It first tries the
existing uniform and endpoint-capable tail certificates, then the reporting
checks. The determinant supplies an additional interval certificate when
those routes fail and the interval lies strictly inside (-1, 1).

Before returning a surviving explanation, the search also evaluates the
determinant bound at that proposed split and correlation. If the bound
rejects the point but no whole-interval certificate succeeds yet, the search
subdivides and continues. A rejected midpoint alone never certifies an
interval.

The work limits remain 4096 visited cells and depth 16. Exhausting either
limit returns an unresolved non-alert. Disabling the determinant reproduces
the inherited search; tests compare its results and certificates exactly
after removing the additional metadata.

The two lost detections
------------------------

Both histories contain a generated 8% slowdown with 100 readings and noise
standard deviation 0.2 around a baseline of 10. The actual generating
correlation is -0.5. The alternative positive correlations below are proposed
explanations evaluated by the detector.

.. list-table:: Surviving explanations after integration
   :header-rows: 1

   * - Actual change location
     - Old surviving rho
     - New surviving rho
     - New determinant probability bound
     - Final decision
   * - After reading 25
     - 0.875
     - 0.8125
     - 0.032029
     - No alert
   * - After reading 50
     - 0.875
     - 0.875
     - 0.082641
     - No alert

The displayed bounds are rounded upward. Both surviving explanations place
the change after reading 1 and pass the size and shape checks. The verifier
recomputes those checks with a separate dense GLS implementation.

For the early-change history, the original search stops after seven cells.
The new search rejects that stopping point and visits an eighth cell, where
correlation 13/16 survives. Numerical integration estimates its directional
probability at 0.0026718, with a quadrature error estimate of about 1.2e-10.
That estimate lies below the 0.01 cutoff, while the certified Chernoff bound
is above it. The numerical estimate guides the next derivation; it does
not authorize rejection by this implementation.

The distinction is now more specific than in the previous checkpoint.
The old bound discarded most covariance information. The new calculation
retains that information, but the Chernoff probability inequality itself
can still be loose. More accurate determinant arithmetic alone will not
close this remaining gap at the fixed tilt.

The five other cases
--------------------

The archived independent and positively correlated 8% examples remain
certified alerts, as do deterministic 20% changes with 24 and 100 readings.
Their interval counts and routes are unchanged. The deterministic 4% example
retains its explanation at the middle split and zero correlation.

No result exhausts a work limit. None of these seven completed runs records
a whole-interval certificate through the new determinant route: the changed
early case finds another survivor before reaching such a certificate, and
the existing routes already cover the four alerts. Tests separately check
a valid determinant interval around the removed explanation, its proof
record, and rejection of damaged records.

Verification and saved evidence
--------------------------------

The `verifier <determinant_verification.py>`_ checks the union of old and new
interval certificates. Every alert must cover the full correlation range
without gaps or overlaps for every split. It reconstructs residual states
by centering the observations independently of the reporting fits and
recomputes each determinant certificate. The existing verifier checks the
inherited routes and surviving explanations.

The detector and verifier share the determinant interval-arithmetic kernel.
The earlier mathematical tests compare that kernel against independent
exact dense determinants. The new 23 integration tests cover full coverage,
unit conversion, exact recovery of the previous search when disabled,
survivors, degenerate data, work limits, calibration validation, and
tampering with proofs, budgets, states, and coverage.

The `summary <data/determinant_reporting_v1_summary.json>`_ records all seven
decisions, witnesses, cell counts, and verification results. The compressed
`replay archive <data/determinant_reporting_v1_diagnosis.json.gz>`_ preserves
the observations, previous and new full results, certificates, source hashes,
and numerical diagnosis. Every result was verified again before archiving,
and the compressed bytes passed a round-trip check. Reproduction commands
are in the `running guide <running.rst>`_.

The next mathematical step
--------------------------

Use the new 13/16 witness to examine a sharper probability inequality.
Two available ingredients are choosing a more effective exponential tilt
and restoring the tilted-density refinement while retaining the full
covariance information. Valid bounds at different tilts can be minimized
because they bound the same p-value; they do not introduce separate tests.
Any proposed refinement still needs a uniform interval proof.

First determine whether those ingredients can certify this new witness.
Then replay the complete search again to discover whether other explanations
remain. Freeze a revised method and evaluate fresh histories only after
that development checkpoint. The middle-change loss remains a separate
calibration question, since its estimated ideal directional probability
already exceeds the existing confidence cutoff.
