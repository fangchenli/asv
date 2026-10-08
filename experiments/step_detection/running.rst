Running the step detection experiments
======================================

The harness compares fitting and penalty-selection methods while retaining
ASV's production scoring rule. Start with the `baseline findings
<baseline_report.rst>`_ for the results of the first run.

Setup and a short run
---------------------

From the repository root, create an environment and build the extension::

    uv venv .venv --python python3
    uv pip install --python .venv/bin/python -e '.[test,dev]' ruff==0.13.3

Check the reporting mathematics
-------------------------------

The `reporting derivation <reporting_uncertainty.rst>`_ has a standalone,
standard-library reference for independent observations with positive
medians and at most one true change. Its checks enumerate sign patterns and
level constraints; they do not run a new statistical experiment::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_reference.py -q

Inspect the known-boundary and unknown-boundary calculations directly::

    .venv/bin/python - <<'PY'
    from experiments.step_detection.reporting_reference import (
        known_boundary_lower,
        single_change_evidence,
    )

    for count in (20, 40):
        values = [10] * count + [12] * count
        print("readings on each side:", count)
        print("known split, lower excess:", known_boundary_lower(values, count))
        print("unknown split:", single_change_evidence(values))
    PY

Both fixed-split lower bounds are 1.5. The unknown-boundary reference returns
``insufficient`` for twenty readings per side and ``evidence`` for forty.
``feasible_splits`` uses retained-observation indices and represents uncertainty
about a single boundary. An empty confidence set returns
``incompatible_model`` and never becomes an alert. Inputs are limited to
200 retained observations; remove missing readings explicitly and do not
treat weighted or correlated readings as independent replicates.

Evaluate the jointly calibrated reporting rule
----------------------------------------------

The `reporting protocol <reporting_protocol.rst>`_ calibrates interval
constraints jointly using fair signs, then freezes them before evaluating
new Gaussian and Laplace histories. Reproduce the held-out evaluation::

    .venv/bin/python -m experiments.step_detection.reporting_study run \
        --frozen experiments/step_detection/data/reporting_v1_frozen.json \
        --output experiments/step_detection/runs/reporting_replay

The runner verifies source, protocol, and native-extension hashes. On a
different build, regenerate the calibration to a new file after verifying
the inherited component-study provenance; do not overwrite the recorded
freeze or call changed-code results an exact replay. The original freeze
command was::

    .venv/bin/python -m experiments.step_detection.reporting_study freeze \
        --inherited experiments/step_detection/data/ablation_v1_frozen.json \
        --output experiments/step_detection/data/reporting_v1_frozen.json

That command refuses to overwrite an existing file. It uses sign sequences
only, without generating evaluation histories. The frozen JSON preserves
all calibration maxima, selected ranks, cutoffs, and failure bounds. The
run saves numerical inputs, fits, existing reports, both confidence-set
results, summaries, and condition breakdowns. Infinite lower bounds are
serialized as the explicit string ``-infinity``.

Run the analytical and pipeline checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_reference.py \
        experiments/step_detection/test_reporting_signs.py \
        experiments/step_detection/test_reporting_study.py -q

Evaluate the direct Gaussian threshold test
-------------------------------------------

The `direct-test protocol <direct_protocol.rst>`_ derives fixed critical
values from Gaussian t and F distributions, without a calibration simulation.
It reuses the previous detector and reporting references on new histories::

    .venv/bin/python -m experiments.step_detection.direct_study run \
        --frozen experiments/step_detection/data/direct_v1_frozen.json \
        --output experiments/step_detection/runs/direct_replay

The original freeze command was::

    .venv/bin/python -m experiments.step_detection.direct_study freeze \
        --inherited experiments/step_detection/data/reporting_v1_frozen.json \
        --output experiments/step_detection/data/direct_v1_frozen.json

Existing freeze files and output directories are never overwritten. The
runner checks the inherited hashes, new source/protocol hashes, native
extension, and SciPy version. Reproduction requires the matching environment;
a changed build needs separately recorded provenance and calibration.

Each result includes the existing fit and all previous evidence decisions.
In ``direct``, array entry i refers to retained-observation split i+1:
``size_t`` tests the relative increase, ``lack_of_fit_f`` tests the largest
improvement from an extra boundary, and ``best_extra_split`` identifies that
boundary. ``null_splits`` contains locations not rejected by either route.
It is not a confidence interval for the change location. Infinite statistics
in deterministic cases are serialized as explicit strings.

Run the new analytical and pipeline checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_direct.py \
        experiments/step_detection/test_direct_study.py -q

Stress variance and independence
--------------------------------

The `stress protocol <direct_stress_protocol.rst>`_ derives covariance
diagnostics and reuses the entire frozen direct-test pipeline. Reproduce
the six paired noise conditions without changing any cutoff::

    .venv/bin/python -m experiments.step_detection.direct_stress run \
        --frozen experiments/step_detection/data/direct_stress_v1_frozen.json \
        --output experiments/step_detection/runs/direct_stress_replay

The original freeze command, committed before generating histories, was::

    .venv/bin/python -m experiments.step_detection.direct_stress freeze \
        --inherited experiments/step_detection/data/direct_v1_frozen.json \
        --output experiments/step_detection/data/direct_stress_v1_frozen.json

The same source/build/SciPy checks apply. Existing freeze files and output
directories are refused. A new build needs separately recorded provenance.
The run saves numerical inputs, complete pipeline outputs, condition
breakdowns, paired gains/losses against the independent constant-variance
control, and rejection routes at the generating split. That split is used
only by the evaluator. ``neither`` must never occur on a direct alert.

Run the mathematical and pipeline checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_direct_stress.py \
        experiments/step_detection/test_reporting_direct.py \
        experiments/step_detection/test_direct_study.py -q

The `stress report <direct_stress_report.rst>`_ links the committed result
index and compressed archives. The index records filenames, counts, and
SHA-256 hashes; decompress each ``*.jsonl.gz`` file with Python's ``gzip``
module. Join inputs by ``id`` to results by ``case.id``. ``pair_id`` connects
the six versions of each base history, so aggregate versions are not
independent statistical trials.

Check the known-covariance reference
--------------------------------------

The `derivation <known_covariance.rst>`_ extends both direct tests to a
supplied Gaussian covariance shape. The checks compare deterministic GLS
calculations and projection identities; they do not generate an evaluation
sample or estimate covariance::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_covariance.py \
        experiments/step_detection/test_reporting_direct.py -q

Use the reference from Python at the repository root::

    import numpy as np
    from experiments.step_detection.reporting_covariance import evidence
    from experiments.step_detection.reporting_direct import critical_values

    n = 16
    positions = np.arange(n)
    values = 10 + 0.001*np.sin(2.1*positions) + 0.8*(positions >= 8)
    amplitudes = np.where(positions >= 8, 3.0, 1.0)
    covariance = (np.outer(amplitudes, amplitudes)
                  * 0.7**np.abs(positions[:, None] - positions[None, :]))
    result = evidence(values, covariance, critical_values(n))
    print(result['status'], result['has_alert'])

This deterministic example returns ``ok True``. The covariance matrix is a
declared input, used unchanged for every candidate mean split. The example
demonstrates the API, not an estimate of detection power. A positive scalar
multiple of the covariance gives the same statistics. The overall noise
scale is estimated from weighted residuals separately under each model.

``insufficient_precision`` means the helper abstained, with ``has_alert=False``
and diagnostic arrays set to None. It is distinct from an ordinary non-alert
with surviving ``null_splits``. Malformed, non-positive-definite, or excessively
ill-conditioned covariance raises ValueError. Histories are limited to 200
readings; this is a mathematical reference with no production integration.

Evaluate the covariance oracle
------------------------------

The `oracle protocol <covariance_protocol.rst>`_ compares the known-covariance
test with the entire preceding pipeline on fresh paired Gaussian histories.
Replay the frozen evaluation with the matching Python, NumPy, SciPy, source,
and native build::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.covariance_study run \
        --frozen experiments/step_detection/data/covariance_v1_frozen.json \
        --output experiments/step_detection/runs/covariance_replay --workers 4

The original freeze command was::

    .venv/bin/python -m experiments.step_detection.covariance_study freeze \
        --inherited experiments/step_detection/data/direct_stress_v1_frozen.json \
        --output experiments/step_detection/data/covariance_v1_frozen.json

The freeze hashes the 24 covariance matrices without generating observations.
Existing freeze files and run directories are refused. A changed environment
requires separately recorded provenance; do not overwrite the recorded freeze.
Workers evaluate independent histories and write results in input order.

Each input's ``covariance_id`` identifies its matrix in
``covariance_shapes.json``. Matrices are saved once, without the unknown
overall scale. The manifest and freeze contain their canonical JSON hashes.
Records contain the old pipeline plus ``oracle`` statistics, numerical status,
and ``oracle_generating_split_route``. ``oracle_alone`` uses the new evidence;
``oracle_gate`` also requires the existing shared-floor report.

``comparisons.json`` counts new detections and false alerts gained/lost versus
the old direct test on the same histories. ``paired.json`` compares each
condition to its paired independent constant-variance control, while
``diagnostics.json`` records numerical abstentions and false-alert rejection
routes. The runner enforces agreement with the old test in the identity
control. No covariance is estimated from the observations in this study.
The `oracle report <covariance_report.rst>`_ links the committed result index
and numerical archives, including all 24 covariance matrices.

Run the study and analytical checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_covariance_study.py \
        experiments/step_detection/test_reporting_covariance.py -q

Check the unknown-correlation reference
---------------------------------------

The `unknown-correlation derivation <unknown_correlation.rst>`_ covers a
stationary Gaussian AR(1) model with constant marginal variance and one
possible mean change. Run its deterministic mathematical checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_ar1.py \
        experiments/step_detection/test_reporting_covariance.py -q

Inspect a fully certified example without generating simulation data::

    from experiments.step_detection.reporting_ar1 import evidence

    values = [10 + 0.1*(-1)**i + 0.8*(i >= 6) for i in range(12)]
    result = evidence(values, max_cells=256, max_depth=8)
    print(result['status'], result['has_alert'])
    print(result['certificate'])

This returns ``certified_alert True``. Every certificate entry identifies a
split, exact rational interval endpoints, a rejection route, and an extra
boundary when the shape test provides the certificate. Intervals together
cover -1<rho<1 at every possible mean split. A point grid is never used to
justify rejection of an entire interval.

``confidence_polynomials`` stores the quadratic residual sum for every split
as exact fraction strings. Compare its value at a candidate rho with
``confidence.residual_bound`` to determine confidence-set membership.
Both use values divided by ``confidence.normalization_scale``, chosen from
the training prefix alone. The fitted predictor in ``confidence`` is used
to build a proper conditional density; it is not a plug-in covariance
estimate for the reporting tests.

``surviving_explanation`` saves a pair that prevents an alert. ``unresolved``
means the declared cell/depth budget was insufficient, and
``insufficient_variation`` identifies an exactly fitted candidate.
All three return ``has_alert=False``. Use ``calibration(n)`` to inspect or
explicitly choose the error allocation; edited cutoff dictionaries are
rejected. The reference accepts 8 to 200 finite observations. The
`frozen study report <ar1_report.rst>`_ records its sensitivity and elapsed time.

Evaluate unknown correlation
----------------------------

The `AR(1) study protocol <ar1_protocol.rst>`_ fixes the predictor, critical
values, and certification budget before fresh observations are generated.
The `report <ar1_report.rst>`_ explains the result; the
`result index <data/ar1_v1_results.json>`_ lists summaries and archive hashes.
Replay it with the matching frozen environment::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.ar1_study run \
        --frozen experiments/step_detection/data/ar1_v1_frozen.json \
        --output experiments/step_detection/runs/ar1_replay --workers 6

The original freeze command was::

    .venv/bin/python -m experiments.step_detection.ar1_study freeze \
        --inherited experiments/step_detection/data/covariance_v1_frozen.json \
        --output experiments/step_detection/data/ar1_v1_frozen.json

Existing freezes and output directories are refused. The runner checks the
entire inherited source/build/environment chain, the new code and protocol,
analytical calibration, work budgets, and covariance hashes. Inputs reference
one of six matrices in ``covariance_shapes.json``; only the oracle uses them.

Records retain the entire inherited pipeline plus ``oracle_reporting`` at
alpha=0.04 and ``ar1`` with its 0.01 confidence and 0.04 reporting allocation.
Both are saved standalone and gated. ``ar1_diagnostics`` contains the
confidence intervals by split, their union, widths, membership of the
generating pair, and its rejection route. Decimal roots describe confidence
regions; they do not replace the exact certificate-based reporting decision.

``diagnostics.json`` summarizes statuses, confidence exclusions, near-one
correlation admissibility, cell counts, and elapsed reference time.
``comparisons.json`` records paired gains/losses versus the old direct test
and both oracle budgets. Unresolved searches count as non-alerts and remain
visible in the status breakdown. Histories sharing a pair identifier across
correlation conditions are dependent.

Run the evaluator and mathematical checks::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_ar1_study.py \
        experiments/step_detection/test_reporting_ar1.py -q

Reproduce the information-loss diagnosis
----------------------------------------

The `mathematical follow-up <ar1_information.rst>`_ reuses one archived missed
slowdown. It generates no new histories and makes no new reporting decisions::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.ar1_information

The command prints JSON with forward and reverse likelihood calculations,
exact versus relaxed conditional fits, averaged evidence, and the full
stationary likelihood near correlation one. It checks the input archive hash
against the frozen result index. The saved output is
``data/ar1_information_v1.json``. Check the algebra independently with::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_ar1_information.py \
        experiments/step_detection/test_reporting_ar1.py -q

Check the full-history confidence reference
-------------------------------------------

The `full-history derivation <full_history_confidence.rst>`_ specifies the
proper predictive mixture and continuous interval checks. The API requires
an external reference unit, expressed in the same units as the observations::

    from experiments.step_detection.reporting_ar1_full import evidence

    result = evidence(values, unit=1)
    print(result['status'], result['witness'])

Reproduce one archived missed slowdown and two deterministic fixtures::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.full_history_diagnosis

The saved output is ``data/full_history_v1_diagnosis.json``. This is a
development diagnosis, not a fresh statistical evaluation. It preserves
the full certificates and witnesses. Verify the density, likelihoods,
and continuous bounds with::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_ar1_full.py \
        experiments/step_detection/test_ar1_information.py \
        experiments/step_detection/test_reporting_ar1.py -q

Check the candidate-specific baseline/jump reference
----------------------------------------------------

The `baseline/jump derivation <baseline_jump_prior.rst>`_ proves coverage for
a separate proper predictor at each proposed change location::

    from experiments.step_detection.reporting_ar1_jump import evidence

    result = evidence(values, unit=1)
    print(result['status'], result['witness'])

``confidence.by_split`` records each visited location's density and confidence
coefficient. ``original_location_evidence`` retains the original level prior
as a control. ``global_evidence`` uses the new jump prior but the old global
location mixture. Reproduce their comparison on the archived history::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.jump_prior_diagnosis

The saved output is ``data/jump_prior_v1_diagnosis.json``. It includes the
exact prediction-cost accounting, all four reference results, independent
GLS witness checks, and three deterministic fixtures. This is development
evidence, not a new statistical study. Run the analytical checks with::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_ar1_jump.py \
        experiments/step_detection/test_reporting_ar1_full.py \
        experiments/step_detection/test_ar1_information.py \
        experiments/step_detection/test_reporting_ar1.py -q

Evaluate candidate-specific prediction
--------------------------------------

The `indexed-study protocol <indexed_protocol.rst>`_ fixes a fresh 360-history
comparison of the conditional reference, the global full-history reference,
both indexed variants, and the matching oracle. The `results
<indexed_report.rst>`_ and `archive index <data/indexed_v1_results.json>`_
preserve the completed comparison. Replay the frozen study::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.indexed_study run \
        --frozen experiments/step_detection/data/indexed_v1_frozen.json \
        --output experiments/step_detection/runs/indexed_replay --workers 6

The original freeze command was::

    .venv/bin/python -m experiments.step_detection.indexed_study freeze \
        --inherited experiments/step_detection/data/ar1_v1_frozen.json \
        --output experiments/step_detection/data/indexed_v1_frozen.json

The runner rejects changed source, settings, environment, or covariance
hashes, and refuses to overwrite an existing freeze or run directory.
The manifest records worker and numerical-library thread settings. Inputs
share six covariance matrices, supplied only to the oracle.

Records retain the inherited pipeline and add ``full_history``,
``indexed_original``, and ``indexed_jump`` with standalone and gated decisions.
The new ``*_diagnostics`` fields check the true pair even when an indexed
search stops before its generating location. ``diagnostics.json`` separates
certificates, surviving witnesses, unresolved searches, insufficient variation,
true-pair exclusions, and confidence fallbacks. ``comparisons.json`` records
paired gains and losses. Unresolved positive histories remain misses.

Run the evaluator checks without generating evaluation-seed observations::

    .venv/bin/python -m pytest experiments/step_detection/test_indexed_study.py -q

After a complete run, reconstruct its certificates and witnesses, regenerate
its observations, reproduce its summaries, and write hashed compressed
archives to a new destination::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.indexed_artifacts \
        --run experiments/step_detection/runs/indexed_replay \
        --frozen experiments/step_detection/data/indexed_v1_frozen.json \
        --destination experiments/step_detection/runs/indexed_replay_archives \
        --workers 6

The archiver refuses to overwrite archive files. Its index contains the
summaries, paired comparisons, diagnostic counts, and SHA-256 hashes.

Check the residual-direction derivation
-----------------------------------------

The `mathematical analysis <residual_direction.rst>`_ removes the two plateau
levels and common noise scale before constructing correlation confidence.
Its helper checks exact density identities and two existing archived
histories; it does not run a new detector or generate evaluation data::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.residual_direction \
        --output experiments/step_detection/runs/residual_direction_replay.json

The output path must not already exist. The checked-in `diagnosis
<data/residual_direction_v1_diagnosis.json>`_ preserves archive/source hashes,
the fixed five-component predictor, density bounds, and endpoint checks.
Run the analytical tests with::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest experiments/step_detection/test_residual_direction.py -q

Check the directional tail certificate
----------------------------------------

The `tail derivation <directional_tail.rst>`_ proves an upper probability
bound across a correlation interval. The helper uses rational arithmetic for
the certificate and separately records floating tail estimates for comparison::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.directional_tail \
        --output experiments/step_detection/runs/directional_tail_replay.json

The output path must be new. The `saved diagnosis
<data/directional_tail_v1_diagnosis.json>`_ includes both archived inputs'
hashes, exact certificate quantities, and the numerical estimates. No full
reporting search or new detection study is performed. Run the checks with::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest experiments/step_detection/test_directional_tail.py \
        experiments/step_detection/test_residual_direction.py -q

Check the complete directional reporting search
------------------------------------------------

The `complete reference <directional_reporting.rst>`_ combines uniform
directional and tilted tail confidence bounds with size and shape checks.
Reproduce the ten verified development runs on two saved histories and
three deterministic fixtures::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.directional_diagnosis \
        --output experiments/step_detection/runs/directional_reporting_replay.json \
        --workers 3

The output path must be new. Every result is verified before saving. The
`compressed diagnosis <data/directional_reporting_v1_diagnosis.json.gz>`_
contains the full certificates and source hashes; open it with Python's
``gzip.open`` and ``json.load``. These runs use development examples and
do not estimate fresh detection or false-alert rates.

Run the mathematical helpers and complete-search checks::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_directional.py \
        experiments/step_detection/test_directional_tail.py \
        experiments/step_detection/test_residual_direction.py -q

Evaluate directional confidence on fresh histories
--------------------------------------------------

The `directional protocol <directional_protocol.rst>`_ freezes the two
directional variants and all inherited comparisons before evaluation.
The `result report <directional_report.rst>`_ describes its 58/144 detections
and zero false alerts among 216 null histories. Reproduce the 360 histories
with::

    PYTHONINTMAXSTRDIGITS=0 \
        VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.directional_study run \
        --frozen experiments/step_detection/data/directional_v1_frozen.json \
        --output experiments/step_detection/runs/directional_replay --workers 6

The original freeze command was::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.directional_study freeze \
        --inherited experiments/step_detection/data/indexed_v1_frozen.json \
        --output experiments/step_detection/data/directional_v1_frozen.json

The freeze is already saved. Both commands refuse to overwrite their output.
The runner rejects changed sources, settings, environment, or covariance
hashes. Records retain every inherited result and add ``directional_uniform``
and ``directional_tail``, including standalone and gated decisions, complete
certificates, witnesses, and generating-pair diagnostics. Test the evaluator
using its separate test stream with::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest experiments/step_detection/test_directional_study.py -q

After completion, verify and archive the complete run to a new destination::

    PYTHONINTMAXSTRDIGITS=0 \
        VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.directional_artifacts \
        --run experiments/step_detection/runs/directional_replay \
        --frozen experiments/step_detection/data/directional_v1_frozen.json \
        --destination experiments/step_detection/runs/directional_replay_archives \
        --workers 6 \
        --cache experiments/step_detection/runs/directional_replay/verification_cache.json

Adding ``--verify-available`` verifies only complete saved records while a
run is in progress; it does not create archives. The optional cache is reused
only when source hashes and the complete input/result pair match. Final
archiving also regenerates all inputs and summary tables and checks every
compressed file against its uncompressed content.

``PYTHONINTMAXSTRDIGITS=0`` permits serialization and parsing of the large
exact fractions produced by these trusted generated inputs. It does not
change numerical values or detector decisions. The initial run reached
Python's default 4300-digit conversion limit after saving 91 histories.
The recovery wrapper preserved those results and resumed the same frozen
evaluation, saving each subsequent completed history immediately::

    PYTHONINTMAXSTRDIGITS=0 \
        VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.directional_recovery \
        --run experiments/step_detection/runs/directional_v1 \
        --frozen experiments/step_detection/data/directional_v1_frozen.json \
        --workers 6

Use this only for an interrupted run. It verifies the saved prefix, frozen
sources, original manifest, worker count, and numerical-library thread
settings before resuming. ``recovery.json`` records the preserved prefix
hashes and wrapper hash; the original manifest is retained. A fresh replay
can use the regular runner with the conversion limit set as shown above.

Reproduce the numerical diagnosis of detections lost to indexed-jump from
the saved, hash-checked archives::

    PYTHONINTMAXSTRDIGITS=0 \
        VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.directional_losses \
        --index experiments/step_detection/data/directional_v1_results.json \
        --output experiments/step_detection/runs/directional_losses_replay.json

This is post-evaluation development analysis. Its floating tail estimates
do not change any frozen decision and are not interval certificates. The
saved output is `the loss diagnosis <data/directional_v1_loss_diagnosis.json>`_.

Full-determinant mathematical checkpoint
-----------------------------------------

The `full-determinant derivation <directional_determinant.rst>`_ adds a
standalone interval certificate and investigates the two saved lost
detections. It uses existing archived inputs and does not run a fresh study::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.determinant_diagnosis \
        --output experiments/step_detection/runs/determinant_replay.json

The command checks input and record hashes and refuses to overwrite output.
The reference output is `the determinant diagnosis
<data/determinant_v1_diagnosis.json>`_. Certificates use exact outward
rational arithmetic. The accompanying GLS/product p-values are numerical
diagnostics and use the previous saved quadrature estimates.

Run the algebra, interval, and archived-witness checks::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest \
        experiments/step_detection/test_directional_determinant.py -q

Complete determinant reporting replay
---------------------------------------

The `integration report <determinant_reporting.rst>`_ records a replay of
seven saved cases. This runs the complete new search and verifies every
result, using the previous saved results as the comparison::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.determinant_replay \
        --output experiments/step_detection/runs/determinant_reporting_replay \
        --workers 2

The output directory must be new. Each completed, verified result is saved
immediately; ``diagnosis.json`` assembles them in input order. Sources and
input archives are hashed. The integer-string setting permits serialization
of long exact rational values from these trusted calculations.

Validate again and create a compressed archive plus a readable JSON summary::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.determinant_replay_artifacts \
        --run experiments/step_detection/runs/determinant_reporting_replay \
        --output experiments/step_detection/runs/determinant_replay.json.gz \
        --summary experiments/step_detection/runs/determinant_replay_summary.json

Both output files must be new. The builder verifies the complete certificates
and checks the source hashes before archiving. It also estimates the new
surviving explanation's directional probability by numerical integration;
that estimate does not determine any reporting decision.

Run the integration and certificate-tampering tests::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_determinant.py -q

Spectral refinement and complete replay
----------------------------------------

The `spectral derivation and results <spectral_refinement.rst>`_ extend the
density bound to several direction counts and two fixed tilts. Mathematical
and integration checks are separate from any fresh simulation study::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest \
        experiments/step_detection/test_directional_spectral.py \
        experiments/step_detection/test_reporting_spectral.py -q

Replay the seven saved cases with the determinant replay as the comparison::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.spectral_replay \
        --output experiments/step_detection/runs/spectral_replay --workers 2

The output directory must be new. Each verified result is saved immediately.
The output also records point and interval certificates at previous surviving
explanations and numerical tail estimates at current survivors.

Verify again and archive the complete replay::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.spectral_replay_artifacts \
        --run experiments/step_detection/runs/spectral_replay \
        --output experiments/step_detection/runs/spectral_replay.json.gz \
        --summary experiments/step_detection/runs/spectral_replay_summary.json

Both outputs must be new. Archived reference results are linked from the
derivation. Only rigorous bounds determine reporting decisions; numerical
quadrature remains a development diagnostic.

Sturm-count mathematical checkpoint
------------------------------------

The `split-aware derivation <sturm_refinement.rst>`_ tightens only the
eigenvalue input to the existing probability formula. Reproduce its saved
witness calculations without running a fresh study or a complete search::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.sturm_diagnosis \
        --output experiments/step_detection/runs/sturm_diagnosis_replay.json

The output must be new. The command checks the previous replay archive and
source hashes and records exact point/interval certificates together with
numerical comparisons that separate approximation losses. The reference is
`the saved Sturm diagnosis <data/sturm_v1_diagnosis.json>`_.

Run the exact-count and interval checks::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest \
        experiments/step_detection/test_directional_sturm.py -q

Complete Sturm reporting replay
--------------------------------

The `Sturm integration report <sturm_reporting.rst>`_ records the complete
search on seven saved cases. Reproduce it using a new output directory::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.sturm_replay \
        --output experiments/step_detection/runs/sturm_replay --workers 2

Each verified case is saved as it finishes. ``diagnosis.json`` combines the
results in input order; the previous spectral replay supplies the comparison.

Re-verify and archive the saved results, including a numerical decomposition
at a changed surviving explanation::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.sturm_replay_artifacts \
        --run experiments/step_detection/runs/sturm_replay \
        --output experiments/step_detection/runs/sturm_replay.json.gz \
        --summary experiments/step_detection/runs/sturm_replay_summary.json

Both output files must be new. The builder checks source hashes and the
original inputs/results, reconstructs every reporting and saved-witness
certificate, and verifies the compressed bytes by a round trip. Numerical
decompositions never determine reporting decisions.

The latest grouped-Fourier result is saved as
`the v2 summary <data/sturm_reporting_v2_heterogeneous_summary.json>`_ and
`the v2 archive <data/sturm_reporting_v2_heterogeneous_diagnosis.json.gz>`_.
To reproduce those names, use a fresh run directory and pass them to the
artifact builder, for example ``runs/sturm_replay_v2`` for the run directory.

Run the integration and proof-integrity checks::

    PYTHONINTMAXSTRDIGITS=0 VECLIB_MAXIMUM_THREADS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m pytest \
        experiments/step_detection/test_reporting_sturm.py -q

Check the Lean pilot
---------------------

The `Lean pilot guide <lean/README.rst>`_ explains the proof boundary and
dependency setup. After fetching its pinned Mathlib cache, run from the
repository root::

    .venv/bin/python experiments/step_detection/lean/verify.py

Pass ``--lake "$HOME/.elan/bin/lake"`` if Lake is not on PATH. This checks
the proof build, archived inputs, full determinant trace and matrix identity,
permitted axioms, and rejection of corrupted probability bounds, intermediate
pivot bounds, and covariance entries.
Use the repository Python environment because the trace exporter imports
the frozen calculation's NumPy/SciPy dependencies. It runs no new statistical
experiments.

Run the original harness
-------------------------

Run a small comparison::

    .venv/bin/python experiments/step_detection/harness.py \
        --output experiments/step_detection/runs/quick \
        --scenarios flat step weak --sizes 100 --seeds 0 1

Use a new output directory for each run. The harness refuses to overwrite
existing results. Generated run directories are ignored by Git.

For a Python-only run, use ``--backends python``. That backend loads the local
source without importing the ASV package and requires only the Python standard
library. The native backend loads the locally built C++ extension. Both use
isolated detector modules, so search instrumentation does not modify an ASV
module imported elsewhere in the process.

Reproduce the initial experiments
---------------------------------

The recorded baseline uses the old correlation search at revision ``d526c2e``.
To reproduce it, check out that revision in a separate worktree and follow
the setup above there. Running these commands on the current branch instead
measures the corrected weighted-median correlation fit.

The baseline uses 15 scenario families, 100 observations, four seeds, both
backends, and all five methods::

    .venv/bin/python experiments/step_detection/harness.py \
        --output experiments/step_detection/runs/baseline

The scaling sample uses flat and single-step histories at 1000 and 10000
observations::

    .venv/bin/python experiments/step_detection/harness.py \
        --output experiments/step_detection/runs/scaling \
        --scenarios flat step --sizes 1000 10000 --seeds 0 \
        --methods approximate current --timeout 120

``--max-exact-size`` defaults to 200; the unrestricted exact method is skipped
above that size. ``--timeout`` is a per-worker limit in seconds. Failed workers
or timeouts stop the run; already completed records remain in ``records.jsonl``.
A final report is written only after all requested measurements succeed.

Compare the original and corrected detector
-------------------------------------------

The version comparison replays identical inputs through revisions ``d526c2e``
and ``123405a``. Both use the same locally built extension. Run from the
repository root with Git history available::

    .venv/bin/python -m experiments.step_detection.compare_versions \
        --output experiments/step_detection/runs/version_comparison

Defaults are 15 synthetic families with seeds 0 through 11 at length 100,
six saved public NumPy histories, and flat/step histories at lengths 1000
and 10000. Each version/backend/case gets three fresh processes. Version
order alternates; timings include input preparation, step detection, graph
coordinate mapping, and the 5 percent regression-reporting rule. Imports
and diagnostic tracing are excluded. Accuracy counts each history once.

``--seeds``, ``--backends``, and ``--repeats`` narrow the run;
``--no-real`` and ``--no-scaling`` omit those corpora. ``--before`` and
``--after`` select other Git revisions. A run saves source snapshots,
input arrays, source and extension hashes, every timing sample, paired
summaries, and untimed candidate traces for changed fits.

``data/numpy_snapshot.json`` is the frozen real-data input, so ordinary
runs require no downloads. ``fetch_numpy.py`` records the six named
operations, environment, source URLs, hashes, and commit/date mappings.
The published graphs contain no uncertainty weights; these replays use
unit weights and have no labeled change-point truth. Refreshing the
snapshot changes the experiment and requires explicitly moving the old
file before running::

    .venv/bin/python -m experiments.step_detection.fetch_numpy

Validate the comparison with::

    .venv/bin/python -m pytest experiments/step_detection/test_compare_versions.py -q

Inspect the five weak-change boundary errors and representative dip/outlier
cases against the exact reference, scoring its candidates with the actual
production closure from each revision::

    .venv/bin/python -m experiments.step_detection.inspect_comparison \
        experiments/step_detection/runs/version_comparison

``--cases`` selects other saved cases of at most 200 retained observations;
``--backend`` chooses the interval backend. This is an untimed diagnostic.
It saves ``inspection-native.json`` by default. Only one independent-error
minimizer is retained at each segment count, so rescoring that frontier does
not establish a global optimum for the correlated score.

See `the version comparison <version_comparison.rst>`_ for measured results
and their interpretation.

Run the near-threshold study
-----------------------------

The `protocol <threshold_protocol.rst>`_ specifies the simulation and setting
selection before evaluation. The `results <threshold_report.rst>`_ report
the observed sensitivity/false-alert tradeoff. Run development first::

    .venv/bin/python -m experiments.step_detection.threshold_study development \
        --output experiments/step_detection/runs/threshold_development

This evaluates 360 histories and writes ``frozen.json`` with one independent
and one bounded-persistence configuration. Evaluate the 600 fresh histories
using that file::

    .venv/bin/python -m experiments.step_detection.threshold_study heldout \
        --frozen experiments/step_detection/runs/threshold_development/frozen.json \
        --output experiments/step_detection/runs/threshold_heldout

Every output directory must be new. Each phase saves all inputs, per-history
results, aggregate results, and breakdowns. Held-out execution checks the
frozen protocol, source, design, and native-extension hashes before generating
cases. In another environment, run both phases to create a freeze matching
that build. The archived original configuration deliberately pins its original
extension binary. Use ``--backend python`` in both phases to avoid the extension.

Validate the runner and its exact reference::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_threshold_study.py \
        experiments/step_detection/test_exact_reference.py -q

Compressed input and result snapshots live in ``data/threshold_v1_*.jsonl.gz``.
For example, inspect the saved held-out results without rerunning fitting::

    import gzip
    import json

    with gzip.open(
        'experiments/step_detection/data/threshold_v1_heldout_records.jsonl.gz', 'rt'
    ) as stream:
        rows = [json.loads(line) for line in stream]

Run the component and calibration study
---------------------------------------

The `component protocol <ablation_protocol.rst>`_ separates the segment
penalty, correlation cap, and floor formulation. It calibrates each scoring
family to the same development false-alert budgets, then measures the actual
rates on fresh histories. Run both phases in the same environment::

    .venv/bin/python -m experiments.step_detection.ablation_study development \
        --output experiments/step_detection/runs/ablation_development

    .venv/bin/python -m experiments.step_detection.ablation_study heldout \
        --frozen experiments/step_detection/runs/ablation_development/frozen.json \
        --output experiments/step_detection/runs/ablation_heldout

Development uses 1200 histories and evaluation uses 2400, with random seeds
disjoint from the earlier study. The held-out run evaluates the eight fixed
factorial settings and the frozen calibrated choices. It rejects changed
protocol, code, extension, or calibration-grid hashes. ``--backend python``
works when supplied to both phases. Each output directory must be new.

Validate the scoring controls and calibration with::

    .venv/bin/python -m pytest experiments/step_detection/test_ablation_study.py -q

See `the component results <ablation_report.rst>`_ for interpretation. Archived
inputs and results use ``data/ablation_v1_*.jsonl.gz``, readable in the same
way as the preceding threshold-study snapshots.

Run the frozen robustness study
--------------------------------

The `robustness protocol <robustness_protocol.rst>`_ reuses the component
study's full-correlation-cap settings at all three calibration budgets.
There is no new configuration search. Freeze settings from a completed
component-development run in the same environment::

    .venv/bin/python -m experiments.step_detection.robustness_study freeze \
        --calibration experiments/step_detection/runs/ablation_development/frozen.json \
        --output experiments/step_detection/runs/robustness_frozen.json

Then run the paired control and stress cases::

    .venv/bin/python -m experiments.step_detection.robustness_study run \
        --frozen experiments/step_detection/runs/robustness_frozen.json \
        --output experiments/step_detection/runs/robustness

The original study used ``data/robustness_v1_frozen.json``; it pins its original
sources and extension build. As with previous studies, a different build
requires a matching calibration/freeze rather than changing the archived
hashes. Supply ``--backend python`` consistently when using that backend.

The study saves all 7680 inputs and per-history results, condition breakdowns,
and paired gains/losses against the Gaussian control. Missing readings are
filtered through production, and fits are mapped back to original revision
coordinates before reporting. Validate these mechanisms with::

    .venv/bin/python -m pytest experiments/step_detection/test_robustness_study.py -q

See `the robustness results <robustness_report.rst>`_ for findings. Archived
``data/robustness_v1_*.jsonl.gz`` files preserve the input and result snapshots.

What the methods compare
------------------------

* ``exact`` and ``approximate`` use the same fixed gamma,
  ``3 * single_segment_error * log(n) / n``. This is the approximate solver's
  default rule, evaluated on the prepared observations.
* ``current`` runs the actual ``solve_potts_autogamma`` implementation, recording
  the penalties, fits, and scores it evaluates.
* ``grid`` replaces only the outer search, using exactly as many candidate
  fits as ``current`` used for that input. It covers the expanded bracket from
  the production golden-section search, then refines between sampled penalties
  that yield different partitions, favoring lower-scoring regions.
* ``hybrid`` runs current search first, then adds a grid/refinement budget equal
  to its candidate count. It keeps the best score encountered across both stages.

Grid and hybrid retain the production objective closure, including its noise
floor, complexity penalty, and residual correlation calculation. Repeated fits
reuse the fitted correlation parameter. The final residual sum still comes
from production code. The tests verify each recorded grid score against a
fresh production evaluation forced to use the same gamma.

Grid and hybrid require ``current`` in the requested methods to establish the
per-input budget. In the initial baseline, current used ten candidate fits,
so grid used ten and hybrid used twenty. Their runtime differences include
both candidate selection and reuse of correlation fits.

Correctness and detection metrics
---------------------------------

Each backend first passes 200 small-input checks against exhaustive partition
enumeration. The oracle evaluates possible segment levels directly from the
observed values, independently of ASV's median and dynamic-programming code.
It checks weights, ties, interval constraints, and subranges. The harness stops
before benchmarking if these checks disagree.

Boundary matching maximizes the number of one-to-one matches within a fixed
tolerance, then minimizes their total location error. The default tolerance is
two input positions, including positions with missing values. For missing-data
cases, truth boundaries are projected to the next retained observation.
Regression alerts are scored after mapping segments into revision coordinates
and applying ASV's existing 5 percent reporting rule.

Synthetic truth includes both physical step positions and explicit expected
regression positions. For example, a brief dip has two true boundaries but
no expected regression alert. Drift has no discrete boundary truth and is
excluded from detection accuracy counts. Even seeds form the development
split and odd seeds form the held-out split. No parameters were tuned during
the initial comparison.

Timing and artifacts
--------------------

Each method/input/backend combination runs in a fresh subprocess. Timings
cover the fit and search, including instrumentation; they exclude interpreter
startup, input preparation, and the separate fixed-gamma scale calculation.
Peak resident memory covers the entire worker, including interpreter overhead
and native allocations. It is unavailable on platforms without ``resource``.

Phase times cover dynamic programming, merging, candidate fitting, and the
correlation fit. Candidate fitting contains the dynamic-programming and
merge phases, so those times overlap. Internal C++ interval queries and cache
hits are not instrumented in this first version.

The corrected helper is recorded as ``rho_fits`` and ``rho_fit`` time.
Historical artifacts instead have ``rho_objective_evaluations`` and
``rho_search`` time from the numerical search.

Each run writes:

* ``manifest.json``: command, settings, interpreter, platform, Git state, and
  SHA256 hashes of detector and harness sources.
* ``inputs.json``: complete generated values, weights, revisions, seeds, and
  expected changes and alerts.
* ``oracle.json``: correctness-check counts and any disagreements.
* ``records.jsonl``: one result per measured method, including every sampled
  gamma, candidate fit, score, timing, memory, and detection metrics.
* ``summary.json``: aggregate metrics and paired method comparisons, including
  descriptive Wilson intervals for the rate of histories with false detections.
* ``report.rst``: a readable summary of the results.

Keep the full artifacts when investigating a changed boundary. They distinguish
an optimizer finding a lower score from a method finding a more accurate change.

Checks
------

The `exact independent reference <exact_reference.rst>`_ has a Python API
and exhaustive correctness checks::

    .venv/bin/python -m pytest experiments/step_detection/test_exact_reference.py -q

It returns the best fit at each feasible segment count and the global
shared-floor optimum. Inputs are already prepared values and weights; the
floor and beta are explicit. Histories are limited to 200 observations.

The shared-floor prototype has analytical tests that do not run the
benchmark harness::

    .venv/bin/python -m pytest experiments/step_detection/test_noise_model.py -q

See `the shared-floor design <shared_noise_floor.rst>`_ for its required
noise-scale and complexity-penalty inputs. Production uses its existing score.

The same test file covers the `bounded-persistence scorer <noise_persistence.rst>`_.
For a fixed candidate, it can be used from Python at the repository root::

    from experiments.step_detection.noise_model import (
        correlation_cap, fit_correlated_score,
    )

    result = fit_correlated_score(
        residuals=[-1, -1, -1, -1, 1, 1, 1, 1],
        weights=[1] * 8,
        segments=1,
        noise_floor=0.1,
        beta=0.2,
        rho_max=correlation_cap(4),
    )

These values illustrate the API; the half-life, floor, and beta are caller
inputs rather than calibrated defaults. The scorer imports ASV's production
correlation helper, so use the environment from the setup section.

Run the experiment tests and the production step-detection tests::

    .venv/bin/python -m pytest \
        experiments/step_detection/test_harness.py test/test_step_detect.py -q

The experiment tests live outside the default production test directory and
must be named explicitly. The native cases are skipped if the extension is
not built; a benchmark run requesting native measurements fails in that case.
