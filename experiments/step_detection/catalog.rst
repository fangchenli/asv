Step-detection experiment map
=============================

This catalog groups the documentation and Python files in this folder by
what they are for. Start at the `current status <current_status.rst>`_ for
the latest result, or at the `plain-language introduction <README.rst>`_ if
you want the broader explanation of ASV's step detector.

Folder layout
-------------

The Python modules stay together at this level because they are imported as
``experiments.step_detection.<module>`` and many experiment commands use
``python -m experiments.step_detection.<module>``. Moving them into separate
folders would change those import paths. This catalog groups them by role
without breaking those commands.

* ``data/`` contains frozen inputs, summaries, and verified archives.
* ``runs/`` contains generated intermediate outputs. These are reproducible
  run products rather than source files.
* ``lean/`` contains the isolated Lean proof pilot.
* ``test_*.py`` files validate the corresponding references, studies, and
  saved certificates.

Start here
----------

* `current_status.rst <current_status.rst>`_: latest result, production
  status, limits, and next target.
* `README.rst <README.rst>`_: plain-language background and chronological
  history of the research.
* `running.rst <running.rst>`_: commands to reproduce studies, tests, and
  verification.
* `improvement_plan.rst <improvement_plan.rst>`_: longer-term implementation
  and evaluation plan.

Fitting, scoring, and noise models
----------------------------------

* `implementation_details.rst <implementation_details.rst>`_: ASV's current
  fitting implementation.
* `mathematical_analysis.rst <mathematical_analysis.rst>`_: mathematical
  structure of the existing fit and exact reference.
* `correlation_design.rst <correlation_design.rst>`_: exact fitting of
  residual correlation.
* `shared_noise_floor.rst <shared_noise_floor.rst>`_ and
  `noise_persistence.rst <noise_persistence.rst>`_: shared noise floor and
  bounded persistence models.
* `exact_reference.rst <exact_reference.rst>`_: exact small-history reference
  for the shared-floor score.
* `baseline_jump_prior.rst <baseline_jump_prior.rst>`_: separates baseline,
  jump, and location-prior costs.
* `baseline_report.rst <baseline_report.rst>`_: initial findings.
* `version_comparison.rst <version_comparison.rst>`_: paired comparison of
  detector versions.

Confidence and reporting references
-----------------------------------

These documents build increasingly careful ways to decide whether a fitted
change is statistically convincing.

* `reporting_uncertainty.rst <reporting_uncertainty.rst>`_: uncertainty when
  asking whether a slowdown exceeds a threshold.
* `direct_protocol.rst <direct_protocol.rst>`_ and
  `direct_report.rst <direct_report.rst>`_: direct test of below-threshold
  explanations.
* `direct_stress_protocol.rst <direct_stress_protocol.rst>`_ and
  `direct_stress_report.rst <direct_stress_report.rst>`_: sensitivity to
  changing variance and serial correlation.
* `known_covariance.rst <known_covariance.rst>`_: generalized least-squares
  reference when covariance is supplied.
* `covariance_protocol.rst <covariance_protocol.rst>`_ and
  `covariance_report.rst <covariance_report.rst>`_: evaluation of that
  known-covariance reference.
* `unknown_correlation.rst <unknown_correlation.rst>`_: confidence set when
  correlation is not known.
* `ar1_protocol.rst <ar1_protocol.rst>`_ and `ar1_report.rst <ar1_report.rst>`_:
  fresh evaluation with unknown AR(1) correlation.
* `ar1_information.rst <ar1_information.rst>`_: diagnosis of information lost
  by the correlation confidence set.
* `full_history_confidence.rst <full_history_confidence.rst>`_: prediction
  using the full history.
* `indexed_protocol.rst <indexed_protocol.rst>`_ and
  `indexed_report.rst <indexed_report.rst>`_: candidate-specific prediction.
* `reporting_protocol.rst <reporting_protocol.rst>`_ and
  `reporting_report.rst <reporting_report.rst>`_: joint calibration of
  reporting checks.

Fresh-data experiments
----------------------

Each protocol states the frozen design; its report records the resulting
measurements.

* Threshold detection: `threshold_protocol.rst <threshold_protocol.rst>`_,
  `threshold_report.rst <threshold_report.rst>`_.
* Component comparison: `ablation_protocol.rst <ablation_protocol.rst>`_,
  `ablation_report.rst <ablation_report.rst>`_.
* Robustness under irregular data: `robustness_protocol.rst <robustness_protocol.rst>`_,
  `robustness_report.rst <robustness_report.rst>`_.

Directional confidence and the latest certificate chain
--------------------------------------------------------

These references focus on the mathematical certificate used to reject a
candidate explanation. Read in order for the derivation, implementation, and
complete saved-case search.

* `residual_direction.rst <residual_direction.rst>`_: remove fitted plateau
  levels and scale before comparing residual patterns.
* `directional_tail.rst <directional_tail.rst>`_: bound the probability of a
  residual direction, including near correlation one.
* `directional_reporting.rst <directional_reporting.rst>`_: turn point and
  interval bounds into a full reporting decision.
* `directional_protocol.rst <directional_protocol.rst>`_ and
  `directional_report.rst <directional_report.rst>`_: evaluate directional
  confidence on fresh histories.
* `directional_determinant.rst <directional_determinant.rst>`_ and
  `determinant_reporting.rst <determinant_reporting.rst>`_: improve the
  determinant part of the certificate and replay saved cases.
* `spectral_refinement.rst <spectral_refinement.rst>`_: use more covariance
  directions and a second tilt.
* `sturm_refinement.rst <sturm_refinement.rst>`_: derive split-aware
  eigenvalue bounds and heterogeneous Fourier density bounds.
* `sturm_reporting.rst <sturm_reporting.rst>`_: complete depth-17 replay and
  remaining saved explanation.
* `data/sturm_reporting_v3_depth17_summary.json <data/sturm_reporting_v3_depth17_summary.json>`_
  and `data/sturm_reporting_v3_depth17_diagnosis.json.gz <data/sturm_reporting_v3_depth17_diagnosis.json.gz>`_:
  verified current replay results and full archive.

Python modules by role
----------------------

Experiment setup and scoring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `harness.py`_: reproducible experiments using ASV's production scoring
  code.
* `fetch_numpy.py`_: save public NumPy timing histories for examples.
* `noise_model.py`_: shared-floor scoring reference.
* `rational_polynomial.py`_: exact polynomial and Bernstein interval tools.
* `exact_reference.py`_: exact small-history independent reference.
* `compare_versions.py`_: paired before/after comparison.
* `inspect_comparison.py`_: explain selected histories with exact fits.

Reporting and confidence references
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `reporting_ar1.py`_, `reporting_ar1_full.py`_, and
  `reporting_ar1_jump.py`_: correlation and full-history predictive
  references.
* `reporting_direct.py`_ and `reporting_covariance.py`_: direct threshold
  tests with independent or supplied covariance.
* `reporting_reference.py`_ and `reporting_signs.py`_: conservative median
  evidence and sign-based calibration.
* `reporting_determinant.py`_, `reporting_directional.py`_,
  `reporting_spectral.py`_, and `reporting_sturm.py`_: successive complete
  reporting references.
* `residual_direction.py`_, `directional_tail.py`_,
  `directional_determinant.py`_, `directional_spectral.py`_, and
  `directional_sturm.py`_: mathematical certificate components.

Study drivers
~~~~~~~~~~~~~

* `threshold_study.py`_, `ablation_study.py`_, and `robustness_study.py`_:
  threshold, component, and robustness experiments.
* `direct_study.py`_ and `direct_stress.py`_: direct-test evaluation and
  assumption stress test.
* `covariance_study.py`_: known-covariance evaluation.
* `ar1_study.py`_ and `indexed_study.py`_: unknown-correlation and
  candidate-specific prediction studies.
* `reporting_study.py`_: reporting calibration study.
* `directional_study.py`_: fresh comparison of directional and predictive
  confidence.
* `sturm_fresh_study.py`_ and `sturm_fresh_protocol.rst`_: frozen independent
  comparison of the rank-group point-certificate refinement.

Diagnostics, replay, and verification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `ar1_information.py`_, `full_history_diagnosis.py`_,
  `jump_prior_diagnosis.py`_, and `directional_losses.py`_: diagnose saved
  cases and explain where confidence bounds lose information.
* `determinant_diagnosis.py`_, `sturm_diagnosis.py`_, and
  `directional_diagnosis.py`_: construct saved mathematical diagnoses.
* `spectral_replay.py`_, `sturm_replay.py`_, and
  `determinant_replay.py`_: replay saved cases.
* `spectral_replay_artifacts.py`_, `sturm_replay_artifacts.py`_,
  `determinant_replay_artifacts.py`_, `directional_artifacts.py`_, and
  `indexed_artifacts.py`_: verify and archive results.
* `spectral_verification.py`_, `sturm_verification.py`_, and
  `determinant_verification.py`_: independently check reporting routes.
* `directional_recovery.py`_: resume a saved evaluation without changing its
  calculations.

Tests
~~~~~

Tests are named by the module or study they validate. Each can be run with
``pytest experiments/step_detection/test_<name>.py``. The current certificate
chain is covered by `test_residual_direction.py`_,
`test_directional_tail.py`_, `test_directional_determinant.py`_,
`test_directional_spectral.py`_, `test_directional_sturm.py`_,
`test_reporting_determinant.py`_, `test_reporting_directional.py`_,
`test_reporting_spectral.py`_, and `test_reporting_sturm.py`_. The remaining
tests cover their corresponding study drivers, fitting references, and
statistical derivations:

* Study and harness checks: `test_harness.py`_, `test_compare_versions.py`_,
  `test_exact_reference.py`_, `test_noise_model.py`_,
  `test_threshold_study.py`_, `test_ablation_study.py`_,
  `test_robustness_study.py`_, `test_direct_study.py`_,
  `test_direct_stress.py`_, `test_covariance_study.py`_,
  `test_ar1_study.py`_, `test_indexed_study.py`_,
  `test_reporting_study.py`_, `test_directional_study.py`_, and
  `test_directional_recovery.py`_, and `test_sturm_fresh_study.py`_.
* Statistical-reference checks: `test_ar1_information.py`_,
  `test_reporting_ar1.py`_, `test_reporting_ar1_full.py`_,
  `test_reporting_ar1_jump.py`_, `test_reporting_covariance.py`_,
  `test_reporting_direct.py`_, `test_reporting_reference.py`_, and
  `test_reporting_signs.py`_.

Older or supporting documents
-----------------------------

* `implementation_details.rst <implementation_details.rst>`_,
  `improvement_plan.rst <improvement_plan.rst>`_,
  `baseline_report.rst <baseline_report.rst>`_, and
  `running.rst <running.rst>`_ are general references.
* `correlation_design.rst <correlation_design.rst>`_,
  `shared_noise_floor.rst <shared_noise_floor.rst>`_,
  `noise_persistence.rst <noise_persistence.rst>`_,
  `baseline_jump_prior.rst <baseline_jump_prior.rst>`_,
  `reporting_uncertainty.rst <reporting_uncertainty.rst>`_, and
  `known_covariance.rst <known_covariance.rst>`_ contain supporting
  derivations.
* `lean/README.rst <lean/README.rst>`_ describes the isolated Lean proof
  pilot.

.. _harness.py: harness.py
.. _fetch_numpy.py: fetch_numpy.py
.. _noise_model.py: noise_model.py
.. _rational_polynomial.py: rational_polynomial.py
.. _exact_reference.py: exact_reference.py
.. _compare_versions.py: compare_versions.py
.. _inspect_comparison.py: inspect_comparison.py
.. _reporting_ar1.py: reporting_ar1.py
.. _reporting_ar1_full.py: reporting_ar1_full.py
.. _reporting_ar1_jump.py: reporting_ar1_jump.py
.. _reporting_direct.py: reporting_direct.py
.. _reporting_covariance.py: reporting_covariance.py
.. _reporting_reference.py: reporting_reference.py
.. _reporting_signs.py: reporting_signs.py
.. _reporting_determinant.py: reporting_determinant.py
.. _reporting_directional.py: reporting_directional.py
.. _reporting_spectral.py: reporting_spectral.py
.. _reporting_sturm.py: reporting_sturm.py
.. _residual_direction.py: residual_direction.py
.. _directional_tail.py: directional_tail.py
.. _directional_determinant.py: directional_determinant.py
.. _directional_spectral.py: directional_spectral.py
.. _directional_sturm.py: directional_sturm.py
.. _threshold_study.py: threshold_study.py
.. _ablation_study.py: ablation_study.py
.. _robustness_study.py: robustness_study.py
.. _direct_study.py: direct_study.py
.. _direct_stress.py: direct_stress.py
.. _covariance_study.py: covariance_study.py
.. _ar1_study.py: ar1_study.py
.. _indexed_study.py: indexed_study.py
.. _reporting_study.py: reporting_study.py
.. _directional_study.py: directional_study.py
.. _sturm_fresh_study.py: sturm_fresh_study.py
.. _sturm_fresh_protocol.rst: sturm_fresh_protocol.rst
.. _ar1_information.py: ar1_information.py
.. _full_history_diagnosis.py: full_history_diagnosis.py
.. _jump_prior_diagnosis.py: jump_prior_diagnosis.py
.. _directional_losses.py: directional_losses.py
.. _determinant_diagnosis.py: determinant_diagnosis.py
.. _sturm_diagnosis.py: sturm_diagnosis.py
.. _directional_diagnosis.py: directional_diagnosis.py
.. _spectral_replay.py: spectral_replay.py
.. _sturm_replay.py: sturm_replay.py
.. _determinant_replay.py: determinant_replay.py
.. _spectral_replay_artifacts.py: spectral_replay_artifacts.py
.. _sturm_replay_artifacts.py: sturm_replay_artifacts.py
.. _determinant_replay_artifacts.py: determinant_replay_artifacts.py
.. _directional_artifacts.py: directional_artifacts.py
.. _indexed_artifacts.py: indexed_artifacts.py
.. _spectral_verification.py: spectral_verification.py
.. _sturm_verification.py: sturm_verification.py
.. _determinant_verification.py: determinant_verification.py
.. _directional_recovery.py: directional_recovery.py
.. _test_residual_direction.py: test_residual_direction.py
.. _test_directional_tail.py: test_directional_tail.py
.. _test_directional_determinant.py: test_directional_determinant.py
.. _test_directional_spectral.py: test_directional_spectral.py
.. _test_directional_sturm.py: test_directional_sturm.py
.. _test_reporting_determinant.py: test_reporting_determinant.py
.. _test_reporting_directional.py: test_reporting_directional.py
.. _test_reporting_spectral.py: test_reporting_spectral.py
.. _test_reporting_sturm.py: test_reporting_sturm.py
.. _test_harness.py: test_harness.py
.. _test_compare_versions.py: test_compare_versions.py
.. _test_exact_reference.py: test_exact_reference.py
.. _test_noise_model.py: test_noise_model.py
.. _test_threshold_study.py: test_threshold_study.py
.. _test_ablation_study.py: test_ablation_study.py
.. _test_robustness_study.py: test_robustness_study.py
.. _test_direct_study.py: test_direct_study.py
.. _test_direct_stress.py: test_direct_stress.py
.. _test_covariance_study.py: test_covariance_study.py
.. _test_ar1_study.py: test_ar1_study.py
.. _test_indexed_study.py: test_indexed_study.py
.. _test_reporting_study.py: test_reporting_study.py
.. _test_directional_study.py: test_directional_study.py
.. _test_directional_recovery.py: test_directional_recovery.py
.. _test_sturm_fresh_study.py: test_sturm_fresh_study.py
.. _test_ar1_information.py: test_ar1_information.py
.. _test_reporting_ar1.py: test_reporting_ar1.py
.. _test_reporting_ar1_full.py: test_reporting_ar1_full.py
.. _test_reporting_ar1_jump.py: test_reporting_ar1_jump.py
.. _test_reporting_covariance.py: test_reporting_covariance.py
.. _test_reporting_direct.py: test_reporting_direct.py
.. _test_reporting_reference.py: test_reporting_reference.py
.. _test_reporting_signs.py: test_reporting_signs.py
