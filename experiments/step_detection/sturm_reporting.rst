The complete search with split-aware eigenvalue bounds
======================================================

The grouped Fourier bound was replayed on all seven saved cases and all seven
results independently verify. Four previous alerts remain alerts; no
previously non-alert case becomes an alert. The early negative-correlation
case now returns a verified surviving explanation at split 2 and rho=99/128.
Its final decision remains a non-alert.

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
Work limits are 4096 visited cells and depth 17. These calculations
refine upper bounds on the same p-value and require no extra confidence
allocation. Reporting's t/F cutoffs remain numerical, as in the previous
references, under the stationary Gaussian common-variance assumptions.

The first lost detection
------------------------

This history has a generated 8% slowdown after reading 25. The detector
must also consider explanations placing the change elsewhere. Its previous
survivor proposed a change after reading 1 and correlation 51/64.

The earlier 51/64 explanation is rejected over its saved neighborhood. The
next point, 101/128, is also excluded by the heterogeneous Fourier bound: its
point certificate is about **0.0087203**, and the bound stays below
**0.0087751** on [51711/65536, 51713/65536]. Both use the 1/8 tilt and all
sixteen directions.

At depth 16, the full search could not certify the narrow interval
[16185/16384, 32371/32768] near rho=0.988. Splitting it once produces two
intervals that the determinant certificate can exclude. The default depth
17 search then finds a different surviving explanation at split 2 and
rho=99/128. Its probability upper bound is 0.012290, above the 0.01 cutoff.
This is a valid reason the history still does not alert: excluding one
explanation exposes another. The early case visits 1,219 cells at depth 17.

The rank-group point refinement computes a lower bound for each covariance
rank, groups four adjacent ranks, and uses their geometric mean in the same
Fourier density integral. This tightens the 99/128 point bound to
**0.0108498**, an 11.7% reduction, while leaving it above 0.01. The floating
tail estimate remains about 0.007726, with an estimated quadrature error of
about 4.5e-10; it is diagnostic only. Thus this refinement does not turn the
case into an alert. It applies only to candidate points; interval coverage
continues to use the existing certified interval route.

The updated `rank-group replay summary
<data/sturm_reporting_v4_rank_groups_summary.json>`_ records all seven saved
cases. The compressed `rank-group replay archive
<data/sturm_reporting_v4_rank_groups_diagnosis.json.gz>`_ preserves the full
proofs and verification details. All seven outcomes verify, with no alert
decisions gained or lost.

The other six cases
--------------------

The second lost detection retains split 1 and rho=7/8. Its new heterogeneous
Fourier bound is about **0.032034**, improved from about 0.04603 but still
above the 0.01 cutoff. The saved numerical directional probability estimate
is about 0.01682, so this remains a separate statistical-bound gap.

The archived independent and positively correlated 8% examples remain
certified alerts. Deterministic 20% changes with 24 and 100 readings also
remain certified alerts, and the deterministic 4% example keeps its middle
split with zero correlation. All four alert partitions and their routes
are identical to the preceding replay.

The Sturm route certifies 588 interior intervals in the early case. The four
alerts still have complete correlation-domain coverage, including the
inherited routes. Separate saved proofs cover the 51/64 and 101/128 intervals
removed from the first lost case.

Checks and reproducibility
---------------------------

The `verifier <sturm_verification.py>`_ reconstructs residual states
independently of the reporting fits, checks new point and interval proofs,
and delegates inherited routes to the previous verifier. Every alert must
cover every split's full correlation domain without gaps or overlaps.
Surviving explanations are independently checked with dense GLS.

The 24 integration tests cover complete coverage, survivors, unit
conversion, calibration, degenerate data, work limits, and damaged proofs
or metadata, including changed eigenvalue precision. The 38 mathematical
tests cover exact Sturm roots, the integer recurrence against its rational
reference, eigenvalue bounds, heterogeneous Fourier integrals, individual
rank groups, and saved point and interval certificates. Disabling the Sturm route recovers
the original directional search exactly after removing the additional
metadata.

The `depth-16 summary <data/sturm_reporting_v2_heterogeneous_summary.json>`_
preserves the earlier unresolved run. The latest `depth-17 summary
<data/sturm_reporting_v3_depth17_summary.json>`_ contains all seven current
decisions. Its compressed `full archive
<data/sturm_reporting_v3_depth17_diagnosis.json.gz>`_ retains the
inputs, previous and current results, certificates, source hashes, and
numerical diagnoses. The archive builder checks each result and saved
witness interval against the preceding archive. Commands are in the
`running guide <running.rst>`_.

This remains development evidence on saved cases. The earlier fresh-study
detection count remains 58/144; the new rank-group method is being evaluated
on a separately frozen fresh set. ASV's production detector is unchanged.
