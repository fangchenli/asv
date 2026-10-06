Frozen evaluation with the true covariance supplied
===================================================

The `known-covariance derivation <known_covariance.rst>`_ corrects the size
and shape tests together. This study measures the resulting sensitivity
and false alerts on fresh histories, with the actual generating covariance
shape supplied. It evaluates an oracle reference, before estimating noise
structure from data.

What is fixed
-------------

Inherit the complete `variance stress study <direct_stress_protocol.rst>`_
configuration, including the shared-floor fit, existing reporting, sign
references, and the old direct test. Append the GLS evidence alone and
as a gate on existing shared-floor alerts. Use the inherited analytical
cutoffs unchanged: total alpha 0.05, split 0.04 for size and 0.01 for shape,
with a 5% slowdown threshold. No calibration simulation or parameter search
is needed for the known-covariance model.

Keep the six conditions: independent or AR(1) correlation 0.7, crossed with
constant, threefold-rising, or threefold-falling noise amplitude. Lengths
are 40 and 100, true increases are 0%, 4%, 5%, 6%, and 8%, base standard
deviations are 0.05 and 0.2, and changes occur in the middle or leave only
four/ten later readings. Use fresh seed labels 1200 through 1209 with the
unchanged stress generator. The seed labels are disjoint from all preceding
evaluation sets using that generator.

There are 400 base histories, with six versions sharing their independent
Gaussian initial draw and innovations. Each condition contains 160 positive
and 240 null histories, of which 80 are exactly 5%. Exactly 5% stays in the
null because the reporting claim is that the true change exceeds 5%.
Across conditions the versions are dependent; aggregate counts are descriptive.

The oracle covariance
---------------------

For amplitude factors s_i and correlation rho, supply::

    R_ij = s_i*s_j*rho**abs(i-j)

The generator uses the stationary initial distribution and noise sigma*s_i*z_i,
so its covariance is sigma squared times this R. Do not supply sigma or
estimate R. Its absolute scale is immaterial to the test. The matrix is
used unchanged for every proposed mean split and extra boundary. The
test still searches all n-1 mean splits; the generating split is used only
by the evaluator to identify rejection routes.

The amplitude transition is tied to the nominal mean boundary in this
design. Supplying R therefore reveals a noise-transition location in the
unequal-variance conditions. That is additional oracle information even
though it does not constrain the mean-location search. Report this limitation
alongside the results. A flat history has the same possible noise transition
without a mean change. All six conditions satisfy the new test's declared
Gaussian one-change model, while only the independent constant-variance
control satisfies the old direct test's model.

Save the 24 condition/length/location matrices once and link each numerical
input to its matrix identifier. Include their canonical JSON hashes in
the freeze. These matrices are determined by the design without generating
evaluation observations. Hash the new runner, protocol, covariance reference,
derivation, and all inherited frozen sources. Record Python, NumPy, and
SciPy versions. Commit this freeze before generating the evaluation data.

Evaluation and checks
---------------------

For every history, retain the complete old pipeline and append covariance
statistics, numerical status, rejection lists, and the route at the generating
split. Numerical abstentions count as non-alerts and are reported separately.
They are not reclassified as surviving statistical explanations.

Report detection and false-alert counts by condition, length, location,
noise, and true change. For each condition, count paired true detections
gained/lost and false alerts added/removed relative to the old direct test
on the same histories, both standalone and gated. Also preserve comparisons
of each condition to its paired independent constant-variance control.
Count size-only, shape-only, and joint rejection routes at the generating
split among false alerts for both methods.

In the independent constant-variance control, require the GLS and old direct
decisions to agree, and their numerical statistics to agree within relative
tolerance 1e-7 and absolute tolerance 1e-8. An unexpected abstention there
stops the run. On all other conditions, save abstentions explicitly. Every
oracle global alert must reject the generating split, and gated decisions
must be a subset of existing shared-floor reports.

Before freezing, test the oracle covariance against the independent innovation
matrix construction, generator pairing with test-only seed labels, complete
reuse of the old pipeline, identity-control agreement, paired counting, and
abstention handling. Do not inspect evaluation seeds during these checks.

Parallel workers may evaluate independent histories; preserve input order
when writing outputs and record the worker count and BLAS thread settings.
Execution parallelism does not change the design, cutoffs, or decision rules.
No cutoff or implementation changes are made after seeing results. The next
decision is whether the known-noise results justify tackling covariance
estimation, with sensitivity losses reported alongside false-alert reductions.
