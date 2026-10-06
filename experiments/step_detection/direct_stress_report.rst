What breaks when benchmark noise changes?
=========================================

The direct test needs its noise assumptions. On fresh histories, its gated
false alerts rise from 3/240 under independent constant-variance noise to
22/240 with rising variance, 22/240 with correlation, and 33/240 with both.
The last condition also reduces true detections from 124/160 to 97/160.
These results support extending the noise model before production use.

The failure has two parts. The test can underestimate uncertainty about the
slowdown's size. It can also mistake noisy excursions for an additional step,
discarding a valid below-threshold explanation. Correcting the first part
alone leaves the second route to false alerts open.

What was held fixed
-------------------

The `mathematical protocol <direct_stress_protocol.rst>`_, code, and inherited
settings were frozen in commit ``1f0efbf`` before generating the histories.
The shared-floor fit, existing reporting, sign references, and direct-test
cutoffs are unchanged from the `preceding study <direct_report.rst>`_.

There are 400 base histories, each transformed into six noise conditions.
All observations have Gaussian noise; correlation is either zero or 0.7.
Noise amplitude is constant, rises threefold, or falls threefold at the
nominal change position. Threefold amplitude means ninefold variance.
The six versions share the same initial draw and innovation sequence.

Each condition includes 160 positive histories with true increases of 6%
or 8%, and 240 null histories with increases of 0%, 4%, or exactly 5%.
The target claim is that the true slowdown exceeds 5%, so exactly 5% belongs
in the null. Histories have 40 or 100 readings, with changes halfway through
or leaving only four or ten later readings. Base noise standard deviations
are 0.05 and 0.2 around a baseline of 10.

The direct test checks every possible change location. At each location it
can reject a below-threshold explanation through a size test or an
extra-boundary shape test. It alerts only after rejecting every location.
The direct gate further requires the existing shared-floor reporter to alert.
That gate certifies a history-level claim, not a particular reported revision.

The direct test's 5% error bound applies only to the independent,
constant-variance control here. The sign references also allow independent
unequal variance, but their independence assumption excludes the correlated
conditions. These are per-history guarantees under the respective models,
not guarantees that every finite sample's observed rate is at most 5%.

Results by noise condition
--------------------------

Each cell below is ``true detections / false alerts``, with denominators
160 and 240 respectively. Production uses ASV's existing defaults;
shared existing uses the frozen experimental fit with existing reporting.

.. list-table:: Same histories, unchanged reporting rules
   :header-rows: 1

   * - Noise
     - Production
     - Shared existing
     - Joint-sign gate
     - Direct gate
   * - Independent, constant
     - 136 / 33
     - 159 / 48
     - 24 / 0
     - 124 / 3
   * - Independent, rises
     - 99 / 23
     - 144 / 55
     - 10 / 0
     - 103 / 22
   * - Independent, falls
     - 94 / 21
     - 134 / 42
     - 11 / 0
     - 67 / 1
   * - Correlated, constant
     - 82 / 20
     - 132 / 38
     - 29 / 2
     - 123 / 22
   * - Correlated, rises
     - 61 / 14
     - 99 / 44
     - 16 / 1
     - 97 / 33
   * - Correlated, falls
     - 33 / 7
     - 95 / 36
     - 22 / 1
     - 77 / 15

Direct-gate false-alert rates are 1.25%, 9.17%, 0.42%, 9.17%, 13.75%,
and 6.25% in table order. Falling variance under independence trades away
sensitivity: it loses 57 of the control's true detections and gains none.
Rising variance loses 27 and gains six. Adding correlation at constant
variance loses nine and gains eight, while adding 19 false alerts and
removing none. Similar detection counts can therefore hide worse reliability.

The rates mix flat, 4%, and 5% null histories. At exactly 5%, direct-gate
false alerts are 3, 16, 1, 20, 22, and 10 out of 80 per condition.
There are also six 4% false alerts with independent rising variance;
correlated rising variance produces ten at 4% and one on a flat history.
The failure is not confined to rounding around the threshold.

Across all versions the direct gate detects 591/960 positive histories
and alerts on 96/1440 null histories. The corresponding production counts
are 505 and 118, but production has fewer false alerts in all three correlated
conditions. An aggregate improvement would conceal that distinction.
Versions of the same base history are dependent, so these totals are
descriptive counts, not 2400 independent trials. This study also omits the
preceding study's largest base noise level; their aggregate rates are not
a before/after comparison.

Why a noisy short plateau causes trouble
----------------------------------------

At a proposed split, the size test estimates
``D = mean_after - 1.05*mean_before``. It uses residuals on both sides to
estimate a common noise variance. That estimate is useful when the variance
really is shared: a long earlier plateau helps when later readings are few.

With 36 quiet readings followed by four readings whose noise amplitude is
three times larger, the common estimate is dominated by the quiet side.
The true variance of D is 4.981 times the expected variance reported by
the frozen test. Reversing the amplitude change gives a ratio of 0.224,
making the variance estimate too large on average. These ratios were derived
before evaluation; they are not fitted to the observed alert rates.

The location breakdown agrees with that mechanism. Under independent rising
variance, the direct gate makes 19/120 false alerts on recent-change null
histories versus 3/120 with a middle split. In the control those counts are
1/120 and 2/120. The standalone direct test makes 37/120 recent-change
false alerts under rising variance; existing reporting filters 18 of them.

For example, ``independent_up-direct-stress-n40-d0.04-s0.02-recent-seed1004``
has a true change from 10 to 10.4 at split 36. Its fitted medians happen to
be approximately 10.060 and 11.034. At the true split, the direct size
statistic is 3.230, exceeding the frozen cutoff 1.799. The shape statistic
is only 4.566, below its cutoff 16.280. This valid 4% explanation is rejected
through the size route, and the full test rejects every other location too.

Why fixing the size test is insufficient
----------------------------------------

The shape test asks whether adding a third plateau explains too much of the
remaining variation. Under the original model, Gaussian residual projections
give the calibrated F distribution. Changing variance or correlating noise
breaks that argument: persistent or unusually variable noise can reward the
extra plateau without an extra change in the underlying benchmark level.

For each false alert, we inspected the test at the true generating split.
For a flat history we used its designated nominal split. This location is
known only to the evaluator; the reporting algorithm still checks all splits.

.. list-table:: Rejection route at that split, among direct-gate false alerts
   :header-rows: 1

   * - Noise
     - Size only
     - Shape only
     - Both
   * - Independent, constant
     - 3
     - 0
     - 0
   * - Independent, rises
     - 10
     - 10
     - 2
   * - Independent, falls
     - 1
     - 0
     - 0
   * - Correlated, constant
     - 8
     - 8
     - 6
   * - Correlated, rises
     - 7
     - 14
     - 12
   * - Correlated, falls
     - 4
     - 9
     - 2

These counts identify which check rejects the valid explanation. They do
not assign each global alert uniquely to one route across every location.

Consider ``independent_up-direct-stress-n40-d0.04-s0.005-recent-seed1000``.
The true split is 36 and the true increase is 4%. Its size statistic there
is -3.956, offering no evidence that the increase exceeds 5%. But adding a
boundary at 39 isolates one noisy final reading and produces F=41.691,
above 16.280. The shape route rejects the valid explanation. The shared
fit itself selects a final one-reading plateau, and the gated rule alerts.

Correlation creates the same issue without unequal variances. In
``correlated_constant-direct-stress-n100-d0.04-s0.02-middle-seed1009``,
the valid split is 50. Its size statistic is -3.664, but adding a boundary
at 76 gives F=61.984, above the length-100 cutoff 16.427. Persistent noise
makes an extra plateau look convincing under the independent-noise test.

Without the existing-report gate, direct false-alert counts are 3, 41, 1,
54, 82, and 51 per condition. Existing reporting removes many, but it is
not a replacement for a calibrated evidence test.

The next mathematical step
--------------------------

Both tests should use the same declared covariance model. Start with
Gaussian noise whose covariance is ``sigma**2 * R``, where R is known and
positive definite while the overall scale sigma is unknown. This permits
a clean reference before estimating correlation or variance ratios from data.

At every proposed split, fit the two plateau levels by generalized least
squares. Equivalently, transform both the observations and the plateau
design using a matrix W satisfying ``W*R*W' = I``. The transformed errors
are independent with common variance. The design must be transformed too;
whitened plateau indicators are generally not ordinary constant segments.

For the size test, use the contrast of those fitted levels and its covariance
from ``(X' * R_inverse * X)_inverse``. Estimate the overall scale from the
weighted residual sum with n-2 degrees of freedom. For the shape scan,
compare nested two- and three-plateau models in that same transformed space.
The t and F arguments can then be derived together. The all-locations
rejection logic and error-budget allocation can be retained.

A supplied R would be an oracle reference, not yet a deployable solution.
The practical question is how to account for uncertainty in R when only a
short benchmark history is available. Simply inserting an estimated R does
not inherit the known-R guarantee. Derive a nuisance-parameter treatment or
an explicitly limited calibration scheme before evaluating that extension
on fresh data. The completed stress cases should remain a diagnostic set.

Artifacts and checks
--------------------

The `frozen configuration <data/direct_stress_v1_frozen.json>`_ records the
inherited cutoffs, source hashes, covariance calculations, and design.
The `result index <data/direct_stress_v1_results.json>`_ records summaries,
paired gains/losses, rejection routes, examples, and SHA-256 hashes for 48
compressed input/result archives. Together they preserve all 2400 numerical
histories and their complete reporting outputs.

All 176 focused tests passed, including 19 new checks of the covariance
calculation, paired noise construction, rejection routes, and unchanged
Python/native pipelines. Archive verification reproduces the summaries and
checks that every global direct alert rejects its generating split.
See `reproduction commands <running.rst>`_. Production behavior is unchanged.
