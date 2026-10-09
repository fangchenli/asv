Independent evaluation of the rank-group Fourier bound
=======================================================

This study evaluates the exact Fourier refinement that bounds four adjacent
covariance directions separately and combines their lower bounds through a
geometric mean. The comparison control runs the same complete reporting
search with the split-aware Sturm route disabled. The candidate enables that
route. Both methods use the same history, calibration, confidence allocation,
size/shape checks, and work limits.

The study freezes source hashes, numerical-library versions, calibration,
correlation conditions, and work limits before generating its observations.
It uses a new stream prefix and seed labels 1602 and 1603, separate from the
earlier directional evaluation. The design has 360 correlated histories
from 120 base histories: 144 slowdowns above the five-percent threshold and
216 null histories. The three correlation versions of each base history are
dependent, so counts by correlation are descriptive rather than independent
trials.

The frozen settings retain the 0.01 confidence allocation, 0.04 reporting
allocation, 4096-cell limit, and depth 17. An unresolved search is a
non-alert. No threshold or setting may be changed after looking at the
evaluation outcomes. A later mathematical change requires another freeze
and another independent evaluation.

The run saves compressed inputs and full per-method records, including
surviving explanations, incomplete searches, cell counts, and elapsed time.
The summary reports alerts and statuses by correlation, length, change,
noise, location, and positive/null label, along with paired gained and lost
alerts. Results are a comparison of experimental references; they are not a
production benchmark or a claim about the C++ detector.

Create and commit the freeze before running the evaluation::

    .venv/bin/python -m experiments.step_detection.sturm_fresh_study freeze \
        --output experiments/step_detection/data/sturm_fresh_v2_frozen.json

Then run the frozen comparison::

    VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        .venv/bin/python -m experiments.step_detection.sturm_fresh_study run \
        --frozen experiments/step_detection/data/sturm_fresh_v2_frozen.json \
        --output experiments/step_detection/runs/sturm_fresh_v2 --workers 6
