"""Design checks for the independently frozen Sturm comparison."""

from collections import Counter

from experiments.step_detection import sturm_fresh_study as study


def test_fresh_design_is_separate_and_has_expected_denominators():
    rows = list(study.cases())
    assert len(rows) == 360
    assert sum(row['positive'] for row in rows) == 144
    assert len({row['id'] for row in rows}) == 360
    assert {row['seed'] for row in rows} == {1602, 1603}
    assert len({row['pair_id'] for row in rows}) == 120
    assert set(Counter(row['pair_id'] for row in rows).values()) == {3}
    assert all(row['pair_id'].startswith('sturm-ar1-fresh-v1-') for row in rows)


def test_freeze_is_stable_and_covers_the_new_certificate_sources():
    frozen = study.freeze()
    assert frozen == study.freeze()
    assert frozen['work'] == {'max_cells': 4096, 'max_depth': 17}
    assert frozen['source_hashes']['directional_sturm.py']
    assert frozen['source_hashes']['sturm_fresh_protocol.rst']
