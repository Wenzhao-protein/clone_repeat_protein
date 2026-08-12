import importlib.util
from pathlib import Path

import pandas as pd

from hurdler.periodicity import scan_repeat_periods
from hurdler.secondary_structure import (
    _optional_dssp_float,
    collapse_dssp_state,
    map_author_template_to_full_sequence,
    score_candidates_with_secondary_structure,
    score_secondary_structure_period,
    select_joint_sequence_structure_candidate,
)


def _load_evidence_script():
    path = (
        Path(__file__).parents[1]
        / "studies/hurdler_validation/scripts/build_secondary_structure_evidence.py"
    )
    spec = importlib.util.spec_from_file_location("build_secondary_structure_evidence", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dssp_states_are_collapsed_to_three_state_alphabet():
    assert "".join(collapse_dssp_state(value) for value in "HGIEBTSC-") == "HHHEECCCC"


def test_dssp_explicit_na_numeric_fields_are_preserved_as_nan():
    assert _optional_dssp_float("1.25") == 1.25
    assert pd.isna(_optional_dssp_float("NA"))


def test_secondary_structure_period_uses_states_transitions_and_boundaries():
    pattern = "HHHHCCCEECC"
    support = score_secondary_structure_period(
        pattern * 5,
        period=len(pattern),
        local_start=0,
        local_end=len(pattern) * 5,
    )
    assert support.informative
    assert support.state_agreement == 1.0
    assert support.transition_agreement == 1.0
    assert support.phase_conservation == 1.0
    assert support.boundary_transition_fraction == 1.0
    assert support.score > 0.95


def test_author_template_maps_around_absolute_anchor():
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 3
    ss, residues = map_author_template_to_full_sequence(
        sequence,
        repeat_start=11,
        full_sequence_origin=1,
        ss3_template="HHHCC",
    )
    assert len(ss) == len(sequence)
    assert "?" not in ss
    assert residues.full_sequence_position.min() == 1
    assert residues.full_sequence_position.max() == len(sequence)
    assert ss[10:15] == "HHHCC"


def test_joint_selector_accepts_smallest_sequence_and_structure_harmonic():
    unit = "ACDEFGHIK"
    full_sequence = unit * 6
    candidates = scan_repeat_periods(
        full_sequence,
        prior_period=len(unit) * 2,
        prior_unit_start=1,
    )
    secondary_structure = "HHHCCCEEC" * 6
    supports = score_candidates_with_secondary_structure(candidates, secondary_structure)
    selected, support, reason = select_joint_sequence_structure_candidate(
        candidates, supports, prior_period=len(unit) * 2
    )
    assert selected is not None
    assert support is not None
    assert selected.period == len(unit)
    assert "jointly supported" in reason


def test_joint_selector_rejects_uninformative_all_helix_annotation():
    unit = "ACDEFGHIK"
    full_sequence = unit * 6
    candidates = scan_repeat_periods(
        full_sequence,
        prior_period=len(unit) * 2,
        prior_unit_start=1,
    )
    supports = score_candidates_with_secondary_structure(candidates, "H" * len(full_sequence))
    selected, support, reason = select_joint_sequence_structure_candidate(
        candidates, supports, prior_period=len(unit) * 2
    )
    assert selected is None
    assert support is None
    assert "source-annotated unit retained" in reason


def test_obsolete_candidate_scan_is_not_mixed_with_current_sequence():
    script = _load_evidence_script()
    candidates = pd.DataFrame(
        [
            {"module_id": "m1", "candidate_rank": 1, "period": 8, "local_start": 4, "local_end": 36},
            {"module_id": "m1", "candidate_rank": 2, "period": 4, "local_start": 8, "local_end": 24},
            {"module_id": "m1", "candidate_rank": 1, "period": 8, "local_start": 14, "local_end": 46},
            {"module_id": "m1", "candidate_rank": 2, "period": 5, "local_start": 20, "local_end": 55},
        ]
    )
    row = {
        "module_id": "m1",
        "full_sequence": "A" * 40,
        "full_sequence_origin": 1,
        "period": 8,
        "repeat_region_start": 5,
        "repeat_region_end": 36,
    }
    selected = script._candidate_records_for_row(row, candidates)
    assert len(selected) == 2
    assert selected[0]["local_start"] == 4
    assert max(record["local_end"] for record in selected) <= 40
