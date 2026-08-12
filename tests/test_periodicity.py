import json

import pytest

from hurdler.periodicity import infer_repeat_boundaries, scan_repeat_periods


def test_exact_tandem_with_flanks_selects_middle_primitive_unit():
    unit = "ACDEFGHIK"
    result, candidates = infer_repeat_boundaries(
        "MNPQ" + unit * 4 + "WYST",
        prior_unit_sequence=unit * 2,
        prior_unit_start=5,
    )
    assert result.period == len(unit)
    assert result.first_module_sequence == unit
    assert result.first_module_start == 5
    assert result.selected_module_index == 2
    assert result.selected_module_start == 14
    assert result.selected_module_sequence == unit
    assert result.repeat_count == 4
    assert result.left_flank_sequence == "MNPQ"
    assert result.right_flank_sequence == "WYST"
    assert result.harmonic_ratio == 2.0
    assert candidates[0].score >= result.score


def test_variable_positions_are_reported_from_all_units():
    units = ["ACDEFGHIK", "ACNEFGHIK", "ACQEFGHIK", "ACDEFGHIK"]
    result, _ = infer_repeat_boundaries(
        "MNPQ" + "".join(units) + "WYST",
        prior_unit_sequence=units[0],
        prior_unit_start=5,
    )
    assert result.period == 9
    assert result.variable_positions == (3,)
    assert result.variable_ranges == ((3, 3),)
    assert result.fixed_mask == "FFVFFFFFF"
    assert result.selected_module_index == 2
    assert result.selected_module_sequence == units[1]
    assert json.loads(result.to_dict()["unit_sequences_json"]) == units


def test_absolute_origin_is_preserved():
    unit = "ACDEFGHIK"
    result, _ = infer_repeat_boundaries(
        unit * 3,
        full_sequence_origin=103,
        prior_unit_sequence=unit,
        prior_unit_start=103,
    )
    assert result.repeat_region_start == 103
    assert result.first_module_end == 111
    assert result.selected_module_start == 112
    assert result.selected_module_end == 120
    assert result.repeat_region_end == 129


def test_homopolymer_is_not_a_supported_repeat():
    with pytest.raises(ValueError, match="No supported repeat period"):
        infer_repeat_boundaries("A" * 60, minimum_period=3)


def test_scan_exposes_longer_source_unit_and_primitive_harmonic():
    unit = "ACDEFGHIKLMNP"
    sequence = "WY" + unit * 4 + "QRST"
    candidates = scan_repeat_periods(
        sequence,
        prior_period=len(unit) * 2,
        prior_unit_start=3,
    )
    by_period = {candidate.period: candidate for candidate in candidates}
    assert len(unit) in by_period
    assert len(unit) * 2 in by_period
    assert by_period[len(unit)].repeat_count >= 4
