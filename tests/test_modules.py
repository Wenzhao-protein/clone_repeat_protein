import json
from pathlib import Path

import pandas as pd

from hurdler.modules import (
    _parse_fasta,
    longest_adjacent_exact_repeat,
    merge_module_catalogs,
    parse_fasta_modules,
    refine_module_boundaries,
    reselect_module_boundaries,
    _selected_middle_unit_from_row,
)
from hurdler.periodicity import scan_repeat_periods


def test_exact_repeat_extraction_prefers_longest_unit():
    assert longest_adjacent_exact_repeat("MNNNACDEFGACDEFGQQQ", minimum=6) == ("ACDEFG", 4)


def test_fasta_module_parser(tmp_path: Path):
    fasta = tmp_path / "designs.fasta"
    fasta.write_text(">DHR_TEST\nMNNNACDEFGACDEFGQQQ\n")
    frame = parse_fasta_modules([fasta], family="DHR", evidence_tier="B", source_url="https://example.test")
    assert len(frame) == 1
    assert frame.iloc[0]["unit_sequence"] == "ACDEFG"
    assert frame.iloc[0]["module_id"] == "designed_DHR_TEST"
    assert frame.iloc[0]["full_sequence"] == "MNNNACDEFGACDEFGQQQ"
    assert frame.iloc[0]["full_sequence_origin"] == 1


def test_repeatsdb_partial_fasta_retains_absolute_origin():
    sequence, origin = _parse_fasta(">PDB 1abc, chain A, range 103:114\nACDEFGHIKLMN\n")
    assert sequence == "ACDEFGHIKLMN"
    assert origin == 103


def test_catalog_deduplicates_within_but_not_across_collections(tmp_path: Path):
    base = {
        "family": "test",
        "unit_sequence": "ACDEFG",
        "evidence_tier": "A",
        "source_name": "source",
        "source_url": "https://example.test",
        "source_accession": "accession",
    }
    frame = pd.DataFrame(
        [
            {**base, "module_id": "natural", "collection": "natural100"},
            {**base, "module_id": "natural_duplicate", "collection": "natural100"},
            {**base, "module_id": "designed", "collection": "designed_all"},
        ]
    )
    catalog = merge_module_catalogs([frame], tmp_path / "catalog.parquet")
    assert len(catalog) == 2
    assert set(catalog["collection"]) == {"natural100", "designed_all"}
    assert catalog["download_date"].astype(str).str.len().gt(0).all()
    manifest = json.loads((tmp_path / "catalog.manifest.json").read_text())
    assert manifest["rows"] == 2
    assert manifest["source_mapping_rows"] == 3
    assert manifest["duplicate_source_rows_collapsed"] == 1
    assert manifest["collection_counts"] == {"designed_all": 1, "natural100": 1}


def test_refinement_keeps_source_unit_and_writes_fixed_variable_tables(tmp_path: Path):
    primitive_units = ["ACDEFGHIK", "ACNEFGHIK", "ACQEFGHIK", "ACDEFGHIK"]
    full_sequence = "MNPQ" + "".join(primitive_units) + "WYST"
    frame = pd.DataFrame(
        [
            {
                "module_id": "designed_double",
                "collection": "designed_all",
                "family": "test",
                "unit_sequence": primitive_units[0] * 2,
                "unit_start": 5,
                "unit_end": 22,
                "full_sequence": full_sequence,
                "full_sequence_origin": 1,
                "evidence_tier": "B",
                "source_name": "test",
                "source_url": "https://example.test",
                "source_accession": "test",
            }
        ]
    )
    units_path = tmp_path / "units.parquet"
    positions_path = tmp_path / "positions.parquet"
    refined = refine_module_boundaries(
        frame,
        unit_alignment_path=units_path,
        position_variability_path=positions_path,
    )
    row = refined.iloc[0]
    assert row["prior_unit_length"] == 18
    assert row["unit_length"] == 9
    assert row["unit_sequence"] == primitive_units[1]
    assert row["selected_module_index"] == 2
    assert row["variable_positions_json"] == "[3]"
    units = pd.read_parquet(units_path)
    positions = pd.read_parquet(positions_path)
    assert len(units) == 4
    assert units.loc[units.is_selected_module, "repeat_index"].tolist() == [2]
    assert positions.loc[~positions.fixed, "module_position"].tolist() == [3]


def test_reselection_uses_frozen_candidate_scores_without_rescanning():
    primitive = "ACDEFGHIK"
    full_sequence = "MNPQ" + primitive * 5 + "WYST"
    source = pd.DataFrame(
        [
            {
                "module_id": "double",
                "collection": "designed_all",
                "family": "test",
                "unit_sequence": primitive,
                "prior_unit_sequence": primitive * 2,
                "prior_unit_start": 5,
                "prior_unit_end": 22,
                "full_sequence": full_sequence,
                "full_sequence_origin": 1,
                "evidence_tier": "B",
                "source_name": "test",
                "source_url": "https://example.test",
                "source_accession": "test",
            }
        ]
    )
    candidates = pd.DataFrame(
        [
            {"module_id": "double", "candidate_rank": rank, **candidate.to_dict()}
            for rank, candidate in enumerate(
                scan_repeat_periods(
                    full_sequence,
                    prior_period=18,
                    prior_unit_start=5,
                ),
                start=1,
            )
        ]
    )
    refined = reselect_module_boundaries(source, candidates)
    assert refined.iloc[0].primitive_period == 9
    assert refined.iloc[0].unit_sequence == primitive
    assert refined.iloc[0].selected_module_index == 3
    assert refined.iloc[0].boundary_reselected_from_candidates


def test_reselection_retains_source_unit_when_scan_has_no_candidates():
    source = pd.DataFrame(
        [
            {
                "module_id": "source_only",
                "collection": "natural100",
                "family": "test",
                "unit_sequence": "ACDEFG",
                "prior_unit_sequence": "ACDEFG",
                "prior_unit_start": 10,
                "prior_unit_end": 15,
                "full_sequence": "M" * 9 + "ACDEFG" + "Q" * 9,
                "full_sequence_origin": 1,
                "evidence_tier": "A",
                "source_name": "test",
                "source_url": "https://example.test",
                "source_accession": "test",
            }
        ]
    )
    refined = reselect_module_boundaries(
        source,
        pd.DataFrame(columns=["module_id", "candidate_rank"]),
        secondary_structure_support=pd.DataFrame(columns=["module_id"]),
    )
    row = refined.iloc[0]
    assert row.boundary_refinement_status == "source_prior_fallback"
    assert row.unit_sequence == "ACDEFG"
    assert row.primitive_period == 6


def test_middle_unit_uses_region_midpoint_and_earlier_tie():
    sequence, start, end, index, count = _selected_middle_unit_from_row(
        {
            "unit_sequences_json": '["AAAA", "CCCCC", "DDDD", "EEEE"]',
            "source_unit_coordinates_json": "[[10, 13], [20, 24], [30, 33], [40, 43]]",
            "repeat_region_start": 10,
            "repeat_region_end": 43,
        },
        default_sequence="AAAA",
        default_start=10,
        default_end=13,
    )
    assert (sequence, start, end, index, count) == ("CCCCC", 20, 24, 2, 4)
