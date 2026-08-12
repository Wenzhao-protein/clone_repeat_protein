import json
from pathlib import Path

import pandas as pd
import pytest

import hurdler.module_experiments as module_experiments
from hurdler.designed_inventory import validate_af3_outputs
from hurdler.module_experiments import (
    finalize_adaptive_copy_results,
    finalize_module_compatibility,
    run_module_compatibility,
    summarize_compatibility,
)
from hurdler.repeatsdb import (
    annotation_repeat_regions,
    select_longest_region_per_protein,
)
from hurdler.structural_repeats import (
    _biotite_chain_and_dssp,
    foldseek_descriptor,
    infer_designed_boundary,
    map_structure_positions_to_full_sequence,
    project_foldseek_descriptor,
    scan_dual_evidence_periods,
)


def test_repeatsdb_pdb_content_loci_and_earlier_middle_unit():
    annotation = {
        "uuid": "pdb-fixture",
        "content": {
            "chain": {"structure": "1abc", "id": "A", "source": "RCSB/PDB"},
            "reviewed": True,
            "features_uniprot": ["P12345"],
            "loci": [
                {"type": "region", "start": 10, "end": 89},
                {"type": "unit", "parent": 0, "start": 10, "end": 29},
                {"type": "unit", "parent": 0, "start": 30, "end": 49},
                {"type": "insertion", "parent": 0, "start": 50, "end": 52},
                {"type": "unit", "parent": 0, "start": 53, "end": 70},
                {"type": "unit", "parent": 0, "start": 71, "end": 89},
            ],
        },
    }
    regions = annotation_repeat_regions(annotation)
    assert len(regions) == 1
    assert regions[0]["annotation_schema"] == "content_loci"
    assert json.loads(regions[0]["unit_coordinates_json"]) == [
        [10, 29],
        [30, 49],
        [53, 70],
        [71, 89],
    ]
    assert json.loads(regions[0]["insertion_coordinates_json"]) == [[50, 52]]
    coordinates = json.loads(regions[0]["unit_coordinates_json"])
    assert coordinates[(len(coordinates) - 1) // 2] == [30, 49]


def test_foldseek_descriptor_projects_terminal_insertions_without_shifting():
    projected, audit = project_foldseek_descriptor(
        "ACDEFG", "XXACDEFGY", "123456789"
    )
    assert projected == "345678"
    assert audit["foldseek_biotite_alignment_coverage"] == 1.0
    assert audit["foldseek_biotite_alignment_identity"] == 1.0
    assert audit["foldseek_biotite_unmapped_residues"] == 0


def test_foldseek_descriptor_rejects_low_coverage_mapping():
    with pytest.raises(ValueError, match="inadequate"):
        project_foldseek_descriptor("ACDEFGHIK", "ACD", "123")


def test_structure_positions_map_across_terminal_tags_and_linker_insertions():
    mapping, audit = map_structure_positions_to_full_sequence(
        "ACDEFGHIK", "MMACDEFGXXHIKHH"
    )
    assert mapping == [2, 3, 4, 5, 6, 7, 10, 11, 12]
    assert audit["structure_full_alignment_coverage"] == 1.0
    assert audit["structure_full_alignment_identity"] == 1.0


def test_af3_validation_materializes_missing_structure_status(tmp_path):
    inventory = tmp_path / "inventory.parquet"
    pd.DataFrame(
        [
            {
                "module_id": "designed_missing",
                "full_sequence": "ACDEFG",
                "structure_inventory_status": "missing_structure_af3_requested",
                "af3_structure_path": str(tmp_path / "missing.cif"),
                "af3_seed": 42,
                "af3_diffusion_samples": 1,
            }
        ]
    ).to_parquet(inventory, index=False)
    result = validate_af3_outputs(inventory, tmp_path / "validation.parquet")
    assert result.iloc[0].af3_validation_status == "missing_structure"
    assert result.iloc[0].exact_sequence_match == False  # noqa: E712


def test_repeatsdb_alphafold_mapped_feature_schema():
    annotation = {
        "uuid": "af-fixture",
        "content": {
            "chain": {"structure": "A0A010", "id": "A", "source": "AlphaFoldDB"},
            "loci": [],
            "features_uniprot": ["A0A010"],
            "features": {
                "RepeatsDB-5b0j.B": {
                    "chain": {"structure": "5b0j", "id": "B"},
                    "loci": [
                        {"type": "unit", "start": 10, "end": 61, "parent": 0},
                        {"type": "unit", "start": 62, "end": 119, "parent": 0},
                        {"type": "unit", "start": 120, "end": 179, "parent": 0},
                        {"type": "insertion", "start": 180, "end": 223, "parent": 0},
                    ],
                }
            },
        },
    }
    regions = annotation_repeat_regions(annotation)
    assert len(regions) == 1
    assert regions[0]["annotation_schema"] == "mapped_repeatsdb_feature"
    assert regions[0]["region_start"] == 10
    assert regions[0]["region_end"] == 223
    assert regions[0]["unit_count"] == 3


def test_longest_region_selection_ties_are_deterministic():
    rows = pd.DataFrame(
        [
            {"protein_key": "uniprot:X", "region_start": 10, "region_end": 100, "unit_count": 4, "reviewed": False, "structure_source": "AlphaFoldDB", "annotation_uuid": "b", "region_locator": "z"},
            {"protein_key": "uniprot:X", "region_start": 20, "region_end": 120, "unit_count": 3, "reviewed": False, "structure_source": "RCSB/PDB", "annotation_uuid": "c", "region_locator": "z"},
            {"protein_key": "uniprot:X", "region_start": 5, "region_end": 105, "unit_count": 3, "reviewed": True, "structure_source": "RCSB/PDB", "annotation_uuid": "a", "region_locator": "z"},
            {"protein_key": "sha256:Y", "region_start": 1, "region_end": 50, "unit_count": 2, "reviewed": True, "structure_source": "RCSB/PDB", "annotation_uuid": "d", "region_locator": "z"},
        ]
    )
    selected = select_longest_region_per_protein(rows)
    assert len(selected) == 2
    # Span is primary; reviewed status resolves the equal 101-residue spans.
    assert selected.loc[selected.protein_key.eq("uniprot:X"), "annotation_uuid"].item() == "a"


def test_dual_evidence_period_scan_requires_both_peaks():
    aa_unit = "ACDEFGHIKLMN"
    dssp_unit = "CCHHHHTTEECC"
    di_unit = "abcdefghijkl"
    frame = scan_dual_evidence_periods(
        aa_unit * 5, dssp_unit * 5, di_unit * 5, minimum_period=6
    )
    period = frame.loc[frame.period.eq(len(aa_unit))].iloc[0]
    assert period.dual_peak_matched
    assert period.lag_thresholds_passed


THR29_PATH = Path(
    "/net/scratch/wendai/projects/hurdler/clone_repeat_protein/"
    "studies/hurdler_validation/step03_module_corpus/runs/run04_extract/raw/"
    "thr/41586_2024_7188_MOESM4_ESM/THR29_design.pdb"
)


@pytest.mark.skipif(not THR29_PATH.exists(), reason="shared THR29 structure is unavailable")
def test_thr29_golden_period_is_68_not_trivial_self_hit():
    mkdssp = Path("/home/wendai/.conda/envs/hurdler/bin/mkdssp")
    mafft = Path("/home/wendai/.conda/envs/hurdler/bin/mafft")
    chain, aa, dssp, chain_id = _biotite_chain_and_dssp(
        THR29_PATH, chain_id=None, dssp_executable=mkdssp
    )
    foldseek_aa, three_di, _ = foldseek_descriptor(THR29_PATH)
    assert aa == foldseek_aa
    candidate = scan_dual_evidence_periods(aa, dssp, three_di).loc[
        lambda frame: frame.period.eq(68)
    ].iloc[0]
    assert candidate.dssp_state_agreement == pytest.approx(0.994117647, abs=1e-6)
    assert candidate.foldseek_3di_identity == pytest.approx(0.911764706, abs=1e-6)
    result, _, units, _ = infer_designed_boundary(
        {
            "module_id": "designed_THR29",
            "structure_path": str(THR29_PATH),
            "full_sequence": "",
            "evidence_tier": "A",
            "family": "THR",
            "source_accession": "THR29",
        },
        dssp_executable=mkdssp,
        mafft_binary=mafft,
    )
    assert result["strict_dual_evidence_passed"] is True
    assert result["period"] == 68
    assert result["selected_module_index"] == 3
    assert len(units) == 6


def test_stage1_cardinality_bins_and_figure_without_rerunning_lookup(tmp_path, monkeypatch):
    catalog = pd.DataFrame(
        [
            {"module_id": "n1", "collection": "natural_all", "module_type": "Natural", "unit_sequence": "ACDEFG", "unit_length": 6, "family": "x", "evidence_tier": "A"},
            {"module_id": "d1", "collection": "designed_all", "module_type": "Designed", "unit_sequence": "HIKLMNPQRST", "unit_length": 11, "family": "y", "evidence_tier": "A"},
        ]
    )
    catalog_path = tmp_path / "catalog.parquet"
    catalog.to_parquet(catalog_path, index=False)
    monkeypatch.setattr(module_experiments.PatternIndex, "load", lambda path: object())

    def solutions(sequence, index):
        if sequence == "ACDEFG":
            return [
                {
                    "plasmid": "pGEX-4T-1",
                    "site_i_enzyme": "BsaI",
                    "site_ii_enzyme": "BsmBI",
                    "site_i_codon_usage_freq": 1.0,
                    "site_ii_codon_usage_freq": 1.0,
                    "orthogonality": 2.0,
                    "site_i_position": 0,
                    "site_ii_position": 3,
                    "candidate_pair_id": 1,
                    "direction": "right",
                }
            ]
        return []

    monkeypatch.setattr(module_experiments, "enumerate_module_solutions", solutions)
    summary, candidates = run_module_compatibility(
        catalog_path, tmp_path / "index", tmp_path / "shards"
    )
    assert len(summary) == 2
    assert summary.hurdler_compatible.tolist() == [False, True] or summary.hurdler_compatible.tolist() == [True, False]
    assert "all_candidate_solutions_json" not in summary
    assert len(candidates) == 1
    summary_path = next((tmp_path / "shards").glob("module_compatibility_shard-*.parquet"))
    candidate_path = next((tmp_path / "shards").glob("module_compatibility_candidates_shard-*.parquet"))
    per_module, all_candidates, binned = finalize_module_compatibility(
        [summary_path], [candidate_path], tmp_path / "final"
    )
    assert len(per_module) == 2 and len(all_candidates) == 1
    assert (binned.compatible_count + binned.incompatible_count).eq(
        binned.total_count
    ).all()
    assert (tmp_path / "final/module_compatibility_by_length.png").stat().st_size > 0
    assert (tmp_path / "final/module_compatibility_by_length.pdf").stat().st_size > 0


def test_summary_keeps_empty_intervening_bins():
    frame = pd.DataFrame(
        {
            "collection": ["Natural", "Designed"],
            "unit_length": [5, 25],
            "hurdler_compatible": [True, False],
        }
    )
    summary = summarize_compatibility(frame)
    assert summary.length_bin.drop_duplicates().tolist() == ["1–10", "11–20", "21–30"]
    assert len(summary) == 6


def test_stage1_large_candidates_stream_to_parquet_without_csv(tmp_path, monkeypatch):
    summary_path = tmp_path / "summary.parquet"
    candidate_path = tmp_path / "candidates.parquet"
    pd.DataFrame(
        {
            "module_id": ["n1"],
            "collection": ["Natural"],
            "unit_sequence": ["ACDEFG"],
            "unit_length": [6],
            "hurdler_compatible": [True],
        }
    ).to_parquet(summary_path, index=False)
    pd.DataFrame(
        {
            "module_id": ["n1", "n1"],
            "candidate_rank": [1, 2],
            "selected_candidate": [True, False],
        }
    ).to_parquet(candidate_path, index=False)
    monkeypatch.setattr(module_experiments, "_CANDIDATE_PANDAS_ROW_LIMIT", 0)

    _, proxy, _ = finalize_module_compatibility(
        [summary_path], [candidate_path], tmp_path / "final"
    )

    output = tmp_path / "final/module_compatibility_candidates.parquet"
    assert len(pd.read_parquet(output)) == 2
    assert proxy.empty
    assert proxy.attrs["candidate_row_count"] == 2
    assert not output.with_suffix(".csv").exists()


def test_stage2_preparation_is_compact_and_keeps_both_capacities(monkeypatch):
    source = pd.DataFrame(
        [
            {
                "module_id": "natural_fixture",
                "collection": "Natural",
                "family": "fixture",
                "unit_sequence": "ACDEFG",
                "unit_length": 6,
                "hurdler_compatible": True,
                "full_sequence": "ACDEFG" * 100,
                "selected_solution_json": json.dumps(
                    {
                        "plasmid": "pUC18",
                        "site_i_position": 0,
                        "site_ii_position": 3,
                        "site_i_recognition_site": "GAATTC",
                        "site_ii_recognition_site": "GGTCTC",
                    }
                ),
            }
        ]
    )
    monkeypatch.setattr(
        module_experiments,
        "_construct_metrics",
        lambda unit, copies, solution, codons, validate_hard_constraints: {
            "dna_sequence": "GCT" * len(unit) * copies,
            "dna_length": 3 * len(unit) * copies,
        },
    )

    prepared = module_experiments.prepare_adaptive_copy_frame(
        source, codon_weights={}, fragment_limits=(1800, 3000)
    )

    assert prepared.fragment_limit_bp.tolist() == [1800, 3000]
    assert prepared.mathematical_max_copies.tolist() == [100, 166]
    assert "full_sequence" not in prepared
    assert prepared.stage2_preparation_status.eq("prepared").all()


def test_adaptive_finalizer_validates_translation_sites_idt_and_proof(tmp_path):
    dna = "GAATTCGAGTTC"  # EFEF with exactly one selected Site-I occurrence
    response_sha = "a" * 64
    sequence_sha = __import__("hashlib").sha256(dna.encode()).hexdigest()
    result = pd.DataFrame(
        [
            {
                "module_id": "natural_fixture",
                "collection": "Natural",
                "unit_sequence": "EF",
                "unit_length": 2,
                "fragment_limit_bp": 1800,
                "external_deduction_bp": 0,
                "available_coding_bp": 1800,
                "verified_max_copies": 2,
                "mathematical_max_copies": 300,
                "final_passed": True,
                "optimization_status": "passed_adaptive_reduced",
                "adaptive_maximum_proof_status": "next_copy_failed_at_100",
                "adaptive_stop_reason": "copy_3_failed_at_100",
                "adaptive_search_trace_json": json.dumps(
                    [
                        {
                            "phase": "linear_escalation",
                            "copies": 3,
                            "generations": 100,
                            "passed": False,
                        }
                    ]
                ),
                "site_i_recognition_site": "GAATTC",
                "site_ii_recognition_site": "GGTCTC",
                "selected_pair_re_site_excess": 0,
                "dna_sequence": dna,
                "dna_length": len(dna),
                "final_ga_weights_json": "{}",
                "idt_status": "passed",
                "idt_complexity_score": 5.0,
                "idt_score_policy": "idt-rule-score-sum-lt10-v1",
                "idt_positive_score_names_json": '["repeat"]',
                "idt_rule_details_json": "[]",
                "idt_scored_sequence_sha256": sequence_sha,
                "idt_response_sha256": response_sha,
            }
        ]
    )
    cap3000 = result.copy()
    cap3000["fragment_limit_bp"] = 3000
    cap3000["available_coding_bp"] = 3000
    cap3000["mathematical_max_copies"] = 500
    result = pd.concat([result, cap3000], ignore_index=True)
    result_path = tmp_path / "result.parquet"
    result.to_parquet(result_path, index=False)
    compatibility = pd.DataFrame(
        [
            {
                "module_id": "natural_fixture",
                "collection": "Natural",
                "unit_sequence": "EF",
                "unit_length": 2,
                "hurdler_compatible": True,
                "selected_plasmid": "pGEX-4T-1",
                "selected_site_i_enzyme": "EcoRI",
                "selected_site_ii_enzyme": "BsaI",
            },
            {
                "module_id": "designed_incompatible",
                "collection": "Designed",
                "unit_sequence": "ACDEFG",
                "unit_length": 6,
                "hurdler_compatible": False,
                "selected_plasmid": None,
                "selected_site_i_enzyme": None,
                "selected_site_ii_enzyme": None,
            },
        ]
    )
    compatibility_path = tmp_path / "compatibility.parquet"
    compatibility.to_parquet(compatibility_path, index=False)
    audit = tmp_path / "audit.jsonl"
    audit.write_text(
        json.dumps(
            {
                "response_sha256": response_sha,
                "summary": {"idt_response_sha256": response_sha},
                "response": [],
            }
        )
        + "\n"
    )
    finalized, trace, summary = finalize_adaptive_copy_results(
        [result_path],
        compatibility_path,
        tmp_path / "final",
        idt_audit_paths=[audit],
    )
    assert finalized.validation_passed.all()
    assert trace.generations.eq(100).all()
    assert len(finalized) == 2
    assert len(summary) == 2
    assert (tmp_path / "final/optimized_constructs.fasta").stat().st_size > 0
    assert (tmp_path / "final/maximum_verified_repeat_copies.png").stat().st_size > 0


def test_adaptive_result_requires_trace_level_maximum_proof():
    from hurdler.module_experiments import _validate_adaptive_result

    dna = "GAATTCGAGTTC"
    response_sha = "b" * 64
    base = {
        "unit_sequence": "EF",
        "verified_max_copies": 2,
        "final_passed": True,
        "available_coding_bp": 1800,
        "fragment_limit_bp": 1800,
        "external_deduction_bp": 0,
        "dna_sequence": dna,
        "site_i_recognition_site": "GAATTC",
        "site_ii_recognition_site": "GGTCTC",
        "selected_pair_re_site_excess": 0,
        "idt_complexity_score": 1.0,
        "idt_score_policy": "idt-rule-score-sum-lt10-v1",
        "idt_scored_sequence_sha256": __import__("hashlib").sha256(dna.encode()).hexdigest(),
        "idt_response_sha256": response_sha,
        "adaptive_maximum_proof_status": "next_copy_failed_at_100",
        "adaptive_search_trace_json": json.dumps(
            [{"copies": 3, "generations": 80, "passed": False}]
        ),
    }
    accepted, reasons = _validate_adaptive_result(base, {response_sha})
    assert accepted
    assert "next_copy_100_generation_failure_missing" in reasons

    capacity = {
        **base,
        "adaptive_maximum_proof_status": "capacity_limit_reached",
        "adaptive_search_upper_bound_copies": 3,
        "adaptive_search_trace_json": "[]",
    }
    accepted, reasons = _validate_adaptive_result(capacity, {response_sha})
    assert accepted
    assert "capacity_limit_proof_mismatch" in reasons


def test_adaptive_finalizer_rejects_a_missing_capacity(tmp_path):
    result = pd.DataFrame(
        [
            {
                "module_id": "natural_fixture",
                "collection": "Natural",
                "unit_sequence": "EF",
                "unit_length": 2,
                "fragment_limit_bp": 1800,
                "verified_max_copies": 0,
                "adaptive_search_trace_json": "[]",
            }
        ]
    )
    compatibility = pd.DataFrame(
        [
            {
                "module_id": "natural_fixture",
                "collection": "Natural",
                "unit_sequence": "EF",
                "hurdler_compatible": True,
            }
        ]
    )
    result_path = tmp_path / "result.parquet"
    compatibility_path = tmp_path / "compatibility.parquet"
    result.to_parquet(result_path, index=False)
    compatibility.to_parquet(compatibility_path, index=False)
    with pytest.raises(ValueError, match="exactly the 1800-bp and 3000-bp"):
        finalize_adaptive_copy_results(
            [result_path],
            compatibility_path,
            tmp_path / "final",
            idt_audit_paths=[],
        )
