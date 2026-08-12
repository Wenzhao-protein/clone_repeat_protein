from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from hurdler.design import (
    DesignRequest,
    analyze_repeat_sequence,
    boundary_confirmation_token,
    bundled_index_dir,
    confirm_repeat_boundary,
    design_construct,
    enumerate_design_candidates,
    plan_purchase_fragments,
    role_enzyme_options,
    simulate_vector_assembly,
    write_design_outputs,
)
from hurdler.constants import PLASMIDS
from hurdler.index import PatternIndex
from hurdler.notebook_ui import PassingMockIDTScorer, RPB1_CTD_FASTA
from hurdler.optimization import recognition_site_count, translate_dna


MODULE = "HIKLMNPQRST"
FULL_PROTEIN = "M" + MODULE * 4 + "K"
REGION_START = 2
REGION_END = 1 + len(MODULE) * 4


@pytest.fixture(scope="module")
def full_index() -> PatternIndex:
    return PatternIndex.load(bundled_index_dir())


def _request(**updates) -> DesignRequest:
    values = {
        "sequence_id": "fixture",
        "full_protein_sequence": FULL_PROTEIN,
        "target_repeat_copies": 4,
        "plasmid": "pGEX-4T-1",
        "confirmed_repeat_start": REGION_START,
        "confirmed_repeat_end": REGION_END,
        "confirmed_period": len(MODULE),
        "confirmation_token": boundary_confirmation_token(
            FULL_PROTEIN, REGION_START, REGION_END, len(MODULE)
        ),
        "optimize": False,
        "population_size": 4,
    }
    values.update(updates)
    return DesignRequest(**values)


def test_sequence_parser_infers_exact_one_and_two_residue_repeats():
    one = analyze_repeat_sequence("M" + "A" * 8 + "K")
    assert one.proposed_period == 1
    assert (one.proposed_start, one.proposed_end) == (2, 9)

    two = analyze_repeat_sequence("M" + "AC" * 5 + "K")
    assert two.proposed_period == 2
    assert (two.proposed_start, two.proposed_end) == (2, 11)


def test_sequence_parser_handles_variable_units_and_prefers_primitive_harmonic():
    variable = analyze_repeat_sequence(
        "M"
        + "ACDEFGHIKL"
        + "ACDEFGHIKM"
        + "ACDEFGHIKN"
        + "ACDEFGHIKQ"
        + "K"
    )
    assert variable.proposed_period == 10
    assert variable.candidates[0].repeat_count == 4

    harmonic = analyze_repeat_sequence("M" + "ACDEACDE" * 4 + "K")
    assert harmonic.proposed_period == 4
    assert harmonic.candidates[0].repeat_count == 8


def test_s_cerevisiae_rpb1_ctd_golden_input_uses_the_earlier_middle_heptad(full_index):
    analysis = analyze_repeat_sequence(RPB1_CTD_FASTA)
    assert analysis.sequence_id == "S_cerevisiae_Rpb1_CTD_26_repeats_plus_C_terminal_tip"
    assert len(analysis.full_protein_sequence) == 192
    assert analysis.proposed_period == 7
    # The fast inference proposes the longest exact core. The user-confirmed
    # biological CTD annotation spans 26 heptads (1..182), including variable
    # terminal repeats, followed by the 10-AA C-terminal tip.
    boundary = confirm_repeat_boundary(
        analysis.full_protein_sequence,
        start=1,
        end=182,
        period=7,
    )
    assert boundary.repeat_count == 26
    assert boundary.middle_unit_index == 13
    assert boundary.middle_unit_start == 85
    assert boundary.middle_module == "YSPTSPS"
    assert boundary.c_terminal_flank == "QKHNENENSR"
    assert boundary.fixed_positions_1based == (1, 2, 3, 4, 5, 6)

    result = design_construct(
        DesignRequest(
            sequence_id=analysis.sequence_id,
            full_protein_sequence=analysis.full_protein_sequence,
            target_repeat_copies=26,
            plasmid="pGEX-4T-1",
            confirmed_repeat_start=1,
            confirmed_repeat_end=182,
            confirmed_period=7,
            confirmation_token=boundary_confirmation_token(
                analysis.full_protein_sequence, 1, 182, 7
            ),
        ),
        index=full_index,
    )
    assert result.status == "hurdler_incompatible"
    assert result.target_protein_sequence == "YSPTSPS" * 26 + "QKHNENENSR"
    assert result.candidates == []
    assert all(
        not enumerate_design_candidates("YSPTSPS", plasmid, full_index)
        for plasmid in PLASMIDS
    )


def test_sequence_parser_accepts_one_fasta_record_and_rejects_multiple():
    result = analyze_repeat_sequence(">example description\n" + MODULE * 3)
    assert result.sequence_id == "example"
    with pytest.raises(ValueError, match="Exactly one"):
        analyze_repeat_sequence(">one\nAAAA\n>two\nAAAA")


def test_confirmation_selects_earlier_middle_variable_unit_and_preserves_flanks():
    units = ("ACDE", "ACDF", "ACDG", "ACDH")
    sequence = "MNP" + "".join(units) + "QRS"
    boundary = confirm_repeat_boundary(
        sequence,
        start=4,
        end=19,
        period=4,
    )
    assert boundary.middle_unit_index == 2
    assert boundary.middle_module == "ACDF"
    assert boundary.n_terminal_flank == "MNP"
    assert boundary.c_terminal_flank == "QRS"
    assert boundary.consensus_module == "ACDE"
    assert boundary.fixed_positions_1based == (1, 2, 3)
    assert boundary.variable_ranges_1based == ((4, 4),)


def test_stale_confirmation_is_rejected():
    token = boundary_confirmation_token(FULL_PROTEIN, REGION_START, REGION_END, len(MODULE))
    with pytest.raises(ValueError, match="stale"):
        confirm_repeat_boundary(
            FULL_PROTEIN + "A",
            start=REGION_START,
            end=REGION_END,
            period=len(MODULE),
            expected_token=token,
        )


def test_request_schema_rejects_legacy_or_unknown_fields():
    with pytest.raises(ValueError, match="Unknown DesignRequest"):
        DesignRequest.from_dict(
            {
                "full_protein_sequence": FULL_PROTEIN,
                "target_repeat_copies": 4,
                "plasmid": "pGEX-4T-1",
                "df2_plasmid": "legacy",
            }
        )


def test_all_role_defaults_and_role_mismatch_error(full_index):
    options = role_enzyme_options(full_index)
    assert options["site_i"] and options["site_ii"] and options["site_iii"]
    candidates = enumerate_design_candidates(MODULE, "pGEX-4T-1", full_index)
    assert candidates
    assert any(int(row["site_ii_position"]) >= len(MODULE) for row in candidates)
    with pytest.raises(ValueError, match="Invalid Site I"):
        enumerate_design_candidates(
            MODULE,
            "pGEX-4T-1",
            full_index,
            site_i_allowlist=("BsaI",),
        )


def test_exact_one_and_two_aa_modules_reach_the_frozen_index(full_index):
    one = enumerate_design_candidates("A", "pGEX-4T-1", full_index)
    two = enumerate_design_candidates("AA", "pGEX-4T-1", full_index)
    assert one and two
    assert all(row["scan_copy_count"] == 2 for row in [one[0], two[0]])


def test_allowlist_filters_each_role_and_keeps_deterministic_ranking(full_index):
    all_rows = enumerate_design_candidates(MODULE, "pGEX-4T-1", full_index)
    chosen = all_rows[0]
    filtered = enumerate_design_candidates(
        MODULE,
        "pGEX-4T-1",
        full_index,
        site_i_allowlist=(chosen["site_i_enzyme"],),
        site_ii_allowlist=(chosen["site_ii_enzyme"],),
        site_iii_allowlist=(chosen["site_iii_enzyme"],),
    )
    assert filtered
    assert all(row["site_i_enzyme"] == chosen["site_i_enzyme"] for row in filtered)
    assert all(row["site_ii_enzyme"] == chosen["site_ii_enzyme"] for row in filtered)
    assert all(row["site_iii_enzyme"] == chosen["site_iii_enzyme"] for row in filtered)
    assert [row["rank"] for row in filtered] == list(range(1, len(filtered) + 1))


def test_missing_confirmation_stops_before_hurdler(full_index):
    request = DesignRequest(
        full_protein_sequence=FULL_PROTEIN,
        target_repeat_copies=4,
        plasmid="pGEX-4T-1",
    )
    result = design_construct(request, index=full_index)
    assert result.status == "needs_boundary_confirmation"
    assert result.candidates == []


def test_no_optimization_is_explicitly_not_orderable(full_index, tmp_path):
    result = design_construct(_request(), index=full_index)
    assert result.status == "not_orderable_not_for_purchase"
    assert result.target_protein_sequence == "M" + MODULE * 4 + "K"
    files = write_design_outputs(result, tmp_path)
    assert "fragment_topology_csv" in files
    assert "optimized_construct_fasta" not in files
    assert not (tmp_path / "purchase_fragments.fasta").exists()


def test_mock_optimized_design_preserves_translation_sites_and_assembly(full_index, tmp_path):
    result = design_construct(
        _request(optimize=True),
        index=full_index,
        idt_scorer=PassingMockIDTScorer(),
    )
    assert result.orderable
    dna = result.optimization["optimized_dna"]
    assert translate_dna(dna) == result.target_protein_sequence
    selected = result.selected_candidate
    assert recognition_site_count(dna, selected["site_i_recognition_site"]) == 1
    assert recognition_site_count(dna, selected["site_ii_recognition_site"]) == 0
    assert all(row["idt_explicit_pass"] for row in result.purchase_fragments)
    assert all(
        row["purchase_sha256"] == hashlib.sha256(row["purchase_sequence"].encode()).hexdigest()
        == row["idt_scored_sequence_sha256"]
        for row in result.purchase_fragments
    )
    assert result.final_plasmid["translation_exact"] is True
    assert result.final_plasmid["diagnostic_digest_fragments_bp"]
    assert "expected_fragment_sizes_bp_json" in result.cloning_steps[-1]
    files = write_design_outputs(result, tmp_path)
    assert Path(files["optimized_construct_fasta"]).stat().st_size > 0
    assert Path(files["purchase_fragments_fasta"]).stat().st_size > 0
    assert Path(files["final_plasmid_fasta"]).stat().st_size > 0


def test_exact_score_ten_reweights_ga_before_retry(full_index):
    class RejectOnce:
        calls = 0

        def score(self, name, sequence):
            self.calls += 1
            sha = hashlib.sha256(sequence.encode()).hexdigest()
            if self.calls == 1:
                return {
                    "idt_status": "failed",
                    # The design layer independently enforces strict <10 even
                    # if a malformed client claims the exact-threshold row passed.
                    "idt_explicit_pass": True,
                    "idt_complexity_score": 10.0,
                    "idt_rule_details_json": json.dumps(
                        [{"name": "High GC content", "score": 10.0, "actual_value": 0.8}]
                    ),
                    "idt_positive_score_names_json": '["High GC content"]',
                    "idt_violation_names_json": '["High GC content"]',
                    "idt_scored_sequence_sha256": sha,
                    "idt_response_sha256": "first",
                }
            return {
                "idt_status": "passed",
                "idt_explicit_pass": True,
                "idt_complexity_score": 0.0,
                "idt_rule_details_json": "[]",
                "idt_positive_score_names_json": "[]",
                "idt_violation_names_json": "[]",
                "idt_scored_sequence_sha256": sha,
                "idt_response_sha256": "later",
            }

    result = design_construct(_request(optimize=True), index=full_index, idt_scorer=RejectOnce())
    assert result.orderable
    attempts = result.optimization["attempts"]
    assert len(attempts) >= 2
    first_after = json.loads(attempts[0]["ga_weights_after_json"])
    first_before = json.loads(attempts[0]["ga_weights_before_json"])
    assert first_after["gc_window_soft_violation"] == 2 * first_before["gc_window_soft_violation"]
    assert attempts[1]["all_fragments_idt_passed"] is True


def test_api_failure_is_a_nonpass_not_a_credential_leaking_exception(full_index):
    candidate = enumerate_design_candidates(MODULE, "pGEX-4T-1", full_index)[0]

    class BrokenScorer:
        def score(self, name, sequence):
            raise RuntimeError("upstream unavailable")

    result = design_construct(
        _request(
            optimize=True,
            selected_candidate_id=candidate["candidate_id"],
        ),
        index=full_index,
        idt_scorer=BrokenScorer(),
    )
    assert result.status == "optimization_failed_not_orderable"
    assert len(result.optimization["attempts"]) == 6
    assert all(not row["all_fragments_idt_passed"] for row in result.optimization["attempts"])
    assert all(row["idt_status"] == "api_failure" for row in result.idt_audit)
    assert "upstream unavailable" not in json.dumps(result.to_dict())


def test_long_construct_is_split_and_reconstructable():
    dna = "GCT" * 1400
    fragments = plan_purchase_fragments("long", dna, "pUC18")
    assert len(fragments) >= 2
    assert all(row["purchase_length_bp"] <= 3000 for row in fragments)
    reconstructed = ""
    for row in fragments:
        core = dna[row["core_start_0based"] : row["core_end_exclusive"]]
        reconstructed += core[row["overlap_with_previous_bp"] :]
    assert reconstructed == dna


@pytest.mark.parametrize(
    "plasmid",
    [
        "pGEX-4T-1",
        "pMAL-c5X",
        "pET-21a(+)",
        "pET-28a(+)",
        "pET-28a(+)_start_codon",
        "pCold_I",
        "pUC18",
        "pQE-3",
    ],
)
def test_all_maintained_plasmid_mcs_simulations_preserve_the_cds(plasmid):
    result = simulate_vector_assembly(plasmid, "ATG" + "GCT" * 4, "MAAAA")
    assert result["assembly_sequence_exact"] is True
    assert result["translation_exact"] is True
    assert result["final_plasmid_length_bp"] > 1000


def test_output_manifest_and_artifacts_contain_no_credential_material(full_index, tmp_path):
    result = design_construct(_request(), index=full_index)
    write_design_outputs(result, tmp_path)
    text = "\n".join(path.read_text(errors="ignore") for path in tmp_path.iterdir() if path.is_file())
    assert "IDT_CLIENT_SECRET" not in text
    assert "IDT_PASSWORD" not in text
    assert ".config/hurdler/idt.env" not in text
    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    assert manifest["credential_material_persisted"] is False
