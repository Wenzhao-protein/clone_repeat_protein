from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import pytest
import pandas as pd

from hurdler.dna_assembly import (
    enumerate_active_latent_pairs,
    load_enzyme_catalog,
    scan_re_sites,
)
from hurdler.cli import main
from hurdler.exact_dna_design import (
    EXACT_DNA_SCHEMA_VERSION,
    ExactDNAQuery,
    ExactDNASelection,
    confirm_exact_dna_route,
    load_exact_dna_enzyme_catalog,
    parse_exact_dna_input,
    query_exact_dna,
    write_exact_dna_outputs,
)
from hurdler.paths import ProjectPaths


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = json.loads(
    (ROOT / "data/examples/rf00059_tpp_riboswitch_array.json").read_text()
)


def _rf00059_query(**overrides) -> ExactDNAQuery:
    values = {
        "schema_version": EXACT_DNA_SCHEMA_VERSION,
        "input_mode": "array",
        "sequence_id": EXAMPLE["example_id"],
        "repeat_unit": EXAMPLE["repeat_unit"],
        "repeat_copies": 4,
        "site_i_allowlist": ("AflII",),
        "site_ii_allowlist": ("ApaI",),
        "plasmid_allowlist": ("pQE-3", "pUC18"),
        "timeout_seconds": 30,
    }
    values.update(overrides)
    return ExactDNAQuery(**values)


@pytest.fixture(scope="module")
def rf00059_result():
    return query_exact_dna(_rf00059_query())


class FakeScorer:
    def __init__(self, score: float = 0.0, *, malformed: bool = False):
        self.score_value = score
        self.malformed = malformed
        self.calls: list[tuple[str, str]] = []

    def score(self, name: str, sequence: str):
        self.calls.append((name, sequence))
        if self.malformed:
            return {
                "idt_complexity_score": float("nan"),
                "idt_score_complete": False,
                "idt_explicit_pass": False,
                "idt_invalid_score_names_json": '["Repeat"]',
            }
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        return {
            "idt_complexity_score": self.score_value,
            "idt_score_complete": True,
            "idt_explicit_pass": self.score_value < 10,
            "idt_status": "pass" if self.score_value < 10 else "fail",
            "idt_response_sha256": hashlib.sha256(("response|" + digest).encode()).hexdigest(),
            "idt_scored_sequence_sha256": digest,
            "idt_positive_score_names_json": "[]" if self.score_value == 0 else '["Repeat"]',
            "idt_rule_details_json": "[]",
            "idt_invalid_score_names_json": "[]",
        }


class BombScorer:
    def score(self, _name: str, _sequence: str):
        raise AssertionError("Bulk mode must not invoke an IDT scorer")


class DiagnosticFailureThenPass(FakeScorer):
    def score(self, name: str, sequence: str):
        if name.endswith("|whole_target_diagnostic"):
            raise RuntimeError("diagnostic-only failure")
        return super().score(name, sequence)


def test_exact_dna_parser_preserves_one_fasta_and_rejects_non_acgt():
    sequence_id, sequence = parse_exact_dna_input(">member details\nac gt\n")
    assert sequence_id == "member"
    assert sequence == "ACGT"
    with pytest.raises(ValueError, match="exactly one FASTA"):
        parse_exact_dna_input(">one\nACGT\n>two\nACGT\n")
    for invalid in ("ACGU", "ACGN"):
        with pytest.raises(ValueError, match="A/C/G/T only"):
            parse_exact_dna_input(invalid)

    array = ExactDNAQuery(
        schema_version=EXACT_DNA_SCHEMA_VERSION,
        input_mode="array",
        repeat_unit="ACGT",
        spacer="TT",
        repeat_copies=3,
    )
    assert array.target_sequence == "ACGTTTACGTTTACGT"
    assert not array.target_sequence.endswith(array.spacer)


def test_exact_catalog_combines_broad_re_roles_with_type_iis_adapters():
    geometries = load_exact_dna_enzyme_catalog()
    assert geometries["ApaI"].site_i_eligible
    assert geometries["ApaI"].site_ii_eligible
    assert any(item.is_type_iis and item.site_iii_eligible for item in geometries.values())


def test_active_latent_enumeration_excludes_active_active_and_two_mismatches():
    paths = ProjectPaths.discover()
    geometries, _legacy = load_enzyme_catalog(
        paths.reference_output, artifact_dir=paths.root / "data/artifacts"
    )
    pair_geometries = {name: geometries[name] for name in ("EcoRI", "HindIII")}
    active = "ACGT" * 8 + "GAATTC" + "ATGCGT" * 20 + "AAGCTT" + "TGCA" * 8
    hits = scan_re_sites("active", active, pair_geometries)
    assert not enumerate_active_latent_pairs(hits, None)

    latent = active.replace("AAGCTT", "GAGCTT")
    schemes = enumerate_active_latent_pairs(
        scan_re_sites("active_latent", latent, pair_geometries), None
    )
    assert schemes
    assert {row["mechanism"] for row in schemes} == {"active+latent"}
    reverse_pair_geometries = {
        name: geometries[name] for name in ("EcoRI", "AciI")
    }
    reverse_target = "A" * 32 + "GCGA" + "A" * 120 + "GAATTC" + "T" * 32
    reverse_schemes = enumerate_active_latent_pairs(
        scan_re_sites("reverse", reverse_target, reverse_pair_geometries), None
    )
    assert reverse_schemes
    assert any(
        row["site_i"].orientation == "-" or row["site_ii"].orientation == "-"
        for row in reverse_schemes
    )

    hind_only = {"HindIII": geometries["HindIII"]}
    assert not scan_re_sites("two_mismatch", "GAGCTC", hind_only)


def test_rf00059_golden_complete_route_is_exact_and_latent_latent(rf00059_result):
    result = rf00059_result
    assert result.status == "hurdler_compatible_molecular"
    assert result.target_length_bp == 432
    assert result.target_sequence == EXAMPLE["repeat_unit"] * 4
    assert result.target_sequence_sha256 == EXAMPLE["target_sha256"]
    assert {row["profile_id"] for row in result.route_candidates} == {"pQE-3", "pUC18"}
    assert {
        (pair["site_i_enzyme"], pair["site_ii_enzyme"])
        for route in result.route_candidates
        for pair in route["pairs"]
    } == {("AflII", "ApaI")}

    confirmed = confirm_exact_dna_route(
        result, ExactDNASelection(result.route_candidates[0]["route_id"], "none")
    )
    assert confirmed.seed["seed_sequence"] == EXAMPLE["repeat_unit"]
    assert confirmed.final_insert_sequence == result.target_sequence
    assert confirmed.final_plasmid["insert_exact"] is True
    assert confirmed.final_plasmid["protected_feature_sequences_preserved"] is True
    assert set(confirmed.final_plasmid["selected_pair_excess_sites"].values()) == {0}
    assert {row["mechanism"] for row in confirmed.transitions} == {"latent+latent"}
    assert confirmed.latent_transitions
    assert any(row["donor_derived"] for row in confirmed.latent_transitions)
    assert all(row["final_target_restored"] for row in confirmed.latent_transitions)
    assert {
        (row["enzyme"], row["orientation"])
        for row in confirmed.latent_transitions
    } >= {("AflII", "+"), ("ApaI", "+")}
    assert all(row["target_base"] for row in confirmed.latent_transitions)
    assert all(row["temporary_active_base"] for row in confirmed.latent_transitions)
    with pytest.raises(ValueError, match="plasmid_profile"):
        confirm_exact_dna_route(
            result,
            ExactDNASelection(
                result.route_candidates[0]["route_id"],
                plasmid_profile="pGEX-4T-1",
            ),
        )


def test_rf00059_regresses_complete_route_v2_production(rf00059_result):
    table_path = ROOT / (
        "studies/hurdler_validation/step06_repetitive_dna_assembly/tables/"
        "arbitrary-dna-complete-route-v2/production/production_target_analysis.parquet"
    )
    if not table_path.is_file():
        pytest.skip("Production regression artifact is not present in this checkout")
    rows = pd.read_parquet(
        table_path,
        filters=[("unit_sequence", "==", EXAMPLE["repeat_unit"]), ("target_copy_count", "==", 4)],
    )
    assert len(rows) == 1
    legacy = rows.iloc[0]
    assert bool(legacy.complete_route_verified)
    assert int(legacy.seed_copy_count) == 1
    assert legacy.plasmid == "pQE-3"
    assert int(legacy.hurdler_step_count) == 6
    assert int(legacy.pair_change_count) == 0

    matching = [
        row
        for row in rf00059_result.route_candidates
        if row["profile_id"] == "pQE-3"
        and row["hurdler_step_count"] == int(legacy.hurdler_step_count)
        and row["pair_change_count"] == int(legacy.pair_change_count)
        and row["seed"]["seed_length_bp"] == len(EXAMPLE["repeat_unit"])
    ]
    assert matching


def test_arbitrary_exact_dna_requires_a_shorter_seed_to_target_route():
    target = "ACGT" * 8 + "GAATTC" + "ATGCGT" * 20 + "GAGCTT" + "TGCA" * 8
    query = ExactDNAQuery(
        schema_version=EXACT_DNA_SCHEMA_VERSION,
        input_mode="exact",
        exact_dna=f">arbitrary_element\n{target}\n",
        site_i_allowlist=("EcoRI",),
        site_ii_allowlist=("HindIII",),
        plasmid_allowlist=("pUC18",),
        timeout_seconds=30,
        max_states=100,
    )
    result = query_exact_dna(query)
    assert result.status == "hurdler_compatible_molecular"
    confirmed = confirm_exact_dna_route(
        result, ExactDNASelection(result.route_candidates[0]["route_id"], "batch")
    )
    assert confirmed.status == "bulk_export_unvalidated"
    assert confirmed.seed["seed_length_bp"] < len(target)
    assert confirmed.transitions
    assert confirmed.final_insert_sequence == target
    assert confirmed.target_sequence_sha256 == hashlib.sha256(target.encode()).hexdigest()


def test_budget_exhaustion_is_unclassified_not_incompatible():
    result = query_exact_dna(_rf00059_query(max_states=1))
    assert result.status == "search_incomplete"
    assert result.termination_reason == "search_budget_exhausted"


def test_batch_mode_makes_no_idt_call_and_exports_exact_files(rf00059_result, tmp_path):
    confirmed = confirm_exact_dna_route(
        rf00059_result,
        ExactDNASelection(rf00059_result.route_candidates[0]["route_id"], "batch"),
        idt_scorer=BombScorer(),
    )
    assert confirmed.status == "bulk_export_unvalidated"
    files = write_exact_dna_outputs(confirmed, tmp_path)
    assert Path(files["target_fasta"]).read_text().splitlines()[-1] == confirmed.target_sequence
    assert (tmp_path / "idt_bulk_input.csv").is_file()
    bulk = pd.read_csv(tmp_path / "idt_bulk_input.csv")
    assert bulk.Sequence.is_unique
    assert not (tmp_path / "idt_raw_audit.jsonl").exists()
    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    assert manifest["credentials_persisted"] is False
    assert manifest["ordering_performed"] is False


def test_live_idt_strict_threshold_invalid_score_and_api_error(rf00059_result):
    # RF00059 uses short primer-pair donors. They remain a valid molecular
    # route but cannot be mislabeled as IDT accepted.
    passing = FakeScorer(0)
    unscored = confirm_exact_dna_route(
        rf00059_result,
        ExactDNASelection(rf00059_result.route_candidates[0]["route_id"], "api"),
        idt_scorer=passing,
    )
    assert unscored.status == "idt_unscored_primer_route"
    assert passing.calls

    # This fixture has a 150-bp exact seed and one 514-bp gBlock donor, so
    # every actual purchase is eligible for the score-sum decision.
    target = "A" * 70 + "GAATTC" + "ATGCGT" * 80 + "GAGCTT" + "T" * 70
    scored_query = ExactDNAQuery(
        schema_version=EXACT_DNA_SCHEMA_VERSION,
        input_mode="exact",
        exact_dna=target,
        site_i_allowlist=("EcoRI",),
        site_ii_allowlist=("HindIII",),
        plasmid_allowlist=("pUC18",),
        timeout_seconds=30,
        max_states=100,
    )
    scored_result = query_exact_dna(scored_query)
    route_id = scored_result.route_candidates[0]["route_id"]
    accepted = confirm_exact_dna_route(
        scored_result, ExactDNASelection(route_id, "api"), idt_scorer=FakeScorer(0)
    )
    assert accepted.status == "idt_accepted_route"
    assert all(row["score"] < 10 for row in accepted.idt_audit)

    boundary = confirm_exact_dna_route(
        scored_result, ExactDNASelection(route_id, "api"), idt_scorer=FakeScorer(10)
    )
    assert boundary.status == "idt_rejected_route"
    assert any(row["score"] == 10 for row in boundary.idt_audit)

    malformed = confirm_exact_dna_route(
        scored_result,
        ExactDNASelection(route_id, "api"),
        idt_scorer=FakeScorer(malformed=True),
    )
    assert malformed.status == "idt_score_error"

    failed = confirm_exact_dna_route(
        scored_result, ExactDNASelection(route_id, "api"), idt_scorer=BombScorer()
    )
    assert failed.status == "idt_api_error"

    diagnostic_failure = confirm_exact_dna_route(
        scored_result,
        ExactDNASelection(route_id, "api"),
        idt_scorer=DiagnosticFailureThenPass(0),
    )
    assert diagnostic_failure.status == "idt_accepted_route"
    assert diagnostic_failure.whole_target_idt["status"] == "diagnostic_api_error"


def test_example_hashes_are_self_consistent():
    unit = EXAMPLE["repeat_unit"]
    target = unit * EXAMPLE["repeat_copies"]
    assert len(unit) == EXAMPLE["repeat_unit_length_bp"] == 108
    assert len(target) == EXAMPLE["target_length_bp"] == 432
    assert hashlib.sha256(unit.encode()).hexdigest() == EXAMPLE["repeat_unit_sha256"]
    assert hashlib.sha256(target.encode()).hexdigest() == EXAMPLE["target_sha256"]


def test_exact_dna_cli_confirms_batch_route_without_http(rf00059_result, tmp_path):
    query = _rf00059_query()
    request = tmp_path / "request.json"
    request.write_text(
        json.dumps(
            {
                "query": asdict(query),
                "selection": {
                    "route_id": rf00059_result.route_candidates[0]["route_id"],
                    "validation_mode": "batch",
                },
            }
        )
    )
    output = tmp_path / "output"
    assert main(
        [
            "dna-assembly",
            "interactive-design",
            "--request",
            str(request),
            "--output-dir",
            str(output),
        ]
    ) == 0
    summary = json.loads((output / "exact_dna_design_summary.json").read_text())
    assert summary["status"] == "bulk_export_unvalidated"
    assert summary["final_insert_sequence"] == query.target_sequence
    assert (output / "idt_bulk_input.csv").is_file()
    assert not (output / "idt_raw_audit.jsonl").exists()
