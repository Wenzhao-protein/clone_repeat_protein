from __future__ import annotations

import hashlib
import json
import copy
import tempfile
from dataclasses import asdict
from pathlib import Path

import pytest
import pandas as pd
from Bio import SeqIO

from hurdler.constants import PLASMIDS
from hurdler.dna_assembly import (
    enumerate_active_latent_pairs,
    load_enzyme_catalog,
    scan_re_sites,
)
from hurdler.cli import main
from hurdler.exact_dna_design import (
    EXACT_DNA_SCHEMA_VERSION,
    IDT_GBLOCK_ONLY_PURCHASE_POLICY,
    ExactDNAQuery,
    ExactDNASelection,
    confirm_best_exact_dna_route,
    confirm_exact_dna_route,
    load_exact_dna_enzyme_catalog,
    parse_exact_dna_input,
    query_exact_dna,
    write_exact_dna_outputs,
    write_exact_dna_minimal_outputs,
    _retain_annotated_route_groups,
)
from hurdler.exact_dna_verification import verify_exact_dna_assembly
from hurdler.plasmid_reference import load_plasmid_reference
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


@pytest.fixture(scope="module")
def rf00059_gblock_result():
    geometries = load_exact_dna_enzyme_catalog()
    return query_exact_dna(
        ExactDNAQuery(
            schema_version=EXACT_DNA_SCHEMA_VERSION,
            input_mode="array",
            sequence_id=EXAMPLE["example_id"],
            repeat_unit=EXAMPLE["repeat_unit"],
            repeat_copies=4,
            site_i_allowlist=tuple(
                sorted(name for name, item in geometries.items() if item.site_i_eligible)
            ),
            site_ii_allowlist=tuple(
                sorted(name for name, item in geometries.items() if item.site_ii_eligible)
            ),
            site_iii_allowlist=tuple(
                sorted(name for name, item in geometries.items() if item.site_iii_eligible)
            ),
            plasmid_allowlist=tuple(sorted(PLASMIDS)),
            purchase_policy=IDT_GBLOCK_ONLY_PURCHASE_POLICY,
            timeout_seconds=60,
        )
    )


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


class RejectUnsplitRF00059Donor(FakeScorer):
    def score(self, name: str, sequence: str):
        self.calls.append((name, sequence))
        score = 200.5 if len(sequence) > 300 else 0.0
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        return {
            "idt_complexity_score": score,
            "idt_score_complete": True,
            "idt_explicit_pass": score < 10,
            "idt_status": "pass" if score < 10 else "fail",
            "idt_response_sha256": hashlib.sha256(("response|" + digest).encode()).hexdigest(),
            "idt_scored_sequence_sha256": digest,
            "idt_positive_score_names_json": "[]" if score == 0 else '["Repeat"]',
            "idt_rule_details_json": "[]",
            "idt_invalid_score_names_json": "[]",
        }


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


def test_route_group_retention_does_not_let_one_pair_hide_other_pairs():
    routes = []
    for pair_index in range(4):
        for rank in range(30 if pair_index == 0 else 1):
            routes.append(
                {
                    "route_id": f"r{pair_index}_{rank}",
                    "profile_id": "pUC18",
                    "scheme_id": "pUC18:inside/inside",
                    "pair_mode": "fixed",
                    "pairs": [
                        {
                            "site_i_enzyme": f"I{pair_index}",
                            "site_ii_enzyme": f"II{pair_index}",
                        }
                    ],
                }
            )
    retained = _retain_annotated_route_groups(routes, routes_per_group=1)
    assert len(retained) == 4
    assert {row["pairs"][0]["site_i_enzyme"] for row in retained} == {
        "I0", "I1", "I2", "I3"
    }


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
    assert Path(files["final_insert_fasta"]).read_text().splitlines()[-1] == confirmed.target_sequence
    assert (tmp_path / "idt_bulk_input.csv").is_file()
    bulk = pd.read_csv(tmp_path / "idt_bulk_input.csv")
    assert bulk.Sequence.is_unique
    assert not (tmp_path / "idt_raw_audit.jsonl").exists()
    assert not (tmp_path / "run_manifest.json").exists()
    assert Path(files["technical_audit_zip"]).is_file()
    assert (tmp_path / "cloning_summary.csv").is_file()
    assert len(confirmed.cloning_steps) == 7
    assert (tmp_path / "step00_plasmid.gb").is_file()
    assert (tmp_path / "step01_insert.gb").is_file()
    assert (tmp_path / "step07_plasmid.gb").is_file()
    for path in sorted(tmp_path.glob("step*.gb")):
        record = SeqIO.read(path, "genbank")
        assert len(record) > 0
        assert record.features
        assert all(
            "hash" not in key.lower() and "sha" not in key.lower()
            for feature in record.features
            for key in feature.qualifiers
        )


def test_independent_verifier_rejects_a_tampered_purchase(rf00059_result):
    route_id = rf00059_result.route_candidates[0]["route_id"]
    confirmed = confirm_exact_dna_route(
        rf00059_result, ExactDNASelection(route_id, "batch")
    )
    route = rf00059_result._route_bodies[route_id]
    primary = copy.deepcopy(confirmed.purchase_fragments[0])
    primary["purchase_sequence"] = "A" + primary["purchase_sequence"][1:]
    verification = verify_exact_dna_assembly(
        query=confirmed.query,
        target_sequence=confirmed.target_sequence,
        route=route,
        primary_fragment=primary,
        database=load_plasmid_reference(),
        geometries=load_exact_dna_enzyme_catalog(),
    )
    assert verification["passed"] is False
    assert any("purchased block" in error for error in verification["errors"])


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

    with tempfile.TemporaryDirectory() as accepted_dir:
        accepted_files = write_exact_dna_outputs(accepted, accepted_dir)
        assert (Path(accepted_dir) / "order_ready_fragments.csv").is_file()
        assert not (Path(accepted_dir) / "idt_bulk_input.csv").exists()
        assert Path(accepted_files["technical_audit_zip"]).is_file()
    with tempfile.TemporaryDirectory() as rejected_dir:
        write_exact_dna_outputs(boundary, rejected_dir)
        assert (Path(rejected_dir) / "idt_bulk_input.csv").is_file()
        assert not (Path(rejected_dir) / "order_ready_fragments.csv").exists()


def test_gblock_only_route_scores_every_purchase_and_writes_two_file_package(
    rf00059_gblock_result, tmp_path
):
    result = confirm_best_exact_dna_route(
        rf00059_gblock_result,
        ExactDNASelection(
            rf00059_gblock_result.route_candidates[0]["route_id"], "api"
        ),
        idt_scorer=FakeScorer(0),
    )
    assert result.status == "idt_accepted_route"
    assert result.final_insert_sequence == EXAMPLE["repeat_unit"] * 4
    assert result.purchase_fragments
    assert all(
        row["product_type"] == "gblock"
        and 125 <= row["purchase_length_bp"] <= 3000
        and row["idt_accepted"] is True
        for row in result.purchase_fragments
    )
    output = tmp_path / "minimal"
    paths = write_exact_dna_minimal_outputs(result, output)
    assert {path.name for path in output.iterdir()} == {
        "cloning_steps.csv",
        "purchase_inserts.fasta",
    }
    rows = pd.read_csv(output / "cloning_steps.csv")
    assert list(rows.columns) == [
        "step",
        "purchase_insert",
        "purchase_length_bp",
        "prepare_insert_with_RE",
        "clone_with_RE",
        "reused_from_step",
        "IDT_accepted",
    ]
    assert rows.IDT_accepted.all()
    assert Path(paths["validation_details_zip"]).is_file()


def test_gblock_only_retries_padding_then_fragmented_route(rf00059_gblock_result):
    scorer = RejectUnsplitRF00059Donor()
    result = confirm_best_exact_dna_route(
        rf00059_gblock_result,
        ExactDNASelection(
            rf00059_gblock_result.route_candidates[0]["route_id"], "api"
        ),
        idt_scorer=scorer,
    )
    assert result.status == "idt_accepted_route"
    assert len(result.cloning_steps) > 2
    assert any(len(sequence) > 300 for _name, sequence in scorer.calls)
    assert all(row["idt_accepted"] is True for row in result.purchase_fragments)
    assert len(result.route_attempts) >= 17


def test_gblock_only_rejects_capacity_below_product_minimum():
    with pytest.raises(ValueError, match="requires max_purchase_bp"):
        _rf00059_query(
            purchase_policy=IDT_GBLOCK_ONLY_PURCHASE_POLICY,
            max_purchase_bp=124,
        )


def test_unknown_purchase_policy_is_rejected():
    with pytest.raises(ValueError, match="purchase_policy"):
        _rf00059_query(purchase_policy="unknown")


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
    summary = pd.read_csv(output / "cloning_summary.csv").iloc[0]
    assert summary["status"] == "bulk_export_unvalidated"
    assert int(summary["target_length_bp"]) == len(query.target_sequence)
    assert (output / "idt_bulk_input.csv").is_file()
    assert not (output / "idt_raw_audit.jsonl").exists()
