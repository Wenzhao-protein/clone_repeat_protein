from __future__ import annotations

import hashlib
import json

import pandas as pd

from hurdler.complete_route import (
    TransitionEvidence,
    build_element_matrix,
    finalize_complete_route_shards,
    find_shortest_purchasable_seed,
    search_complete_repeat_routes,
)
from hurdler.dna_assembly import TargetRecord


def _edge(
    start: int,
    end: int,
    *,
    plasmid: str = "pUC18",
    pair=("EcoRI", "HindIII"),
    product_type: str = "annealed_sticky_end_primer_pair",
    ceiling: int = 60,
):
    target = "ACGTTA" * end
    digest = "ACGTTA" * (end - start)
    purchase_sha = hashlib.sha256(digest.encode()).hexdigest()
    route_id = (
        f"route_{start}_{end}_{plasmid}_{pair[0]}_{pair[1]}_"
        f"{product_type}_{ceiling}"
    )
    return TransitionEvidence(
        recipient_copy_count=start,
        result_copy_count=end,
        donor_copy_count=end - start,
        plasmid=plasmid,
        site_i_enzyme=pair[0],
        site_ii_enzyme=pair[1],
        route_id=route_id,
        hurdle_steps=1,
        purchase_sha256s=(purchase_sha,),
        total_purchase_bp=len(digest),
        maximum_idt_score=None,
        whole_target_idt_status="failed",
        whole_target_idt_score=20.0,
        whole_target_idt_response_sha256="fixture-response",
        final_target_exact=True,
        all_purchase_fragments_accepted=True,
        route_row={
            "route_id": route_id,
            "fragment_purchase_ceiling_bp": ceiling,
            "local_constraints_passed": True,
        },
        step_rows=(
            {
                "route_id": route_id,
                "unintended_cut_count": 0,
                "double_strand_source_verified": True,
            },
        ),
        fragment_rows=(
            {
                "route_id": route_id,
                "fragment_id": f"fragment_{start}_{end}",
                "purchase_sha256": purchase_sha,
                "product_type": product_type,
                "purchase_sequence": digest if product_type.startswith("annealed") else "A" * 100,
                "purchase_length_bp": len(digest),
                "core_length_bp": min(80, len(digest)),
                "idt_response_sha256": "",
            },
        ),
    )


def _base():
    return TargetRecord(
        target_id="fixture-element",
        sequence="ACGTTA",
        cohort="real_element_derived",
        architecture="exact_tandem",
        unit_sequence="ACGTTA",
        copy_count=1,
        source_database="fixture",
        element_id="fixture-element",
    )


def test_shortest_seed_uses_unscored_primer_pair_under_90_bp():
    seed = find_shortest_purchasable_seed("ACGTTA", require_idt=True)
    assert seed is not None
    assert seed.copy_count == 1
    assert seed.product_type == "duplexed_seed_oligo_pair"
    assert seed.idt_status == "not_applicable_primer_pair_under_90bp"
    assert seed.accepted


def test_complete_route_uses_one_plasmid_and_minimizes_steps():
    transitions = {
        (1, 2): [_edge(1, 2)],
        (2, 4): [_edge(2, 4, pair=("EcoRI", "SpeI"))],
        # Direct 1->4 is one HURDLER cycle and must beat 1->2->4.
        (1, 4): [_edge(1, 4, pair=("MluI", "SpeI"))],
    }

    result = search_complete_repeat_routes(
        _base(),
        {},
        pd.DataFrame({"pUC18": [True]}, index=["EcoRI"]),
        require_idt=False,
        target_copy_counts=(2, 4),
        max_copy_count=4,
        transition_provider=lambda start, end: transitions.get((start, end), []),
    )
    targets = result["targets"].set_index("target_copy_count")
    assert targets.loc[2, "complete_route_verified"]
    assert targets.loc[4, "complete_route_verified"]
    assert targets.loc[4, "hurdler_step_count"] == 1
    route = result["selected_routes"].query("target_copy_count == 4").iloc[0]
    assert route.plasmid == "pUC18"
    assert route.transition_count == 1
    assert route.final_target_exact


def test_plasmid_cannot_change_inside_a_complete_route():
    transitions = {
        (1, 2): [_edge(1, 2, plasmid="pUC18")],
        (2, 4): [_edge(2, 4, plasmid="pCold_I")],
    }
    plasmids = pd.DataFrame(
        {"pUC18": [True], "pCold_I": [True]}, index=["EcoRI"]
    )
    result = search_complete_repeat_routes(
        _base(),
        {},
        plasmids,
        require_idt=False,
        target_copy_counts=(2, 4),
        max_copy_count=4,
        transition_provider=lambda start, end: transitions.get((start, end), []),
    )
    targets = result["targets"].set_index("target_copy_count")
    assert targets.loc[2, "complete_route_verified"]
    assert not targets.loc[4, "complete_route_verified"]


def test_element_matrix_retains_all_copy_results_without_any_pass_collapse():
    rows = []
    for copies in (2, 4, 8, 16, 32):
        passed = copies in {2, 4, 8, 16}
        rows.append(
            {
                "source_database": "Rfam",
                "element_id": "RFfixture",
                "unit_sequence": "ACGTTA",
                "unit_length_bp": 6,
                "target_copy_count": copies,
                "seed_copy_count": 1,
                "complete_route_verified": passed,
                "final_target_exact": passed,
                "failure_reason": "" if passed else "no_complete_route",
                "complete_route_id": f"route-{copies}" if passed else "",
            }
        )
    matrix = build_element_matrix(pd.DataFrame(rows))
    assert len(matrix) == 1
    assert matrix.iloc[0].successful_target_count == 4
    assert not matrix.iloc[0].all_five_complete
    assert matrix.iloc[0].maximum_verified_copy_count == 16


def test_live_idt_is_applied_only_to_candidate_purchases_and_score_10_fails():
    class ThresholdScorer:
        def __init__(self):
            self.calls = []

        def score(self, name, sequence):
            self.calls.append((name, sequence))
            return {
                "idt_status": "failed",
                "idt_explicit_pass": False,
                "idt_complexity_score": 10.0,
                "idt_response_sha256": "response-at-threshold",
                "idt_scored_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            }

    rejected = _edge(1, 2, product_type="gblock", ceiling=3000)
    accepted = _edge(1, 2, product_type="annealed_sticky_end_primer_pair", ceiling=60)
    scorer = ThresholdScorer()
    result = search_complete_repeat_routes(
        _base(),
        {},
        pd.DataFrame({"pUC18": [True]}, index=["EcoRI"]),
        idt_scorer=scorer,
        require_idt=True,
        target_copy_counts=(2,),
        max_copy_count=2,
        transition_provider=lambda start, end: (
            [rejected, accepted] if (start, end) == (1, 2) else []
        ),
    )
    target = result["targets"].iloc[0]
    assert target.complete_route_verified
    assert result["selected_routes"].iloc[0].transition_route_ids_json == json.dumps(
        [accepted.route_id]
    )
    assert len(scorer.calls) == 1
    audit = result["candidate_routes"]
    assert set(audit.all_purchase_fragments_accepted) == {False, True}


def test_complete_finalizer_enforces_five_targets_and_writes_headline(tmp_path):
    transitions = {
        (1, copies): [_edge(1, copies)] for copies in (2, 4, 8, 16, 32)
    }
    result = search_complete_repeat_routes(
        _base(),
        {},
        pd.DataFrame({"pUC18": [True]}, index=["EcoRI"]),
        require_idt=False,
        transition_provider=lambda start, end: transitions.get((start, end), []),
    )
    shard = tmp_path / "shard_00000"
    shard.mkdir()
    for name, frame in result.items():
        frame.to_parquet(shard / f"complete_route_{name}.parquet", index=False)
    (shard / "idt_audit.jsonl").touch()
    (shard / "complete_route_manifest.json").write_text(
        json.dumps(
            {
                "version": "arbitrary-dna-complete-route-v2",
                "shard_index": 0,
                "shard_count": 1,
                "element_rows": 1,
                "idt_required": True,
                "output_rows": {name: len(frame) for name, frame in result.items()},
            }
        )
    )
    output = tmp_path / "final"
    finalized = finalize_complete_route_shards(
        [shard], output, expected_public_elements=1, expected_real_targets=5
    )
    assert len(finalized["element_matrix"]) == 1
    assert (output / "production_target_analysis.parquet").is_file()
    headline = json.loads((output / "production_headline_summary.json").read_text())
    assert headline["real_exact_targets"] == 5
    assert headline["real_exact_targets_complete"] == 5
