from __future__ import annotations

import hashlib

from hurdler.index import PatternIndex
from hurdler.optimization import translate_dna
from hurdler.plasmid_reference import (
    CutterSite,
    decide_cutter_silencing,
    load_plasmid_reference,
    validate_plasmid_reference,
)
from hurdler.protein_index import ProteinPatternIndex
from hurdler.vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    DesignSelection,
    design_construct_v2,
    design_query,
)


def test_bundled_plasmid_database_has_seven_references_eight_profiles_and_four_schemes():
    database = load_plasmid_reference()
    summary = validate_plasmid_reference(database)
    assert summary["physical_reference_count"] == 7
    assert summary["profile_count"] == 8
    assert summary["scheme_count"] == 32
    assert all(len(database.schemes_for(profile.profile_id)) == 4 for profile in database.profiles)
    assert database.profile("pET-28a(+)").reference_id == database.profile("pET-28a(+)_start_codon").reference_id


def test_pet28_ncoi_start_codon_is_never_silenceable():
    import hurdler.plasmid_reference as reference_module

    database = load_plasmid_reference()
    profile = database.profile("pET-28a(+)_start_codon")
    reference = database.reference(profile.reference_id)
    ncoi = next(
        cutter
        for cutter in reference_module._candidate_cutter_sites(reference, profile)
        if "NcoI" in cutter.enzyme_aliases
    )
    decision = decide_cutter_silencing(database, profile.profile_id, ncoi)
    assert decision.allowed is False
    assert decision.status == "protected_feature"
    assert any(value.startswith("start_codon:") for value in decision.overlapping_features)


def test_stop_codon_is_the_only_protected_silencing_exception():
    database = load_plasmid_reference()
    stop = next(
        feature
        for feature in database.reference("pGEX-4T-1").features
        if feature.feature_class == "stop_codon"
    )
    start, end = stop.intervals[0]
    cutter = CutterSite(
        side="right",
        location="outside",
        canonical_enzyme="synthetic_stop_fixture",
        enzyme_aliases=("synthetic_stop_fixture",),
        recognition_site="TAATAG",
        oriented_start=start,
        oriented_end=min(end, start + 6),
        physical_start=start,
        physical_end=min(end, start + 6),
        top_cut_oriented=start,
        bottom_cut_oriented=start + 4,
        overhang_length=-4,
        outside_mcs_occurrences=1,
    )
    decision = decide_cutter_silencing(database, "pGEX-4T-1", cutter)
    assert decision.allowed is True
    assert decision.status == "stop_rescue_then_silence"
    assert decision.stop_rescue_sequence == "TAA"


def test_protein_only_index_restores_all_776_pairs_without_losing_legacy_pairs():
    protein = ProteinPatternIndex.load("data/artifacts/vector-aware-hurdler-v2")
    legacy = PatternIndex.load("data/artifacts/legacy-optimized-v1")
    assert len(protein.pair_table) == 776
    assert len(legacy.pair_table) == 512
    columns = ["site_i_enzyme", "site_ii_enzyme"]
    new_pairs = set(map(tuple, protein.pair_table[columns].itertuples(index=False, name=None)))
    old_pairs = set(map(tuple, legacy.pair_table[columns].itertuples(index=False, name=None)))
    assert old_pairs <= new_pairs
    assert len(new_pairs - old_pairs) == 264


def _compatible_query() -> CompatibilityQuery:
    return CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode="split",
        sequence_id="v2_smoke",
        repeat_module="ACDEFGHIKLMNPQRSTVWY",
        repeat_copies=3,
    )


def test_query_is_protein_first_and_returns_pair_to_profile_to_scheme_routes():
    result = design_query(_compatible_query())
    assert result.status == "compatible_unoptimized"
    assert result.protein_candidates
    assert result.vector_routes
    row = result.vector_routes[0]
    assert row["site_i_enzyme"] and row["site_ii_enzyme"]
    assert row["profile_id"] and row["scheme_id"]
    assert row["protein_index_version"] == "vector-aware-hurdler-v2"
    assert all(not route["cutter_reuse"] for route in result.vector_routes)


def test_batch_mode_never_calls_idt_and_preserves_translation():
    query = _compatible_query()
    queried = design_query(query)
    route = queried.vector_routes[0]

    class MustNotBeCalled:
        def score(self, name: str, sequence: str):
            raise AssertionError("batch mode made an HTTP/scorer call")

    request = DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection(
            route["candidate_id"], route["profile_id"], route["scheme_id"], route["site_iii_options"][0]
        ),
        validation_mode="batch",
        population_size=4,
        generation_schedule=(10, 100),
    )
    result = design_construct_v2(request, idt_scorer=MustNotBeCalled())
    assert result.status == "optimized_unvalidated_batch"
    assert translate_dna(result.final_dna_sequence) == result.final_protein_sequence
    assert result.idt_audit == []
    assert result.final_plasmid["cds_translation_exact"] is True
    assert result.final_plasmid["protected_feature_sequences_preserved"] is True
    assert result.final_plasmid["restoration_exact"] is True


def test_api_mode_scores_exact_purchase_fragment_hashes():
    query = _compatible_query()
    queried = design_query(query)
    route = queried.vector_routes[0]

    class PassingScorer:
        def score(self, name: str, sequence: str):
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            return {
                "idt_status": "passed",
                "idt_explicit_pass": True,
                "idt_complexity_score": 0.0,
                "idt_score_complete": True,
                "idt_rule_details_json": "[]",
                "idt_positive_score_names_json": "[]",
                "idt_violation_names_json": "[]",
                "idt_scored_sequence_sha256": digest,
                "idt_response_sha256": hashlib.sha256(("response:" + digest).encode()).hexdigest(),
            }

    request = DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection(
            route["candidate_id"], route["profile_id"], route["scheme_id"], route["site_iii_options"][0]
        ),
        validation_mode="api",
        population_size=4,
        generation_schedule=(10, 100),
    )
    result = design_construct_v2(request, idt_scorer=PassingScorer())
    assert result.status == "idt_accepted"
    assert result.idt_audit
    for fragment, audit in zip(result.primary_fragments, result.idt_audit, strict=True):
        assert fragment["purchase_sha256"] == audit["dna_sha256"]
        assert float(fragment["idt_complexity_score"]) < 10


def test_split_adaptive_search_reaches_hard_upper_bound_with_machine_readable_proof():
    query = _compatible_query()
    queried = design_query(query)
    route = queried.vector_routes[0]
    request = DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection(
            route["candidate_id"], route["profile_id"], route["scheme_id"], route["site_iii_options"][0]
        ),
        validation_mode="batch",
        max_repeat_copies=3,
        population_size=4,
        generation_schedule=(10, 100),
    )
    result = design_construct_v2(request)
    assert result.status == "optimized_unvalidated_batch"
    assert result.selected_route["maximum_verified_repeat_copies"] == 3
    assert result.termination_reason == "reached_local_upper_bound"
    assert "binary_short" in {row["phase"] for row in result.optimization_attempts}
    assert any(row["copies"] == 3 and row["passed"] for row in result.optimization_attempts)


def test_full_sequence_mode_preserves_every_variant_residue():
    unit_a = "ACDEFGHIKLMNPQRSTVWY"
    unit_b = "ACDEYGHIKLMNPQRSTVWY"
    protein = "M" + unit_a + unit_b + unit_a + "G"
    query = CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode="full",
        sequence_id="heterogeneous_repeat",
        full_protein_sequence=protein,
        repeat_region_start=2,
        repeat_region_end=61,
        repeat_period=20,
    )
    result = design_query(query)
    assert result.final_protein_sequence == protein
    assert result.confirmed_boundary["unit_sequences"] == [unit_a, unit_b, unit_a]
