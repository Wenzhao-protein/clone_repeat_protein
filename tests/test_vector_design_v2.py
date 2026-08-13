from __future__ import annotations

import hashlib
from dataclasses import asdict, replace
from pathlib import Path

import pytest

from hurdler.index import PatternIndex
from hurdler.optimization import translate_dna
from hurdler.plasmid_reference import (
    CutterSite,
    decide_cutter_silencing,
    load_plasmid_reference,
    validate_plasmid_reference,
)
from hurdler.protein_index import ProteinPatternIndex
from hurdler.progress import DesignProgressEvent
from hurdler.vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRouteUniverse,
    DesignRequestV2,
    DesignSelection,
    _exact_split_boundary,
    _adapt_ga_parameters_from_idt,
    _idt_feedback_guidance,
    _merge_idt_feedback_guidance,
    filter_route_universe,
    design_construct_v2,
    design_query,
    write_design_outputs_v2,
)


def test_idt_feedback_coordinates_are_clipped_to_mutable_coding_core():
    fragment = {
        "fragment_id": "secondary",
        "purchase_sequence": "AAAACCGGTTAACCTTTT",
        "core_start_bp": 4,
        "core_end_bp": 14,
    }
    audit = {
        "fragment_id": "secondary",
        "idt_rule_details_json": (
            '[{"name":"Overall Repeat","score":4.0,"threshold_value":40},'
            '{"name":"Repeat Length (Fragment)","score":3.0,'
            '"repeated_segment":"CCGG","forward_locations":[0,4,10]}]'
        ),
    }
    guidance = _idt_feedback_guidance([audit], [fragment])
    assert guidance["repeat_coverage_threshold"] == 40
    assert [0, 10] in guidance["target_ranges"]
    assert [0, 4] in guidance["target_ranges"]
    assert [6, 10] in guidance["target_ranges"]
    assert guidance["repeat_aware_steps"] == 40_000


def test_idt_feedback_accumulates_distinct_worst_windows_across_rounds():
    first = {
        "repeat_coverage_threshold": 40.0,
        "repeat_windows": [{"start": 10, "end": 100, "threshold": 88.0}],
        "target_ranges": [[10, 100]],
        "avoid_segments": ["AACCGGTTAACC"],
        "repeat_aware_steps": 40_000,
        "hotspot_mutation_rate": 0.65,
    }
    second = {
        "repeat_coverage_threshold": 40.0,
        "repeat_windows": [{"start": 300, "end": 390, "threshold": 88.0}],
        "target_ranges": [[300, 390]],
        "avoid_segments": ["TTGGCCAATTGG"],
        "repeat_aware_steps": 40_000,
        "hotspot_mutation_rate": 0.65,
    }
    merged = _merge_idt_feedback_guidance(first, second)
    assert len(merged["repeat_windows"]) == 2
    assert merged["target_ranges"] == [[10, 100], [300, 390]]
    assert merged["avoid_segments"] == ["AACCGGTTAACC", "TTGGCCAATTGG"]


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


def _restoration_filter_universe() -> DesignRouteUniverse:
    candidate = {
        "candidate_id": "candidate-1",
        "site_i_enzyme": "SiteI",
        "site_ii_enzyme": "SiteII",
        "site_iii_options": ["SiteIII"],
    }
    route_common = {
        "candidate_id": "candidate-1",
        "profile_id": "plasmid-1",
        "site_i_enzyme": "SiteI",
        "site_ii_enzyme": "SiteII",
        "site_iii_options": ["SiteIII"],
    }
    return DesignRouteUniverse(
        confirmed_boundary={"period": 6},
        protein_candidates=[candidate],
        vector_routes=[
            {
                **route_common,
                "scheme_id": "unmodified-101",
                "restoration_length_bp": 101,
                "cutter_reuse": False,
            },
            {
                **route_common,
                "scheme_id": "fallback-100",
                "restoration_length_bp": 100,
                "cutter_reuse": True,
            },
        ],
        final_protein_sequence="ACDEFGACDEFG",
    )


def test_restoration_limit_is_inclusive_and_precedes_cutter_reuse_fallback():
    query = replace(_compatible_query(), max_restoration_length_bp=100)
    result = filter_route_universe(_restoration_filter_universe(), query)
    assert result.status == "compatible_unoptimized"
    assert [row["scheme_id"] for row in result.vector_routes] == ["fallback-100"]
    assert result.request["max_restoration_length_bp"] == 100

    unbounded = filter_route_universe(
        _restoration_filter_universe(), replace(query, max_restoration_length_bp=None)
    )
    assert [row["scheme_id"] for row in unbounded.vector_routes] == ["unmodified-101"]


def test_restoration_limit_zero_and_all_routes_filtered_message():
    query = replace(_compatible_query(), max_restoration_length_bp=0)
    result = filter_route_universe(_restoration_filter_universe(), query)
    assert result.status == "no_vector_route"
    assert result.vector_routes == []
    assert "all require restoration longer than 0 bp" in result.message


@pytest.mark.parametrize("value", [-1, 1.5, True, "100"])
def test_restoration_limit_rejects_negative_or_non_integer_values(value):
    with pytest.raises(ValueError, match="non-negative integer"):
        replace(_compatible_query(), max_restoration_length_bp=value)


def test_compatibility_query_missing_restoration_limit_remains_unbounded():
    payload = {
        "schema_version": DESIGN_SCHEMA_VERSION_V2,
        "input_mode": "split",
        "sequence_id": "legacy_request",
        "repeat_module": "ACDEFG",
        "repeat_copies": 2,
    }
    assert CompatibilityQuery.from_dict(payload).max_restoration_length_bp is None


def test_secondary_search_bounds_validate_and_legacy_payload_is_unbounded():
    query = _compatible_query()
    route = design_query(query).vector_routes[0]
    common = dict(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection(
            route["candidate_id"], route["profile_id"], route["scheme_id"],
            route["site_iii_options"][0],
        ),
    )
    request = DesignRequestV2(**common, minimum_secondary_copies=2)
    assert request.maximum_secondary_copies is None
    assert DesignRequestV2.from_dict(asdict(request)).maximum_secondary_copies is None
    with pytest.raises(ValueError, match="cannot be smaller"):
        DesignRequestV2(**common, minimum_secondary_copies=3, maximum_secondary_copies=2)
    with pytest.raises(ValueError, match="positive integer"):
        DesignRequestV2(**common, maximum_secondary_copies=True)


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
        minimum_secondary_copies=100,
    )
    result = design_construct_v2(request, idt_scorer=MustNotBeCalled())
    assert result.status == "optimized_unvalidated_batch"
    assert translate_dna(result.final_dna_sequence) == result.final_protein_sequence
    assert result.idt_audit == []
    assert result.final_plasmid["cds_translation_exact"] is True
    assert result.final_plasmid["protected_feature_sequences_preserved"] is True
    assert result.final_plasmid["restoration_exact"] is True
    assert result.rdl_plan["minimum_secondary_bypassed_by_single_purchase"] is True


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


def test_live_idt_missing_numeric_score_is_a_fatal_visible_result():
    query = _compatible_query()
    route = design_query(query).vector_routes[0]

    class MissingScoreScorer:
        def score(self, name: str, sequence: str):
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            return {
                "idt_status": "scored_unclassified",
                "idt_explicit_pass": None,
                "idt_complexity_score": None,
                "idt_score_complete": False,
                "idt_invalid_score_names_json": '["Overall Repeat"]',
                "idt_rule_details_json": '[{"name":"Overall Repeat","score":null}]',
                "idt_positive_score_names_json": "[]",
                "idt_violation_names_json": "[]",
                "idt_scored_sequence_sha256": digest,
                "idt_response_sha256": "response-hash",
            }

    result = design_construct_v2(
        DesignRequestV2(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            query=query,
            selection=DesignSelection(
                route["candidate_id"], route["profile_id"], route["scheme_id"],
                route["site_iii_options"][0],
            ),
            validation_mode="api",
            population_size=4,
            max_idt_feedback_rounds=2,
            generation_schedule=(10, 100),
        ),
        idt_scorer=MissingScoreScorer(),
    )
    assert result.status == "idt_score_error"
    assert result.termination_reason == "idt_score_error"
    assert "Overall Repeat" in result.message
    assert not result.final_dna_sequence
    assert result.idt_feedback_history[-1]["response_sha256"] == "response-hash"


def test_live_idt_transport_failure_is_not_misclassified_as_ga_rejection():
    query = _compatible_query()
    route = design_query(query).vector_routes[0]

    class BrokenScorer:
        def score(self, name: str, sequence: str):
            raise TimeoutError("retry budget exhausted")

    result = design_construct_v2(
        DesignRequestV2(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            query=query,
            selection=DesignSelection(
                route["candidate_id"], route["profile_id"], route["scheme_id"],
                route["site_iii_options"][0],
            ),
            validation_mode="api",
            population_size=4,
            generation_schedule=(10, 100),
        ),
        idt_scorer=BrokenScorer(),
    )
    assert result.status == "idt_api_error"
    assert "TimeoutError" in result.message
    assert not result.final_dna_sequence


def test_rejected_idt_feedback_is_bounded_adaptive_and_never_resubmits_same_dna():
    query = _compatible_query()
    route = design_query(query).vector_routes[0]

    class RejectingScorer:
        def __init__(self):
            self.hashes: list[str] = []

        def score(self, name: str, sequence: str):
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            self.hashes.append(digest)
            return {
                "idt_status": "failed",
                "idt_explicit_pass": False,
                "idt_complexity_score": 12.0,
                "idt_score_complete": True,
                "idt_rule_details_json": '[{"name":"Overall Repeat","score":12.0}]',
                "idt_positive_score_names_json": '["Overall Repeat"]',
                "idt_violation_names_json": "[]",
                "idt_scored_sequence_sha256": digest,
                "idt_response_sha256": hashlib.sha256(("response:" + digest).encode()).hexdigest(),
            }

    scorer = RejectingScorer()
    result = design_construct_v2(
        DesignRequestV2(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            query=query,
            selection=DesignSelection(
                route["candidate_id"], route["profile_id"], route["scheme_id"],
                route["site_iii_options"][0],
            ),
            validation_mode="api",
            population_size=16,
            minimum_secondary_copies=100,
            max_idt_feedback_rounds=2,
            generations_per_feedback_round=1,
            generation_schedule=(10, 100),
        ),
        idt_scorer=scorer,
    )
    direct = [
        row for row in result.optimization_attempts
        if row.get("component") == "direct_primary"
    ]
    assert len(direct) == 2
    assert len(scorer.hashes) == len(set(scorer.hashes))
    assert result.termination_reason == "minimum_secondary_exceeds_capacity"
    assert result.ga_parameter_history
    assert result.ga_parameter_history[0]["parameter_tier"] == "10_to_20"
    assert result.ga_parameter_history[0]["new_population_size"] == 20
    assert result.ga_parameter_history[0]["new_mutation_rate"] == pytest.approx(0.088)


@pytest.mark.parametrize(
    ("score", "expected_population", "expected_mutation", "expected_crossover", "tier"),
    [
        (12.0, 20, 0.088, 0.77, "10_to_20"),
        (25.0, 24, 0.10, 0.80, "20_to_50"),
        (60.0, 32, 0.12, 0.85, "50_or_more"),
    ],
)
def test_idt_score_tiers_adapt_ga_parameters(
    score, expected_population, expected_mutation, expected_crossover, tier
):
    query = _compatible_query()
    route = design_query(query).vector_routes[0]
    request = DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection(
            route["candidate_id"], route["profile_id"], route["scheme_id"],
            route["site_iii_options"][0],
        ),
        population_size=16,
    )
    population, mutation, crossover, policy = _adapt_ga_parameters_from_idt(
        request,
        score,
        population_size=16,
        mutation_rate=0.08,
        crossover_rate=0.75,
    )
    assert population == expected_population
    assert mutation == pytest.approx(expected_mutation)
    assert crossover == pytest.approx(expected_crossover)
    assert policy["tier"] == tier


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
        assembly_strategy="legacy_adaptive_max",
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


def test_exact_target_rdl_reuses_one_secondary_and_emits_progress(tmp_path):
    n_cap = "MGSHHHHHHSSGIEGRSSGYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGGGGSGGGGSLEVLFQGPDLPKLVKLLKSSNEEILLKALRALAEIASGG"
    module = "NEQIQAVIDAGALPALVQLLSSPNEQILQEALWALSNIASGG"
    c_cap = "NEQIQAVIDAGALPALVQLLSSPNEQILQEALWALSNIASGGNEQKQAVKEAGALEKLEQLQSHENEKIQKEAQEALEKLQSHGGGLEVLFQGPSSGEFGGGGSMVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYKGGGSSGHHHHHH"
    query = CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode="split",
        sequence_id="exact_25_rdl",
        n_cap=n_cap,
        repeat_module=module,
        repeat_copies=25,
        c_cap=c_cap,
    )
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

    events: list[DesignProgressEvent] = []
    checkpoints: list[dict] = []
    result = design_construct_v2(
        DesignRequestV2(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            query=query,
            selection=DesignSelection(
                route["candidate_id"], route["profile_id"], route["scheme_id"],
                route["site_iii_options"][0],
            ),
            validation_mode="api",
            population_size=4,
            generation_schedule=(10, 100),
            minimum_secondary_copies=12,
            max_idt_feedback_rounds=3,
        ),
        idt_scorer=PassingScorer(),
        progress_callback=events.append,
        checkpoint_callback=checkpoints.append,
    )
    assert result.status == "idt_accepted"
    plan = result.rdl_plan
    assert plan["primary_repeat_copies"] + plan["secondary_reuse_count"] * plan["secondary_repeat_copies"] == 25
    assert plan["primary_repeat_copies"] >= 1
    assert plan["secondary_reuse_count"] == 1
    assert plan["secondary_repeat_copies"] == result.maximum_secondary_evidence["maximum_verified_copies"]
    assert plan["secondary_repeat_copies"] >= 12
    assert plan["minimum_secondary_satisfied"] is True
    assert result.ga_elite_candidates
    assert {
        row["fragment_kind"] for row in result.ga_elite_candidates
    } == {"primary", "secondary"}
    assert all(row["dna_sequence"] for row in result.ga_elite_candidates)
    assert result.idt_feedback_history
    assert plan["same_secondary_reused_every_round"] is True
    assert result.final_protein_sequence == n_cap + module * 25 + c_cap
    assert len(result.final_protein_sequence) == 1_524
    assert len(result.final_dna_sequence) == 4_572
    assert translate_dna(result.final_dna_sequence) == result.final_protein_sequence
    assert all(row["purchase_length_bp"] <= 3000 for row in [*result.primary_fragments, *result.secondary_fragments])
    assert all(row["passed"] for row in result.intermediate_validations)
    assert events and events[0].stage == "design"
    assert any(event.stage == "ga" and event.status == "running" for event in events)
    scored_events = [
        event for event in events
        if event.stage == "idt" and event.status == "fragment_scored"
    ]
    assert scored_events
    assert all(event.event_name == "idt_fragment_scored" for event in scored_events)
    assert [event.idt_evaluation_index for event in scored_events] == list(
        range(1, len(scored_events) + 1)
    )
    assert all(event.idt_fragment_name for event in scored_events)
    assert all(event.idt_classification == "passed" for event in scored_events)
    assert all(event.idt_response_sha256 for event in scored_events)
    assert events[-1].stage == "rdl_plan" and events[-1].status == "completed"
    assert checkpoints
    assert checkpoints[-1]["repeat_copies"] == plan["secondary_repeat_copies"]
    assert checkpoints[-1]["idt_complexity_score"] < 10
    assert hashlib.sha256(checkpoints[-1]["purchase_sequence"].encode()).hexdigest() == checkpoints[-1]["purchase_sha256"]
    paths = write_design_outputs_v2(result, tmp_path)
    expected_steps = plan["secondary_reuse_count"] + 1
    plasmid_files = sorted(tmp_path.glob("step*_plasmid.gb"))
    insert_files = sorted(tmp_path.glob("step*_insert.gb"))
    assert len(plasmid_files) == expected_steps + 1  # includes step00
    assert len(insert_files) == expected_steps
    assert paths["final_plasmid_genbank"].endswith("final_plasmid.gb")
    assert Path(paths["idt_score_history_csv"]).stat().st_size > 0
    for suffix in ("png", "pdf", "svg"):
        assert Path(paths[f"idt_score_trajectory_{suffix}"]).stat().st_size > 0
    from Bio import SeqIO

    final_record = SeqIO.read(plasmid_files[-1], "genbank")
    assert str(final_record.seq) == result.final_plasmid["final_plasmid_sequence"]
    target_cds = next(
        feature for feature in final_record.features
        if feature.type == "CDS" and "repeat-protein CDS" in feature.qualifiers.get("label", [""])[0]
    )
    assert target_cds.qualifiers["translation"][0] == result.final_protein_sequence
    secondary_records = [SeqIO.read(path, "genbank") for path in insert_files[1:]]
    assert {
        "".join(next(feature for feature in record.features if feature.type == "source").qualifiers["purchase_sha256"])
        for record in secondary_records
    } == {result.secondary_fragments[0]["purchase_sha256"]}


def test_exact_rdl_internal_primary_boundary_supports_one_physical_copy():
    query = CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode="split",
        sequence_id="one_copy_primary",
        n_cap="M",
        repeat_module="ACDEFG",
        repeat_copies=25,
        c_cap="G",
    )
    boundary = _exact_split_boundary(query, 1)
    assert boundary.repeat_count == 1
    assert boundary.unit_sequences == ("ACDEFG",)
    assert boundary.n_terminal_flank == "M"
    assert boundary.c_terminal_flank == "G"


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
