import json

import pandas as pd
import pytest

import hurdler.ga_optimization as ga_optimization
from hurdler.ga_optimization import (
    GAPopulationState,
    GA_SCORE_PROFILE,
    GA_RE_SITE_POLICY,
    adaptive_copy_search,
    adjust_ga_score_profile_from_idt,
    ga_sequence_metrics,
    genetic_refine_dna,
    load_restriction_sites,
    repeat_aware_synonymous_refine,
    repeated_kmer_coverage,
    repeated_re_site_excess,
)
from hurdler.idt import summarize_complexity_response, write_cached_response
from hurdler.optimization import translate_dna


def test_repeated_re_site_is_an_explicit_ga_score_term():
    repeated = "GAATTCGAATTC"  # EFEF, with two EcoRI sites
    single = "GAGTTCGAATTT"  # same translation, without a repeated EcoRI site
    assert translate_dna(repeated) == translate_dna(single) == "EFEF"
    assert repeated_re_site_excess(repeated, ("GAATTC",)) == 1
    weights = {codon: 1.0 for codon in ("GAA", "GAG", "TTC", "TTT")}
    repeated_metrics = ga_sequence_metrics(repeated, weights, ("GAATTC",), {})
    single_metrics = ga_sequence_metrics(single, weights, ("GAATTC",), {})
    assert repeated_metrics["repeated_re_site_excess"] == 1
    assert repeated_metrics["ga_score"] > single_metrics["ga_score"]


def test_ga_refinement_preserves_translation_and_selected_site_limits():
    original = "GCA" * 4
    weights = {"GCG": 1.0, "GCC": 1.0, "GCA": 1.0, "GCT": 1.0}
    refined, metrics = genetic_refine_dna(
        original,
        locked_positions=set(),
        selected_site_limits={"GCAGCA": 0},
        recognition_sites=("GCAGCA",),
        codon_weights=weights,
        seed=42,
        population_size=4,
        generations=2,
    )
    assert translate_dna(refined) == "AAAA"
    assert "GCAGCA" not in refined
    assert metrics["selected_re_site_excess"] == 0
    assert "ga_initial_repeated_re_site_excess" in metrics
    assert metrics["ga_repeated_re_site_excess_removed"] >= 0


def test_warm_start_ga_retains_ranked_elites_and_continues_generation_count():
    original = "GCA" * 12
    weights = {"GCG": 1.0, "GCC": 1.0, "GCA": 1.0, "GCT": 1.0}
    first, first_metrics = genetic_refine_dna(
        original,
        locked_positions=set(),
        selected_site_limits={},
        recognition_sites=(),
        codon_weights=weights,
        seed=7,
        population_size=8,
        generations=2,
        elite_seed_count=4,
        capture_population_state=True,
    )
    first_state = first_metrics["ga_population_state"]
    assert isinstance(first_state, GAPopulationState)
    assert first_state.total_generations == 2
    assert len(first_state.elite_sequences) == 4
    second, second_metrics = genetic_refine_dna(
        original,
        locked_positions=set(),
        selected_site_limits={},
        recognition_sites=(),
        codon_weights=weights,
        seed=999,
        population_size=8,
        generations=3,
        elite_seed_count=4,
        capture_population_state=True,
        population_state=first_state,
    )
    second_state = second_metrics["ga_population_state"]
    assert second_state.total_generations == 5
    assert set(first_state.elite_sequences) & set(second_state.elite_sequences)
    assert translate_dna(first) == translate_dna(second) == "A" * 12


def test_repeat_aware_seed_reduces_coverage_and_long_repeats_without_changing_protein():
    protein = "NEQIQAVIDAGALPALVQLLSSP" * 5
    original = "".join(ga_optimization.GENETIC_CODE[aa][0] for aa in protein)
    refined, metrics = repeat_aware_synonymous_refine(
        original,
        locked_positions={0, 7},
        seed=42,
        steps=10_000,
    )
    assert translate_dna(refined) == protein
    assert refined[:3] == original[:3]
    assert refined[21:24] == original[21:24]
    assert repeated_kmer_coverage(refined, 8) < repeated_kmer_coverage(original, 8)
    assert metrics["repeat_aware_final_repeated_13mer"] <= metrics[
        "repeat_aware_initial_repeated_13mer"
    ]
    assert metrics["repeat_aware_final_repeated_14mer"] <= metrics[
        "repeat_aware_initial_repeated_14mer"
    ]


def test_local_gate_is_selected_pair_only_and_gc_is_left_to_idt():
    metrics = ga_sequence_metrics(
        "GCG" * 20,
        {"GCG": 1.0},
        recognition_sites=(),
        selected_site_limits={},
    )
    assert metrics["ga_gc_bounds_passed"] is False
    assert metrics["ga_local_constraints_passed"] is True


def test_restriction_site_catalog_canonicalizes_reverse_complements(tmp_path):
    catalog = tmp_path / "restriction_sites.csv"
    catalog.write_text("site\nGGTCTC\nGAGACC\nGAATTC\n")
    assert load_restriction_sites(catalog) == ("GAATTC", "GAGACC")


def test_adaptive_copy_search_uses_binary_then_single_copy_escalation():
    calls = []

    def evaluate(copies, generations):
        calls.append((copies, generations))
        passed = copies <= 3 or (copies == 4 and generations >= 40)
        return {
            "passed": passed,
            "dna_sequence": "ATG" * copies,
            "ga_score": float(copies),
            "repeated_re_site_excess": 0 if passed else 1,
            "selected_re_site_excess": 0,
        }

    copies, result, trace, reason = adaptive_copy_search(
        1,
        6,
        short_generations=10,
        generation_schedule=(10, 20, 40, 60, 80, 100),
        evaluate=evaluate,
    )
    assert copies == 4
    assert result["dna_sequence"] == "ATG" * 4
    assert reason == "copy_5_failed_at_100"
    assert (4, 10) in calls and (4, 20) in calls and (4, 40) in calls
    assert (5, 100) in calls
    assert all(item["phase"] in {"binary_short", "linear_escalation"} for item in trace)


def test_terminal_construction_failure_is_still_probed_through_100_generations():
    calls = []

    def evaluate(copies, generations):
        calls.append((copies, generations))
        return {"passed": False, "terminal": True, "error": "cannot construct"}

    copies, result, trace, reason = adaptive_copy_search(
        2,
        3,
        short_generations=10,
        generation_schedule=(10, 20, 40, 60, 80, 100),
        evaluate=evaluate,
    )
    assert copies == 0 and result is None
    assert reason == "copy_2_failed_at_100"
    assert (2, 100) in calls
    assert any(item["copies"] == 2 and item["generations"] == 100 for item in trace)


def test_idt_api_failure_is_an_audited_nonpass_not_a_shard_exception():
    class BrokenScorer:
        def score(self, name, sequence):
            raise TimeoutError("temporary outage")

    result = ga_optimization._score_idt_candidate(
        BrokenScorer(), "module|copies=2", "ATGATG"
    )
    assert result["idt_status"] == "api_failure"
    assert result["idt_explicit_pass"] is None
    assert result["idt_request_attempted"] is True
    assert result["idt_score_complete"] is False
    assert result["idt_scored_sequence_sha256"]


def test_adaptive_refinement_searches_mathematical_not_legacy_bound(monkeypatch):
    def construct_metrics(unit, copies, payload, codon_weights, **kwargs):
        assert kwargs["validate_hard_constraints"] is False
        return {"dna_sequence": "ATG" * copies}

    def refine(dna, **kwargs):
        copies = len(dna) // 3
        passed = copies <= 5
        return dna, {
            "ga_local_constraints_passed": passed,
            "repeated_re_site_excess": 0 if passed else 1,
            "selected_re_site_excess": 0,
            "ga_score": float(copies),
            "ga_generations": kwargs["generations"],
        }

    monkeypatch.setattr(ga_optimization, "_construct_metrics", construct_metrics)
    monkeypatch.setattr(ga_optimization, "genetic_refine_dna", refine)
    payload = {
        "module_id": "test",
        "fragment_limit_bp": 1800,
        "unit_sequence": "M",
        "verified_max_copies": 3,
        "mathematical_max_copies": 7,
        "dna_sequence": "ATG" * 3,
        "site_i_position": 0,
        "site_ii_position": 0,
        "site_i_recognition_site": "",
        "site_ii_recognition_site": "",
        "site_iii_sites": "",
    }
    result = ga_optimization._adaptive_refine_payload(
        payload,
        codon_weights={"ATG": 1.0},
        recognition_sites=(),
        seed=42,
        population_size=4,
        short_generations=10,
        generation_schedule=(10, 20, 40, 60, 80, 100),
    )
    assert result["pre_adaptive_verified_max_copies"] == 3
    assert result["adaptive_search_upper_bound_copies"] == 7
    assert result["adaptive_verified_max_copies"] == 5
    assert result["verified_max_copies"] == 5
    assert result["dna_sequence"] == "ATG" * 5
    assert result["adaptive_stop_reason"] == "copy_6_failed_at_100"


def test_adaptive_refinement_preserves_a_known_orderable_lower_bound(monkeypatch):
    def construct_metrics(unit, copies, payload, codon_weights, **kwargs):
        return {"dna_sequence": "ATG" * len(unit) * copies}

    def refine(dna, **kwargs):
        copies = len(dna) // 3
        passed = copies <= 3
        return dna, {
            "ga_local_constraints_passed": passed,
            "repeated_re_site_excess": 0 if passed else 1,
            "selected_re_site_excess": 0,
            "ga_score": float(copies),
            "ga_generations": kwargs["generations"],
        }

    class PassingScorer:
        def __init__(self):
            self.calls = []

        def score(self, name, sequence):
            self.calls.append((name, sequence))
            return {
                "idt_api_called": True,
                "idt_cache_hit": False,
                "idt_status": "passed",
                "idt_explicit_pass": True,
                "idt_violation_count": 0,
                "idt_violation_names_json": "[]",
                "idt_rule_scores_json": "{}",
                "idt_rule_details_json": "[]",
                "idt_scored_sequence_sha256": "known",
                "idt_scored_sequence_unchanged": True,
                "idt_result_selected_by_index": True,
            }

    monkeypatch.setattr(ga_optimization, "_construct_metrics", construct_metrics)
    monkeypatch.setattr(ga_optimization, "genetic_refine_dna", refine)
    scorer = PassingScorer()
    result = ga_optimization._adaptive_refine_payload(
        {
            "module_id": "known",
            "fragment_limit_bp": 3000,
            "unit_sequence": "M",
            "verified_max_copies": 3,
            "mathematical_max_copies": 7,
            "dna_sequence": "ATG" * 3,
            "known_orderable_copies": 3,
            "known_orderable_dna_sequence": "ATG" * 3,
            "known_orderable_dna_pre_ga": "ATG" * 3,
            "site_i_position": 0,
            "site_ii_position": 0,
            "site_i_recognition_site": "",
            "site_ii_recognition_site": "",
            "site_iii_sites": "",
        },
        codon_weights={"ATG": 1.0},
        recognition_sites=(),
        seed=42,
        population_size=4,
        short_generations=10,
        generation_schedule=(10, 20, 40, 60, 80, 100),
        idt_scorer=scorer,
        require_idt_orderable=True,
    )
    assert result["adaptive_search_minimum_copies"] == 3
    assert result["adaptive_known_orderable_copies"] == 3
    assert result["verified_max_copies"] == 3
    assert result["adaptive_stop_reason"] == "copy_4_failed_at_100"
    assert result["adaptive_orderable_passed"] is True
    assert scorer.calls and scorer.calls[0][1] == "ATG" * 3


def test_capacity_rows_reuse_smaller_cap_as_monotonic_lower_bound(
    monkeypatch, tmp_path
):
    inputs = tmp_path / "inputs.parquet"
    pd.DataFrame(
        [
            {
                "module_id": "same",
                "fragment_limit_bp": 1800,
                "unit_sequence": "M",
                "mathematical_max_copies": 600,
                "verified_max_copies": 2,
                "dna_sequence": "ATGATG",
            },
            {
                "module_id": "same",
                "fragment_limit_bp": 3000,
                "unit_sequence": "M",
                "mathematical_max_copies": 1000,
                "verified_max_copies": 2,
                "dna_sequence": "ATGATG",
            },
        ]
    ).to_parquet(inputs, index=False)
    observed = []

    def fake_adaptive(payload, **kwargs):
        observed.append(dict(payload))
        copies = 3 if payload["fragment_limit_bp"] == 1800 else 5
        return {
            **payload,
            "verified_max_copies": copies,
            "dna_sequence": "ATG" * copies,
            "dna_sequence_pre_ga": "ATG" * copies,
            "ga_status": "passed",
            "idt_status": "passed",
        }

    monkeypatch.setattr(ga_optimization, "credentials_available", lambda: True)
    monkeypatch.setattr(ga_optimization, "load_codon_weights", lambda path: {})
    monkeypatch.setattr(ga_optimization, "load_restriction_sites", lambda path: ())
    monkeypatch.setattr(ga_optimization, "_adaptive_refine_payload", fake_adaptive)
    result = ga_optimization.refine_construct_table(
        inputs,
        tmp_path / "codons.csv",
        tmp_path / "sites.csv",
        tmp_path / "out",
        use_idt=True,
        adaptive_copy_search_enabled=True,
    )
    assert len(result) == 2
    assert "known_orderable_copies" not in observed[0]
    assert observed[1]["known_orderable_copies"] == 3
    assert observed[1]["known_orderable_dna_sequence"] == "ATG" * 3


def test_idt_response_is_never_called_passed_without_explicit_evidence():
    assert summarize_complexity_response({"score": 0.2})["idt_status"] == "scored_unclassified"
    assert summarize_complexity_response({"passed": True})["idt_status"] == "passed"
    assert summarize_complexity_response({"isComplex": True})["idt_status"] == "failed"


def test_idt_response_matches_named_sequence_before_classification():
    response = {"results": [{"Name": "a", "passed": True}, {"Name": "b", "passed": False}]}
    summary = summarize_complexity_response(response, name="b")
    assert summary["idt_status"] == "failed"
    assert summary["idt_result_matched_by_name"] is True


def test_idt_gblocks_response_uses_ordered_rule_lists():
    response = [
        [
            {"Name": "Overall Repeat", "IsViolated": False, "Score": 1.2},
            {"Name": "GC content", "IsViolated": False, "Score": 0.0},
        ],
        [
            {"Name": "Overall Repeat", "IsViolated": True, "Score": 23.5},
            {"Name": "GC content", "IsViolated": False, "Score": 0.0},
        ],
    ]
    passed = summarize_complexity_response(response, name="first", sequence_index=0)
    failed = summarize_complexity_response(response, name="second", sequence_index=1)
    assert passed["idt_status"] == "passed"
    assert passed["idt_violation_count"] == 0
    assert passed["idt_result_selected_by_index"] is True
    assert failed["idt_status"] == "failed"
    assert failed["idt_violation_count"] == 1
    assert '"Overall Repeat"' in failed["idt_violation_names_json"]
    assert '"actual_value"' in failed["idt_rule_details_json"]


def test_idt_rule_geometry_is_retained_for_actionable_ga_feedback():
    summary = summarize_complexity_response(
        [[{
            "Name": "Repeat Length (Fragment)",
            "IsViolated": True,
            "Score": 7.0,
            "ActualValue": 23.0,
            "DisplayText": "repeat at two locations",
            "RepeatedSegment": "AACCGGTTAACCGGTTAACCGGT",
            "ForwardLocations": [100, 352],
            "ReverseLocations": [710],
            "StartIndex": 100,
            "TerminalEnd": 5,
            "ThresholdOutput": {
                "ThresholdType": 1,
                "Value": 13,
                "WindowLength": 90,
            },
        }]],
        sequence_index=0,
    )
    details = json.loads(summary["idt_rule_details_json"])
    assert details[0]["display_text"] == "repeat at two locations"
    assert details[0]["repeated_segment"] == "AACCGGTTAACCGGTTAACCGGT"
    assert details[0]["forward_locations"] == [100, 352]
    assert details[0]["reverse_locations"] == [710]
    assert details[0]["terminal_end"] == 5
    assert details[0]["threshold_value"] == 13
    assert details[0]["threshold_window_length"] == 90


def test_idt_empty_position_matched_rule_list_is_orderable():
    summary = summarize_complexity_response([[]], name="orderable", sequence_index=0)
    assert summary["idt_status"] == "passed"
    assert summary["idt_explicit_pass"] is True
    assert summary["idt_violation_count"] == 0
    assert summary["idt_result_selected_by_index"] is True


def test_idt_score_sum_policy_ignores_violation_boolean_for_acceptance():
    summary = summarize_complexity_response(
        [[{"Name": "diagnostic", "IsViolated": True, "Score": 9.5}]],
        sequence_index=0,
    )
    assert summary["idt_status"] == "passed"
    assert summary["idt_violation_count"] == 1
    assert summary["idt_complexity_score"] == pytest.approx(9.5)


def test_idt_score_sum_fails_strictly_at_ten_and_missing_score_is_unclassified():
    exact_ten = summarize_complexity_response(
        [[
            {"Name": "repeat", "IsViolated": False, "Score": 6.0},
            {"Name": "hairpin", "IsViolated": False, "Score": 4.0},
        ]],
        sequence_index=0,
    )
    missing = summarize_complexity_response(
        [[{"Name": "missing", "IsViolated": False}]], sequence_index=0
    )
    assert exact_ten["idt_status"] == "failed"
    assert exact_ten["idt_complexity_score"] == 10.0
    assert missing["idt_status"] == "scored_unclassified"
    assert missing["idt_score_complete"] is False


def test_idt_positive_nonviolated_rule_still_reweights_ga():
    summary = {
        "idt_rule_details_json": (
            '[{"name":"Overall Repeat","is_violated":false,"score":2.0}]'
        )
    }
    adjusted, changes = adjust_ga_score_profile_from_idt(GA_SCORE_PROFILE, summary)
    assert adjusted["repeated_14mer"] > GA_SCORE_PROFILE["repeated_14mer"]
    assert changes


def test_idt_rejection_reweights_corresponding_ga_components():
    summary = {
        "idt_rule_details_json": (
            '[{"name":"Hairpin Length","is_violated":true,"score":3.0,'
            '"actual_value":13},{"name":"Overall Repeat","is_violated":true,'
            '"score":23.5,"actual_value":98.0}]'
        )
    }
    adjusted, changes = adjust_ga_score_profile_from_idt(GA_SCORE_PROFILE, summary)
    assert adjusted["hairpin_10mer_proxy"] > GA_SCORE_PROFILE["hairpin_10mer_proxy"]
    assert adjusted["repeated_14mer"] > GA_SCORE_PROFILE["repeated_14mer"]
    assert adjusted["gc_window_soft_violation"] == GA_SCORE_PROFILE["gc_window_soft_violation"]
    assert {item["score_component"] for item in changes} >= {
        "hairpin_10mer_proxy",
        "repeated_14mer",
    }


def test_adaptive_refinement_uses_idt_gate_and_feedback(monkeypatch):
    def construct_metrics(unit, copies, payload, codon_weights, **kwargs):
        return {"dna_sequence": "ATG" * copies}

    def refine(dna, **kwargs):
        copies = len(dna) // 3
        return dna, {
            "ga_local_constraints_passed": True,
            "repeated_re_site_excess": 0,
            "selected_re_site_excess": 0,
            "ga_score": float(copies),
            "ga_generations": kwargs["generations"],
        }

    class FakeScorer:
        def __init__(self):
            self.calls = []

        def score(self, name, sequence):
            fields = dict(item.split("=", 1) for item in name.split("|")[1:])
            copies = int(fields["copies"])
            generations = int(fields["generations"])
            self.calls.append((copies, generations))
            passed = copies <= 3 or (copies == 4 and generations >= 40)
            details = [] if passed else [{
                "name": "Overall Repeat",
                "is_violated": True,
                "score": 20.0,
                "actual_value": 95.0,
            }]
            return {
                "idt_api_called": True,
                "idt_cache_hit": False,
                "idt_status": "passed" if passed else "failed",
                "idt_explicit_pass": passed,
                "idt_violation_count": 0 if passed else 1,
                "idt_violation_names_json": "[]" if passed else '["Overall Repeat"]',
                "idt_rule_scores_json": "{}",
                "idt_rule_details_json": __import__("json").dumps(details),
                "idt_scored_sequence_sha256": "test-sha",
            }

    monkeypatch.setattr(ga_optimization, "_construct_metrics", construct_metrics)
    monkeypatch.setattr(ga_optimization, "genetic_refine_dna", refine)
    scorer = FakeScorer()
    payload = {
        "module_id": "test",
        "fragment_limit_bp": 1800,
        "unit_sequence": "M",
        "verified_max_copies": 3,
        "mathematical_max_copies": 6,
        "dna_sequence": "ATG" * 3,
        "site_i_position": 0,
        "site_ii_position": 0,
        "site_i_recognition_site": "",
        "site_ii_recognition_site": "",
        "site_iii_sites": "",
    }
    result = ga_optimization._adaptive_refine_payload(
        payload,
        codon_weights={"ATG": 1.0},
        recognition_sites=(),
        seed=42,
        population_size=4,
        short_generations=10,
        generation_schedule=(10, 20, 40, 60, 80, 100),
        idt_scorer=scorer,
        require_idt_orderable=True,
    )
    assert result["verified_max_copies"] == 4
    assert result["idt_status"] == "passed"
    assert result["adaptive_orderable_passed"] is True
    assert (4, 40) in scorer.calls
    assert (5, 100) in scorer.calls
    trace = __import__("json").loads(result["adaptive_search_trace_json"])
    rejected = [item for item in trace if item["idt_status"] == "failed"]
    assert rejected
    assert any(item["idt_feedback_adjustments_json"] != "[]" for item in rejected)


def test_repeated_nonselected_re_sites_are_soft_before_idt(monkeypatch):
    def construct_metrics(unit, copies, payload, codon_weights, **kwargs):
        return {"dna_sequence": "ATG" * copies}

    def refine(dna, **kwargs):
        return dna, {
            "ga_local_constraints_passed": True,
            "repeated_re_site_excess": 3,
            "selected_re_site_excess": 0,
            "ga_score": 30_000.0,
            "ga_generations": kwargs["generations"],
        }

    class PassingScorer:
        def __init__(self):
            self.calls = []

        def score(self, name, sequence):
            self.calls.append((name, sequence))
            return {
                "idt_api_called": True,
                "idt_cache_hit": False,
                "idt_status": "passed",
                "idt_explicit_pass": True,
                "idt_violation_count": 0,
                "idt_violation_names_json": "[]",
                "idt_rule_scores_json": "{}",
                "idt_rule_details_json": "[]",
                "idt_scored_sequence_sha256": "test-sha",
            }

    monkeypatch.setattr(ga_optimization, "_construct_metrics", construct_metrics)
    monkeypatch.setattr(ga_optimization, "genetic_refine_dna", refine)
    scorer = PassingScorer()
    result = ga_optimization._adaptive_refine_payload(
        {
            "module_id": "soft-re-sites",
            "fragment_limit_bp": 1800,
            "unit_sequence": "MMM",
            "verified_max_copies": 1,
            "mathematical_max_copies": 2,
            "dna_sequence": "ATG" * 3,
            "site_i_position": 0,
            "site_ii_position": 0,
            "site_i_recognition_site": "",
            "site_ii_recognition_site": "",
            "site_iii_sites": "",
        },
        codon_weights={"ATG": 1.0},
        recognition_sites=("ATG",),
        seed=42,
        population_size=4,
        short_generations=10,
        generation_schedule=(10, 20, 40, 60, 80, 100),
        idt_scorer=scorer,
        require_idt_orderable=True,
    )
    assert result["verified_max_copies"] == 2
    assert result["adaptive_orderable_passed"] is True
    assert result["repeated_re_site_excess"] == 3
    assert result["ga_re_site_policy"] == GA_RE_SITE_POLICY
    assert scorer.calls


def test_idt_cache_records_request_hash_but_not_sequence(tmp_path):
    destination = tmp_path / "response.json"
    write_cached_response(
        {"results": [{"Name": "test", "passed": True}]},
        destination,
        sequences=[{"Name": "test", "Sequence": "ATGGCT"}],
    )
    cached = destination.read_text()
    assert "ATGGCT" not in cached
    assert '"length_bp": 6' in cached
    assert '"sha256"' in cached
