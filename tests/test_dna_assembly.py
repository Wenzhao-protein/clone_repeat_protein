from __future__ import annotations

import hashlib

import pandas as pd
import pytest

from hurdler.dna_assembly import (
    EnzymeGeometry,
    TargetRecord,
    build_synthetic_factorial,
    enumerate_active_latent_pairs,
    expand_element_inventory,
    is_exact_repeat_gain_route,
    plan_target,
    reverse_complement,
    scan_re_sites,
)


def _geometry(
    enzyme: str,
    site: str,
    ovhgseq: str,
    top: int,
    bottom: int,
    *,
    site_i: bool = False,
    site_ii: bool = False,
    site_iii: bool = False,
) -> EnzymeGeometry:
    return EnzymeGeometry(
        enzyme=enzyme,
        canonical_enzyme=enzyme,
        recognition_site=site,
        ovhg=-len(ovhgseq),
        ovhgseq=ovhgseq,
        top_cut_offset=top,
        bottom_cut_offset=bottom,
        elucidate="fixture",
        is_type_iis=site_iii,
        site_i_eligible=site_i,
        site_ii_eligible=site_ii,
        site_iii_eligible=site_iii,
        methylation_compatible=True,
        ligation_ok=True,
        no_star_activity=True,
    )


def _catalog():
    return {
        "EcoRI": _geometry("EcoRI", "GAATTC", "AATT", 1, 5, site_i=True, site_ii=True),
        "HindIII": _geometry("HindIII", "AAGCTT", "AGCT", 1, 5, site_i=True, site_ii=True),
        "BsaI": _geometry("BsaI", "GGTCTC", "NNNN", 7, 11, site_iii=True),
    }


def _plasmids():
    return pd.DataFrame(
        {"pUC18": [True, True]},
        index=["EcoRI", "HindIII"],
    )


def _active_latent_target():
    # The first base of the terminal GAGCTT is the single-base difference
    # from HindIII AAGCTT, and lies between the inward EcoRI/HindIII cuts.
    return "GAATTC" + "ATGCGT" * 20 + "GAGCTT"


def test_active_and_one_base_latent_hits_have_explicit_cut_geometry():
    sequence = _active_latent_target()
    hits = scan_re_sites("target", sequence, _catalog())
    eco = next(hit for hit in hits if hit.enzyme == "EcoRI" and hit.state == "active")
    hind = next(hit for hit in hits if hit.enzyme == "HindIII" and hit.state == "latent")
    assert (eco.top_cut, eco.bottom_cut) == (1, 5)
    assert hind.mismatch_position == len(sequence) - 6
    assert hind.target_base == "G" and hind.active_base == "A"


def test_active_active_is_not_this_method_and_two_mismatches_are_not_latent():
    catalog = _catalog()
    active_sequence = "GAATTC" + "ATGCGT" * 20 + "AAGCTT"
    active_hits = scan_re_sites("active", active_sequence, catalog)
    assert enumerate_active_latent_pairs(active_hits, _plasmids()) == []

    two_mismatch = "GAATTC" + "ATGCGT" * 20 + "GAGCTA"
    hits = scan_re_sites("two", two_mismatch, catalog)
    assert not any(hit.enzyme == "HindIII" and hit.start == len(two_mismatch) - 6 for hit in hits)


def test_latent_mutation_is_donor_derived_between_the_inward_cuts():
    hits = scan_re_sites("target", _active_latent_target(), _catalog())
    pairs = enumerate_active_latent_pairs(hits, _plasmids())
    pair = next(
        item
        for item in pairs
        if item["site_i"].enzyme == "EcoRI" and item["site_ii"].enzyme == "HindIII"
    )
    assert pair["mechanism"] == "active+latent"
    assert pair["mutation_between_overhangs"] is True
    mutation = pair["donor_derived_mutation_positions"][0]
    assert pair["replacement_start"] <= mutation < pair["replacement_end"]


def test_latent_latent_pair_is_supported_but_active_active_is_not():
    sequence = "GAATTA" + "ATGCGT" * 20 + "GAGCTT"
    hits = scan_re_sites("both-latent", sequence, _catalog())
    pairs = enumerate_active_latent_pairs(hits, _plasmids())
    pair = next(
        item
        for item in pairs
        if item["site_i"].enzyme == "EcoRI" and item["site_ii"].enzyme == "HindIII"
    )
    assert pair["mechanism"] == "latent+latent"
    assert len(pair["donor_derived_mutation_positions"]) == 2


def test_nonpalindromic_reverse_strand_hit_is_reported():
    bau = _geometry("BauI", "CACGAG", "ACGA", 1, 5, site_i=True)
    hits = scan_re_sites("reverse", "TTTCTCGTGAAA", {"BauI": bau})
    hit = next(item for item in hits if item.state == "active")
    assert hit.orientation == "-"
    assert hit.observed == "CTCGTG"
    assert hit.top_cut != hit.bottom_cut


def test_extra_selected_enzyme_occurrence_rejects_route():
    sequence = "GAATTC" + "ATGCGT" * 10 + "GAATTC" + "ATGCGT" * 10 + "GAGCTT"
    result = plan_target(
        TargetRecord("extra-site", sequence), _catalog(), _plasmids(), require_idt=False
    )
    assert not result["summary"].iloc[0].hurdler_compatible


def test_planner_preserves_exact_final_target_and_retains_passing_route():
    target = TargetRecord("exact", _active_latent_target(), cohort="real_construct")
    result = plan_target(
        target,
        _catalog(),
        _plasmids(),
        require_idt=False,
    )
    assert result["summary"].iloc[0].hurdler_compatible
    passing = result["routes"].query("passed")
    assert not passing.empty
    assert passing.final_sequence_exact.all()
    expected = hashlib.sha256(target.sequence.encode()).hexdigest()
    assert set(passing.final_sequence_sha256) == {expected}
    assert set(passing.target_sequence_sha256) == {expected}


def test_route_identity_changes_with_fragmentation_ceiling():
    target = TargetRecord("route-identity", _active_latent_target())
    wide = plan_target(
        target, _catalog(), _plasmids(), require_idt=False,
        max_purchase_bp=3000,
    )
    narrow = plan_target(
        target, _catalog(), _plasmids(), require_idt=False,
        max_purchase_bp=60,
    )
    assert set(wide["routes"].route_id).isdisjoint(set(narrow["routes"].route_id))


def test_under_90bp_primer_fragments_do_not_require_type_iis_adapter():
    target = TargetRecord("primer-without-type-iis", _active_latent_target())
    catalog_without_type_iis = {
        name: geometry
        for name, geometry in _catalog().items()
        if not geometry.site_iii_eligible
    }
    result = plan_target(
        target,
        catalog_without_type_iis,
        _plasmids(),
        require_idt=False,
        max_purchase_bp=60,
    )
    passing = result["routes"].query("passed")
    assert not passing.empty
    assert passing.site_iii_left.eq("not_required_primer_pair").all()
    fragments = result["fragments"].loc[
        result["fragments"].route_id.isin(passing.route_id)
    ]
    assert not fragments.empty
    assert fragments.core_length_bp.lt(90).all()
    assert fragments.product_type.eq("annealed_sticky_end_primer_pair").all()
    assert fragments.idt_status.eq(
        "not_applicable_primer_pair_under_90bp"
    ).all()


def test_minimum_replacement_filters_formal_but_nonexpanding_routes():
    target = TargetRecord("minimum-insert", _active_latent_target())
    baseline = plan_target(target, _catalog(), _plasmids(), require_idt=False)
    replacement_bp = int(baseline["routes"].replacement_length_bp.max())
    meaningful = plan_target(
        target,
        _catalog(),
        _plasmids(),
        require_idt=False,
        min_replacement_bp=replacement_bp,
    )
    assert meaningful["summary"].iloc[0].hurdler_compatible
    assert meaningful["routes"].replacement_length_bp.ge(replacement_bp).all()

    impossible = plan_target(
        target,
        _catalog(),
        _plasmids(),
        require_idt=False,
        min_replacement_bp=replacement_bp + 1,
    )
    assert not impossible["summary"].iloc[0].hurdler_compatible
    assert impossible["routes"].empty


def test_exact_repeat_gain_requires_the_exact_shorter_precursor():
    unit = "ACGTTG"
    target = TargetRecord(
        "repeat-gain", unit * 4, unit_sequence=unit, copy_count=4
    )
    assert is_exact_repeat_gain_route(target, 6, 12, repeat_unit_gain=1)
    assert is_exact_repeat_gain_route(target, 3, 9, repeat_unit_gain=1)
    assert not is_exact_repeat_gain_route(target, 3, 8, repeat_unit_gain=1)

    missing_metadata = TargetRecord("missing-unit", unit * 4)
    with pytest.raises(ValueError, match="unit_sequence and copy_count"):
        is_exact_repeat_gain_route(missing_metadata, 6, 12)


def test_idt_rejection_falls_back_to_unscored_under_90bp_primer_pairs():
    class RejectingScorer:
        def score(self, name, sequence):
            return {
                "idt_status": "failed",
                "idt_explicit_pass": False,
                "idt_complexity_score": 10.0,
                "idt_response_sha256": "response",
                "idt_scored_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            }

    result = plan_target(
        TargetRecord("rejected", _active_latent_target()),
        _catalog(),
        _plasmids(),
        idt_scorer=RejectingScorer(),
        require_idt=True,
    )
    assert result["summary"].iloc[0].hurdler_compatible
    assert result["summary"].iloc[0].breakpoint_attempts_json == "[200, 60]"
    accepted = result["fragments"].loc[
        result["fragments"].route_id.isin(result["routes"].query("passed").route_id)
    ]
    assert accepted.product_type.eq("annealed_sticky_end_primer_pair").all()
    assert accepted.idt_status.eq("not_applicable_primer_pair_under_90bp").all()
    assert accepted.idt_score.isna().all()
    assert accepted.purchase_sequence_count.eq(2).all()
    assert accepted.primer_pair_exposes_sticky_ends.all()
    for row in accepted.itertuples(index=False):
        assert row.primer_forward_5to3.startswith(row.left_digest_overhang)
        assert row.primer_reverse_5to3.startswith(row.right_digest_overhang)
        forward_annealing_region = row.primer_forward_5to3[
            len(row.left_digest_overhang) :
        ]
        reverse_annealing_region = row.primer_reverse_5to3[
            len(row.right_digest_overhang) :
        ]
        assert reverse_annealing_region == reverse_complement(forward_annealing_region)


def test_idt_rejection_retries_shorter_breakpoints_until_fragments_pass():
    class LengthAwareScorer:
        def __init__(self):
            self.calls = []

        def score(self, name, sequence):
            self.calls.append((name, sequence))
            is_whole = name.endswith("|whole_target")
            score = 20.0 if is_whole or len(sequence) > 55 else 0.0
            return {
                "idt_status": "passed" if score < 10 else "failed",
                "idt_explicit_pass": score < 10,
                "idt_complexity_score": score,
                "idt_response_sha256": hashlib.sha256(f"response:{name}:{score}".encode()).hexdigest(),
                "idt_scored_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            }

    target = TargetRecord("rescued", _active_latent_target())
    scorer = LengthAwareScorer()
    result = plan_target(
        target,
        _catalog(),
        _plasmids(),
        idt_scorer=scorer,
        require_idt=True,
    )
    summary = result["summary"].iloc[0]
    assert summary.hurdler_compatible
    assert summary.fragment_rescued_by_hurdler
    assert summary.breakpoint_attempts_json == "[200, 60]"
    assert result["routes"].query("passed").maximum_idt_score.isna().all()
    assert result["routes"].query("passed").final_sequence_exact.all()
    accepted = result["fragments"].loc[
        result["fragments"].route_id.isin(result["routes"].query("passed").route_id)
    ]
    assert accepted.idt_status.eq("not_applicable_primer_pair_under_90bp").all()
    assert all(len(sequence) >= 90 for name, sequence in scorer.calls if not name.endswith("whole_target"))


def test_synthetic_factorial_is_exactly_900_seeded_cases():
    first = build_synthetic_factorial(seed=42)
    second = build_synthetic_factorial(seed=42)
    assert len(first) == first.target_id.nunique() == 900
    assert first.target_id.tolist() == second.target_id.tolist()
    assert first.sequence.tolist() == second.sequence.tolist()
    assert set(first.architecture) == {
        "exact_tandem",
        "fixed_spacer",
        "alternating_ab",
        "nonrepetitive_control",
    }


def test_whole_target_under_90bp_is_not_sent_to_idt():
    class MustNotBeCalled:
        def score(self, name, sequence):
            raise AssertionError("under-90-bp target must not be sent to IDT")

    short = TargetRecord("short", "GAATTC" + "ATGC" * 15 + "GAGCTT")
    assert len(short.sequence) < 90
    result = plan_target(
        short, _catalog(), _plasmids(), idt_scorer=MustNotBeCalled(), require_idt=True
    )
    summary = result["summary"].iloc[0]
    assert summary.whole_target_idt_status == "not_applicable_primer_pair_under_90bp"
    assert pd.isna(summary.whole_target_idt_score)


def test_primer_pair_cutoff_is_strictly_below_90bp():
    class PassingScorer:
        def __init__(self):
            self.calls = []

        def score(self, name, sequence):
            self.calls.append((name, sequence))
            return {
                "idt_status": "passed",
                "idt_explicit_pass": True,
                "idt_complexity_score": 0.0,
                "idt_response_sha256": "response",
                "idt_scored_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            }

    # With these fixture cut offsets, the replaceable donor core is ten bases
    # shorter than the complete target.  Thus the two targets exercise 89 bp
    # (primer pair, no fragment score) and exactly 90 bp (synthesis fragment,
    # score required), respectively.
    scorer_89 = PassingScorer()
    result_89 = plan_target(
        TargetRecord("core-89", "GAATTC" + "A" * 87 + "GAGCTT"),
        _catalog(),
        _plasmids(),
        idt_scorer=scorer_89,
        require_idt=True,
    )
    accepted_89 = result_89["fragments"].loc[
        result_89["fragments"].route_id.isin(result_89["routes"].query("passed").route_id)
    ]
    assert accepted_89.core_length_bp.eq(89).all()
    assert accepted_89.product_type.eq("annealed_sticky_end_primer_pair").all()
    assert all(name.endswith("|whole_target") for name, _ in scorer_89.calls)

    scorer_90 = PassingScorer()
    result_90 = plan_target(
        TargetRecord("core-90", "GAATTC" + "A" * 88 + "GAGCTT"),
        _catalog(),
        _plasmids(),
        idt_scorer=scorer_90,
        require_idt=True,
    )
    accepted_90 = result_90["fragments"].loc[
        result_90["fragments"].route_id.isin(result_90["routes"].query("passed").route_id)
    ]
    assert accepted_90.core_length_bp.eq(90).all()
    assert accepted_90.product_type.ne("annealed_sticky_end_primer_pair").all()
    assert any(not name.endswith("|whole_target") for name, _ in scorer_90.calls)


def test_public_elements_expand_to_five_exact_copy_counts():
    elements = pd.DataFrame(
        [
            {
                "element_id": "regulatory-one",
                "element_sequence": "AUGC",
                "source_database": "fixture",
                "source_url": "https://example.test/source",
            }
        ]
    )
    expanded = expand_element_inventory(elements)
    assert expanded.copy_count.tolist() == [2, 4, 8, 16, 32]
    assert expanded.sequence.tolist() == ["ATGC" * count for count in (2, 4, 8, 16, 32)]
    assert expanded.cohort.eq("real_element_derived").all()
