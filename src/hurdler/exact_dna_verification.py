"""Independent molecular verification for exact-DNA HURDLER routes.

The route planner searches precursor graphs in :mod:`hurdler.dna_assembly`.
This module deliberately does not import or call that planner.  It rebuilds
every insert from the declared purchased molecules, independently simulates
the circular vector cut arc, and verifies each intermediate sequence from the
physical fragment order.
"""

from __future__ import annotations

import json
from typing import Any, Mapping, Sequence

from Bio import Restriction

from .dna_assembly import EnzymeGeometry, reverse_complement
from .plasmid_reference import PlasmidReferenceDatabase, VectorCutScheme


def _circular_slice(sequence: str, start: int, end: int) -> str:
    length = len(sequence)
    start %= length
    end %= length
    if start == end:
        return sequence
    return sequence[start:end] if start < end else sequence[start:] + sequence[:end]


def _occurrences(sequence: str, motif: str, *, circular: bool) -> set[tuple[int, str]]:
    motif = motif.upper()
    extended = sequence + sequence[: max(0, len(motif) - 1)] if circular else sequence
    result: set[tuple[int, str]] = set()
    for oriented, strand in {(motif, "+"), (reverse_complement(motif), "-")}:
        for start in range(len(sequence)):
            if extended.startswith(oriented, start):
                result.add((start, strand))
    return result


def _vector_overhang(oriented: str, top_cut: int, bottom_cut: int) -> str:
    lower, upper = sorted((int(top_cut), int(bottom_cut)))
    return _circular_slice(oriented, lower, upper)


def _verify_fragment(
    fragment: Mapping[str, Any],
    *,
    expected_left_overhang: str,
    expected_right_overhang: str,
    geometries: Mapping[str, EnzymeGeometry],
) -> list[str]:
    errors: list[str] = []
    core = str(fragment.get("core_sequence", ""))
    left = str(fragment.get("left_digest_overhang", ""))
    right = str(fragment.get("right_digest_overhang", ""))
    if left != expected_left_overhang or right != expected_right_overhang:
        errors.append("insert cohesive ends do not match the selected recipient sites")
    product = str(fragment.get("product_type", ""))
    if product in {"annealed_sticky_end_primer_pair", "duplexed_seed_oligo_pair"}:
        forward = str(fragment.get("primer_forward_5to3") or fragment.get("purchase_sequence", ""))
        reverse = str(
            fragment.get("primer_reverse_5to3")
            or fragment.get("secondary_purchase_sequence", "")
        )
        if forward != left + core:
            errors.append("forward insert oligo does not encode its declared cohesive end and core")
        if reverse != right + reverse_complement(core):
            errors.append("reverse insert oligo does not encode its declared cohesive end and core")
    else:
        purchase = str(fragment.get("purchase_sequence", ""))
        left_adapter = str(fragment.get("left_adapter_sequence", ""))
        right_adapter = str(fragment.get("right_adapter_sequence", ""))
        if purchase != left_adapter + core + right_adapter:
            errors.append("purchased block does not equal left adapter + core + right adapter")
        for side in ("left", "right"):
            enzyme = str(fragment.get(f"{side}_adapter_enzyme", ""))
            adapter = str(fragment.get(f"{side}_adapter_sequence", ""))
            geometry = geometries.get(enzyme)
            if geometry is None or not geometry.site_iii_eligible:
                errors.append(f"{side} adapter lacks a supported Site-III enzyme")
            elif not _occurrences(adapter, geometry.recognition_site, circular=False):
                errors.append(f"{side} adapter lacks the declared {enzyme} recognition site")
    return errors


def _mapped_position(
    position: int,
    *,
    replacement_start: int,
    replacement_end: int,
    assembled_length: int,
    direction: str,
) -> int | None:
    if position < replacement_start:
        return position
    if position >= replacement_end:
        return position - (replacement_end - replacement_start) + assembled_length
    if direction == "left_to_right":
        return position if position < replacement_start + assembled_length else None
    source_start = replacement_end - assembled_length
    return replacement_start + position - source_start if position >= source_start else None


def _activate_declared_latent_sites(
    exact_intermediate: str,
    exact_result: str,
    edge: Mapping[str, Any],
    *,
    assembled_length: int,
    geometries: Mapping[str, EnzymeGeometry],
) -> tuple[str, list[dict[str, Any]], list[str]]:
    """Activate one-base latent sites from geometry, not planner step hashes."""
    mutable = list(exact_intermediate)
    evidence: list[dict[str, Any]] = []
    errors: list[str] = []
    replacement_start = int(edge["replacement_start"])
    replacement_end = int(edge["replacement_end"])
    direction = str(edge["direction"])
    for role in ("site_i", "site_ii"):
        enzyme = str(edge[f"{role}_enzyme"])
        state = str(edge[f"{role}_state"])
        geometry = geometries[enzyme]
        orientation = str(edge[f"{role}_orientation"])
        motif = (
            geometry.recognition_site
            if orientation == "+"
            else reverse_complement(geometry.recognition_site)
        )
        start = int(edge[f"{role}_position"])
        observed = exact_result[start : start + len(motif)]
        differences = [
            index for index, (actual, active) in enumerate(zip(observed, motif))
            if actual != active
        ]
        expected = 0 if state == "active" else 1
        if len(differences) != expected:
            errors.append(f"{enzyme} {state} declaration disagrees with the exact result")
            continue
        if not differences:
            continue
        mismatch = start + differences[0]
        mapped = _mapped_position(
            mismatch,
            replacement_start=replacement_start,
            replacement_end=replacement_end,
            assembled_length=assembled_length,
            direction=direction,
        )
        if mapped is None:
            continue
        mutable[mapped] = motif[differences[0]]
        evidence.append(
            {
                "enzyme": enzyme,
                "role": "Site I" if role == "site_i" else "Site II",
                "orientation": orientation,
                "state": "latent_temporarily_activated",
                "mismatch_position_in_result_0based": mismatch,
                "mismatch_position_in_intermediate_0based": mapped,
                "target_base": observed[differences[0]],
                "temporary_active_base": motif[differences[0]],
            }
        )
    return "".join(mutable), evidence, errors


def _protected_feature_errors(
    database: PlasmidReferenceDatabase,
    scheme: VectorCutScheme,
    final_sequence: str,
) -> list[str]:
    profile = database.profile(scheme.profile_id)
    reference = database.reference(profile.reference_id)
    circular = final_sequence + final_sequence
    errors: list[str] = []
    for feature in reference.features:
        if not feature.protected:
            continue
        if any(
            start < profile.mcs_end and end > profile.mcs_start
            for start, end in feature.intervals
        ):
            continue
        for start, end in feature.intervals:
            segment = reference.sequence[start:end]
            if segment and segment not in circular and reverse_complement(segment) not in circular:
                errors.append(f"missing protected feature {feature.feature_class}:{feature.label}")
                break
    return errors


def verify_exact_dna_assembly(
    *,
    query: Mapping[str, Any],
    target_sequence: str,
    route: Mapping[str, Any],
    primary_fragment: Mapping[str, Any],
    database: PlasmidReferenceDatabase,
    geometries: Mapping[str, EnzymeGeometry],
) -> dict[str, Any]:
    """Reconstruct and verify every physical cloning step independently."""
    errors: list[str] = []
    scheme = next(
        (item for item in database.schemes if item.scheme_id == route["scheme_id"]),
        None,
    )
    if scheme is None or scheme.left_cutter is None or scheme.right_cutter is None:
        return {"passed": False, "errors": ["selected vector cut scheme is unavailable"], "steps": []}
    profile = database.profile(str(route["profile_id"]))
    reference = database.reference(profile.reference_id)
    oriented = (
        reference.sequence
        if profile.expression_strand == 1
        else reverse_complement(reference.sequence)
    )
    left_cut = int(scheme.left_cutter.top_cut_oriented)
    right_cut = int(scheme.right_cutter.top_cut_oriented)
    retained = _circular_slice(oriented, right_cut, left_cut)
    removed = _circular_slice(oriented, left_cut, right_cut)
    if retained != scheme.retained_backbone_sequence:
        errors.append("independent circular digest disagrees with the retained backbone")
    if removed != scheme.removed_arc_sequence:
        errors.append("independent circular digest disagrees with the removed MCS arc")
    for cutter in (scheme.left_cutter, scheme.right_cutter):
        observed = _circular_slice(oriented, cutter.oriented_start, cutter.oriented_end)
        if observed not in {
            cutter.recognition_site,
            reverse_complement(cutter.recognition_site),
        }:
            errors.append(f"{cutter.canonical_enzyme} is absent at its declared vector cut site")

    left_vector_overhang = _vector_overhang(
        oriented, scheme.left_cutter.top_cut_oriented, scheme.left_cutter.bottom_cut_oriented
    )
    right_vector_overhang = _vector_overhang(
        oriented, scheme.right_cutter.top_cut_oriented, scheme.right_cutter.bottom_cut_oriented
    )
    expected_primary = (
        scheme.left_restoration_sequence
        + str(route["seed"]["seed_sequence"])
        + scheme.right_restoration_sequence
    )
    if str(primary_fragment.get("core_sequence", "")) != expected_primary:
        errors.append("primary fragment does not restore the vector and install the exact seed")
    errors.extend(
        _verify_fragment(
            primary_fragment,
            expected_left_overhang=left_vector_overhang,
            expected_right_overhang=right_vector_overhang,
            geometries=geometries,
        )
    )

    current_insert = str(route["seed"]["seed_sequence"])
    steps: list[dict[str, Any]] = []
    primary_plasmid = (
        retained + scheme.left_restoration_sequence + current_insert + scheme.right_restoration_sequence
    )
    steps.append(
        {
            "step": 1,
            "stage": "primary_seed_installation",
            "input_plasmid": "step00_plasmid.gb",
            "insert_name": "step01_insert",
            "insert_sequence": str(primary_fragment.get("digest_fragment_sequence", expected_primary)),
            "insert_core_sequence": expected_primary,
            "purchase_fragment_ids": [str(primary_fragment.get("fragment_id", "primary_seed"))],
            "enzymes": [
                scheme.left_cutter.canonical_enzyme,
                scheme.right_cutter.canonical_enzyme,
            ],
            "result_insert_sequence": current_insert,
            "plasmid_sequence": primary_plasmid,
            "output_plasmid": "step01_plasmid.gb",
            "latent_sites": [],
        }
    )

    step_number = 1
    for edge_index, edge in enumerate(route["edges"], start=1):
        start = int(edge["replacement_start"])
        end = int(edge["replacement_end"])
        if str(edge.get("precursor_sequence", "")) != current_insert:
            errors.append(f"edge {edge_index} precursor differs from independently reconstructed insert")
        fragments = sorted(edge["fragments"], key=lambda row: int(row["target_start"]))
        if not fragments:
            errors.append(f"edge {edge_index} has no physical donor fragments")
            continue
        cursor = start
        for fragment in fragments:
            if int(fragment["target_start"]) != cursor:
                errors.append(f"edge {edge_index} donor fragments do not cover one contiguous interval")
            cursor = int(fragment["target_end"])
        if cursor != end:
            errors.append(f"edge {edge_index} donor fragments do not reach the replacement end")
        full_core = "".join(str(row["core_sequence"]) for row in fragments)
        exact_result = current_insert[:start] + full_core + current_insert[start:]
        ordered = fragments if edge["direction"] == "left_to_right" else list(reversed(fragments))
        assembled_parts: list[str] = []
        left_role = "site_i" if int(edge["site_i_position"]) <= int(edge["site_ii_position"]) else "site_ii"
        right_role = "site_ii" if left_role == "site_i" else "site_i"
        for local_index, fragment in enumerate(ordered, start=1):
            core = str(fragment["core_sequence"])
            if edge["direction"] == "left_to_right":
                assembled_parts.append(core)
            else:
                assembled_parts.insert(0, core)
            assembled = "".join(assembled_parts)
            exact_intermediate = current_insert[:start] + assembled + current_insert[start:]
            is_edge_final = local_index == len(ordered)
            latent_evidence: list[dict[str, Any]] = []
            if is_edge_final:
                intermediate = exact_result
            else:
                intermediate, latent_evidence, activation_errors = _activate_declared_latent_sites(
                    exact_intermediate,
                    exact_result,
                    edge,
                    assembled_length=len(assembled),
                    geometries=geometries,
                )
                errors.extend(f"edge {edge_index}: {item}" for item in activation_errors)
            expected_left = str(edge[f"{left_role}_overhang"])
            expected_right = str(edge[f"{right_role}_overhang"])
            errors.extend(
                f"edge {edge_index} fragment {local_index}: {item}"
                for item in _verify_fragment(
                    fragment,
                    expected_left_overhang=expected_left,
                    expected_right_overhang=expected_right,
                    geometries=geometries,
                )
            )
            step_number += 1
            plasmid = (
                retained
                + scheme.left_restoration_sequence
                + intermediate
                + scheme.right_restoration_sequence
            )
            steps.append(
                {
                    "step": step_number,
                    "stage": "hurdler_growth",
                    "edge_index": edge_index,
                    "edge_local_step": local_index,
                    "input_plasmid": f"step{step_number - 1:02d}_plasmid.gb",
                    "insert_name": f"step{step_number:02d}_insert",
                    "insert_sequence": str(fragment.get("digest_fragment_sequence", "")),
                    "insert_core_sequence": core,
                    "purchase_fragment_ids": [str(fragment["fragment_id"])],
                    "enzymes": [str(edge["site_i_enzyme"]), str(edge["site_ii_enzyme"])],
                    "site_iii_enzymes": [
                        value for value in (
                            str(fragment.get("left_adapter_enzyme", "")),
                            str(fragment.get("right_adapter_enzyme", "")),
                        ) if value and value != "not_required"
                    ],
                    "result_insert_sequence": intermediate,
                    "plasmid_sequence": plasmid,
                    "output_plasmid": f"step{step_number:02d}_plasmid.gb",
                    "latent_sites": latent_evidence,
                }
            )
        current_insert = exact_result

    if current_insert != target_sequence:
        errors.append("final independently assembled insert differs from the exact target")
    final_plasmid = (
        retained + scheme.left_restoration_sequence + current_insert + scheme.right_restoration_sequence
    )
    errors.extend(_protected_feature_errors(database, scheme, final_plasmid))
    selected_sites: dict[str, int] = {}
    for pair in route["pairs"]:
        for role in ("site_i", "site_ii"):
            enzyme = str(pair[f"{role}_enzyme"])
            motif = str(pair[f"{role}_recognition_site"])
            excess = len(_occurrences(final_plasmid, motif, circular=True)) - len(
                _occurrences(target_sequence, motif, circular=False)
            )
            selected_sites[enzyme] = max(selected_sites.get(enzyme, 0), excess)
    if any(value > 0 for value in selected_sites.values()):
        errors.append("final plasmid contains an excess selected-pair restriction site")
    return {
        "passed": not errors,
        "errors": errors,
        "verifier": "independent-circular-digest-ligation-v1",
        "initial_plasmid_sequence": oriented[right_cut:] + oriented[:right_cut],
        "rotation_origin_0based": right_cut,
        "retained_backbone_sequence": retained,
        "removed_arc_sequence": removed,
        "left_vector_overhang": left_vector_overhang,
        "right_vector_overhang": right_vector_overhang,
        "steps": steps,
        "final_insert_sequence": current_insert,
        "final_plasmid_sequence": final_plasmid,
        "selected_pair_excess_sites": selected_sites,
        "protected_features_preserved": not any(
            item.startswith("missing protected feature") for item in errors
        ),
        "query_input_mode": str(query.get("input_mode", "")),
    }
