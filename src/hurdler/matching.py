"""Fast sequence matching against a :class:`~hurdler.index.PatternIndex`."""

from __future__ import annotations

from math import ceil

import pandas as pd

from .constants import PLASMIDS, validate_protein_sequence
from .index import PatternIndex, encode_pattern
from .rules import RuleProfile, LEGACY_OPTIMIZED_V1
from .schemas import MatchResult


def expand_short_module(module: str, minimum_length: int = 6) -> tuple[str, int]:
    normalized = validate_protein_sequence(module)
    copies = max(1, ceil(minimum_length / len(normalized)))
    return normalized * copies, copies


def match_module(
    module: str,
    plasmid: str,
    index: PatternIndex,
    *,
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
    expand_short: bool = True,
) -> MatchResult:
    """Return the first deterministic legacy-compatible match for a module."""
    normalized = validate_protein_sequence(module)
    if plasmid not in PLASMIDS:
        raise ValueError(f"Unknown plasmid: {plasmid!r}")
    if len(normalized) < 6 and expand_short:
        effective, expansion_copies = expand_short_module(normalized)
    else:
        effective, expansion_copies = normalized, 1
    module_length = len(effective)
    sequence = effective * 2 if rules.double_module_before_matching else effective
    plasmid_bit = PLASMIDS.index(plasmid)

    for left_position in range(len(sequence) - 2):
        left = sequence[left_position : left_position + 3]
        maximum_right = min(len(sequence) - 3, left_position + module_length - 1)
        for right_position in range(left_position + rules.minimum_start_distance, maximum_right + 1):
            distance = right_position - left_position
            if not rules.distance_is_valid(distance, module_length):
                continue
            right = sequence[right_position : right_position + 3]
            candidates = (
                (left, right, "right", left_position, right_position),
                (right, left, "left", right_position, left_position),
            )
            for site_i, site_ii, direction, site_i_position, site_ii_position in candidates:
                key = encode_pattern(site_i, site_ii, direction)
                index_position = index.locate(key)
                if index_position is None:
                    continue
                if not int(index.plasmid_masks[index_position]) & (1 << plasmid_bit):
                    continue
                pair_id = int(index.best_pair_ids[index_position, plasmid_bit])
                if pair_id == 65535:
                    pair_id = None
                return MatchResult(
                    module=normalized,
                    effective_module=effective,
                    original_length=len(normalized),
                    effective_length=module_length,
                    expansion_copies=expansion_copies,
                    plasmid=plasmid,
                    success=True,
                    site_i_3mer=site_i,
                    site_ii_3mer=site_ii,
                    direction=direction,
                    site_i_position=site_i_position,
                    site_ii_position=site_ii_position,
                    pattern_key=key,
                    solution_count=int(index.solution_counts[index_position, plasmid_bit]),
                    best_pair_id=pair_id,
                )
    return MatchResult(
        module=normalized,
        effective_module=effective,
        original_length=len(normalized),
        effective_length=module_length,
        expansion_copies=expansion_copies,
        plasmid=plasmid,
        success=False,
    )


def query_all_plasmids(
    module: str,
    index: PatternIndex,
    *,
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
    expand_short: bool = True,
) -> list[MatchResult]:
    """Query all plasmids in one sequence scan while preserving legacy order.

    Each plasmid receives the same first match it would receive from an
    independent :func:`match_module` call, but candidate 3-mers are encoded
    only once.  This is the hot path for the 3.37-million-motif screen.
    """
    normalized = validate_protein_sequence(module)
    if len(normalized) < 6 and expand_short:
        effective, expansion_copies = expand_short_module(normalized)
    else:
        effective, expansion_copies = normalized, 1
    module_length = len(effective)
    sequence = effective * 2 if rules.double_module_before_matching else effective
    resolved: dict[int, MatchResult] = {}

    for left_position in range(len(sequence) - 2):
        left = sequence[left_position : left_position + 3]
        maximum_right = min(len(sequence) - 3, left_position + module_length - 1)
        for right_position in range(left_position + rules.minimum_start_distance, maximum_right + 1):
            distance = right_position - left_position
            if not rules.distance_is_valid(distance, module_length):
                continue
            right = sequence[right_position : right_position + 3]
            candidates = (
                (left, right, "right", left_position, right_position),
                (right, left, "left", right_position, left_position),
            )
            for site_i, site_ii, direction, site_i_position, site_ii_position in candidates:
                index_position = index.locate(encode_pattern(site_i, site_ii, direction))
                if index_position is None:
                    continue
                compatible_mask = int(index.plasmid_masks[index_position])
                for bit, plasmid in enumerate(PLASMIDS):
                    if bit in resolved or not compatible_mask & (1 << bit):
                        continue
                    pair_id = int(index.best_pair_ids[index_position, bit])
                    resolved[bit] = MatchResult(
                        module=normalized,
                        effective_module=effective,
                        original_length=len(normalized),
                        effective_length=module_length,
                        expansion_copies=expansion_copies,
                        plasmid=plasmid,
                        success=True,
                        site_i_3mer=site_i,
                        site_ii_3mer=site_ii,
                        direction=direction,
                        site_i_position=site_i_position,
                        site_ii_position=site_ii_position,
                        pattern_key=encode_pattern(site_i, site_ii, direction),
                        solution_count=int(index.solution_counts[index_position, bit]),
                        best_pair_id=None if pair_id == 65535 else pair_id,
                    )
                if len(resolved) == len(PLASMIDS):
                    break
            if len(resolved) == len(PLASMIDS):
                break
        if len(resolved) == len(PLASMIDS):
            break

    return [
        resolved.get(
            bit,
            MatchResult(
                module=normalized,
                effective_module=effective,
                original_length=len(normalized),
                effective_length=module_length,
                expansion_copies=expansion_copies,
                plasmid=plasmid,
                success=False,
            ),
        )
        for bit, plasmid in enumerate(PLASMIDS)
    ]


def materialize_best_solution(result: MatchResult, index: PatternIndex) -> dict[str, object]:
    """Attach enzyme and DNA fields to a successful match."""
    payload = result.to_dict()
    if not result.success or result.best_pair_id is None:
        return payload
    pair = index.pair_table.loc[index.pair_table["pair_id"] == result.best_pair_id]
    if pair.empty:
        raise ValueError(f"Missing enzyme pair {result.best_pair_id}")
    pair_row = pair.iloc[0]
    site_i = index.site_i_table[
        (index.site_i_table["site_i_enzyme"] == pair_row["site_i_enzyme"])
        & (index.site_i_table["site_i_3mer_aa"] == result.site_i_3mer)
    ]
    site_ii = index.site_ii_table[
        (index.site_ii_table["site_ii_enzyme"] == pair_row["site_ii_enzyme"])
        & (index.site_ii_table["site_ii_3mer_aa"] == result.site_ii_3mer)
    ]
    if site_i.empty or site_ii.empty:
        raise ValueError("Pattern index points to a missing site variant")
    payload.update(
        {
            "site_i_enzyme": pair_row["site_i_enzyme"],
            "site_ii_enzyme": pair_row["site_ii_enzyme"],
            "site_iii_enzymes": pair_row["site_iii_enzymes"],
            "site_iii_sites": pair_row.get("site_iii_sites", ""),
            "site_i_ovhg": int(pair_row["site_i_ovhg"]),
            "site_ii_ovhg": int(pair_row["site_ii_ovhg"]),
            "orthogonality": float(pair_row["orthogonality"]),
            "site_i_9mer_bp": site_i.iloc[0]["site_i_9mer_bp"],
            "site_i_recognition_site": site_i.iloc[0].get("site_i_recognition_site", ""),
            "site_i_codon_usage_freq": float(site_i.iloc[0].get("site_i_codon_usage_freq", 0.0)),
            "site_ii_9mer_bp_original": site_ii.iloc[0]["site_ii_9mer_bp_original"],
            "site_ii_9mer_bp_mutated": site_ii.iloc[0]["site_ii_9mer_bp_mutated"],
            "site_ii_recognition_site": site_ii.iloc[0].get("site_ii_recognition_site", ""),
            "site_ii_codon_usage_freq": float(site_ii.iloc[0].get("site_ii_codon_usage_freq", 0.0)),
        }
    )
    return payload


def enumerate_module_solutions(
    module: str,
    index: PatternIndex,
    *,
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
    expand_short: bool = True,
) -> list[dict[str, object]]:
    """Materialize every retained pattern/position/enzyme-pair candidate."""
    if index.directory is None:
        return [
            materialize_best_solution(result, index)
            for result in query_all_plasmids(module, index, rules=rules, expand_short=expand_short)
            if result.success
        ]
    catalog_root = index.directory / str(index.metadata.get("solution_catalog", "pattern_solutions.parquet"))
    if not catalog_root.exists():
        return [
            materialize_best_solution(result, index)
            for result in query_all_plasmids(module, index, rules=rules, expand_short=expand_short)
            if result.success
        ]

    normalized = validate_protein_sequence(module)
    if len(normalized) < 6 and expand_short:
        effective, expansion_copies = expand_short_module(normalized)
    else:
        effective, expansion_copies = normalized, 1
    module_length = len(effective)
    sequence = effective * 2 if rules.double_module_before_matching else effective
    occurrences: dict[int, list[tuple[int, str, str, str, int, int]]] = {bit: [] for bit in range(len(PLASMIDS))}
    for left_position in range(len(sequence) - 2):
        left = sequence[left_position : left_position + 3]
        maximum_right = min(len(sequence) - 3, left_position + module_length - 1)
        for right_position in range(left_position + rules.minimum_start_distance, maximum_right + 1):
            if not rules.distance_is_valid(right_position - left_position, module_length):
                continue
            right = sequence[right_position : right_position + 3]
            candidates = (
                (left, right, "right", left_position, right_position),
                (right, left, "left", right_position, left_position),
            )
            for site_i, site_ii, direction, site_i_position, site_ii_position in candidates:
                key = encode_pattern(site_i, site_ii, direction)
                index_position = index.locate(key)
                if index_position is None:
                    continue
                mask = int(index.plasmid_masks[index_position])
                for bit in range(len(PLASMIDS)):
                    if mask & (1 << bit):
                        occurrences[bit].append((key, site_i, site_ii, direction, site_i_position, site_ii_position))

    rows: list[dict[str, object]] = []
    pair_metadata_cache: dict[tuple[int, str, str], dict[str, object]] = {}
    for bit, plasmid in enumerate(PLASMIDS):
        if not occurrences[bit]:
            continue
        keys = sorted({item[0] for item in occurrences[bit]})
        catalog_path = catalog_root / f"plasmid={plasmid}" / "part-00000.parquet"
        catalog = pd.read_parquet(catalog_path, filters=[("pattern_key", "in", keys)])
        by_key = {int(key): group for key, group in catalog.groupby("pattern_key", sort=False)}
        for key, site_i, site_ii, direction, site_i_position, site_ii_position in occurrences[bit]:
            for candidate in by_key.get(key, pd.DataFrame()).itertuples(index=False):
                result = MatchResult(
                    module=normalized,
                    effective_module=effective,
                    original_length=len(normalized),
                    effective_length=module_length,
                    expansion_copies=expansion_copies,
                    plasmid=plasmid,
                    success=True,
                    site_i_3mer=site_i,
                    site_ii_3mer=site_ii,
                    direction=direction,
                    site_i_position=site_i_position,
                    site_ii_position=site_ii_position,
                    pattern_key=key,
                    solution_count=int(candidate.codon_variant_count),
                    best_pair_id=int(candidate.pair_id),
                )
                payload = result.to_dict()
                metadata_key = (int(candidate.pair_id), site_i, site_ii)
                if metadata_key not in pair_metadata_cache:
                    materialized = materialize_best_solution(result, index)
                    pair_metadata_cache[metadata_key] = {
                        field: value for field, value in materialized.items() if field not in payload
                    }
                payload.update(pair_metadata_cache[metadata_key])
                payload["candidate_pair_id"] = int(candidate.pair_id)
                payload["codon_variant_count"] = int(candidate.codon_variant_count)
                rows.append(payload)
    return rows
