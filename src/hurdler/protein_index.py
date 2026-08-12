"""Plasmid-independent HURDLER protein-pattern index.

Unlike the frozen legacy index, this artifact does not apply any whole-vector
restriction-site mask while it is built.  Vector compatibility is evaluated
later against an explicit annotated cut scheme.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .constants import THREE_MER_SPACE, validate_protein_sequence
from .index import _as_bool, encode_3mer, encode_pattern
from .io import require_columns, sha256_file, utc_now, write_json_atomic
from .matching import expand_short_module
from .rules import LEGACY_OPTIMIZED_V1, RuleProfile


PROTEIN_INDEX_SCHEMA_VERSION = "protein-pattern-index-v2"
PROTEIN_INDEX_VERSION = "vector-aware-hurdler-v2"


@dataclass(frozen=True)
class ProteinPatternIndex:
    keys: np.ndarray
    offsets: np.ndarray
    pair_ids: np.ndarray
    variant_counts: np.ndarray
    pair_table: pd.DataFrame
    site_i_table: pd.DataFrame
    site_ii_table: pd.DataFrame
    metadata: dict[str, Any]
    directory: Path | None = None

    @classmethod
    def load(cls, directory: str | Path) -> "ProteinPatternIndex":
        root = Path(directory)
        arrays = np.load(root / "protein_pattern_index.npz", allow_pickle=False)
        metadata = json.loads((root / "metadata.json").read_text())
        if metadata.get("schema_version") != PROTEIN_INDEX_SCHEMA_VERSION:
            raise ValueError(f"Unsupported protein index schema: {metadata.get('schema_version')!r}")
        index = cls(
            keys=arrays["keys"],
            offsets=arrays["offsets"],
            pair_ids=arrays["pair_ids"],
            variant_counts=arrays["variant_counts"],
            pair_table=pd.read_parquet(root / "enzyme_pairs.parquet"),
            site_i_table=pd.read_parquet(root / "site_i_variants.parquet"),
            site_ii_table=pd.read_parquet(root / "site_ii_variants.parquet"),
            metadata=metadata,
            directory=root,
        )
        index.validate()
        return index

    def validate(self) -> None:
        if len(self.offsets) != len(self.keys) + 1:
            raise ValueError("CSR offsets must have pattern_count + 1 entries")
        if int(self.offsets[-1]) != len(self.pair_ids) or len(self.pair_ids) != len(self.variant_counts):
            raise ValueError("CSR protein-index arrays have inconsistent lengths")
        if len(self.keys) and np.any(self.keys[1:] <= self.keys[:-1]):
            raise ValueError("Protein-pattern keys must be strictly increasing")
        if len(self.pair_table) != int(self.metadata.get("enzyme_pair_count", -1)):
            raise ValueError("Protein enzyme-pair cardinality differs from metadata")

    def pairs_for_key(self, key: int) -> tuple[np.ndarray, np.ndarray]:
        position = int(np.searchsorted(self.keys, np.uint32(key)))
        if position >= len(self.keys) or int(self.keys[position]) != int(key):
            empty = np.empty(0, dtype=np.uint16)
            return empty, empty
        start, end = int(self.offsets[position]), int(self.offsets[position + 1])
        return self.pair_ids[start:end], self.variant_counts[start:end]


def build_protein_pattern_index(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    orthogonality_path: str | Path | None = None,
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
) -> dict[str, Any]:
    source = Path(input_dir)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    site_i = pd.read_csv(source / "hurdler_site_i_dataframe.csv")
    site_ii = pd.read_csv(source / "hurdler_site_ii_dataframe.csv")
    selected_iii = pd.read_csv(source / "selected_site_iii_enzymes.csv")
    matrix = pd.read_csv(source / "site_i_site_ii_pairing_matrix.csv", index_col=0).apply(_as_bool)
    require_columns(
        site_i,
        {"site_i_enzyme", "site_i_ovhg", "site_i_3mer_aa", "site_i_9mer_bp", "site_i_codon_usage_freq"},
        source=str(source / "hurdler_site_i_dataframe.csv"),
    )
    require_columns(
        site_ii,
        {"site_ii_enzyme", "site_ii_ovhg", "site_ii_3mer_aa", "site_ii_9mer_bp_original", "site_ii_9mer_bp_mutated", "site_ii_search_direction", "site_ii_codon_usage_freq"},
        source=str(source / "hurdler_site_ii_dataframe.csv"),
    )

    fidelity: dict[tuple[str, str], float] = {}
    if orthogonality_path is not None and Path(orthogonality_path).exists():
        for row in pd.read_csv(orthogonality_path).itertuples(index=False):
            fidelity[(str(row.re1), str(row.re2))] = float(row.orthogonality)
            fidelity[(str(row.re2), str(row.re1))] = float(row.orthogonality)
    iii_by_overhang = {
        int(overhang): tuple(sorted(group["enzyme"].astype(str).unique()))
        for overhang, group in selected_iii.groupby("ovhg")
    }
    iii_site_by_enzyme = dict(zip(selected_iii["enzyme"].astype(str), selected_iii["site"].astype(str)))
    i_groups = {name: group.reset_index(drop=True) for name, group in site_i.groupby("site_i_enzyme")}
    ii_groups = {name: group.reset_index(drop=True) for name, group in site_ii.groupby("site_ii_enzyme")}
    pair_rows: list[dict[str, Any]] = []
    raw_keys: list[np.ndarray] = []
    raw_pairs: list[np.ndarray] = []
    for enzyme_i in matrix.index.astype(str):
        if enzyme_i not in i_groups:
            continue
        for enzyme_ii in matrix.columns.astype(str):
            if not bool(matrix.loc[enzyme_i, enzyme_ii]) or enzyme_ii not in ii_groups:
                continue
            rows_i, rows_ii = i_groups[enzyme_i], ii_groups[enzyme_ii]
            ovhg_i, ovhg_ii = int(rows_i.iloc[0].site_i_ovhg), int(rows_ii.iloc[0].site_ii_ovhg)
            site_iii = iii_by_overhang.get(ovhg_ii, ())
            if not site_iii:
                continue
            orthogonality = fidelity.get((enzyme_i, enzyme_ii))
            if orthogonality is None:
                orthogonality = 4.0 if ovhg_i != ovhg_ii or rules.missing_fidelity_is_compatible else 0.0
            pair_id = len(pair_rows)
            pair_rows.append({
                "pair_id": pair_id,
                "site_i_enzyme": enzyme_i,
                "site_ii_enzyme": enzyme_ii,
                "site_i_recognition_site": str(rows_i.iloc[0].get("site_i_recognition_site", "")),
                "site_ii_recognition_site": str(rows_ii.iloc[0].get("site_ii_recognition_site", "")),
                "site_i_ovhg": ovhg_i,
                "site_ii_ovhg": ovhg_ii,
                "site_iii_enzymes": ",".join(site_iii),
                "site_iii_sites": ",".join(iii_site_by_enzyme[name] for name in site_iii),
                "orthogonality": float(orthogonality),
            })
            i_codes = np.fromiter((encode_3mer(x) for x in rows_i.site_i_3mer_aa), dtype=np.uint16)
            ii_codes = np.fromiter((encode_3mer(x) for x in rows_ii.site_ii_3mer_aa), dtype=np.uint16)
            directions = (rows_ii.site_ii_search_direction.astype(str).to_numpy() == "left").astype(np.uint32)
            keys = (((i_codes.astype(np.uint32)[:, None] * THREE_MER_SPACE + ii_codes.astype(np.uint32)[None, :]) * 2) + directions[None, :]).ravel()
            raw_keys.append(keys.astype(np.uint32, copy=False))
            raw_pairs.append(np.full(keys.size, pair_id, dtype=np.uint16))
    if not raw_keys:
        raise ValueError("No protein-level enzyme pairs were generated")

    keys_all = np.concatenate(raw_keys)
    pairs_all = np.concatenate(raw_pairs)
    order = np.lexsort((pairs_all, keys_all))
    keys_sorted, pairs_sorted = keys_all[order], pairs_all[order]
    starts = np.flatnonzero(np.r_[True, (keys_sorted[1:] != keys_sorted[:-1]) | (pairs_sorted[1:] != pairs_sorted[:-1])])
    distinct_keys = keys_sorted[starts]
    distinct_pairs = pairs_sorted[starts]
    counts = np.diff(np.r_[starts, len(keys_sorted)])
    key_starts = np.flatnonzero(np.r_[True, distinct_keys[1:] != distinct_keys[:-1]])
    unique_keys = distinct_keys[key_starts]
    offsets = np.r_[key_starts, len(distinct_keys)].astype(np.uint32)
    variant_counts = np.minimum(counts, np.iinfo(np.uint16).max).astype(np.uint16)
    np.savez_compressed(
        destination / "protein_pattern_index.npz",
        keys=unique_keys.astype(np.uint32),
        offsets=offsets,
        pair_ids=distinct_pairs.astype(np.uint16),
        variant_counts=variant_counts,
    )
    pair_table = pd.DataFrame(pair_rows)
    pair_table.to_parquet(destination / "enzyme_pairs.parquet", index=False)
    site_i.to_parquet(destination / "site_i_variants.parquet", index=False)
    site_ii.to_parquet(destination / "site_ii_variants.parquet", index=False)
    metadata = {
        "schema_version": PROTEIN_INDEX_SCHEMA_VERSION,
        "artifact_version": PROTEIN_INDEX_VERSION,
        "created_at": utc_now(),
        "rule_profile": rules.to_dict(),
        "protein_only": True,
        "whole_plasmid_mask_applied": False,
        "pattern_count": int(len(unique_keys)),
        "pattern_pair_count": int(len(distinct_pairs)),
        "candidate_variant_count": int(len(keys_all)),
        "enzyme_pair_count": int(len(pair_table)),
        "source_hashes": {
            name: sha256_file(source / name)
            for name in ("hurdler_site_i_dataframe.csv", "hurdler_site_ii_dataframe.csv", "selected_site_iii_enzymes.csv", "site_i_site_ii_pairing_matrix.csv")
        },
    }
    write_json_atomic(metadata, destination / "metadata.json")
    return metadata


def enumerate_protein_solutions(
    module: str,
    index: ProteinPatternIndex,
    *,
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
    expand_short: bool = True,
) -> list[dict[str, Any]]:
    normalized = validate_protein_sequence(module)
    effective, expansion_copies = expand_short_module(normalized) if len(normalized) < 6 and expand_short else (normalized, 1)
    module_length = len(effective)
    sequence = effective * 2 if rules.double_module_before_matching else effective
    pair_lookup = index.pair_table.set_index("pair_id", drop=False)
    solutions: list[dict[str, Any]] = []
    seen: set[tuple[int, int, int, str]] = set()
    for left_position in range(len(sequence) - 2):
        left = sequence[left_position:left_position + 3]
        maximum_right = min(len(sequence) - 3, left_position + module_length - 1)
        for right_position in range(left_position + rules.minimum_start_distance, maximum_right + 1):
            distance = right_position - left_position
            if not rules.distance_is_valid(distance, module_length):
                continue
            right = sequence[right_position:right_position + 3]
            for site_i, site_ii, direction, site_i_position, site_ii_position in (
                (left, right, "right", left_position, right_position),
                (right, left, "left", right_position, left_position),
            ):
                pattern_key = encode_pattern(site_i, site_ii, direction)
                pair_ids, counts = index.pairs_for_key(pattern_key)
                for pair_id, count in zip(pair_ids, counts, strict=True):
                    marker = (int(pair_id), site_i_position, site_ii_position, direction)
                    if marker in seen:
                        continue
                    seen.add(marker)
                    pair = pair_lookup.loc[int(pair_id)]
                    i_rows = index.site_i_table[(index.site_i_table.site_i_enzyme == pair.site_i_enzyme) & (index.site_i_table.site_i_3mer_aa == site_i)]
                    ii_rows = index.site_ii_table[(index.site_ii_table.site_ii_enzyme == pair.site_ii_enzyme) & (index.site_ii_table.site_ii_3mer_aa == site_ii) & (index.site_ii_table.site_ii_search_direction == direction)]
                    if i_rows.empty or ii_rows.empty:
                        continue
                    i_row = i_rows.sort_values("site_i_codon_usage_freq", ascending=False).iloc[0]
                    ii_row = ii_rows.sort_values("site_ii_codon_usage_freq", ascending=False).iloc[0]
                    solutions.append({
                        "module": normalized,
                        "effective_module": effective,
                        "original_length": len(normalized),
                        "effective_length": module_length,
                        "expansion_copies": expansion_copies,
                        "pair_id": int(pair_id),
                        "pattern_key": int(pattern_key),
                        "direction": direction,
                        "site_i_position": site_i_position,
                        "site_ii_position": site_ii_position,
                        "site_i_3mer": site_i,
                        "site_ii_3mer": site_ii,
                        "site_i_enzyme": str(pair.site_i_enzyme),
                        "site_ii_enzyme": str(pair.site_ii_enzyme),
                        "site_iii_enzymes": str(pair.site_iii_enzymes),
                        "site_iii_sites": str(pair.site_iii_sites),
                        "site_i_recognition_site": str(pair.site_i_recognition_site),
                        "site_ii_recognition_site": str(pair.site_ii_recognition_site),
                        "site_i_9mer_bp": str(i_row.site_i_9mer_bp),
                        "site_ii_9mer_bp_original": str(ii_row.site_ii_9mer_bp_original),
                        "site_ii_9mer_bp_mutated": str(ii_row.site_ii_9mer_bp_mutated),
                        "site_i_ovhg": int(pair.site_i_ovhg),
                        "site_ii_ovhg": int(pair.site_ii_ovhg),
                        "orthogonality": float(pair.orthogonality),
                        "codon_variant_count": int(count),
                        "codon_usage_score": float(i_row.site_i_codon_usage_freq) + float(ii_row.site_ii_codon_usage_freq),
                    })
    return sorted(
        solutions,
        key=lambda row: (
            -float(row["codon_usage_score"]),
            -float(row["orthogonality"]),
            str(row["site_i_enzyme"]),
            str(row["site_ii_enzyme"]),
            int(row["site_i_position"]),
            int(row["site_ii_position"]),
        ),
    )
