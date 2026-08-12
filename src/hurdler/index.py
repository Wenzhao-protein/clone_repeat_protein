"""Build and load the compact, versioned HURDLER pattern index."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import AA_TO_INT, AMINO_ACIDS, PLASMIDS, SCHEMA_VERSION, THREE_MER_SPACE
from .io import require_columns, sha256_file, utc_now, write_json_atomic
from .rules import RuleProfile, LEGACY_OPTIMIZED_V1


def encode_3mer(sequence: str) -> int:
    if len(sequence) != 3:
        raise ValueError(f"Expected a 3mer, got {sequence!r}")
    try:
        return (AA_TO_INT[sequence[0]] * 20 + AA_TO_INT[sequence[1]]) * 20 + AA_TO_INT[sequence[2]]
    except KeyError as exc:
        raise ValueError(f"Invalid 3mer: {sequence!r}") from exc


def decode_3mer(code: int) -> str:
    if not 0 <= int(code) < THREE_MER_SPACE:
        raise ValueError(f"3mer code out of range: {code}")
    first, remainder = divmod(int(code), 400)
    second, third = divmod(remainder, 20)
    return AMINO_ACIDS[first] + AMINO_ACIDS[second] + AMINO_ACIDS[third]


def encode_pattern(site_i: str, site_ii: str, direction: str) -> int:
    direction_bit = 0 if direction == "right" else 1 if direction == "left" else None
    if direction_bit is None:
        raise ValueError(f"Invalid direction: {direction!r}")
    return ((encode_3mer(site_i) * THREE_MER_SPACE + encode_3mer(site_ii)) * 2) + direction_bit


def decode_pattern(key: int) -> tuple[str, str, str]:
    pair, direction_bit = divmod(int(key), 2)
    site_i_code, site_ii_code = divmod(pair, THREE_MER_SPACE)
    return decode_3mer(site_i_code), decode_3mer(site_ii_code), "right" if direction_bit == 0 else "left"


def _as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().map({"true": True, "false": False}).fillna(False).astype(bool)


def plasmid_mask_from_row(row: pd.Series, prefix: str) -> int:
    mask = 0
    for bit, plasmid in enumerate(PLASMIDS):
        value = row[f"{prefix}{plasmid}"]
        if bool(value):
            mask |= 1 << bit
    return mask


@dataclass(frozen=True)
class PatternIndex:
    keys: np.ndarray
    plasmid_masks: np.ndarray
    solution_counts: np.ndarray
    best_pair_ids: np.ndarray
    pair_table: pd.DataFrame
    site_i_table: pd.DataFrame
    site_ii_table: pd.DataFrame
    metadata: dict[str, object]
    directory: Path | None = None

    @classmethod
    def load(cls, directory: str | Path) -> "PatternIndex":
        root = Path(directory)
        arrays = np.load(root / "pattern_index.npz", allow_pickle=False)
        metadata = json.loads((root / "metadata.json").read_text())
        if metadata.get("schema_version") != SCHEMA_VERSION:
            raise ValueError(f"Unsupported index schema: {metadata.get('schema_version')}")
        result = cls(
            keys=arrays["keys"],
            plasmid_masks=arrays["plasmid_masks"],
            solution_counts=arrays["solution_counts"],
            best_pair_ids=arrays["best_pair_ids"],
            pair_table=pd.read_parquet(root / "enzyme_pairs.parquet"),
            site_i_table=pd.read_parquet(root / "site_i_variants.parquet"),
            site_ii_table=pd.read_parquet(root / "site_ii_variants.parquet"),
            metadata=metadata,
            directory=root,
        )
        result.validate()
        return result

    def validate(self) -> None:
        n = len(self.keys)
        if not (len(self.plasmid_masks) == n == len(self.solution_counts) == len(self.best_pair_ids)):
            raise ValueError("Pattern-index arrays have inconsistent lengths")
        if self.best_pair_ids.shape != (n, len(PLASMIDS)):
            raise ValueError(f"Unexpected best_pair_ids shape: {self.best_pair_ids.shape}")
        if self.solution_counts.shape != (n, len(PLASMIDS)):
            raise ValueError(f"Unexpected solution_counts shape: {self.solution_counts.shape}")
        if n and (np.any(self.keys[1:] <= self.keys[:-1])):
            raise ValueError("Pattern keys must be strictly increasing")

    def locate(self, key: int) -> int | None:
        position = int(np.searchsorted(self.keys, np.uint32(key)))
        if position < len(self.keys) and int(self.keys[position]) == int(key):
            return position
        return None


def build_pattern_index(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    orthogonality_path: str | Path | None = None,
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
) -> dict[str, object]:
    """Build a sparse index from the latest optimized-notebook tables."""
    source = Path(input_dir)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)

    site_i = pd.read_csv(source / "hurdler_site_i_dataframe.csv")
    site_ii = pd.read_csv(source / "hurdler_site_ii_dataframe.csv")
    selected_iii = pd.read_csv(source / "selected_site_iii_enzymes.csv")
    matrix = pd.read_csv(source / "site_i_site_ii_pairing_matrix.csv", index_col=0)

    required_i = {"site_i_enzyme", "site_i_ovhg", "site_i_3mer_aa", "site_i_9mer_bp", "site_i_codon_usage_freq"}
    required_ii = {"site_ii_enzyme", "site_ii_ovhg", "site_ii_3mer_aa", "site_ii_9mer_bp_original", "site_ii_9mer_bp_mutated", "site_ii_search_direction", "site_ii_codon_usage_freq"}
    required_i.update(f"site_i_{p}" for p in PLASMIDS)
    required_ii.update(f"site_ii_{p}" for p in PLASMIDS)
    require_columns(site_i, required_i, source=str(source / "hurdler_site_i_dataframe.csv"))
    require_columns(site_ii, required_ii, source=str(source / "hurdler_site_ii_dataframe.csv"))

    for column in [c for c in site_i.columns if c.startswith("site_i_p")]:
        site_i[column] = _as_bool(site_i[column])
    for column in [c for c in site_ii.columns if c.startswith("site_ii_p")]:
        site_ii[column] = _as_bool(site_ii[column])
    matrix = matrix.apply(_as_bool)

    fidelity: dict[tuple[str, str], float] = {}
    if orthogonality_path is not None and Path(orthogonality_path).exists():
        fidelity_frame = pd.read_csv(orthogonality_path)
        for row in fidelity_frame.itertuples(index=False):
            fidelity[(str(row.re1), str(row.re2))] = float(row.orthogonality)
            fidelity[(str(row.re2), str(row.re1))] = float(row.orthogonality)

    iii_by_overhang = {
        int(overhang): tuple(sorted(group["enzyme"].astype(str).unique()))
        for overhang, group in selected_iii.groupby("ovhg")
    }
    iii_site_by_enzyme = dict(zip(selected_iii["enzyme"].astype(str), selected_iii["site"].astype(str)))
    i_groups = {name: group.reset_index(drop=True) for name, group in site_i.groupby("site_i_enzyme")}
    ii_groups = {name: group.reset_index(drop=True) for name, group in site_ii.groupby("site_ii_enzyme")}

    pair_rows: list[dict[str, object]] = []
    raw_keys: list[np.ndarray] = []
    raw_masks: list[np.ndarray] = []
    raw_scores: list[np.ndarray] = []
    raw_pair_ids: list[np.ndarray] = []

    for enzyme_i in matrix.index.astype(str):
        if enzyme_i not in i_groups:
            continue
        for enzyme_ii in matrix.columns.astype(str):
            if not bool(matrix.loc[enzyme_i, enzyme_ii]) or enzyme_ii not in ii_groups:
                continue
            rows_i = i_groups[enzyme_i]
            rows_ii = ii_groups[enzyme_ii]
            ovhg_i = int(rows_i.iloc[0]["site_i_ovhg"])
            ovhg_ii = int(rows_ii.iloc[0]["site_ii_ovhg"])
            site_iii = iii_by_overhang.get(ovhg_ii, ())
            if not site_iii:
                continue
            mask_i = plasmid_mask_from_row(rows_i.iloc[0], "site_i_")
            mask_ii = plasmid_mask_from_row(rows_ii.iloc[0], "site_ii_")
            plasmid_mask = mask_i & mask_ii
            if not plasmid_mask:
                continue
            orthogonality = fidelity.get((enzyme_i, enzyme_ii))
            if orthogonality is None:
                orthogonality = 4.0 if ovhg_i != ovhg_ii or rules.missing_fidelity_is_compatible else 0.0
            pair_id = len(pair_rows)
            pair_rows.append({
                "pair_id": pair_id,
                "site_i_enzyme": enzyme_i,
                "site_ii_enzyme": enzyme_ii,
                "site_i_ovhg": ovhg_i,
                "site_ii_ovhg": ovhg_ii,
                "site_iii_enzymes": ",".join(site_iii),
                "site_iii_sites": ",".join(iii_site_by_enzyme[name] for name in site_iii),
                "orthogonality": float(orthogonality),
                "plasmid_mask": plasmid_mask,
            })

            i_codes = np.fromiter((encode_3mer(x) for x in rows_i["site_i_3mer_aa"]), dtype=np.uint16)
            ii_codes = np.fromiter((encode_3mer(x) for x in rows_ii["site_ii_3mer_aa"]), dtype=np.uint16)
            directions = (rows_ii["site_ii_search_direction"].astype(str).to_numpy() == "left").astype(np.uint32)
            keys = (((i_codes.astype(np.uint32)[:, None] * THREE_MER_SPACE + ii_codes.astype(np.uint32)[None, :]) * 2) + directions[None, :]).ravel()
            scores = (
                rows_i["site_i_codon_usage_freq"].to_numpy(dtype=np.float32)[:, None]
                + rows_ii["site_ii_codon_usage_freq"].to_numpy(dtype=np.float32)[None, :]
                + np.float32(orthogonality / 1000.0)
            ).ravel()
            raw_keys.append(keys.astype(np.uint32, copy=False))
            raw_masks.append(np.full(keys.size, plasmid_mask, dtype=np.uint8))
            raw_scores.append(scores.astype(np.float32, copy=False))
            raw_pair_ids.append(np.full(keys.size, pair_id, dtype=np.uint16))

    if not raw_keys:
        raise ValueError("No valid enzyme pairs were generated")
    keys_all = np.concatenate(raw_keys)
    masks_all = np.concatenate(raw_masks)
    scores_all = np.concatenate(raw_scores)
    pairs_all = np.concatenate(raw_pair_ids)

    key_order = np.argsort(keys_all, kind="stable")
    sorted_keys = keys_all[key_order]
    group_starts = np.flatnonzero(np.r_[True, sorted_keys[1:] != sorted_keys[:-1]])
    unique_keys = sorted_keys[group_starts]
    plasmid_masks = np.bitwise_or.reduceat(masks_all[key_order], group_starts).astype(np.uint8)
    solution_counts = np.zeros((len(unique_keys), len(PLASMIDS)), dtype=np.uint16)
    best_pair_ids = np.full((len(unique_keys), len(PLASMIDS)), np.iinfo(np.uint16).max, dtype=np.uint16)
    solution_catalog_root = destination / "pattern_solutions.parquet"
    solution_catalog_rows = 0

    for bit, _plasmid in enumerate(PLASMIDS):
        compatible = (masks_all & (1 << bit)) != 0
        compatible_order = np.lexsort((-scores_all[compatible], keys_all[compatible]))
        original_indices = np.flatnonzero(compatible)[compatible_order]
        compatible_keys = keys_all[original_indices]
        first = np.r_[True, compatible_keys[1:] != compatible_keys[:-1]]
        selected_indices = original_indices[first]
        selected_keys = keys_all[selected_indices]
        selected_positions = np.searchsorted(unique_keys, selected_keys)
        best_pair_ids[selected_positions, bit] = pairs_all[selected_indices]

        sorted_compatible_keys = np.sort(keys_all[compatible], kind="stable")
        starts = np.flatnonzero(np.r_[True, sorted_compatible_keys[1:] != sorted_compatible_keys[:-1]])
        counts = np.diff(np.r_[starts, len(sorted_compatible_keys)])
        positions = np.searchsorted(unique_keys, sorted_compatible_keys[starts])
        solution_counts[positions, bit] = np.minimum(counts, np.iinfo(np.uint16).max).astype(np.uint16)

        # Retain every distinct enzyme-pair candidate in a normalized Parquet
        # catalog.  The compact NPZ remains the latency-sensitive query index.
        import pyarrow as pa
        import pyarrow.parquet as pq

        pair_order = np.lexsort((pairs_all[compatible], keys_all[compatible]))
        catalog_keys = keys_all[compatible][pair_order]
        catalog_pairs = pairs_all[compatible][pair_order]
        catalog_starts = np.flatnonzero(
            np.r_[
                True,
                (catalog_keys[1:] != catalog_keys[:-1]) | (catalog_pairs[1:] != catalog_pairs[:-1]),
            ]
        )
        variant_counts = np.diff(np.r_[catalog_starts, len(catalog_keys)])
        partition = solution_catalog_root / f"plasmid={_plasmid}"
        partition.mkdir(parents=True, exist_ok=True)
        pq.write_table(
            pa.table(
                {
                    "pattern_key": catalog_keys[catalog_starts].astype(np.uint32),
                    "pair_id": catalog_pairs[catalog_starts].astype(np.uint16),
                    "codon_variant_count": np.minimum(variant_counts, np.iinfo(np.uint16).max).astype(np.uint16),
                }
            ),
            partition / "part-00000.parquet",
            compression="zstd",
        )
        solution_catalog_rows += len(catalog_starts)

    np.savez_compressed(
        destination / "pattern_index.npz",
        keys=unique_keys.astype(np.uint32),
        plasmid_masks=plasmid_masks,
        solution_counts=solution_counts,
        best_pair_ids=best_pair_ids,
    )
    pair_table = pd.DataFrame(pair_rows)
    pair_table.to_parquet(destination / "enzyme_pairs.parquet", index=False)
    site_i.to_parquet(destination / "site_i_variants.parquet", index=False)
    site_ii.to_parquet(destination / "site_ii_variants.parquet", index=False)

    metadata = {
        "schema_version": SCHEMA_VERSION,
        "rule_profile": rules.to_dict(),
        "created_at": utc_now(),
        "input_dir": str(source.absolute()),
        "pattern_count": int(len(unique_keys)),
        "candidate_count": int(len(keys_all)),
        "enzyme_pair_count": int(len(pair_table)),
        "normalized_solution_rows": int(solution_catalog_rows),
        "solution_catalog": "pattern_solutions.parquet",
        "plasmids": list(PLASMIDS),
        "source_hashes": {
            name: sha256_file(source / name)
            for name in (
                "hurdler_site_i_dataframe.csv",
                "hurdler_site_ii_dataframe.csv",
                "selected_site_iii_enzymes.csv",
                "site_i_site_ii_pairing_matrix.csv",
            )
        },
    }
    write_json_atomic(metadata, destination / "metadata.json")
    return metadata
