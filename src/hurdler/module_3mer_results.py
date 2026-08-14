"""Export the production repeat-module HURDLER/3-mer result table.

The public module result catalog deliberately omits the bulky selected-solution
JSON.  It does retain the selected Site-I/Site-II start positions, so the exact
3-mer amino-acid windows can be recovered without rerunning HURDLER.  This
module performs that recovery and checks every recovered window against the
frozen lookup artifact before writing a compact spreadsheet-oriented table.
"""

from __future__ import annotations

import argparse
import json
from math import ceil
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from .constants import PLASMIDS, validate_protein_sequence


DEFAULT_SOURCE = Path("data/results/natural_designed_repeat_protein_hurdler_idt.csv")
DEFAULT_OUTPUT = Path("data/results/repeatsdb_designed_hurdler_3mer_results.csv")
DEFAULT_INDEX = Path("data/artifacts/legacy-optimized-v1")
EXPECTED_PRODUCTION_COUNTS = {"Natural": 25_913, "Designed": 182}

REQUIRED_SOURCE_COLUMNS = {
    "module_id",
    "display_name",
    "source_accession",
    "collection",
    "family",
    "middle_module_sequence_aa",
    "middle_module_length_aa",
    "hurdler_compatible",
    "selected_plasmid",
    "selected_site_i_enzyme",
    "selected_site_ii_enzyme",
    "selected_direction",
    "selected_site_i_position",
    "selected_site_ii_position",
    "corpus_version",
    "hurdler_rules_version",
}

OUTPUT_COLUMNS = (
    "sequence_id",
    "source_accession",
    "display_name",
    "collection",
    "family",
    "selected_module_sequence_aa",
    "module_length_aa",
    "hurdler_compatible",
    "selected_plasmid",
    "selected_re_pair",
    "site_i_enzyme",
    "site_i_3mer_aa",
    "site_ii_enzyme",
    "site_ii_3mer_aa",
    "three_mer_aa_pair",
    "direction",
    "site_i_start_aa_0based",
    "site_ii_start_aa_0based",
    "corpus_version",
    "hurdler_rules_version",
)


def _required_text(value: Any, field: str, sequence_id: str) -> str:
    if pd.isna(value) or not str(value).strip():
        raise ValueError(f"{sequence_id}: compatible row is missing {field}")
    return str(value).strip()


def _required_position(value: Any, field: str, sequence_id: str) -> int:
    if pd.isna(value):
        raise ValueError(f"{sequence_id}: compatible row is missing {field}")
    numeric = float(value)
    if not numeric.is_integer() or numeric < 0:
        raise ValueError(f"{sequence_id}: {field} must be a non-negative integer")
    return int(numeric)


def _read_lookup(index_dir: str | Path, expected_rules: set[str]) -> tuple[set, set, set]:
    root = Path(index_dir)
    metadata = json.loads((root / "metadata.json").read_text())
    artifact_rule = str(metadata.get("rule_profile", {}).get("name", ""))
    if expected_rules != {artifact_rule}:
        raise ValueError(
            "Source and lookup rule profiles differ: "
            f"source={sorted(expected_rules)!r}, lookup={artifact_rule!r}"
        )

    site_i = pd.read_parquet(
        root / "site_i_variants.parquet",
        columns=["site_i_enzyme", "site_i_3mer_aa"],
    )
    site_ii = pd.read_parquet(
        root / "site_ii_variants.parquet",
        columns=["site_ii_enzyme", "site_ii_3mer_aa"],
    )
    pairs = pd.read_parquet(
        root / "enzyme_pairs.parquet",
        columns=["site_i_enzyme", "site_ii_enzyme", "plasmid_mask"],
    )
    valid_site_i = set(
        site_i[["site_i_enzyme", "site_i_3mer_aa"]].itertuples(
            index=False, name=None
        )
    )
    valid_site_ii = set(
        site_ii[["site_ii_enzyme", "site_ii_3mer_aa"]].itertuples(
            index=False, name=None
        )
    )
    valid_pairs = {
        (str(row.site_i_enzyme), str(row.site_ii_enzyme), plasmid)
        for row in pairs.itertuples(index=False)
        for bit, plasmid in enumerate(PLASMIDS)
        if int(row.plasmid_mask) & (1 << bit)
    }
    return valid_site_i, valid_site_ii, valid_pairs


def build_module_3mer_results(
    source: pd.DataFrame,
    index_dir: str | Path,
    *,
    expected_counts: dict[str, int] | None = None,
) -> pd.DataFrame:
    """Recover and validate the selected Site-I/Site-II 3-mer AA windows."""
    missing = sorted(REQUIRED_SOURCE_COLUMNS - set(source.columns))
    if missing:
        raise ValueError(f"Source result table is missing columns: {missing}")
    if source["module_id"].isna().any() or not source["module_id"].is_unique:
        raise ValueError("module_id must be populated and globally unique")

    observed_counts = source["collection"].value_counts().to_dict()
    if expected_counts is not None and observed_counts != expected_counts:
        raise ValueError(
            f"Unexpected production collection counts: {observed_counts}; "
            f"expected {expected_counts}"
        )

    expected_rules = set(source["hurdler_rules_version"].dropna().astype(str))
    valid_site_i, valid_site_ii, valid_pairs = _read_lookup(
        index_dir, expected_rules
    )

    records: list[dict[str, Any]] = []
    for row in source.itertuples(index=False):
        sequence_id = str(row.module_id)
        unit = validate_protein_sequence(str(row.middle_module_sequence_aa))
        if len(unit) != int(row.middle_module_length_aa):
            raise ValueError(f"{sequence_id}: module sequence/length mismatch")
        compatible = bool(row.hurdler_compatible)
        record: dict[str, Any] = {
            "sequence_id": sequence_id,
            "source_accession": str(row.source_accession),
            "display_name": str(row.display_name),
            "collection": str(row.collection),
            "family": "" if pd.isna(row.family) else str(row.family),
            "selected_module_sequence_aa": unit,
            "module_length_aa": len(unit),
            "hurdler_compatible": compatible,
            "selected_plasmid": "",
            "selected_re_pair": "",
            "site_i_enzyme": "",
            "site_i_3mer_aa": "",
            "site_ii_enzyme": "",
            "site_ii_3mer_aa": "",
            "three_mer_aa_pair": "",
            "direction": "",
            "site_i_start_aa_0based": pd.NA,
            "site_ii_start_aa_0based": pd.NA,
            "corpus_version": str(row.corpus_version),
            "hurdler_rules_version": str(row.hurdler_rules_version),
        }
        if compatible:
            plasmid = _required_text(row.selected_plasmid, "selected_plasmid", sequence_id)
            enzyme_i = _required_text(
                row.selected_site_i_enzyme, "selected_site_i_enzyme", sequence_id
            )
            enzyme_ii = _required_text(
                row.selected_site_ii_enzyme, "selected_site_ii_enzyme", sequence_id
            )
            direction = _required_text(
                row.selected_direction, "selected_direction", sequence_id
            )
            position_i = _required_position(
                row.selected_site_i_position, "selected_site_i_position", sequence_id
            )
            position_ii = _required_position(
                row.selected_site_ii_position, "selected_site_ii_position", sequence_id
            )

            expansion_copies = max(1, ceil(6 / len(unit)))
            effective = unit * expansion_copies if len(unit) < 6 else unit
            scanned = effective * 2
            three_mer_i = scanned[position_i : position_i + 3]
            three_mer_ii = scanned[position_ii : position_ii + 3]
            if len(three_mer_i) != 3 or len(three_mer_ii) != 3:
                raise ValueError(f"{sequence_id}: selected 3-mer position is out of range")
            expected_direction = "right" if position_i < position_ii else "left"
            if direction != expected_direction:
                raise ValueError(
                    f"{sequence_id}: direction {direction!r} conflicts with selected positions"
                )
            if (enzyme_i, three_mer_i) not in valid_site_i:
                raise ValueError(
                    f"{sequence_id}: Site-I {enzyme_i}/{three_mer_i} is absent from lookup"
                )
            if (enzyme_ii, three_mer_ii) not in valid_site_ii:
                raise ValueError(
                    f"{sequence_id}: Site-II {enzyme_ii}/{three_mer_ii} is absent from lookup"
                )
            if (enzyme_i, enzyme_ii, plasmid) not in valid_pairs:
                raise ValueError(
                    f"{sequence_id}: selected RE pair is unsupported by {plasmid}"
                )
            record.update(
                {
                    "selected_plasmid": plasmid,
                    "selected_re_pair": f"{enzyme_i} / {enzyme_ii}",
                    "site_i_enzyme": enzyme_i,
                    "site_i_3mer_aa": three_mer_i,
                    "site_ii_enzyme": enzyme_ii,
                    "site_ii_3mer_aa": three_mer_ii,
                    "three_mer_aa_pair": f"{three_mer_i} / {three_mer_ii}",
                    "direction": direction,
                    "site_i_start_aa_0based": position_i,
                    "site_ii_start_aa_0based": position_ii,
                }
            )
        records.append(record)

    result = pd.DataFrame(records, columns=OUTPUT_COLUMNS)
    integer_columns = ["site_i_start_aa_0based", "site_ii_start_aa_0based"]
    result[integer_columns] = result[integer_columns].astype("Int64")
    compatible = result["hurdler_compatible"]
    if result.loc[compatible, "three_mer_aa_pair"].eq("").any():
        raise AssertionError("A compatible module is missing its 3-mer AA pair")
    if result.loc[~compatible, "three_mer_aa_pair"].ne("").any():
        raise AssertionError("An incompatible module unexpectedly has a 3-mer AA pair")
    return result


def export_module_3mer_results(
    source_path: str | Path = DEFAULT_SOURCE,
    output_path: str | Path = DEFAULT_OUTPUT,
    index_dir: str | Path = DEFAULT_INDEX,
    *,
    expected_counts: dict[str, int] | None = EXPECTED_PRODUCTION_COUNTS,
) -> pd.DataFrame:
    source = pd.read_csv(source_path, low_memory=False)
    result = build_module_3mer_results(
        source, index_dir, expected_counts=expected_counts
    )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(destination, index=False, lineterminator="\n")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Export the RepeatsDB/designed HURDLER result with selected 3-mer AA pairs."
    )
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--index-dir", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--expected-natural", type=int, default=25_913)
    parser.add_argument("--expected-designed", type=int, default=182)
    args = parser.parse_args(argv)
    result = export_module_3mer_results(
        args.source,
        args.output,
        args.index_dir,
        expected_counts={
            "Natural": args.expected_natural,
            "Designed": args.expected_designed,
        },
    )
    compatible = int(result["hurdler_compatible"].sum())
    print(
        f"Wrote {len(result):,} modules ({compatible:,} HURDLER-compatible) "
        f"to {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
