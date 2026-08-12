#!/usr/bin/env python3
"""Select the middle module from each complete RepeatsDB repeat region."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from hurdler.modules import apply_natural_middle_unit_annotations, merge_module_catalogs
from hurdler.periodicity import BOUNDARY_METHOD_VERSION


def selected_is_geometric_middle(row: pd.Series) -> bool:
    coordinates = [(int(value[0]), int(value[1])) for value in json.loads(row.source_unit_coordinates_json)]
    midpoint = (int(row.repeat_region_start) + int(row.repeat_region_end)) / 2
    expected = min(
        range(len(coordinates)),
        key=lambda index: (
            abs((coordinates[index][0] + coordinates[index][1]) / 2 - midpoint),
            coordinates[index][0],
        ),
    )
    return (
        int(row.selected_module_index) == expected + 1
        and int(row.selected_module_start) == coordinates[expected][0]
        and int(row.selected_module_end) == coordinates[expected][1]
        and row.unit_sequence == json.loads(row.unit_sequences_json)[expected]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--target-count", type=int, default=100)
    args = parser.parse_args()
    source = pd.read_parquet(args.input)
    pool_audit = args.output.with_name(args.output.stem + "_pool_boundary_audit.parquet")
    pool_units = args.output.with_name(args.output.stem + "_pool_unit_alignment.parquet")
    pool_positions = args.output.with_name(
        args.output.stem + "_pool_position_variability.parquet"
    )
    refined = apply_natural_middle_unit_annotations(
        source,
        output_path=pool_audit,
        unit_alignment_path=pool_units,
        position_variability_path=pool_positions,
        workers=args.workers,
    )
    eligible = refined.loc[
        refined.boundary_refinement_status.eq("source_annotation_middle_unit")
    ].copy()
    eligible = eligible.sort_values(
        ["family", "reviewed", "source_accession", "source_chain", "unit_start", "module_id"],
        ascending=[True, False, True, True, True, True],
    ).drop_duplicates("unit_sequence", keep="first")
    if len(eligible) < args.target_count:
        raise RuntimeError(
            f"Only {len(eligible)} unique middle units are available; target={args.target_count}"
        )
    selected = eligible.head(args.target_count).copy()
    selected.to_parquet(
        args.output.with_name(args.output.stem + "_boundary_audit.parquet"), index=False
    )
    selected.to_csv(
        args.output.with_name(args.output.stem + "_boundary_audit.csv"), index=False
    )
    selected_ids = set(selected.module_id)
    for pool_path, suffix in (
        (pool_units, "_unit_alignment"),
        (pool_positions, "_position_variability"),
    ):
        selected_table = pd.read_parquet(pool_path)
        selected_table = selected_table.loc[selected_table.module_id.isin(selected_ids)].copy()
        destination = args.output.with_name(args.output.stem + suffix + ".parquet")
        selected_table.to_parquet(destination, index=False)
        selected_table.to_csv(destination.with_suffix(".csv"), index=False)
    catalog = merge_module_catalogs([selected], args.output)
    validation = {
        "boundary_method_version": BOUNDARY_METHOD_VERSION,
        "source_pool_rows": len(source),
        "resolved_pool_rows": len(eligible),
        "target_count": args.target_count,
        "catalog_rows": len(catalog),
        "middle_unit_resolved_in_pool": int(
            refined.boundary_refinement_status.eq("source_annotation_middle_unit").sum()
        ),
        "selected_rows": len(selected),
        "selected_unique_sequences": int(selected.unit_sequence.nunique()),
        "all_have_multiple_aligned_units": bool(selected.repeat_count.ge(2).all()),
        "all_middle_sequences_match_length": bool(
            selected.unit_sequence.str.len().eq(selected.unit_length).all()
        ),
        "all_selected_modules_are_middle": bool(
            selected.apply(selected_is_geometric_middle, axis=1).all()
            and selected.module_selection_policy.eq(
                "repeat-region-middle-unit-tie-earlier-v1"
            ).all()
        ),
        "all_auth_label_chains_resolved": bool(
            selected.rcsb_auth_chain_id.astype(str).str.len().gt(0).all()
            and selected.rcsb_label_chain_id.astype(str).str.len().gt(0).all()
        ),
        "output": str(args.output),
    }
    validation["passed"] = bool(
        validation["catalog_rows"] == args.target_count
        and validation["selected_rows"] == args.target_count
        and validation["selected_unique_sequences"] == args.target_count
        and validation["all_have_multiple_aligned_units"]
        and validation["all_middle_sequences_match_length"]
        and validation["all_selected_modules_are_middle"]
        and validation["all_auth_label_chains_resolved"]
    )
    args.output.with_name(args.output.stem + "_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0 if validation["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
