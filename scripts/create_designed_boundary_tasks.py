#!/usr/bin/env python3
"""Create one-row, recoverable designed DSSP/Foldseek CPU taskfiles."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


HURDLER = Path("/home/wendai/.conda/envs/hurdler/bin/hurdler")
MKDSSP = Path("/home/wendai/.conda/envs/hurdler/bin/mkdssp")
MAFFT = Path("/home/wendai/.conda/envs/hurdler/bin/mafft")
FOLDSEEK = Path("/net/software/utils/foldseek")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--scratch-run-dir", type=Path, required=True)
    parser.add_argument("--final-catalog", type=Path, required=True)
    args = parser.parse_args()
    inventory = args.inventory.resolve()
    source = pd.read_parquet(inventory).sort_values(
        "module_id", kind="mergesort"
    ).reset_index(drop=True)
    run_dir = args.run_dir.resolve()
    task_dir = run_dir / "taskfiles"
    raw_dir = args.scratch_run_dir.resolve() / "raw"
    task_dir.mkdir(parents=True, exist_ok=True)
    task_lines: list[str] = []
    mappings: list[Path] = []
    exclusions: list[Path] = []
    candidates: list[Path] = []
    units: list[Path] = []
    positions: list[Path] = []
    index_rows: list[dict[str, object]] = []
    shard_count = len(source)
    for shard_index, row in source.iterrows():
        shard_dir = raw_dir / f"shard_{shard_index:03d}"
        catalog = shard_dir / "designed.parquet"
        candidate = shard_dir / "candidates.parquet"
        unit = shard_dir / "units.parquet"
        position = shard_dir / "positions.parquet"
        exclusion = shard_dir / "exclusions.parquet"
        command = " ".join(
            [
                str(HURDLER),
                "infer-designed-boundaries",
                "--input",
                str(inventory),
                "--output",
                str(catalog),
                "--candidates-output",
                str(candidate),
                "--units-output",
                str(unit),
                "--positions-output",
                str(position),
                "--exclusions-output",
                str(exclusion),
                "--dssp-engine",
                "biotite",
                "--mkdssp",
                str(MKDSSP),
                "--foldseek",
                str(FOLDSEEK),
                "--mafft",
                str(MAFFT),
                "--shard-index",
                str(shard_index),
                "--shard-count",
                str(shard_count),
            ]
        )
        task_lines.append(command)
        mappings.append(catalog.with_name("designed_source_mappings.parquet"))
        exclusions.append(exclusion)
        candidates.append(candidate)
        units.append(unit)
        positions.append(position)
        index_rows.append(
            {
                "task_id": shard_index + 1,
                "shard_index": shard_index,
                "shard_count": shard_count,
                "module_id": row.module_id,
                "structure_inventory_status": row.structure_inventory_status,
                "author_structure_path": row.author_structure_path,
                "af3_structure_path": row.af3_structure_path,
                "output": catalog,
            }
        )
    (task_dir / "tasks.txt").write_text("\n".join(task_lines) + "\n")
    with (task_dir / "task_index.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0]))
        writer.writeheader()
        writer.writerows(index_rows)
    final_command = " ".join(
        [
            str(HURDLER),
            "infer-designed-boundaries",
            "--output",
            str(args.final_catalog.resolve()),
            "--finalize-mappings",
            *(str(path) for path in mappings),
            "--finalize-exclusions",
            *(str(path) for path in exclusions),
            "--finalize-candidate-tables",
            *(str(path) for path in candidates),
            "--finalize-unit-tables",
            *(str(path) for path in units),
            "--finalize-position-tables",
            *(str(path) for path in positions),
        ]
    )
    (task_dir / "finalize_task.txt").write_text(final_command + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
