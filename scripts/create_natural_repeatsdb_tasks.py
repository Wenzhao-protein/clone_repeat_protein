#!/usr/bin/env python3
"""Create resumable RepeatsDB-direct CPU shard and finalization taskfiles."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PYTHON_ENV = Path("/home/wendai/.conda/envs/hurdler/bin/hurdler")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--scratch-run-dir", type=Path, required=True)
    parser.add_argument("--final-catalog", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=128)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.shards < 1 or args.workers < 1:
        raise ValueError("shards and workers must be positive")
    inventory = args.inventory.resolve()
    run_dir = args.run_dir.resolve()
    task_dir = run_dir / "taskfiles"
    raw_dir = args.scratch_run_dir.resolve() / "raw"
    cache_dir = raw_dir / "fasta_cache"
    task_dir.mkdir(parents=True, exist_ok=True)
    task_lines: list[str] = []
    mapping_paths: list[Path] = []
    inventory_paths: list[Path] = []
    exclusion_paths: list[Path] = []
    index_rows: list[dict[str, object]] = []
    for shard_index in range(args.shards):
        shard_dir = raw_dir / f"shard_{shard_index:03d}"
        output = shard_dir / "natural.parquet"
        mapping = shard_dir / "natural_source_mappings.parquet"
        region_inventory = shard_dir / "natural_region_inventory.parquet"
        exclusion = shard_dir / "natural_exclusions.csv"
        command = " ".join(
            [
                str(PYTHON_ENV),
                "curate-modules",
                "--all-repeatsdb",
                "--one-per-protein",
                "--annotation-inventory",
                str(inventory),
                "--natural-output",
                str(output),
                "--natural-mappings-output",
                str(mapping),
                "--natural-exclusions-output",
                str(exclusion),
                "--natural-cache-dir",
                str(cache_dir),
                "--natural-workers",
                str(args.workers),
                "--natural-shard-index",
                str(shard_index),
                "--natural-shard-count",
                str(args.shards),
            ]
        )
        task_lines.append(command)
        mapping_paths.append(mapping)
        inventory_paths.append(region_inventory)
        exclusion_paths.append(exclusion)
        index_rows.append(
            {
                "task_id": shard_index + 1,
                "shard_index": shard_index,
                "shard_count": args.shards,
                "workers": args.workers,
                "output": output,
                "source_mappings": mapping,
                "region_inventory": region_inventory,
                "exclusions": exclusion,
            }
        )
    (task_dir / "tasks.txt").write_text("\n".join(task_lines) + "\n")
    with (task_dir / "task_index.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0]))
        writer.writeheader()
        writer.writerows(index_rows)
    final_command = " ".join(
        [
            str(PYTHON_ENV),
            "curate-modules",
            "--natural-output",
            str(args.final_catalog.resolve()),
            "--annotation-inventory",
            str(inventory),
            "--finalize-natural-mappings",
            *(str(path) for path in mapping_paths),
            "--finalize-natural-inventories",
            *(str(path) for path in inventory_paths),
            "--finalize-natural-exclusions",
            *(str(path) for path in exclusion_paths),
        ]
    )
    (task_dir / "finalize_task.txt").write_text(final_command + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
