#!/usr/bin/env python3
"""Generate missing-shard, one-module rescue tasks for module optimization."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shlex
from pathlib import Path

import pandas as pd


REQUIRED_OUTPUTS = (
    "module_hurdler_results.parquet",
    "module_hurdler_candidates.parquet",
    "optimized_constructs.parquet",
)


def complete(directory: Path) -> bool:
    return all((directory / name).is_file() and (directory / name).stat().st_size for name in REQUIRED_OUTPUTS)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--primary-root", type=Path, required=True)
    parser.add_argument("--rescue-root", type=Path, required=True)
    parser.add_argument(
        "--completed-root",
        type=Path,
        action="append",
        default=[],
        help="Skip modules already complete in any of these earlier rescue roots",
    )
    parser.add_argument("--task-dir", type=Path, required=True)
    parser.add_argument("--primary-shards", type=int, default=32)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    repo = args.repo.absolute()
    catalog_path = repo / "studies/hurdler_validation/step03_module_corpus/tables/module_catalog.parquet"
    catalog = pd.read_parquet(catalog_path).reset_index(drop=True)
    missing_primary = {
        shard
        for shard in range(args.primary_shards)
        if not complete(args.primary_root / f"shard_{shard:03d}")
    }
    cases = [
        (index, str(row.module_id))
        for index, row in catalog.iterrows()
        if index % args.primary_shards in missing_primary
        and not any(complete(root / f"module_{index:03d}") for root in args.completed_root)
    ]
    # Avoid putting all longest/similar modules next to each other in the array.
    cases.sort(key=lambda item: hashlib.sha256(item[1].encode()).hexdigest())
    hurdler = Path("/home/wendai/.conda/envs/hurdler/bin/hurdler")
    index_dir = Path(
        "/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/"
        "step01_reference_lookup/runs/run01_production/raw/legacy-optimized-v1"
    )
    task_rows = []
    for catalog_index, module_id in cases:
        command = [
            str(hurdler), "optimize-modules",
            "--catalog", str(catalog_path),
            "--index-dir", str(index_dir),
            "--output-dir", str(args.rescue_root / f"module_{catalog_index:03d}"),
            "--fragment-limits", "1800", "3000",
            "--codon-usage", str(repo / "data/reference_output/codon_usage.csv"),
            "--shard-index", str(catalog_index),
            "--shard-count", str(len(catalog)),
            "--workers", str(args.workers),
        ]
        task_rows.append((catalog_index, module_id, command))
    args.task_dir.mkdir(parents=True, exist_ok=True)
    with (args.task_dir / "tasks.txt").open("w") as handle:
        for _index, _module_id, command in task_rows:
            rendered = " ".join(shlex.quote(part) for part in command)
            handle.write(f"set -o pipefail; {rendered}; rc=$?; date -Is; exit \"$rc\"\n")
    with (args.task_dir / "task_index.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["task_index", "catalog_index", "module_id", "source_primary_shard", "output_dir"])
        for task_index, (catalog_index, module_id, _command) in enumerate(task_rows, start=1):
            writer.writerow([
                task_index,
                catalog_index,
                module_id,
                catalog_index % args.primary_shards,
                args.rescue_root / f"module_{catalog_index:03d}",
            ])
    print({"missing_primary_shards": sorted(missing_primary), "rescue_tasks": len(task_rows)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
