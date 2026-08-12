#!/usr/bin/env python3
"""Assemble complete primary shards with per-module rescue outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


FILES = (
    "module_hurdler_results.parquet",
    "module_hurdler_candidates.parquet",
    "optimized_constructs.parquet",
)


def complete(directory: Path) -> bool:
    return all((directory / name).is_file() and (directory / name).stat().st_size for name in FILES)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--primary-root", type=Path, required=True)
    parser.add_argument(
        "--rescue-root",
        type=Path,
        action="append",
        required=True,
        help="One or more per-module rescue roots, in preference order",
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--primary-shards", type=int, default=32)
    args = parser.parse_args()
    catalog = pd.read_parquet(args.catalog).reset_index(drop=True)
    source_by_shard: dict[str, str] = {}
    for shard in range(args.primary_shards):
        primary = args.primary_root / f"shard_{shard:03d}"
        module_indices = list(range(shard, len(catalog), args.primary_shards))
        if complete(primary):
            source_dirs = [primary]
            source_by_shard[str(shard)] = "primary"
        else:
            source_dirs = []
            missing = []
            rescue_roots_used: set[str] = set()
            for index in module_indices:
                candidates = [root / f"module_{index:03d}" for root in args.rescue_root]
                selected = next((path for path in candidates if complete(path)), None)
                if selected is None:
                    missing.append(" OR ".join(str(path) for path in candidates))
                else:
                    source_dirs.append(selected)
                    rescue_roots_used.add(str(selected.parent))
            if missing:
                raise FileNotFoundError("Incomplete rescue modules:\n" + "\n".join(missing))
            source_by_shard[str(shard)] = "per_module_rescue:" + ",".join(sorted(rescue_roots_used))
        destination = args.output_root / f"shard_{shard:03d}"
        destination.mkdir(parents=True, exist_ok=True)
        for filename in FILES:
            frames = [pd.read_parquet(directory / filename) for directory in source_dirs]
            combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
            combined.to_parquet(destination / filename, index=False)
    validation = {
        "catalog_modules": len(catalog),
        "primary_shards": args.primary_shards,
        "source_counts": pd.Series(source_by_shard).value_counts().to_dict(),
        "source_by_shard": source_by_shard,
        "passed": True,
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "assembly_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
