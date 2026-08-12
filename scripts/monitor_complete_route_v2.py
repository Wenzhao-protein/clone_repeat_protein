#!/usr/bin/env python3
"""Summarize recoverable complete-route production without loading route hits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=512)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--task-file", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    completed = set()
    for manifest_path in sorted(args.raw_root.glob("shard_*/complete_route_manifest.json")):
        manifest = json.loads(manifest_path.read_text())
        shard = int(manifest["shard_index"])
        completed.add(shard)
        target_path = manifest_path.parent / "complete_route_targets.parquet"
        targets = pd.read_parquet(
            target_path,
            columns=[
                "complete_route_verified", "failure_reason",
                "whole_target_idt_status",
            ],
        )
        rows.append(
            {
                "shard_index": shard,
                "element_rows": int(manifest["element_rows"]),
                "target_rows": len(targets),
                "complete_targets": int(targets.complete_route_verified.sum()),
                "api_unclassified_targets": int(
                    targets.whole_target_idt_status.isin(
                        ["api_failure", "api_unclassified", "scored_unclassified"]
                    ).sum()
                ),
            }
        )
    missing = sorted(set(range(args.shard_count)) - completed)
    commands = args.task_file.read_text().splitlines()
    report = {
        "version": "arbitrary-dna-complete-route-v2",
        "expected_shards": args.shard_count,
        "completed_shards": len(completed),
        "missing_shard_count": len(missing),
        "missing_shards": missing,
        "element_rows": sum(row["element_rows"] for row in rows),
        "target_rows": sum(row["target_rows"] for row in rows),
        "complete_targets": sum(row["complete_targets"] for row in rows),
        "api_unclassified_targets": sum(
            row["api_unclassified_targets"] for row in rows
        ),
        "missing_only_rerun_commands": [commands[index] for index in missing],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value for key, value in report.items()
                if key not in {"missing_shards", "missing_only_rerun_commands"}
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
