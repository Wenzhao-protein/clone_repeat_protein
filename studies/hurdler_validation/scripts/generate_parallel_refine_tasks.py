#!/usr/bin/env python3
"""Generate recoverable Digs tasks for parallel adaptive GA/IDT refinement."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shlex
from pathlib import Path

import pandas as pd


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--constructs", type=Path, required=True)
    parser.add_argument("--groups", type=int, default=32)
    parser.add_argument("--total-shards", type=int, default=128)
    parser.add_argument("--shards-per-group", type=int, default=4)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--population-size", type=int, default=64)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    repo = args.repo.resolve()
    run_dir = args.run_dir.resolve()
    constructs = args.constructs.resolve()
    if args.groups * args.shards_per_group != args.total_shards:
        raise ValueError("groups * shards-per-group must equal total-shards")
    frame = pd.read_parquet(constructs)
    if len(frame) < args.total_shards:
        raise ValueError("Each scientific shard must receive at least one row")
    if frame.duplicated(["module_id", "fragment_limit_bp"]).any():
        raise ValueError("Construct rows must be unique by module and capacity")

    raw = run_dir / "raw"
    taskfiles = run_dir / "taskfiles"
    raw.mkdir(parents=True, exist_ok=True)
    taskfiles.mkdir(parents=True, exist_ok=True)
    wrapper = repo / "studies/hurdler_validation/scripts/run_with_idt_credentials.sh"
    group_runner = repo / "studies/hurdler_validation/scripts/run_parallel_refine_group.py"
    python = Path("/home/wendai/.conda/envs/hurdler/bin/python")
    hurdler = Path("/home/wendai/.conda/envs/hurdler/bin/hurdler")
    codon_usage = repo / "data/reference_output/codon_usage.csv"
    restriction_sites = repo / "data/reference_output/restriction_enzyme.csv"

    commands: list[str] = []
    index_rows: list[dict[str, object]] = []
    for group_index in range(args.groups):
        command = [
            str(wrapper),
            str(python),
            str(group_runner),
            "--hurdler",
            str(hurdler),
            "--constructs",
            str(constructs),
            "--output-root",
            str(raw),
            "--codon-usage",
            str(codon_usage),
            "--restriction-sites",
            str(restriction_sites),
            "--group-index",
            str(group_index),
            "--group-count",
            str(args.groups),
            "--total-shards",
            str(args.total_shards),
            "--shards-per-group",
            str(args.shards_per_group),
            "--max-workers",
            str(args.workers),
            "--population-size",
            str(args.population_size),
            "--short-generations",
            "10",
            "--generation-schedule",
            "10",
            "20",
            "40",
            "60",
            "80",
            "100",
            "--seed",
            str(args.seed),
            "--use-idt",
        ]
        rendered = shlex.join(command)
        commands.append(rendered)
        index_rows.append(
            {
                "task_index": group_index + 1,
                "case_id": f"soft_re_idt_group_{group_index:03d}",
                "expected_output": str(
                    raw / "groups" / f"group_{group_index:03d}.json"
                ),
                "command": rendered,
            }
        )

    (taskfiles / "tasks.txt").write_text("\n".join(commands) + "\n")
    with (taskfiles / "task_index.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0]))
        writer.writeheader()
        writer.writerows(index_rows)
    manifest = {
        "run_id": args.run_id,
        "status": "prepared",
        "scientific_rule_change": (
            "non-selected repeated RE sites are a soft GA score term and do not "
            "block IDT scoring"
        ),
        "ga_re_site_policy": (
            "nonselected-re-sites-soft-score-selected-sites-hard-v2"
        ),
        "constructs": str(constructs),
        "constructs_sha256": file_sha256(constructs),
        "expected_rows": len(frame),
        "expected_module_cap_keys": int(
            frame[["module_id", "fragment_limit_bp"]].drop_duplicates().shape[0]
        ),
        "groups": args.groups,
        "total_scientific_shards": args.total_shards,
        "shards_per_group": args.shards_per_group,
        "processes_per_group": args.workers,
        "population_size": args.population_size,
        "seed": args.seed,
        "credentials": "runtime IDT_CREDENTIAL_FILE (contents and resolved path excluded)",
        "resources": {
            "partition": "cpu",
            "cpu": args.workers,
            "mem": "8G",
            "time": "02:00:00",
        },
    }
    (run_dir / "run.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
