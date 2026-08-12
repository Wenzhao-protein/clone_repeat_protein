#!/usr/bin/env python3
"""Generate absolute, recoverable Digs tasks for arbitrary-DNA planning."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
from pathlib import Path

import pandas as pd


REPO = Path("/home/wendai/projects/hurdler/clone_repeat_protein")
STUDY = REPO / "studies" / "hurdler_validation"
SCRATCH = Path(
    "/net/scratch/wendai/projects/hurdler/clone_repeat_protein/"
    "studies/hurdler_validation"
)
PYTHON = Path("/home/wendai/.conda/envs/hurdler/bin/python")
VERSION = "arbitrary-dna-active-latent-v1"


def shell_line(arguments: list[str]) -> str:
    command = " ".join(shlex.quote(value) for value in arguments)
    return (
        f"set -o pipefail; export PYTHONPATH={shlex.quote(str(REPO / 'src'))}; "
        f"{command}; rc=$?; date -Is; exit \"$rc\""
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--run-name", default="run001_production")
    parser.add_argument("--shards", type=int, default=256)
    parser.add_argument("--without-idt", action="store_true")
    args = parser.parse_args()
    if args.shards < 1:
        raise ValueError("--shards must be positive")

    catalog = args.catalog.absolute()
    row_count = len(pd.read_parquet(catalog, columns=["target_id"]))
    shard_count = min(args.shards, row_count)
    local_run = STUDY / "step06_repetitive_dna_assembly" / "runs" / args.run_name
    scratch_run = SCRATCH / "step06_repetitive_dna_assembly" / "runs" / args.run_name
    task_dir = local_run / "taskfiles"
    task_dir.mkdir(parents=True, exist_ok=True)
    raw_root = scratch_run / "raw"
    commands: list[str] = []
    index_rows: list[dict[str, object]] = []
    shard_dirs: list[Path] = []
    for shard_index in range(shard_count):
        output_dir = raw_root / f"shard_{shard_index:05d}"
        shard_dirs.append(output_dir)
        arguments = [
            str(PYTHON),
            "-m",
            "hurdler",
            "dna-assembly",
            "plan",
            "--catalog",
            str(catalog),
            "--reference-dir",
            str(REPO / "data" / "reference_output"),
            "--artifact-dir",
            str(REPO / "output"),
            "--output-dir",
            str(output_dir),
            "--shard-index",
            str(shard_index),
            "--shard-count",
            str(shard_count),
        ]
        if not args.without_idt:
            arguments.extend(
                [
                    "--use-idt",
                    "--credential-mode",
                    "path",
                    "--credential-path",
                    str(Path.home() / ".config" / "hurdler" / "idt.env"),
                ]
            )
        commands.append(shell_line(arguments))
        index_rows.append(
            {
                "task_index": shard_index + 1,
                "shard_index": shard_index,
                "shard_count": shard_count,
                "expected_input_rows": len(range(shard_index, row_count, shard_count)),
                "output_dir": str(output_dir),
                "summary_output": str(output_dir / "dna_assembly_summary.parquet"),
                "idt_audit": str(output_dir / "idt_audit.jsonl") if not args.without_idt else "",
            }
        )
    (task_dir / "tasks.txt").write_text("\n".join(commands) + "\n")
    (task_dir / "shard_dirs.txt").write_text("\n".join(map(str, shard_dirs)) + "\n")
    with (task_dir / "task_index.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0]))
        writer.writeheader()
        writer.writerows(index_rows)

    smoke_catalog = scratch_run / "inputs" / "smoke_catalog.parquet"
    smoke_catalog.parent.mkdir(parents=True, exist_ok=True)
    pd.read_parquet(catalog).head(1).to_parquet(smoke_catalog, index=False)
    smoke_output = scratch_run / "smoke" / "raw"
    smoke_arguments = [
        str(PYTHON),
        "-m",
        "hurdler",
        "dna-assembly",
        "plan",
        "--catalog",
        str(smoke_catalog),
        "--reference-dir",
        str(REPO / "data" / "reference_output"),
        "--artifact-dir",
        str(REPO / "output"),
        "--output-dir",
        str(smoke_output),
    ]
    if not args.without_idt:
        smoke_arguments.extend(
            [
                "--use-idt",
                "--credential-mode",
                "path",
                "--credential-path",
                str(Path.home() / ".config" / "hurdler" / "idt.env"),
            ]
        )
    (task_dir / "smoke_tasks.txt").write_text(shell_line(smoke_arguments) + "\n")

    final_output = STUDY / "step06_repetitive_dna_assembly" / "tables" / VERSION
    finalize_arguments = [
        str(PYTHON),
        "-m",
        "hurdler",
        "dna-assembly",
        "finalize",
        "--shard-dir-list",
        str(task_dir / "shard_dirs.txt"),
        "--output-dir",
        str(final_output),
        "--figure-dir",
        str(STUDY / "step06_repetitive_dna_assembly" / "figures" / VERSION),
    ]
    (task_dir / "finalize_task.txt").write_text(shell_line(finalize_arguments) + "\n")
    run_manifest = {
        "version": VERSION,
        "catalog": str(catalog),
        "target_rows": row_count,
        "shard_count": shard_count,
        "idt_required": not args.without_idt,
        "resources": {
            "partition": "cpu",
            "cpu_per_task": 1,
            "memory": "8G",
            "walltime": "02:00:00",
            "array_throttle": 16,
            "group_size": 1,
            "container": "/net/software/containers/universal.sif",
        },
        "credential_path": "~/.config/hurdler/idt.env",
        "credential_contents_recorded": False,
    }
    (local_run / "run.json").write_text(json.dumps(run_manifest, indent=2) + "\n")
    print(json.dumps(run_manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
