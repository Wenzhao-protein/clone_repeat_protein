#!/usr/bin/env python3
"""Create recoverable Digs task files for complete-route v2 production."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


VERSION = "arbitrary-dna-complete-route-v2"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--scratch-run", type=Path, required=True)
    parser.add_argument("--local-run", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=512)
    parser.add_argument(
        "--credential-path",
        type=Path,
        default=Path("~/.config/hurdler/idt.env").expanduser(),
    )
    args = parser.parse_args()
    repo = args.repo.absolute()
    catalog = args.catalog.absolute()
    scratch = args.scratch_run.absolute()
    local = args.local_run.absolute()
    task_dir = local / "taskfiles"
    task_dir.mkdir(parents=True, exist_ok=True)
    tasks = []
    shard_dirs = []
    for shard in range(args.shards):
        output = scratch / "raw" / f"shard_{shard:05d}"
        shard_dirs.append(str(output))
        tasks.append(
            "set -o pipefail; "
            f"export PYTHONPATH={repo / 'src'}; "
            "export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; "
            "export OPENBLAS_NUM_THREADS=1; "
            f"{Path.home() / '.conda/envs/hurdler/bin/python'} -m hurdler "
            "dna-assembly plan-complete "
            f"--catalog {catalog} "
            f"--reference-dir {repo / 'data/reference_output'} "
            f"--artifact-dir {repo / 'data/artifacts'} "
            f"--output-dir {output} "
            f"--shard-index {shard} --shard-count {args.shards} "
            "--use-idt --credential-mode path "
            f"--credential-path {args.credential_path.absolute()}; "
            "rc=$?; date -Is; exit \"$rc\""
        )
    (task_dir / "tasks.txt").write_text("\n".join(tasks) + "\n")
    (task_dir / "shard_dirs.txt").write_text(
        "\n".join(shard_dirs) + "\n"
    )
    finalize = (
        "set -o pipefail; "
        f"export PYTHONPATH={repo / 'src'}; "
        f"{Path.home() / '.conda/envs/hurdler/bin/python'} -m hurdler "
        "dna-assembly finalize-complete "
        f"--shard-dir-list {task_dir / 'shard_dirs.txt'} "
        f"--output-dir {repo / 'studies/hurdler_validation/step06_repetitive_dna_assembly/tables' / VERSION / 'production'} "
        "--expected-elements 29042 --expected-real-targets 145210 "
        f"--figure-dir {repo / 'studies/hurdler_validation/step06_repetitive_dna_assembly/figures' / VERSION}; "
        "rc=$?; date -Is; exit \"$rc\""
    )
    (task_dir / "finalize_task.txt").write_text(finalize + "\n")
    (local / "run.json").write_text(
        json.dumps(
            {
                "version": VERSION,
                "catalog": str(catalog),
                "scratch_run": str(scratch),
                "shard_count": args.shards,
                "submission_status": "taskfiles_created",
                "idt_required": True,
                "credential_path": "~/.config/hurdler/idt.env",
                "credential_contents_recorded": False,
                "resources": {
                    "partition": "cpu",
                    "cpu_per_task": 1,
                    "memory": "8G",
                    "walltime": "02:00:00",
                    "array_throttle": 16,
                    "nested_multiprocessing": False,
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(task_dir / "tasks.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
