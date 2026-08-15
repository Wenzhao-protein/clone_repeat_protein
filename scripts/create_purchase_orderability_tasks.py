#!/usr/bin/env python3
"""Create the smoke and production Digs tasks for purchase orderability."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--scratch-run", type=Path, required=True)
    parser.add_argument("--local-run", type=Path, required=True)
    parser.add_argument("--credential-path", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.absolute()
    raw = args.raw_root.absolute()
    scratch = args.scratch_run.absolute()
    local = args.local_run.absolute()
    credential = args.credential_path.absolute()
    task_dir = local / "taskfiles"
    task_dir.mkdir(parents=True, exist_ok=True)
    common = (
        "set -o pipefail; "
        f"export PYTHONPATH={repo / 'src'}; "
        "export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; "
        "export OPENBLAS_NUM_THREADS=1; "
    )
    python = Path.home() / ".conda/envs/hurdler/bin/python"
    smoke = (
        common
        + f"{python} {repo / 'scripts/smoke_purchase_orderability.py'} "
        + f"--raw-root {raw} --output-dir {scratch / 'smoke'} "
        + f"--credential-path {credential}; rc=$?; date -Is; exit \"$rc\""
    )
    production = (
        common
        + f"{python} -m hurdler dna-assembly audit-purchases "
        + f"--raw-root {raw} --output-dir {scratch / 'production'} "
        + "--expected-shards 512 --expected-routes 15535 --expected-elements 3129 "
        + f"--use-idt --credential-mode path --credential-path {credential}; "
        + "rc=$?; date -Is; exit \"$rc\""
    )
    (task_dir / "smoke_tasks.txt").write_text(smoke + "\n")
    (task_dir / "production_tasks.txt").write_text(production + "\n")
    (local / "run.json").write_text(
        json.dumps(
            {
                "version": "complete-route-purchase-orderability-v1",
                "input_run": raw.parent.name,
                "expected_shards": 512,
                "expected_found_routes": 15535,
                "expected_elements_with_routes": 3129,
                "credential_source": "external_owner_only_env",
                "credential_path_recorded": False,
                "resources": {
                    "partition": "cpu",
                    "cpu_per_task": 1,
                    "memory": "8G",
                    "walltime": "00:30:00",
                    "group_size": 1,
                    "array_throttle": 1,
                    "nested_multiprocessing": False,
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(task_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
