#!/usr/bin/env python3
"""Create an exact missing-only recovery taskfile for Stage 2.

The production task index is the authority.  A shard is complete only when
its result Parquet opens and contains both fragment-capacity rows for every
module assigned to that task.  Partial IDT JSONL files are deliberately kept:
the rerun's exact-DNA cache reuses every valid response already obtained.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-index", type=Path, required=True)
    parser.add_argument("--tasks", type=Path, required=True)
    parser.add_argument("--output-taskfile", type=Path, required=True)
    parser.add_argument("--output-index", type=Path, required=True)
    parser.add_argument(
        "--task-id-max",
        type=int,
        help="Restrict recovery to original task IDs at or below this value",
    )
    args = parser.parse_args()

    index = pd.read_csv(args.task_index)
    commands = args.tasks.read_text().splitlines()
    if len(commands) != len(index):
        raise ValueError(
            f"Task/index length mismatch: {len(commands)} != {len(index)}"
        )
    if args.task_id_max is not None:
        if args.task_id_max < 1:
            raise ValueError("--task-id-max must be positive")
        index = index.loc[index.task_id.astype(int).le(args.task_id_max)].copy()

    missing_rows: list[dict[str, object]] = []
    missing_commands: list[str] = []
    for row in index.to_dict(orient="records"):
        position = int(row["task_id"]) - 1
        result = Path(str(row["result_output"]))
        expected_rows = 2 * int(row["module_count"])
        reason = ""
        observed_rows = 0
        if not result.is_file():
            reason = "missing_result_parquet"
        else:
            try:
                observed_rows = int(pq.ParquetFile(result).metadata.num_rows)
            except Exception as exc:
                reason = f"unreadable_result_parquet:{type(exc).__name__}"
            else:
                if observed_rows != expected_rows:
                    reason = (
                        f"result_row_count_mismatch:{observed_rows}!={expected_rows}"
                    )
        if reason:
            missing_commands.append(commands[position])
            missing_rows.append(
                {
                    **row,
                    "original_task_line": position + 1,
                    "expected_result_rows": expected_rows,
                    "observed_result_rows": observed_rows,
                    "recovery_reason": reason,
                }
            )

    args.output_taskfile.parent.mkdir(parents=True, exist_ok=True)
    args.output_taskfile.write_text(
        "\n".join(missing_commands) + ("\n" if missing_commands else "")
    )
    recovery = pd.DataFrame(missing_rows)
    recovery.to_csv(args.output_index, index=False)
    summary = {
        "source_tasks_considered": len(index),
        "complete_tasks": len(index) - len(recovery),
        "missing_tasks": len(recovery),
        "output_taskfile": str(args.output_taskfile.absolute()),
        "output_index": str(args.output_index.absolute()),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
