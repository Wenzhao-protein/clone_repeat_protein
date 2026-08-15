#!/usr/bin/env python3
"""Create a missing-only task file from a V2 production task index."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-index", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.task_index.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    missing = [row for row in rows if not Path(row["expected_output"]).exists()]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite missing-only task file: {args.output}")
    args.output.write_text("".join(row["command"] + "\n" for row in missing))
    print(json.dumps({"source_tasks": len(rows), "missing_tasks": len(missing), "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
