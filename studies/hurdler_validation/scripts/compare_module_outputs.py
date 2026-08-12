#!/usr/bin/env python3
"""Verify that two deterministic module-optimization runs are identical."""

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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    row_counts: dict[str, int] = {}
    for filename in FILES:
        expected = pd.read_parquet(args.expected / filename)
        observed = pd.read_parquet(args.observed / filename)
        pd.testing.assert_frame_equal(expected, observed, check_like=False)
        row_counts[filename] = len(observed)
    payload = {"passed": True, "files": list(FILES), "row_counts": row_counts}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
