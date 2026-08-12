#!/usr/bin/env python3
"""Digs smoke and deterministic multiprocessing check for module boundaries."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import pandas as pd

from hurdler.modules import parse_dhr_supplement, refine_module_boundaries


COMPARISON_COLUMNS = [
    "module_id",
    "unit_sequence",
    "unit_start",
    "unit_end",
    "primitive_period",
    "repeat_count",
    "fixed_mask",
    "variable_positions_json",
    "periodicity_score",
]


def frame_hash(frame: pd.DataFrame) -> str:
    payload = frame[COMPARISON_COLUMNS].sort_values("module_id").to_csv(index=False)
    return hashlib.sha256(payload.encode()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dhr-text", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = parse_dhr_supplement(args.dhr_text)
    smoke = source[source.source_accession.isin(["DHR12", "DHR13", "DHR17", "DHR18"])].copy()
    if len(smoke) != 4:
        raise RuntimeError(f"Expected four smoke modules, found {len(smoke)}")

    start = time.perf_counter()
    serial = refine_module_boundaries(smoke, workers=1)
    serial_seconds = time.perf_counter() - start
    start = time.perf_counter()
    parallel = refine_module_boundaries(smoke, workers=4)
    parallel_seconds = time.perf_counter() - start
    serial_hash = frame_hash(serial)
    parallel_hash = frame_hash(parallel)
    if serial_hash != parallel_hash:
        raise RuntimeError("Serial and multiprocessing boundary outputs differ")
    parallel.to_parquet(args.output_dir / "boundary_smoke.parquet", index=False)
    parallel.to_csv(args.output_dir / "boundary_smoke.csv", index=False)
    report = {
        "rows": len(parallel),
        "serial_seconds": serial_seconds,
        "parallel_seconds": parallel_seconds,
        "speedup": serial_seconds / parallel_seconds,
        "serial_sha256": serial_hash,
        "parallel_sha256": parallel_hash,
        "deterministic": True,
        "all_refined": bool(parallel.boundary_refinement_status.eq("refined").all()),
        "primitive_periods": dict(
            zip(parallel.source_accession, parallel.primitive_period.astype(int), strict=True)
        ),
    }
    report["passed"] = report["deterministic"] and report["all_refined"]
    (args.output_dir / "benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
