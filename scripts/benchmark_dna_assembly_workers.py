#!/usr/bin/env python3
"""Measure serial versus process-parallel local planning on one Digs node."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from hurdler.dna_assembly import TargetRecord, load_enzyme_catalog, plan_target


_GEOMETRIES = None
_PLASMIDS = None


def initialize(reference_dir: str, artifact_dir: str) -> None:
    global _GEOMETRIES, _PLASMIDS
    _GEOMETRIES, _PLASMIDS = load_enzyme_catalog(reference_dir, artifact_dir=artifact_dir)


def evaluate(payload: dict[str, object]) -> tuple[str, str]:
    target = TargetRecord(
        target_id=str(payload["target_id"]),
        sequence=str(payload["sequence"]),
        cohort=str(payload.get("cohort", "benchmark")),
        architecture=str(payload.get("architecture", "benchmark")),
    )
    result = plan_target(target, _GEOMETRIES, _PLASMIDS, require_idt=False)
    summary = result["summary"].iloc[0].to_dict()
    stable = {
        "target_id": target.target_id,
        "hurdler_compatible": bool(summary["hurdler_compatible"]),
        "candidate_pair_count": int(summary["candidate_pair_count"]),
        "best_route_id": str(summary["best_route_id"]),
    }
    digest = hashlib.sha256(json.dumps(stable, sort_keys=True).encode()).hexdigest()
    return target.target_id, digest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=80)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 8, 16])
    args = parser.parse_args()
    frame = pd.read_parquet(args.catalog)
    frame["length_bp"] = frame.sequence.str.len()
    frame = frame.sort_values(["length_bp", "target_id"], kind="mergesort").reset_index(drop=True)
    indexes = pd.Series(range(len(frame))).sample(
        n=min(args.sample_size, len(frame)), random_state=42
    ).sort_values()
    payloads = frame.iloc[indexes][
        ["target_id", "sequence", "cohort", "architecture"]
    ].to_dict("records")
    runs = []
    reference: dict[str, str] | None = None
    for workers in args.workers:
        started = time.perf_counter()
        if workers == 1:
            initialize(str(args.reference_dir), str(args.artifact_dir))
            results = [evaluate(payload) for payload in payloads]
        else:
            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=initialize,
                initargs=(str(args.reference_dir), str(args.artifact_dir)),
            ) as executor:
                results = list(executor.map(evaluate, payloads, chunksize=1))
        elapsed = time.perf_counter() - started
        observed = dict(results)
        if reference is None:
            reference = observed
        if observed != reference:
            raise RuntimeError(f"Scientific outputs differ at worker count {workers}")
        runs.append(
            {
                "workers": workers,
                "runtime_seconds": elapsed,
                "targets_per_second": len(payloads) / elapsed,
                "speedup_vs_serial": runs[0]["runtime_seconds"] / elapsed if runs else 1.0,
            }
        )
    best = max(runs, key=lambda row: row["targets_per_second"])
    payload = {
        "version": "arbitrary-dna-active-latent-v1",
        "host_cpu_count": os.cpu_count(),
        "sample_rows": len(payloads),
        "scientific_equivalence": True,
        "runs": runs,
        "fastest_worker_count": best["workers"],
        "production_decision": (
            "taskrunner_array_parallelism_no_nested_pool"
            if best["workers"] > 1
            else "serial_tasks"
        ),
        "production_array_concurrency": 16,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
