#!/usr/bin/env python3
"""Benchmark v2 molecular routing with 1/2/4/8/16 local processes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from hurdler.complete_route import search_complete_repeat_routes
from hurdler.dna_assembly import (
    DNA_COMPLETE_ROUTE_VERSION,
    TargetRecord,
    load_enzyme_catalog,
    validate_dna,
)


_GEOMETRIES = None
_PLASMIDS = None


def initialize(reference_dir: str, artifact_dir: str) -> None:
    global _GEOMETRIES, _PLASMIDS
    _GEOMETRIES, _PLASMIDS = load_enzyme_catalog(
        reference_dir, artifact_dir=artifact_dir
    )


def evaluate(payload: dict[str, object]) -> tuple[str, str]:
    unit = validate_dna(str(payload["unit_sequence"]))
    target = TargetRecord(
        target_id=str(payload["element_id"]),
        sequence=unit,
        cohort="real_element_derived",
        architecture="exact_tandem",
        source_database=str(payload["source_database"]),
        element_id=str(payload["element_id"]),
        unit_sequence=unit,
        copy_count=1,
    )
    result = search_complete_repeat_routes(
        target,
        _GEOMETRIES,
        _PLASMIDS,
        require_idt=False,
    )
    stable = result["targets"][[
        "target_copy_count", "complete_route_verified", "failure_reason",
        "hurdler_step_count", "plasmid",
    ]].fillna("").to_dict("records")
    key = f"{payload['source_database']}|{payload['element_id']}"
    digest = hashlib.sha256(
        json.dumps(stable, sort_keys=True).encode()
    ).hexdigest()
    return key, digest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=16)
    parser.add_argument(
        "--include-element",
        default="",
        help="Force one element_id into the deterministic benchmark sample.",
    )
    parser.add_argument(
        "--workers", type=int, nargs="+", default=[1, 2, 4, 8, 16]
    )
    args = parser.parse_args()
    frame = pd.read_parquet(args.catalog)
    elements = (
        frame.loc[frame.cohort.eq("real_element_derived")]
        .sort_values(
            ["source_database", "element_id", "copy_count"],
            kind="mergesort",
        )
        .drop_duplicates(["source_database", "element_id"])
        .reset_index(drop=True)
    )
    # Use a deterministic length-stratified sample rather than only the many
    # ~30-bp CRISPR repeats.
    elements["unit_length_bp"] = elements.unit_sequence.str.len()
    ordered = elements.sort_values(
        ["unit_length_bp", "source_database", "element_id"], kind="mergesort"
    ).reset_index(drop=True)
    count = min(args.sample_size, len(ordered))
    positions = sorted(
        set(round(index * (len(ordered) - 1) / max(1, count - 1)) for index in range(count))
    )
    payloads = ordered.iloc[positions][
        ["source_database", "element_id", "unit_sequence"]
    ].to_dict("records")
    if args.include_element and not any(
        str(row["element_id"]) == args.include_element for row in payloads
    ):
        forced = elements.loc[
            elements.element_id.astype(str).eq(args.include_element),
            ["source_database", "element_id", "unit_sequence"],
        ]
        if forced.empty:
            raise ValueError(
                f"Benchmark include-element not found: {args.include_element}"
            )
        payloads[-1] = forced.iloc[0].to_dict()
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
            raise RuntimeError(
                f"Scientific outputs differ at worker count {workers}"
            )
        runs.append(
            {
                "workers": workers,
                "runtime_seconds": elapsed,
                "elements_per_second": len(payloads) / elapsed,
                "speedup_vs_serial": (
                    runs[0]["runtime_seconds"] / elapsed if runs else 1.0
                ),
            }
        )
    fastest = max(runs, key=lambda row: row["elements_per_second"])
    report = {
        "version": DNA_COMPLETE_ROUTE_VERSION,
        "host_cpu_count": os.cpu_count(),
        "sample_elements": len(payloads),
        "forced_element": args.include_element,
        "scientific_equivalence": True,
        "runs": runs,
        "fastest_worker_count": fastest["workers"],
        "production_decision": "16-way taskrunner array; no nested pool",
        "production_array_concurrency": 16,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
