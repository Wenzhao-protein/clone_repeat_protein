#!/usr/bin/env python3
"""Prepare larger-cap recovery rows from proven smaller-cap constructs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from hurdler.optimization import translate_dna


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--optimized", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int)
    args = parser.parse_args()

    frame = pd.read_parquet(args.optimized)
    by_key = {
        (str(row.module_id), int(row.fragment_limit_bp)): row._asdict()
        for row in frame.itertuples(index=False)
    }
    rows: list[dict[str, object]] = []
    for module_id in sorted(frame.module_id.astype(str).unique()):
        smaller = by_key[(module_id, 1800)]
        larger = by_key[(module_id, 3000)]
        smaller_copies = int(smaller.get("verified_max_copies", 0) or 0)
        larger_copies = int(larger.get("verified_max_copies", 0) or 0)
        if not bool(smaller.get("final_passed")) or smaller_copies <= larger_copies:
            continue
        dna = smaller.get("dna_sequence")
        if not isinstance(dna, str) or not dna:
            raise RuntimeError(f"{module_id}: smaller-cap passing DNA is missing")
        payload = dict(smaller)
        payload.update(
            fragment_limit_bp=3000,
            mathematical_max_copies=int(larger["mathematical_max_copies"]),
            verified_max_copies=smaller_copies,
            known_orderable_copies=smaller_copies,
            known_orderable_dna_sequence=dna,
            known_orderable_dna_pre_ga=smaller.get("dna_sequence_pre_ga", dna),
            recovery_original_cap3000_verified_max_copies=larger_copies,
            recovery_source_fragment_limit_bp=1800,
            recovery_policy="larger-cap-carries-proven-smaller-cap-lower-bound-v1",
            final_passed=False,
            final_status="pending_monotonic_recovery",
        )
        rows.append(payload)

    recovery = pd.DataFrame(rows)
    if args.expected_rows is not None and len(recovery) != args.expected_rows:
        raise RuntimeError(
            f"Expected {args.expected_rows} non-monotonic rows, observed {len(recovery)}"
        )
    translation_mismatches = sum(
        translate_dna(str(row.known_orderable_dna_sequence))
        != str(row.unit_sequence) * int(row.known_orderable_copies)
        for row in recovery.itertuples(index=False)
    )
    validation = {
        "source_rows": len(frame),
        "source_modules": int(frame.module_id.nunique()),
        "recovery_rows": len(recovery),
        "unique_recovery_modules": int(recovery.module_id.nunique()),
        "translation_mismatches": int(translation_mismatches),
        "all_target_3000bp": bool(recovery.fragment_limit_bp.eq(3000).all()),
        "all_known_lower_bounds_positive": bool(
            recovery.known_orderable_copies.gt(0).all()
        ),
    }
    validation["passed"] = bool(
        len(frame) == 498
        and frame.module_id.nunique() == 249
        and len(recovery) == recovery.module_id.nunique()
        and translation_mismatches == 0
        and validation["all_target_3000bp"]
        and validation["all_known_lower_bounds_positive"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    recovery.to_parquet(args.output, index=False)
    args.validation.parent.mkdir(parents=True, exist_ok=True)
    args.validation.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    if not validation["passed"]:
        raise RuntimeError(json.dumps(validation, indent=2))
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
