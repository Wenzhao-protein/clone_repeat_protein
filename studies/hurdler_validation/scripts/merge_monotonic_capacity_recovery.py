#!/usr/bin/env python3
"""Merge strict recovery rows and enforce cross-cap maximum monotonicity."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

import pandas as pd

from hurdler.optimization import translate_dna


KEY = ["module_id", "fragment_limit_bp"]


def sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--recovery", type=Path, required=True)
    parser.add_argument("--recovery-validation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--promote-primary", action="store_true")
    args = parser.parse_args()

    original = pd.read_parquet(args.original)
    recovery = pd.read_parquet(args.recovery)
    recovery_validation = json.loads(args.recovery_validation.read_text())
    recovery_keys = set(zip(recovery.module_id.astype(str), recovery.fragment_limit_bp))
    if len(recovery_keys) != len(recovery):
        raise RuntimeError("Recovery rows are not unique by module and fragment cap")
    kept = original.loc[
        ~pd.Series(
            list(zip(original.module_id.astype(str), original.fragment_limit_bp)),
            index=original.index,
        ).isin(recovery_keys)
    ]
    combined = pd.concat([kept, recovery], ignore_index=True).sort_values(KEY)

    pivot = combined.pivot(index="module_id", columns="fragment_limit_bp", values="verified_max_copies")
    decreasing = pivot.index[pivot[3000] < pivot[1800]].astype(str).tolist()
    translation_mismatches = 0
    idt_hash_mismatches = 0
    locked_window_mismatches = 0
    for row in combined.itertuples(index=False):
        dna = getattr(row, "dna_sequence", None)
        if not isinstance(dna, str) or not dna:
            continue
        translation_mismatches += int(
            translate_dna(dna) != str(row.unit_sequence) * int(row.verified_max_copies)
        )
        idt_hash_mismatches += int(
            str(getattr(row, "idt_scored_sequence_sha256", "")) != sha256(dna)
        )
        pre_ga = str(row.dna_sequence_pre_ga)
        for start in (int(row.site_i_position), int(row.site_ii_position)):
            locked_window_mismatches += int(
                dna[start * 3 : (start + 3) * 3]
                != pre_ga[start * 3 : (start + 3) * 3]
            )

    recovery_failed = int((~recovery.final_passed.fillna(False)).sum())
    recovery_boundary_unproven = int(
        (~recovery.adaptive_boundary_proven.fillna(False)).sum()
    )
    recovery_below_known = int(
        (
            pd.to_numeric(recovery.verified_max_copies)
            < pd.to_numeric(recovery.known_orderable_copies)
        ).sum()
    )
    validation = {
        "original_rows": len(original),
        "recovery_rows": len(recovery),
        "combined_rows": len(combined),
        "combined_modules": int(combined.module_id.nunique()),
        "duplicate_module_cap_rows": int(combined.duplicated(KEY).sum()),
        "cross_cap_decreasing_modules": decreasing,
        "recovery_failed_rows": recovery_failed,
        "recovery_boundary_unproven_rows": recovery_boundary_unproven,
        "recovery_below_known_lower_bound_rows": recovery_below_known,
        "translation_mismatches": translation_mismatches,
        "idt_hash_mismatches": idt_hash_mismatches,
        "locked_hurdler_window_mismatches": locked_window_mismatches,
        "strict_recovery_validation_passed": bool(recovery_validation.get("passed")),
    }
    validation["passed"] = bool(
        len(combined) == 498
        and combined.module_id.nunique() == 249
        and validation["duplicate_module_cap_rows"] == 0
        and not decreasing
        and recovery_failed == 0
        and recovery_boundary_unproven == 0
        and recovery_below_known == 0
        and translation_mismatches == 0
        and idt_hash_mismatches == 0
        and locked_window_mismatches == 0
        and recovery_validation.get("passed") is True
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(args.output, index=False)
    combined.to_csv(args.output.with_suffix(".csv"), index=False)
    with args.output.with_suffix(".fasta").open("w") as handle:
        for row in combined.itertuples(index=False):
            dna = getattr(row, "dna_sequence", None)
            if not isinstance(dna, str) or not dna:
                continue
            handle.write(
                f">{row.module_id}|cap={int(row.fragment_limit_bp)}|"
                f"copies={int(row.verified_max_copies)}|idt={row.idt_status}\n"
            )
            for start in range(0, len(dna), 80):
                handle.write(dna[start : start + 80] + "\n")
    args.validation.parent.mkdir(parents=True, exist_ok=True)
    args.validation.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    if not validation["passed"]:
        raise RuntimeError(json.dumps(validation, indent=2))

    if args.promote_primary:
        backup = args.original.with_name("optimized_constructs_independent_caps.parquet")
        if not backup.exists():
            shutil.copyfile(args.original, backup)
        temporary = args.original.with_suffix(".promoting.parquet")
        combined.to_parquet(temporary, index=False)
        os.replace(temporary, args.original)
        combined.to_csv(args.original.with_suffix(".csv"), index=False)
        shutil.copyfile(args.output.with_suffix(".fasta"), args.original.with_suffix(".fasta"))
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
