#!/usr/bin/env python3
"""Validate and atomically promote a complete adaptive construct result table."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

import pandas as pd

from hurdler.ga_optimization import GA_RE_SITE_POLICY
from hurdler.optimization import translate_dna


KEY = ["module_id", "fragment_limit_bp"]


def dna_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-validation", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument(
        "--backup-stem", default="optimized_constructs_pre_soft_re_site_policy_v2"
    )
    args = parser.parse_args()

    frame = pd.read_parquet(args.source).sort_values(KEY).reset_index(drop=True)
    source_validation = json.loads(args.source_validation.read_text())
    duplicate_rows = int(frame.duplicated(KEY).sum())
    policy_mismatches = int(
        frame.ga_re_site_policy.fillna("").ne(GA_RE_SITE_POLICY).sum()
    )
    search_expected = (
        pd.to_numeric(frame.pre_adaptive_verified_max_copies, errors="coerce")
        .fillna(0)
        .gt(0)
    )
    boundary_unproven = int(
        (search_expected & ~frame.adaptive_boundary_proven.fillna(False)).sum()
    )
    final_passed = frame.final_passed.fillna(False).astype(bool)
    passing_missing_dna = 0
    passing_idt_failures = 0
    translation_mismatches = 0
    idt_hash_mismatches = 0
    locked_window_mismatches = 0
    for row in frame.loc[final_passed].itertuples(index=False):
        dna = getattr(row, "dna_sequence", None)
        if not isinstance(dna, str) or not dna:
            passing_missing_dna += 1
            continue
        passing_idt_failures += int(
            str(row.idt_status) != "passed" or not bool(row.idt_explicit_pass)
        )
        translation_mismatches += int(
            translate_dna(dna)
            != str(row.unit_sequence) * int(row.verified_max_copies)
        )
        idt_hash_mismatches += int(
            str(row.idt_scored_sequence_sha256) != dna_sha256(dna)
        )
        pre_ga = str(row.dna_sequence_pre_ga)
        for start in (int(row.site_i_position), int(row.site_ii_position)):
            locked_window_mismatches += int(
                dna[start * 3 : (start + 3) * 3]
                != pre_ga[start * 3 : (start + 3) * 3]
            )
    pivot = frame.pivot(
        index="module_id", columns="fragment_limit_bp", values="verified_max_copies"
    )
    decreasing = pivot.index[pivot[3000] < pivot[1800]].astype(str).tolist()
    validation = {
        "source": str(args.source.resolve()),
        "source_validation": str(args.source_validation.resolve()),
        "source_validation_passed": bool(source_validation.get("passed")),
        "rows": len(frame),
        "modules": int(frame.module_id.nunique()),
        "duplicate_module_cap_rows": duplicate_rows,
        "ga_re_site_policy": GA_RE_SITE_POLICY,
        "ga_re_site_policy_mismatch_rows": policy_mismatches,
        "adaptive_search_expected_rows": int(search_expected.sum()),
        "adaptive_boundary_unproven_rows": boundary_unproven,
        "final_passed_rows": int(final_passed.sum()),
        "passing_missing_dna_rows": passing_missing_dna,
        "passing_idt_failures": passing_idt_failures,
        "translation_mismatches": translation_mismatches,
        "idt_hash_mismatches": idt_hash_mismatches,
        "locked_hurdler_window_mismatches": locked_window_mismatches,
        "cross_cap_decreasing_modules": decreasing,
    }
    validation["passed"] = bool(
        source_validation.get("passed") is True
        and len(frame) == 498
        and frame.module_id.nunique() == 249
        and duplicate_rows == 0
        and policy_mismatches == 0
        and boundary_unproven == 0
        and passing_missing_dna == 0
        and passing_idt_failures == 0
        and translation_mismatches == 0
        and idt_hash_mismatches == 0
        and locked_window_mismatches == 0
        and not decreasing
    )
    args.validation.parent.mkdir(parents=True, exist_ok=True)
    args.validation.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    if not validation["passed"]:
        raise RuntimeError(json.dumps(validation, indent=2))

    args.target.parent.mkdir(parents=True, exist_ok=True)
    backups: dict[str, str] = {}
    for suffix in (".parquet", ".csv", ".fasta"):
        current = args.target.with_suffix(suffix)
        backup = args.target.with_name(args.backup_stem).with_suffix(suffix)
        replacement = args.source.with_suffix(suffix)
        if current.exists() and not backup.exists():
            shutil.copyfile(current, backup)
        if replacement.exists():
            temporary = current.with_suffix(current.suffix + ".promoting")
            shutil.copyfile(replacement, temporary)
            os.replace(temporary, current)
        backups[suffix] = str(backup)
    maximum = frame.loc[final_passed].copy()
    maximum_path = args.target.parent / "maximum_passed_constructs.parquet"
    maximum.to_parquet(maximum_path, index=False)
    maximum.to_csv(maximum_path.with_suffix(".csv"), index=False)
    with maximum_path.with_suffix(".fasta").open("w") as handle:
        for row in maximum.itertuples(index=False):
            dna = str(row.dna_sequence)
            handle.write(
                f">{row.module_id}|cap={int(row.fragment_limit_bp)}|"
                f"copies={int(row.verified_max_copies)}|idt={row.idt_status}\n"
            )
            for start in range(0, len(dna), 80):
                handle.write(dna[start : start + 80] + "\n")
    validation["promoted_target"] = str(args.target.resolve())
    validation["backups"] = backups
    validation["maximum_passed_construct_rows"] = len(maximum)
    validation["maximum_passed_constructs"] = str(maximum_path.resolve())
    args.validation.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
