#!/usr/bin/env python3
"""Merge GA-refined constructs and enforce local hard constraints."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

from hurdler.ga_optimization import GA_RE_SITE_POLICY
from hurdler.optimization import translate_dna


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-root", type=Path, required=True)
    parser.add_argument("--table-dir", type=Path, required=True)
    parser.add_argument("--expected-shards", type=int, required=True)
    parser.add_argument("--expected-rows", type=int, required=True)
    parser.add_argument("--require-adaptive", action="store_true")
    parser.add_argument("--require-idt-orderable", action="store_true")
    args = parser.parse_args()
    frames = []
    missing = []
    for index in range(args.expected_shards):
        path = args.shard_root / f"shard_{index:03d}" / "optimized_constructs_ga.parquet"
        if not path.is_file():
            missing.append(str(path))
            continue
        frame = pd.read_parquet(path)
        frame["ga_source_shard"] = index
        frames.append(frame)
    if missing:
        raise FileNotFoundError("Missing GA shards:\n" + "\n".join(missing))
    result = pd.concat(frames, ignore_index=True)
    duplicate_rows = int(result.duplicated(["module_id", "fragment_limit_bp"]).sum())
    ga_re_site_policy_mismatches = int(
        result.get("ga_re_site_policy", pd.Series("", index=result.index))
        .fillna("")
        .ne(GA_RE_SITE_POLICY)
        .sum()
    )
    translated_mismatches = 0
    locked_mismatches = 0
    for row in result.itertuples(index=False):
        dna = getattr(row, "dna_sequence", None)
        if not isinstance(dna, str) or not dna:
            continue
        expected_protein = str(row.unit_sequence) * int(row.verified_max_copies)
        translated_mismatches += int(translate_dna(dna) != expected_protein)
        original = str(row.dna_sequence_pre_ga)
        for start in (int(row.site_i_position), int(row.site_ii_position)):
            locked_mismatches += int(dna[start * 3 : (start + 3) * 3] != original[start * 3 : (start + 3) * 3])
    applicable = result.dna_sequence.notna() & result.dna_sequence.astype(str).ne("")
    local_passed = result.get(
        "ga_local_constraints_passed", pd.Series(False, index=result.index)
    ).fillna(False).astype(bool)
    local_failures = int((applicable & ~local_passed).sum())
    adaptive_missing = 0
    adaptive_trace_parse_failures = 0
    adaptive_copy_bound_failures = 0
    adaptive_constraint_failures = 0
    adaptive_route_failures = 0
    adaptive_idt_gate_missing_rows = 0
    adaptive_idt_trace_unscored_local_passes = 0
    adaptive_idt_trace_semantic_failures = 0
    adaptive_idt_feedback_missing = 0
    adaptive_boundary_proven = pd.Series(False, index=result.index, dtype=bool)
    adaptive_boundary_evidence = pd.Series(
        "not_applicable_no_hurdler_solution", index=result.index, dtype=object
    )
    search_expected = pd.Series(False, index=result.index, dtype=bool)
    if args.require_adaptive:
        required_columns = {
            "pre_adaptive_verified_max_copies",
            "adaptive_search_upper_bound_copies",
            "adaptive_verified_max_copies",
            "adaptive_search_trace_json",
            "adaptive_stop_reason",
            "ga_adaptive_constraints_passed",
        }
        absent = required_columns.difference(result.columns)
        if absent:
            raise RuntimeError(f"Adaptive result columns are missing: {sorted(absent)}")
        search_expected = (
            pd.to_numeric(result.pre_adaptive_verified_max_copies, errors="coerce")
            .fillna(0)
            .gt(0)
        )
        if args.require_idt_orderable:
            orderability_gate = result.get(
                "adaptive_orderability_gate",
                pd.Series(False, index=result.index),
            ).fillna(False).astype(bool)
            adaptive_idt_gate_missing_rows = int(
                (search_expected & ~orderability_gate).sum()
            )
        trace_present = result.adaptive_search_trace_json.fillna("").astype(str).str.startswith("[")
        stop_present = result.adaptive_stop_reason.fillna("").astype(str).ne("")
        adaptive_missing = int((search_expected & ~(trace_present & stop_present)).sum())
        for value in result.loc[search_expected & trace_present, "adaptive_search_trace_json"]:
            try:
                if not isinstance(json.loads(str(value)), list):
                    adaptive_trace_parse_failures += 1
            except (TypeError, ValueError, json.JSONDecodeError):
                adaptive_trace_parse_failures += 1
        adaptive_copies = pd.to_numeric(result.adaptive_verified_max_copies, errors="coerce").fillna(-1)
        search_upper = pd.to_numeric(
            result.adaptive_search_upper_bound_copies, errors="coerce"
        ).fillna(-1)
        adaptive_copy_bound_failures = int(
            ((adaptive_copies < 0) | (adaptive_copies > search_upper)).sum()
        )
        adaptive_passed = result.ga_adaptive_constraints_passed.fillna(False).astype(bool)
        adaptive_constraint_failures = int((applicable & ~adaptive_passed).sum())
        for row in result.loc[search_expected].itertuples(index=True):
            try:
                row_index = int(row.Index)
                trace = json.loads(str(row.adaptive_search_trace_json))
                minimum = int(
                    getattr(
                        row,
                        "adaptive_search_minimum_copies",
                        max(
                            1,
                            (
                                max(int(row.site_i_position), int(row.site_ii_position))
                                + 3
                                + len(str(row.unit_sequence))
                                - 1
                            )
                            // len(str(row.unit_sequence)),
                        ),
                    )
                )
                upper = int(row.adaptive_search_upper_bound_copies)
                best = int(row.adaptive_verified_max_copies)
                short = int(row.adaptive_short_generations)
                schedule = {
                    int(value)
                    for value in str(row.adaptive_generation_schedule).split(",")
                    if value
                } | {short, 100}
                if not isinstance(trace, list) or not trace:
                    adaptive_route_failures += 1
                    continue
                if args.require_idt_orderable:
                    for item in trace:
                        # Non-selected repeated RE sites are a soft GA score
                        # term.  They must not suppress the required IDT call;
                        # only the selected-site/GC local hard gate applies.
                        local_candidate = bool(
                            item.get("ga_local_constraints_passed")
                        )
                        if local_candidate and not item.get("idt_api_called"):
                            adaptive_idt_trace_unscored_local_passes += 1
                        expected_pass = bool(
                            local_candidate
                            and item.get("idt_explicit_pass") is True
                        )
                        if bool(item.get("passed")) != expected_pass:
                            adaptive_idt_trace_semantic_failures += 1
                        if (
                            item.get("idt_api_called")
                            and item.get("idt_explicit_pass") is False
                            and int(item.get("idt_violation_count") or 0) > 0
                            and str(
                                item.get("idt_feedback_adjustments_json", "[]")
                            ) == "[]"
                        ):
                            adaptive_idt_feedback_missing += 1
                if any(
                    item.get("phase") not in {"binary_short", "linear_escalation"}
                    or (
                        item.get("phase") == "binary_short"
                        and int(item.get("generations", -1)) != short
                    )
                    for item in trace
                ):
                    adaptive_route_failures += 1
                    continue
                passed_copies = [int(item["copies"]) for item in trace if item.get("passed")]
                if max(passed_copies, default=0) != best or not (0 <= best <= upper):
                    adaptive_route_failures += 1
                    continue
                linear_by_copy: dict[int, list[dict[str, object]]] = {}
                for item in trace:
                    if item.get("phase") == "linear_escalation":
                        linear_by_copy.setdefault(int(item["copies"]), []).append(item)
                route_invalid = False
                for evaluations in linear_by_copy.values():
                    generations = [int(item["generations"]) for item in evaluations]
                    if (
                        generations != sorted(set(generations))
                        or any(value not in schedule for value in generations)
                        or any(item.get("passed") for item in evaluations[:-1])
                    ):
                        route_invalid = True
                        break
                stop = str(row.adaptive_stop_reason)
                if stop == "reached_local_upper_bound":
                    route_invalid |= best != upper
                    if not route_invalid:
                        adaptive_boundary_proven.loc[row_index] = True
                        adaptive_boundary_evidence.loc[row_index] = (
                            "mathematical_fragment_cap_reached"
                        )
                elif stop.startswith("copy_") and stop.endswith("_failed_at_100"):
                    failed_copy = int(stop.split("_")[1])
                    expected_failed = best + 1 if best else minimum
                    evaluations = linear_by_copy.get(failed_copy, [])
                    route_invalid |= bool(
                        failed_copy != expected_failed
                        or not evaluations
                        or int(evaluations[-1]["generations"]) != 100
                        or evaluations[-1].get("passed")
                    )
                    if not route_invalid:
                        adaptive_boundary_proven.loc[row_index] = True
                        adaptive_boundary_evidence.loc[row_index] = stop
                elif stop.startswith("copy_") and stop.endswith(
                    "_terminal_construction_failure"
                ):
                    # A builder exception is useful diagnostic evidence, but
                    # it does not satisfy the user's requirement that N+1 be
                    # attempted through 100 GA generations.
                    route_invalid = True
                else:
                    route_invalid = True
                adaptive_route_failures += int(route_invalid)
            except (AttributeError, KeyError, TypeError, ValueError, json.JSONDecodeError):
                adaptive_route_failures += 1
        result["adaptive_boundary_proven"] = adaptive_boundary_proven
        result["adaptive_boundary_evidence"] = adaptive_boundary_evidence
    api_called = result.get("idt_api_called", pd.Series(False, index=result.index)).fillna(False).astype(bool)
    unchanged = result.get(
        "idt_scored_sequence_unchanged", pd.Series(False, index=result.index)
    ).fillna(False).astype(bool)
    matched_response = (
        result.get("idt_result_matched_by_name", pd.Series(False, index=result.index))
        .fillna(False)
        .astype(bool)
        | result.get("idt_result_selected_by_index", pd.Series(False, index=result.index))
        .fillna(False)
        .astype(bool)
    )
    idt_api_not_called = int((applicable & ~api_called).sum())
    idt_scored_sequence_changed = int((applicable & api_called & ~unchanged).sum())
    idt_response_unmatched = int((applicable & api_called & ~matched_response).sum())
    raw_response_missing = 0
    raw_response_files = result.get("idt_raw_response_file", pd.Series("", index=result.index))
    for path in raw_response_files.loc[applicable & api_called].dropna():
        raw_response_missing += int(not Path(str(path)).is_file())
    idt_not_passed = int((applicable & ~result.idt_status.eq("passed")).sum())
    idt_rule_violation_total = int(
        pd.to_numeric(
            result.get("idt_violation_count", pd.Series(0, index=result.index)),
            errors="coerce",
        )
        .fillna(0)
        .sum()
    )
    validation = {
        "expected_shards": args.expected_shards,
        "observed_shards": len(frames),
        "expected_rows": args.expected_rows,
        "observed_rows": len(result),
        "duplicate_module_cap_rows": duplicate_rows,
        "ga_re_site_policy": GA_RE_SITE_POLICY,
        "ga_re_site_policy_mismatch_rows": ga_re_site_policy_mismatches,
        "translated_mismatches": translated_mismatches,
        "locked_hurdler_window_mismatches": locked_mismatches,
        "applicable_constructs": int(applicable.sum()),
        "ga_local_constraint_failures": local_failures,
        "adaptive_required": args.require_adaptive,
        "adaptive_search_missing_rows": adaptive_missing,
        "adaptive_trace_parse_failures": adaptive_trace_parse_failures,
        "adaptive_copy_bound_failures": adaptive_copy_bound_failures,
        "adaptive_constraint_failures": adaptive_constraint_failures,
        "adaptive_route_failures": adaptive_route_failures,
        "adaptive_idt_orderability_required": args.require_idt_orderable,
        "adaptive_idt_gate_missing_rows": adaptive_idt_gate_missing_rows,
        "adaptive_idt_trace_unscored_local_passes": adaptive_idt_trace_unscored_local_passes,
        "adaptive_idt_trace_semantic_failures": adaptive_idt_trace_semantic_failures,
        "adaptive_idt_feedback_missing": adaptive_idt_feedback_missing,
        "adaptive_boundary_proven_rows": int(adaptive_boundary_proven.sum()),
        "adaptive_boundary_unproven_rows": int(
            (search_expected & ~adaptive_boundary_proven).sum()
        ),
        "idt_api_not_called_rows": idt_api_not_called,
        "idt_scored_sequence_changed_rows": idt_scored_sequence_changed,
        "idt_response_unmatched_rows": idt_response_unmatched,
        "idt_raw_response_missing_rows": raw_response_missing,
        "idt_not_passed_rows": idt_not_passed,
        "idt_rule_violation_total": idt_rule_violation_total,
        "idt_status_counts": {
            str(status): int(count)
            for status, count in result.idt_status.value_counts(dropna=False).items()
        },
        "final_passed": int(result.final_passed.sum()),
    }
    validation["passed"] = bool(
        len(result) == args.expected_rows
        and duplicate_rows == 0
        and ga_re_site_policy_mismatches == 0
        and translated_mismatches == 0
        and locked_mismatches == 0
        and local_failures == 0
        and adaptive_missing == 0
        and adaptive_trace_parse_failures == 0
        and adaptive_copy_bound_failures == 0
        and adaptive_constraint_failures == 0
        and adaptive_route_failures == 0
        and adaptive_idt_gate_missing_rows == 0
        and adaptive_idt_trace_unscored_local_passes == 0
        and adaptive_idt_trace_semantic_failures == 0
        and adaptive_idt_feedback_missing == 0
        and int((search_expected & ~adaptive_boundary_proven).sum()) == 0
        and idt_api_not_called == 0
        and idt_scored_sequence_changed == 0
        and idt_response_unmatched == 0
        and raw_response_missing == 0
        and (not args.require_idt_orderable or idt_not_passed == 0)
    )
    if not validation["passed"]:
        raise RuntimeError(json.dumps(validation, indent=2))
    result = result.sort_values(["collection", "module_id", "fragment_limit_bp"]).reset_index(drop=True)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    primary = args.table_dir / "optimized_constructs.parquet"
    if primary.is_file() and not (args.table_dir / "optimized_constructs_pre_ga.parquet").exists():
        shutil.copy2(primary, args.table_dir / "optimized_constructs_pre_ga.parquet")
    result.to_parquet(primary, index=False)
    result.to_csv(args.table_dir / "optimized_constructs.csv", index=False)
    result.to_parquet(args.table_dir / "optimized_constructs_ga.parquet", index=False)
    maximum_mask = (
        result.dna_sequence.notna()
        & result.dna_sequence.astype(str).ne("")
        & result.get(
            "ga_local_constraints_passed", pd.Series(False, index=result.index)
        )
        .fillna(False)
        .astype(bool)
    )
    if args.require_adaptive:
        maximum_mask &= result.adaptive_boundary_proven.fillna(False).astype(bool)
    if args.require_idt_orderable:
        maximum_mask &= (
            result.get(
                "adaptive_orderable_passed",
                pd.Series(False, index=result.index),
            )
            .fillna(False)
            .astype(bool)
            & result.idt_status.eq("passed")
            & result.final_passed.fillna(False).astype(bool)
        )
    maximum = result.loc[maximum_mask].copy()
    maximum.to_parquet(args.table_dir / "maximum_passed_constructs.parquet", index=False)
    maximum.to_csv(args.table_dir / "maximum_passed_constructs.csv", index=False)
    if args.require_adaptive:
        trace_rows: list[dict[str, object]] = []
        for row in result.itertuples(index=False):
            raw_trace = getattr(row, "adaptive_search_trace_json", None)
            if not isinstance(raw_trace, str):
                continue
            for evaluation_index, evaluation in enumerate(json.loads(raw_trace)):
                trace_rows.append(
                    {
                        "module_id": row.module_id,
                        "collection": row.collection,
                        "family": row.family,
                        "in_designed_primary100": row.in_designed_primary100,
                        "fragment_limit_bp": row.fragment_limit_bp,
                        "evaluation_index": evaluation_index,
                        **evaluation,
                    }
                )
        trace_frame = pd.DataFrame(trace_rows)
        trace_frame.to_parquet(args.table_dir / "adaptive_copy_search_trace.parquet", index=False)
        trace_frame.to_csv(args.table_dir / "adaptive_copy_search_trace.csv", index=False)
        validation["adaptive_trace_rows"] = len(trace_frame)
    fasta = args.table_dir / "optimized_constructs.fasta"
    maximum_fasta = args.table_dir / "maximum_passed_constructs.fasta"
    with fasta.open("w") as handle, maximum_fasta.open("w") as maximum_handle:
        for row in result.itertuples(index=False):
            dna = getattr(row, "dna_sequence", None)
            if not isinstance(dna, str) or not dna:
                continue
            handle.write(f">{row.module_id}|cap={row.fragment_limit_bp}|copies={row.verified_max_copies}|idt={row.idt_status}\n")
            maximum_handle.write(f">{row.module_id}|cap={row.fragment_limit_bp}|copies={row.verified_max_copies}|idt={row.idt_status}\n")
            for start in range(0, len(dna), 80):
                handle.write(dna[start : start + 80] + "\n")
                maximum_handle.write(dna[start : start + 80] + "\n")
    validation["output_sha256"] = sha256(primary)
    (args.table_dir / "ga_validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
