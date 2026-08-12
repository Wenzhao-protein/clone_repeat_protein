#!/usr/bin/env python3
"""Create recoverable Digs taskfiles for the two module experiments."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from hurdler.module_experiments import (
    STAGE2_INPUT_COLUMNS,
    prepare_adaptive_copy_frame,
)
from hurdler.optimization import load_codon_weights


HURDLER = Path("/home/wendai/.conda/envs/hurdler/bin/hurdler")
DEFAULT_INDEX = Path(
    "/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/"
    "hurdler_validation/step01_reference_lookup/runs/run01_production/raw/"
    "legacy-optimized-v1"
)
DEFAULT_CODONS = Path(
    "/home/wendai/projects/hurdler/clone_repeat_protein/"
    "data/reference_output/codon_usage.csv"
)
DEFAULT_RE_SITES = Path(
    "/home/wendai/projects/hurdler/clone_repeat_protein/"
    "data/reference_output/restriction_enzyme.csv"
)


def _write_index(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def create_stage1(args: argparse.Namespace) -> None:
    catalog = args.catalog.absolute()
    run_dir = args.run_dir.absolute()
    raw_dir = args.scratch_run_dir.absolute() / "raw"
    task_dir = run_dir / "taskfiles"
    task_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_parquet(catalog)
    shard_count = min(max(1, args.shards), len(frame))
    tasks: list[str] = []
    index_rows: list[dict[str, object]] = []
    summaries: list[Path] = []
    candidates: list[Path] = []
    for shard in range(shard_count):
        output_dir = raw_dir / f"shard_{shard:05d}"
        suffix = f"shard-{shard:05d}-of-{shard_count:05d}"
        summary = output_dir / f"module_compatibility_{suffix}.parquet"
        candidate = output_dir / f"module_compatibility_candidates_{suffix}.parquet"
        tasks.append(
            " ".join(
                [
                    str(HURDLER),
                    "module-compatibility",
                    "--catalog",
                    str(catalog),
                    "--index-dir",
                    str(args.index_dir.absolute()),
                    "--output-dir",
                    str(output_dir),
                    "--shard-index",
                    str(shard),
                    "--shard-count",
                    str(shard_count),
                ]
            )
        )
        summaries.append(summary)
        candidates.append(candidate)
        index_rows.append(
            {
                "task_id": shard + 1,
                "shard_index": shard,
                "shard_count": shard_count,
                "input_rows": len(frame.iloc[shard::shard_count]),
                "summary_output": summary,
                "candidate_output": candidate,
            }
        )
    (task_dir / "tasks.txt").write_text("\n".join(tasks) + "\n")
    (task_dir / "smoke_tasks.txt").write_text(tasks[0] + "\n")
    _write_index(task_dir / "task_index.csv", index_rows)
    finalize = " ".join(
        [
            str(HURDLER),
            "module-compatibility",
            "--output-dir",
            str(args.final_output_dir.absolute()),
            "--finalize-summaries",
            *(str(path) for path in summaries),
            "--finalize-candidates",
            *(str(path) for path in candidates),
        ]
    )
    (task_dir / "finalize_task.txt").write_text(finalize + "\n")


def create_stage2(args: argparse.Namespace) -> None:
    compatibility = args.compatibility.absolute()
    run_dir = args.run_dir.absolute()
    raw_dir = args.scratch_run_dir.absolute() / "raw"
    task_dir = run_dir / "taskfiles"
    task_dir.mkdir(parents=True, exist_ok=True)
    available_columns = set(pq.read_schema(compatibility).names)
    read_columns = [
        column for column in STAGE2_INPUT_COLUMNS if column in available_columns
    ]
    source = pd.read_parquet(compatibility, columns=read_columns)
    source = source.loc[source.hurdler_compatible.astype(bool)].sort_values(
        ["collection", "module_id"], kind="mergesort"
    ).reset_index(drop=True)
    if args.limit_modules is not None:
        if args.limit_modules < 1:
            raise ValueError("--limit-modules must be positive")
        source = source.head(args.limit_modules).copy()
    codon_weights = load_codon_weights(args.codon_usage.absolute())
    if source.empty:
        raise ValueError("Stage 1 has no compatible modules")
    if args.modules_per_task < 1:
        raise ValueError("--modules-per-task must be positive")
    groups = [
        source.iloc[start : start + args.modules_per_task].copy()
        for start in range(0, len(source), args.modules_per_task)
    ]
    shard_count = len(groups)
    tasks: list[str] = []
    index_rows: list[dict[str, object]] = []
    results: list[Path] = []
    audits: list[Path] = []
    prepared_index_rows: list[dict[str, object]] = []
    input_dir = raw_dir / "module_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    for shard, module_source in enumerate(groups):
        output_dir = raw_dir / f"shard_{shard:05d}"
        module_input = input_dir / f"module_{shard:05d}.parquet"
        rows = prepare_adaptive_copy_frame(
            module_source,
            codon_weights=codon_weights,
            fragment_limits=(1800, 3000),
            external_deduction_bp=args.external_deduction_bp,
        ).sort_values("fragment_limit_bp", kind="mergesort")
        rows.to_parquet(module_input, index=False)
        prepared_index_rows.extend(
            {
                "collection": row.collection,
                "module_id": row.module_id,
                "unit_length": int(row.unit_length),
                "fragment_limit_bp": int(row.fragment_limit_bp),
                "mathematical_max_copies": int(row.mathematical_max_copies),
                "stage2_preparation_status": row.stage2_preparation_status,
                "module_input": str(module_input),
            }
            for row in rows.itertuples(index=False)
        )
        result = output_dir / "optimized_constructs_ga.parquet"
        audit = output_dir / "idt_optimization_responses.jsonl"
        tasks.append(
            " ".join(
                [
                    str(HURDLER),
                    "adaptive-copy-search",
                    "--constructs",
                    str(module_input),
                    "--output-dir",
                    str(output_dir),
                    "--codon-usage",
                    str(args.codon_usage.absolute()),
                    "--restriction-sites",
                    str(args.restriction_sites.absolute()),
                    "--shard-index",
                    "0",
                    "--shard-count",
                    "1",
                    "--population-size",
                    str(args.population_size),
                    "--seed",
                    "42",
                    "--short-generations",
                    "10",
                    "--generation-schedule",
                    "10",
                    "20",
                    "40",
                    "60",
                    "80",
                    "100",
                    "--idt-policy",
                    "idt-rule-score-sum-lt10-v1",
                ]
            )
        )
        results.append(result)
        audits.append(audit)
        index_rows.append(
            {
                "task_id": shard + 1,
                "shard_index": shard,
                "shard_count": shard_count,
                "module_count": int(len(module_source)),
                "module_ids": ",".join(module_source.module_id.astype(str)),
                "collections": ",".join(module_source.collection.astype(str)),
                "unit_lengths": ",".join(module_source.unit_length.astype(int).astype(str)),
                "fragment_limits_bp": "1800,3000",
                "module_input": module_input,
                "result_output": result,
                "idt_audit": audit,
            }
        )
    (task_dir / "tasks.txt").write_text("\n".join(tasks) + "\n")
    (task_dir / "smoke_tasks.txt").write_text(tasks[0] + "\n")
    _write_index(task_dir / "task_index.csv", index_rows)
    pd.DataFrame(prepared_index_rows).to_parquet(
        raw_dir / "adaptive_copy_input_index.parquet", index=False
    )
    pd.DataFrame(prepared_index_rows).to_csv(
        raw_dir / "adaptive_copy_input_index.csv", index=False
    )
    result_list = task_dir / "finalize_result_paths.txt"
    audit_list = task_dir / "finalize_idt_audit_paths.txt"
    result_list.write_text("\n".join(str(path) for path in results) + "\n")
    audit_list.write_text("\n".join(str(path) for path in audits) + "\n")
    finalize = " ".join(
        [
            str(HURDLER),
            "adaptive-copy-search",
            "--compatibility",
            str(compatibility),
            "--output-dir",
            str(args.final_output_dir.absolute()),
            "--finalize-result-list",
            str(result_list),
            "--idt-audit-list",
            str(audit_list),
        ]
    )
    (task_dir / "finalize_task.txt").write_text(finalize + "\n")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    subparsers = root.add_subparsers(dest="stage", required=True)
    stage1 = subparsers.add_parser("stage1")
    stage1.add_argument("--catalog", type=Path, required=True)
    stage1.add_argument("--index-dir", type=Path, default=DEFAULT_INDEX)
    stage1.add_argument("--run-dir", type=Path, required=True)
    stage1.add_argument("--scratch-run-dir", type=Path, required=True)
    stage1.add_argument("--final-output-dir", type=Path, required=True)
    stage1.add_argument("--shards", type=int, default=128)
    stage1.set_defaults(function=create_stage1)

    stage2 = subparsers.add_parser("stage2")
    stage2.add_argument("--compatibility", type=Path, required=True)
    stage2.add_argument("--run-dir", type=Path, required=True)
    stage2.add_argument("--scratch-run-dir", type=Path, required=True)
    stage2.add_argument("--final-output-dir", type=Path, required=True)
    stage2.add_argument("--codon-usage", type=Path, default=DEFAULT_CODONS)
    stage2.add_argument("--restriction-sites", type=Path, default=DEFAULT_RE_SITES)
    stage2.add_argument("--external-deduction-bp", type=int, default=0)
    stage2.add_argument("--population-size", type=int, default=16)
    stage2.add_argument(
        "--modules-per-task",
        type=int,
        default=1,
        help="Group this many independent modules in each recoverable CPU task",
    )
    stage2.add_argument(
        "--limit-modules",
        type=int,
        help="Generate only the first N stable module tasks (benchmark/smoke only)",
    )
    stage2.set_defaults(function=create_stage2)
    return root


def main() -> int:
    args = parser().parse_args()
    args.function(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
