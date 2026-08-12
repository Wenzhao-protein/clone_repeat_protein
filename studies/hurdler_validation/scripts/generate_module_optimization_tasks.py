#!/usr/bin/env python3
"""Generate recoverable Digs tasks for a versioned HURDLER/IDT analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shlex
from pathlib import Path

import pandas as pd


def command_line(command: list[str]) -> str:
    rendered = " ".join(shlex.quote(value) for value in command)
    return f'set -o pipefail; {rendered}; rc=$?; date -Is; exit "$rc"'


def write_tasks(run: Path, rows: list[tuple[str, list[str], str]]) -> None:
    taskfiles = run / "taskfiles"
    taskfiles.mkdir(parents=True, exist_ok=True)
    with (taskfiles / "tasks.txt").open("w") as handle:
        for _, command, _ in rows:
            handle.write(command_line(command) + "\n")
    with (taskfiles / "task_index.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["task_index", "case_id", "expected_output", "command"])
        for index, (case_id, command, expected) in enumerate(rows, start=1):
            writer.writerow([index, case_id, expected, " ".join(map(shlex.quote, command))])


def write_run_manifest(run: Path, payload: dict[str, object]) -> None:
    run.mkdir(parents=True, exist_ok=True)
    (run / "run.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo", type=Path, default=Path("/home/wendai/projects/hurdler/clone_repeat_protein")
    )
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument(
        "--optimization-shards",
        type=int,
        default=0,
        help=(
            "Number of HURDLER shards; 0 creates one recoverable shard per "
            "catalog module to minimize array tail latency"
        ),
    )
    parser.add_argument("--adaptive-groups", type=int, default=32)
    parser.add_argument("--processes-per-group", type=int, default=4)
    parser.add_argument("--analysis-tag", default="periodic_v4")
    parser.add_argument("--run-start", type=int, default=42)
    parser.add_argument("--hurdler-workers", type=int, default=8)
    args = parser.parse_args()
    repo = args.repo.absolute()
    study = repo / "studies" / "hurdler_validation"
    scratch_study = Path("/net/scratch") / study.relative_to("/home")
    env = Path("/home/wendai/.conda/envs/hurdler")
    hurdler = env / "bin" / "hurdler"
    catalog = pd.read_parquet(args.catalog)
    module_count = len(catalog)
    optimization_shards = args.optimization_shards or module_count
    if optimization_shards < 1:
        raise ValueError("optimization-shards must be positive, or zero for one shard per module")
    if optimization_shards > module_count:
        raise ValueError("optimization-shards cannot exceed the number of catalog modules")
    expected_rows = module_count * 2
    index_dir = (
        scratch_study
        / "step01_reference_lookup/runs/run01_production/raw/legacy-optimized-v1"
    )
    table_dir = study / "step04_module_optimization/tables" / args.analysis_tag
    table_dir.mkdir(parents=True, exist_ok=True)

    run37 = study / "step04_module_optimization/runs" / f"run{args.run_start}_{args.analysis_tag}_hurdler"
    # HURDLER shard outputs must be visible to the merge job. Some Digs nodes
    # expose node-private /net/scratch views, so authoritative shard tables use
    # shared /home; the large lookup index remains a read-only scratch input.
    raw37 = study / "step04_module_optimization/runs" / run37.name / "raw"
    rows37: list[tuple[str, list[str], str]] = []
    for shard in range(optimization_shards):
        output = raw37 / f"shard_{shard:03d}"
        rows37.append(
            (
                f"{args.analysis_tag}_hurdler_{shard:03d}",
                [
                    str(hurdler),
                    "optimize-modules",
                    "--catalog",
                    str(args.catalog),
                    "--index-dir",
                    str(index_dir),
                    "--output-dir",
                    str(output),
                    "--fragment-limits",
                    "1800",
                    "3000",
                    "--codon-usage",
                    str(repo / "data/reference_output/codon_usage.csv"),
                    "--shard-index",
                    str(shard),
                    "--shard-count",
                    str(optimization_shards),
                    "--workers",
                    str(args.hurdler_workers),
                ],
                str(output / "optimized_constructs.parquet"),
            )
        )
    write_tasks(run37, rows37)
    module_manifest = catalog.copy().reset_index(drop=True)
    module_manifest.insert(0, "catalog_row", range(len(module_manifest)))
    module_manifest.insert(
        1,
        "shard_index",
        [index % optimization_shards for index in range(len(module_manifest))],
    )
    module_manifest.insert(2, "task_index", module_manifest.shard_index + 1)
    module_manifest["unit_sequence_sha256"] = module_manifest.unit_sequence.map(
        lambda value: hashlib.sha256(str(value).encode()).hexdigest()
    )
    module_manifest["expected_output"] = module_manifest.shard_index.map(
        lambda shard: str(
            raw37 / f"shard_{int(shard):03d}" / "optimized_constructs.parquet"
        )
    )
    manifest_columns = [
        "task_index",
        "shard_index",
        "catalog_row",
        "module_id",
        "collection",
        "family",
        "unit_length",
        "unit_sequence_sha256",
        "selected_module_index",
        "selected_module_start",
        "selected_module_end",
        "selected_module_policy",
        "expected_output",
    ]
    module_manifest[manifest_columns].to_csv(
        run37 / "taskfiles/module_shard_manifest.csv", index=False
    )
    write_run_manifest(
        run37,
        {
            "run_id": run37.name,
            "status": "ready",
            "catalog": str(args.catalog),
            "module_count": module_count,
            "shards": optimization_shards,
            "sharding_policy": (
                "one_module_per_shard"
                if optimization_shards == module_count
                else "round_robin_module_shards"
            ),
            "resources": {"partition": "cpu", "cpu": args.hurdler_workers, "mem": "12G", "time": "02:00:00"},
        },
    )

    run38 = study / "step04_module_optimization/runs" / f"run{args.run_start + 1}_{args.analysis_tag}_hurdler_merge"
    candidate_output = table_dir / "module_hurdler_candidates.parquet"
    rows38 = [
        (
            f"{args.analysis_tag}_hurdler_merge",
            [
                str(env / "bin/python"),
                str(study / "scripts/merge_optimization_shards.py"),
                "--shard-root",
                str(raw37),
                "--catalog",
                str(args.catalog),
                "--table-dir",
                str(table_dir),
                "--candidate-output",
                str(candidate_output),
                "--expected-shards",
                str(optimization_shards),
            ],
            str(table_dir / "optimization_validation.json"),
        )
    ]
    write_tasks(run38, rows38)
    write_run_manifest(
        run38,
        {"run_id": run38.name, "status": f"waiting_for_{run37.name}"},
    )

    run39 = study / "step04_module_optimization/runs" / f"run{args.run_start + 2}_{args.analysis_tag}_idt_adaptive"
    # IDT responses are small but remain on shared /home because a prior Digs
    # smoke found node-private scratch visibility for these adaptive shards.
    raw39 = run39 / "raw"
    total_shards = args.adaptive_groups * args.processes_per_group
    rows39: list[tuple[str, list[str], str]] = []
    for group in range(args.adaptive_groups):
        rows39.append(
            (
                f"{args.analysis_tag}_idt_group_{group:03d}",
                [
                    str(study / "scripts/run_with_idt_credentials.sh"),
                    str(env / "bin/python"),
                    str(study / "scripts/run_parallel_refine_group.py"),
                    "--hurdler",
                    str(hurdler),
                    "--constructs",
                    str(table_dir / "optimized_constructs.parquet"),
                    "--output-root",
                    str(raw39),
                    "--codon-usage",
                    str(repo / "data/reference_output/codon_usage.csv"),
                    "--restriction-sites",
                    str(repo / "data/reference_output/restriction_enzyme.csv"),
                    "--group-index",
                    str(group),
                    "--group-count",
                    str(args.adaptive_groups),
                    "--total-shards",
                    str(total_shards),
                    "--shards-per-group",
                    str(args.processes_per_group),
                    "--max-workers",
                    str(args.processes_per_group),
                    "--population-size",
                    "64",
                    "--short-generations",
                    "10",
                    "--generation-schedule",
                    "10",
                    "20",
                    "40",
                    "60",
                    "80",
                    "100",
                    "--seed",
                    "42",
                    "--use-idt",
                ],
                str(raw39 / "groups" / f"group_{group:03d}.json"),
            )
        )
    write_tasks(run39, rows39)
    write_run_manifest(
        run39,
        {
            "run_id": run39.name,
            "status": f"waiting_for_{run38.name}",
            "expected_rows": expected_rows,
            "groups": args.adaptive_groups,
            "processes_per_group": args.processes_per_group,
            "total_scientific_shards": total_shards,
            "resources": {"partition": "cpu", "cpu": 4, "mem": "8G", "time": "02:00:00"},
            "credentials": "runtime IDT_CREDENTIAL_FILE (contents and resolved path excluded)",
        },
    )

    run40 = study / "step04_module_optimization/runs" / f"run{args.run_start + 3}_{args.analysis_tag}_idt_merge"
    rows40 = [
        (
            f"{args.analysis_tag}_idt_merge",
            [
                str(env / "bin/python"),
                str(study / "scripts/merge_ga_shards.py"),
                "--shard-root",
                str(raw39),
                "--table-dir",
                str(table_dir),
                "--expected-shards",
                str(total_shards),
                "--expected-rows",
                str(expected_rows),
                "--require-adaptive",
                "--require-idt-orderable",
            ],
            str(table_dir / "ga_validation.json"),
        )
    ]
    write_tasks(run40, rows40)
    write_run_manifest(
        run40,
        {"run_id": run40.name, "status": f"waiting_for_{run39.name}"},
    )

    run41 = study / "step04_module_optimization/runs" / f"run{args.run_start + 4}_{args.analysis_tag}_summary"
    rows41 = [
        (
            f"{args.analysis_tag}_final_summary",
            [
                str(env / "bin/python"),
                str(study / "scripts/build_module_summary.py"),
                "--repo",
                str(repo),
                "--catalog",
                str(args.catalog),
                "--hurdler-results",
                str(table_dir / "module_hurdler_results.parquet"),
                "--optimized-constructs",
                str(table_dir / "optimized_constructs.parquet"),
                "--output-dir",
                str(table_dir),
                "--figure-dir",
                str(study / "step04_module_optimization/figures" / args.analysis_tag),
            ],
            str(table_dir / "module_final_summary_validation.json"),
        )
    ]
    write_tasks(run41, rows41)
    write_run_manifest(
        run41,
        {"run_id": run41.name, "status": f"waiting_for_{run40.name}"},
    )
    print(
        json.dumps(
            {
                "module_count": module_count,
                "expected_construct_rows": expected_rows,
                "hurdler_tasks": len(rows37),
                "adaptive_groups": len(rows39),
                "adaptive_scientific_shards": total_shards,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
