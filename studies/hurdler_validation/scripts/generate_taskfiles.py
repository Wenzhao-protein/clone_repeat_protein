#!/usr/bin/env python3
"""Generate inspectable, absolute Digs task files for this study."""

from __future__ import annotations

import argparse
import csv
import shlex
from pathlib import Path

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def command_line(command: list[str]) -> str:
    rendered = " ".join(shlex.quote(item) for item in command)
    return f"set -o pipefail; {rendered}; rc=$?; date -Is; exit \"$rc\""


def write_tasks(path: Path, rows: list[tuple[str, list[str]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for _case_id, command in rows:
            handle.write(command_line(command) + "\n")
    with path.with_name("task_index.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["task_index", "case_id", "command"])
        for index, (case_id, command) in enumerate(rows, start=1):
            writer.writerow([index, case_id, " ".join(shlex.quote(item) for item in command)])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path("/home/wendai/projects/hurdler/clone_repeat_protein"))
    parser.add_argument("--env", type=Path, default=Path("/home/wendai/.conda/envs/hurdler"))
    args = parser.parse_args()
    repo = args.repo.absolute()
    study = repo / "studies" / "hurdler_validation"
    scratch = Path("/net/scratch") / study.relative_to("/home")
    hurdler = args.env / "bin" / "hurdler"
    smoke_index_dir = scratch / "step01_reference_lookup" / "runs" / "run01_smoke" / "raw" / "legacy-optimized-v1"
    index_dir = scratch / "step01_reference_lookup" / "runs" / "run01_production" / "raw" / "legacy-optimized-v1"

    smoke_run = study / "step01_reference_lookup" / "runs" / "run01_smoke" / "taskfiles"
    smoke_rows = [
        (
            "pattern_index_smoke",
            [str(hurdler), "lookup", "build", "--rules", "legacy-optimized-v1", "--input-dir", str(repo / "output"), "--orthogonality", str(repo / "data/reference_output/orthogonality.csv"), "--output-dir", str(smoke_index_dir)],
        )
    ]
    write_tasks(smoke_run / "tasks.txt", smoke_rows)

    lookup_run = study / "step01_reference_lookup" / "runs" / "run01_production" / "taskfiles"
    lookup_rows = [
        (
            "reference_manifest",
            [str(hurdler), "reference", "build", "--reference-dir", str(repo / "data/reference_output"), "--output", str(study / "step01_reference_lookup/tables/reference_manifest.json")],
        ),
        (
            "pattern_index",
            [str(hurdler), "lookup", "build", "--rules", "legacy-optimized-v1", "--input-dir", str(repo / "output"), "--orthogonality", str(repo / "data/reference_output/orthogonality.csv"), "--output-dir", str(index_dir)],
        ),
    ]
    write_tasks(lookup_run / "tasks.txt", lookup_rows)

    short_run = study / "step02_success_landscape" / "runs" / "run01_production" / "taskfiles"
    short_raw = scratch / "step02_success_landscape" / "runs" / "run01_production" / "raw" / "short_shards"
    short_rows: list[tuple[str, list[str]]] = []
    for length in range(1, 5):
        short_rows.append((f"k{length}_all", [str(hurdler), "screen-short", "--index-dir", str(index_dir), "--output-dir", str(short_raw), "--length", str(length)]))
    for first in AMINO_ACIDS:
        for second in AMINO_ACIDS:
            prefix = first + second
            short_rows.append((f"k5_{prefix}", [str(hurdler), "screen-short", "--index-dir", str(index_dir), "--output-dir", str(short_raw), "--length", "5", "--prefix", prefix]))
    write_tasks(short_run / "tasks.txt", short_rows)

    finalize_run = study / "step02_success_landscape" / "runs" / "run01_finalize" / "taskfiles"
    finalize_rows = [
        (
            "finalize_short_motifs",
            [str(hurdler), "screen-short", "--index-dir", str(index_dir), "--output-dir", str(short_raw), "--length", "1", "--finalize"],
        )
    ]
    write_tasks(finalize_run / "tasks.txt", finalize_rows)

    rate_run = study / "step02_success_landscape" / "runs" / "run02_monte_carlo" / "taskfiles"
    rate_output = study / "step02_success_landscape" / "tables" / "success_rate_7_60.csv"
    rate_rows = [
        (
            "independent_6aa",
            [str(hurdler), "success-rate", "--index-dir", str(index_dir), "--output", str(study / "step02_success_landscape" / "tables" / "success_rate_6.csv"), "--min-length", "6", "--max-length", "6", "--tests", "1000", "--seed", "420006"],
        ),
        (
            "legacy_7_60",
            [str(hurdler), "success-rate", "--index-dir", str(index_dir), "--output", str(rate_output), "--min-length", "7", "--max-length", "60", "--tests", "1000", "--seed", "42"],
        )
    ]
    write_tasks(rate_run / "tasks.txt", rate_rows)

    optimization_run = study / "step04_module_optimization" / "runs" / "run03_production" / "taskfiles"
    optimization_raw = scratch / "step04_module_optimization" / "runs" / "run03_production" / "raw"
    catalog = study / "step03_module_corpus" / "tables" / "module_catalog.parquet"
    optimization_rows: list[tuple[str, list[str]]] = []
    optimization_shards = 32
    for shard_index in range(optimization_shards):
        optimization_rows.append(
            (
                f"module_optimization_{shard_index:03d}",
                [
                    str(hurdler),
                    "optimize-modules",
                    "--catalog",
                    str(catalog),
                    "--index-dir",
                    str(index_dir),
                    "--output-dir",
                    str(optimization_raw / f"shard_{shard_index:03d}"),
                    "--fragment-limits",
                    "1800",
                    "3000",
                    "--codon-usage",
                    str(repo / "data" / "reference_output" / "codon_usage.csv"),
                    "--shard-index",
                    str(shard_index),
                    "--shard-count",
                    str(optimization_shards),
                ],
            )
        )
    write_tasks(optimization_run / "tasks.txt", optimization_rows)

    legacy_run = study / "step05_reproducibility" / "runs" / "run03_legacy_notebooks" / "taskfiles"
    legacy_script = study / "scripts" / "execute_legacy_case.py"
    legacy_scratch = scratch / "step05_reproducibility" / "runs" / "run03_legacy_notebooks" / "raw"
    legacy_artifacts = study / "step05_reproducibility" / "legacy_notebooks"
    legacy_cases = [
        ("utils_get_codon_usage", "notebooks/utils/get_codon_usage.ipynb", "utils"),
        ("utils_methylation_check", "notebooks/utils/methylation_check.ipynb", "utils"),
        ("utils_neb_buffer", "notebooks/utils/neb_buffer_activity_check.ipynb", "utils"),
        ("utils_plasmid_check", "notebooks/utils/plasmid_check.ipynb", "utils"),
        ("utils_re_pair_fidelity", "notebooks/utils/re_pair_fidelity.ipynb", "utils"),
        ("utils_get_re_sites", "notebooks/utils/get_re_sites.ipynb", "utils"),
        ("enzyme_selection", "notebooks/enzyme_selection/enzyme_selection_analysis.ipynb", "root"),
        ("enzyme_candidates", "notebooks/enzyme_selection/inspect_site_candidates.ipynb", "root"),
        ("enzyme_3mer", "notebooks/enzyme_selection/re_3mer_analysis.ipynb", "root"),
        ("enzyme_plasmid", "notebooks/enzyme_selection/re_plasmid_compatibility.ipynb", "root"),
        ("hurdler_site_combinations", "notebooks/hurdler/hurdler_site_combination_analysis.ipynb", "root"),
        ("codon_benchmark", "notebooks/codon_optimization/codon_opt_benchmark.ipynb", "codon"),
        ("codon_reverse_translate", "notebooks/codon_optimization/reverse_translate.ipynb", "codon"),
        ("sec_analysis", "notebooks/sec/result_analysis_total.ipynb", "sec"),
        ("archive_enzyme_candidates", "archive/notebooks/inspect_site_candidates_executed.ipynb", "root"),
        ("archive_get_re_site", "archive/get_re_dict/get_re_site.ipynb", "archive_get_re"),
        ("archive_get_re_pair", "archive/get_re_dict/get_re_site_pair.ipynb", "archive_get_re"),
        ("archive_sankey", "archive/get_re_dict/sankey_summary.ipynb", "archive_get_re"),
    ]
    legacy_rows: list[tuple[str, list[str]]] = []
    for case_id, source, cwd_kind in legacy_cases:
        legacy_rows.append(
            (
                case_id,
                [
                    str(args.env / "bin" / "python"),
                    str(legacy_script),
                    "--case-id",
                    case_id,
                    "--source",
                    str(repo / source),
                    "--repo",
                    str(repo),
                    "--scratch-root",
                    str(legacy_scratch),
                    "--artifact-root",
                    str(legacy_artifacts),
                    "--cwd-kind",
                    cwd_kind,
                ],
            )
        )
    write_tasks(legacy_run / "tasks.txt", legacy_rows)
    compatibility_case_ids = {
        "utils_neb_buffer",
    }
    compatibility_scratch = (
        scratch
        / "step05_reproducibility"
        / "runs"
        / "run35_legacy_compatibility_final"
        / "raw"
    )
    compatibility_rows = []
    for case_id, command in legacy_rows:
        if case_id not in compatibility_case_ids:
            continue
        compatibility_rows.append(
            (
                case_id,
                [
                    str(compatibility_scratch) if item == str(legacy_scratch) else item
                    for item in command
                ],
            )
        )
    write_tasks(
        study
        / "step05_reproducibility"
        / "runs"
        / "run35_legacy_compatibility_final"
        / "taskfiles"
        / "tasks.txt",
        compatibility_rows,
    )

    ga_run = study / "step04_module_optimization" / "runs" / "run05_ga_refinement" / "taskfiles"
    ga_raw = scratch / "step04_module_optimization" / "runs" / "run05_ga_refinement" / "raw"
    merged_constructs = study / "step04_module_optimization" / "tables" / "optimized_constructs.parquet"
    ga_rows: list[tuple[str, list[str]]] = []
    ga_shards = 32
    idt_wrapper = study / "scripts" / "run_with_idt_credentials.sh"
    for shard_index in range(ga_shards):
        ga_rows.append(
            (
                f"ga_refinement_{shard_index:03d}",
                [
                    str(idt_wrapper),
                    str(hurdler),
                    "refine-ga",
                    "--constructs",
                    str(merged_constructs),
                    "--output-dir",
                    str(ga_raw / f"shard_{shard_index:03d}"),
                    "--codon-usage",
                    str(repo / "data" / "reference_output" / "codon_usage.csv"),
                    "--restriction-sites",
                    str(repo / "data" / "reference_output" / "restriction_enzyme.csv"),
                    "--shard-index",
                    str(shard_index),
                    "--shard-count",
                    str(ga_shards),
                    "--population-size",
                    "64",
                    "--generations",
                    "100",
                    "--seed",
                    "42",
                    "--use-idt",
                ],
            )
        )
    write_tasks(ga_run / "tasks.txt", ga_rows)

    adaptive_run = (
        study
        / "step04_module_optimization"
        / "runs"
        / "run07e_idt_feedback_adaptive_copy_refinement"
        / "taskfiles"
    )
    adaptive_raw = (
        study
        / "step04_module_optimization"
        / "runs"
        / "run07e_idt_feedback_adaptive_copy_refinement"
        / "raw"
    )
    deterministic_constructs = (
        study
        / "step04_module_optimization"
        / "tables"
        / "optimized_constructs_pre_ga.parquet"
    )
    adaptive_group_count = 32
    adaptive_processes = 4
    adaptive_total_shards = adaptive_group_count * adaptive_processes
    adaptive_rows: list[tuple[str, list[str]]] = []
    for group_index in range(adaptive_group_count):
        adaptive_rows.append(
            (
                f"adaptive_copy_refinement_group_{group_index:03d}",
                [
                    str(idt_wrapper),
                    str(args.env / "bin" / "python"),
                    str(study / "scripts" / "run_parallel_refine_group.py"),
                    "--hurdler",
                    str(hurdler),
                    "--constructs",
                    str(deterministic_constructs),
                    "--output-root",
                    str(adaptive_raw),
                    "--codon-usage",
                    str(repo / "data" / "reference_output" / "codon_usage.csv"),
                    "--restriction-sites",
                    str(repo / "data" / "reference_output" / "restriction_enzyme.csv"),
                    "--group-index",
                    str(group_index),
                    "--group-count",
                    str(adaptive_group_count),
                    "--total-shards",
                    str(adaptive_total_shards),
                    "--shards-per-group",
                    str(adaptive_processes),
                    "--max-workers",
                    str(adaptive_processes),
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
            )
        )
    write_tasks(adaptive_run / "tasks.txt", adaptive_rows)

    adaptive_merge_run = (
        study
        / "step04_module_optimization"
        / "runs"
        / "run08e_idt_feedback_adaptive_merge"
        / "taskfiles"
    )
    adaptive_merge_rows = [
        (
            "merge_mathbound_adaptive_copy_results",
            [
                str(args.env / "bin" / "python"),
                str(study / "scripts" / "merge_ga_shards.py"),
                "--shard-root",
                str(adaptive_raw),
                "--table-dir",
                str(study / "step04_module_optimization" / "tables"),
                "--expected-shards",
                str(adaptive_total_shards),
                "--expected-rows",
                "498",
                "--require-adaptive",
                "--require-idt-orderable",
            ],
        )
    ]
    write_tasks(adaptive_merge_run / "tasks.txt", adaptive_merge_rows)

    validation_run = study / "step05_reproducibility" / "runs" / "run36_code_validation_frozen_final" / "taskfiles"
    validation_rows = [
        (
            "code_validation_final",
            [
                str(args.env / "bin" / "python"),
                str(study / "scripts" / "validate_codebase.py"),
                "--repo",
                str(repo),
                "--output",
                str(study / "step05_reproducibility" / "tables" / "code_validation.json"),
            ],
        )
    ]
    write_tasks(validation_run / "tasks.txt", validation_rows)

    final_notebook_run = study / "step05_reproducibility" / "runs" / "run25_notebooks_idt_final" / "taskfiles"
    notebook_output = study / "step05_reproducibility" / "notebooks"
    html_output = study / "step05_reproducibility" / "html"
    notebook_cases = [
        ("reference_manifest", "notebooks/reference/01_reference_manifest.ipynb", "01_reference_manifest"),
        ("lookup_qc", "notebooks/reference/02_lookup_qc.ipynb", "02_lookup_qc"),
        ("hurdler_query", "notebooks/tasks/01_hurdler_query.ipynb", "03_hurdler_query"),
        ("success_rate_1_60", "notebooks/tasks/02_success_rate_1_60.ipynb", "04_success_rate_1_60"),
        ("repeat_module_benchmark", "notebooks/tasks/03_repeat_module_benchmark.ipynb", "05_repeat_module_benchmark"),
    ]
    final_notebook_rows = [
        (
            case_id,
            [
                str(args.env / "bin" / "python"),
                str(study / "scripts" / "execute_notebook.py"),
                str(repo / source),
                str(notebook_output / f"{stem}_executed.ipynb"),
                str(html_output / f"{stem}.html"),
                "--cwd",
                str(repo),
            ],
        )
        for case_id, source, stem in notebook_cases
    ]
    write_tasks(final_notebook_run / "tasks.txt", final_notebook_rows)

    figure_run = study / "step05_reproducibility" / "runs" / "run38_figure_report_frozen" / "taskfiles"
    figure_rows = [
        (
            "figure_report",
            [
                str(args.env / "bin" / "python"),
                str(study / "scripts" / "build_figure_report.py"),
                "--repo",
                str(repo),
                "--output-dir",
                str(study / "step05_reproducibility" / "figures"),
            ],
        )
    ]
    write_tasks(figure_run / "tasks.txt", figure_rows)

    report_run = study / "step05_reproducibility" / "runs" / "run27_reproducibility_report_idt_final" / "taskfiles"
    report_command = [
        str(args.env / "bin" / "python"),
        str(study / "scripts" / "build_reproducibility_report.py"),
        "--repo",
        str(repo),
        "--output-dir",
        str(study / "step05_reproducibility" / "tables"),
    ]
    write_tasks(report_run / "tasks.txt", [("reproducibility_report", report_command)])

    status_run = study / "step05_reproducibility" / "runs" / "run28_status_notebook_idt_final" / "taskfiles"
    status_command = [
        str(args.env / "bin" / "python"),
        str(study / "scripts" / "execute_notebook.py"),
        str(repo / "notebooks" / "tasks" / "04_reproducibility_status.ipynb"),
        str(notebook_output / "06_reproducibility_status_executed.ipynb"),
        str(html_output / "06_reproducibility_status.html"),
        "--cwd",
        str(repo),
    ]
    write_tasks(status_run / "tasks.txt", [("reproducibility_status", status_command)])
    write_tasks(
        study / "step05_reproducibility" / "runs" / "run37_reproducibility_frozen" / "taskfiles" / "tasks.txt",
        [("reproducibility_report_refresh", report_command)],
    )
    write_tasks(
        study / "step05_reproducibility" / "runs" / "run30_status_notebook_frozen" / "taskfiles" / "tasks.txt",
        [("reproducibility_status_final", status_command)],
    )
    print(
        f"smoke_tasks={len(smoke_rows)} lookup_tasks={len(lookup_rows)} "
        f"short_tasks={len(short_rows)} finalize_tasks={len(finalize_rows)} "
        f"rate_tasks={len(rate_rows)} optimization_tasks={len(optimization_rows)} "
        f"legacy_notebook_tasks={len(legacy_rows)} ga_tasks={len(ga_rows)} "
        f"adaptive_tasks={len(adaptive_rows)} adaptive_merge_tasks={len(adaptive_merge_rows)} "
        f"validation_tasks={len(validation_rows)} final_notebook_tasks={len(final_notebook_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
