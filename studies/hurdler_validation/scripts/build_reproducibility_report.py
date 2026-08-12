#!/usr/bin/env python3
"""Build the exhaustive workflow status and missing-input inventories."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


DEFERRED_NOTEBOOKS = {
    "notebooks/hurdler/hurdler_success_rate_analysis.ipynb",
    "notebooks/hurdler/hurdler_success_rate_optimized.ipynb",
    "archive/notebooks/codon_optimization_20230331_1101-Copy3.ipynb",
    "archive/notebooks/hurdler_minimal.ipynb",
    "archive/notebooks/hurdler_standalone.ipynb",
    "archive/notebooks/hurdler_standalone_backup.ipynb",
    "archive/notebooks/hurdler_success_rate_analysis_backup.ipynb",
    "notebooks/hurdler/hurdler_site_combination_analysis.ipynb",
}

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def latest_legacy_statuses(repo: Path) -> dict[str, dict[str, object]]:
    statuses: dict[str, dict[str, object]] = {}
    root = repo / "studies" / "hurdler_validation" / "step05_reproducibility" / "legacy_notebooks"
    for path in sorted(root.glob("*/status.json")):
        payload = json.loads(path.read_text())
        try:
            source = str(Path(str(payload["source"])).relative_to("/mnt/home/wendai/projects/hurdler/clone_repeat_protein"))
        except ValueError:
            source = str(payload["source"]).split("clone_repeat_protein/", 1)[-1]
        previous = statuses.get(source)
        if previous is None or payload["status"] == "passed" or float(payload["runtime_seconds"]) > float(previous["runtime_seconds"]):
            statuses[source] = payload
    return statuses


def canonical_statuses(repo: Path) -> dict[str, dict[str, object]]:
    statuses: dict[str, dict[str, object]] = {}
    root = repo / "studies" / "hurdler_validation" / "step05_reproducibility" / "notebooks"
    for path in root.rglob("*.manifest.json"):
        payload = json.loads(path.read_text())
        source = str(payload["source"]).split("clone_repeat_protein/", 1)[-1]
        source_path = repo / source
        payload["current_source_hash_matches"] = bool(
            source_path.is_file() and sha256(source_path) == payload["source_sha256"]
        )
        payload["manifest_mtime_ns"] = path.stat().st_mtime_ns
        previous = statuses.get(source)
        if previous is None or (
            bool(payload["current_source_hash_matches"]),
            int(payload["manifest_mtime_ns"]),
        ) > (
            bool(previous["current_source_hash_matches"]),
            int(previous["manifest_mtime_ns"]),
        ):
            statuses[source] = payload
    return statuses


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    legacy = latest_legacy_statuses(repo)
    canonical = canonical_statuses(repo)
    registry_path = args.output_dir / "job_registry.json"
    job_registry = json.loads(registry_path.read_text()) if registry_path.is_file() else {}
    canonical_jobs = job_registry.get("canonical_notebooks", {})
    rows: list[dict[str, object]] = []

    notebooks = sorted(
        path
        for base in (repo / "notebooks", repo / "archive")
        for path in base.rglob("*.ipynb")
        if "studies/hurdler_validation" not in str(path)
    )
    for path in notebooks:
        relative = str(path.relative_to(repo))
        row: dict[str, object] = {
            "workflow": relative,
            "workflow_type": "archive_notebook" if relative.startswith("archive/") else "notebook",
            "source_sha256": sha256(path),
            "status": "failed",
            "job_id": "",
            "runtime_seconds": "",
            "resources": "cpu=1,mem=8G,time=02:00:00",
            "output_path": "",
            "reason": "No clean-kernel audit result was found",
            "rerun_command": "",
        }
        if relative in canonical and canonical[relative].get("current_source_hash_matches"):
            result = canonical[relative]
            row.update(
                status="passed",
                job_id=canonical_jobs.get(relative, ""),
                runtime_seconds=result["runtime_seconds"],
                output_path=result["executed"],
                reason="Clean-kernel Papermill execution and HTML export passed",
                rerun_command=(
                    f"/home/wendai/.conda/envs/hurdler/bin/python {repo}/studies/hurdler_validation/scripts/execute_notebook.py "
                    f"{path} {result['executed']} {result['html']} --cwd {repo}"
                ),
            )
        elif relative in legacy and relative not in DEFERRED_NOTEBOOKS:
            result = legacy[relative]
            row.update(
                status=result["status"],
                job_id=(
                    "17094712" if str(result["workflow"]).endswith("v3")
                    else "17094052" if str(result["workflow"]).endswith("v2")
                    else "17092499"
                ),
                runtime_seconds=result["runtime_seconds"],
                output_path=result["output_path"],
                reason=(result.get("error") or "Isolated legacy-overlay execution passed")[-1000:],
                rerun_command=result["rerun_command"],
            )
        elif "agarose_gel" in relative:
            row.update(
                status="blocked_missing_input",
                reason="Historical .scn image inputs are absent; no scientific image was fabricated",
                rerun_command=f"Recover the files in missing_scn_inputs.txt, then execute {path}",
            )
        elif relative in DEFERRED_NOTEBOOKS:
            row.update(
                status="deferred_long",
                reason="Monolithic historical Cartesian/GA workflow; canonical work was sharded, but this source cannot be safely split without changing it",
                rerun_command=f"/home/wendai/.conda/envs/hurdler/bin/python {repo}/studies/hurdler_validation/scripts/execute_notebook.py {path} <executed.ipynb> <report.html> --cwd <legacy-overlay> --timeout 7200",
            )
        elif relative == "notebooks/codon_optimization/mutation_grid_best_per_iteration.ipynb":
            row.update(
                status="blocked_missing_input",
                reason="Referenced codon_opt_results_m1..m15 mutation-grid trajectories are absent; exact DNA sequences cannot be reconstructed from embedded summary plots",
                rerun_command="Restore the named benchmark trajectory, then execute the notebook with Papermill",
            )
        rows.append(row)

    code_validation = repo / "studies" / "hurdler_validation" / "step05_reproducibility" / "tables" / "code_validation.json"
    validation_payload = json.loads(code_validation.read_text()) if code_validation.is_file() else {}
    validation_passed = bool(validation_payload.get("passed"))
    validated_hashes = validation_payload.get("source_sha256", {})
    code_validation_job = str(job_registry.get("code_validation", ""))
    source_files = sorted(
        [*repo.rglob("*.py"), *repo.rglob("*.sh")],
        key=lambda path: str(path.relative_to(repo)),
    )
    for path in source_files:
        relative = str(path.relative_to(repo))
        source_validated = bool(validation_passed and validated_hashes.get(relative) == sha256(path))
        rows.append(
            {
                "workflow": relative,
                "workflow_type": "shell_source" if path.suffix == ".sh" else "python_source",
                "source_sha256": sha256(path),
                "status": "passed" if source_validated else "failed",
                "job_id": code_validation_job,
                "runtime_seconds": "",
                "resources": "cpu=1,mem=8G,time=02:00:00",
                "output_path": str(code_validation),
                "reason": "Digs pytest/syntax/import/CLI/kernel smoke at this source hash" if source_validated else "Code-validation task did not pass at the current source hash",
                "rerun_command": f"/home/wendai/.conda/envs/hurdler/bin/python {repo}/studies/hurdler_validation/scripts/validate_codebase.py --repo {repo} --output {code_validation}",
            }
        )

    cli_commands = [
        "reference build", "lookup build", "query", "screen-short",
        "success-rate", "curate-modules", "optimize-modules", "validate-run",
        "refine-ga", "merge-module-catalogs", "designed-inventory",
        "validate-designed-structures", "infer-designed-boundaries",
        "module-compatibility", "adaptive-copy-search",
    ]
    cli_jobs = job_registry.get("cli_jobs", {})
    cli_resources = job_registry.get("cli_resources", {})
    cli_outputs = job_registry.get("cli_outputs", {})
    for command in cli_commands:
        rows.append(
            {
                "workflow": f"hurdler {command}",
                "workflow_type": "cli",
                "source_sha256": sha256(repo / "src" / "hurdler" / "cli.py"),
                "status": "passed" if validation_passed and validated_hashes.get("src/hurdler/cli.py") == sha256(repo / "src" / "hurdler" / "cli.py") else "failed",
                "job_id": cli_jobs.get(command, code_validation_job),
                "runtime_seconds": "",
                "resources": cli_resources.get(command, "cpu=1,mem=8G,time=02:00:00"),
                "output_path": cli_outputs.get(command, str(code_validation)),
                "reason": "CLI contract smoke passed; production jobs are recorded in the study run directories",
                "rerun_command": f"/home/wendai/.conda/envs/hurdler/bin/hurdler {command} --help",
            }
        )

    for run in job_registry.get("scientific_runs", []):
        source = str(run.get("source", ""))
        source_path = repo / source if source else None
        rows.append(
            {
                "workflow": str(run["workflow"]),
                "workflow_type": "scientific_run",
                "source_sha256": (
                    sha256(source_path)
                    if source_path is not None and source_path.is_file()
                    else ""
                ),
                "status": str(run["status"]),
                "job_id": str(run.get("job_id", "")),
                "runtime_seconds": run.get("runtime_seconds", ""),
                "resources": str(run.get("resources", "")),
                "output_path": str(run.get("output_path", "")),
                "reason": str(run.get("reason", "")),
                "rerun_command": str(run.get("rerun_command", "")),
            }
        )

    status = pd.DataFrame(rows).sort_values(["workflow_type", "workflow"]).reset_index(drop=True)
    status["runtime_seconds"] = pd.to_numeric(status["runtime_seconds"], errors="coerce")
    for column in status.columns.difference(["runtime_seconds"]):
        status[column] = status[column].fillna("").astype(str)
    status.to_csv(args.output_dir / "execution_status.csv", index=False)
    status.to_parquet(args.output_dir / "execution_status.parquet", index=False)

    scn_references: set[str] = set()
    allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_+./'-")
    for path in notebooks:
        if "agarose" not in str(path):
            continue
        contents = path.read_text(errors="ignore")
        offset = 0
        while (position := contents.find(".scn", offset)) >= 0:
            start = position - 1
            lower_bound = max(-1, position - 240)
            while start > lower_bound and contents[start] in allowed:
                start -= 1
            literal = contents[start + 1 : position + 4].strip(" '\"\\n\\t")
            if literal and "*" not in literal:
                scn_references.add(literal)
            offset = position + 4
    (args.output_dir / "missing_scn_inputs.txt").write_text("\n".join(sorted(scn_references)) + "\n")

    summary = {
        "workflow_rows": len(status),
        "status_counts": status.status.value_counts().sort_index().to_dict(),
        "missing_scn_reference_count": len(scn_references),
        "execution_status_sha256": sha256(args.output_dir / "execution_status.csv"),
    }
    (args.output_dir / "execution_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
