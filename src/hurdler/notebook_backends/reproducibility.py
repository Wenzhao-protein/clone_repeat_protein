"""Backend for V2 notebook 11: repository/workspace reproducibility audit."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from ..artifacts import ArtifactRegistry
from ..io import sha256_file
from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "11_reproducibility",
    "Reproducibility dashboard",
    "Validate historical path preservation, registered artifacts and all manifests in the active workspace.",
    production_workflows=("reports",),
)

DEFERRED_LONG = (
    "notebooks/hurdler/hurdler_success_rate_analysis.ipynb",
    "notebooks/hurdler/hurdler_success_rate_optimized.ipynb",
    "notebooks/hurdler/hurdler_site_combination_analysis.ipynb",
    "archive/notebooks/codon_optimization_20230331_1101-Copy3.ipynb",
    "archive/notebooks/hurdler_minimal.ipynb",
    "archive/notebooks/hurdler_standalone.ipynb",
    "archive/notebooks/hurdler_standalone_backup.ipynb",
    "archive/notebooks/hurdler_success_rate_analysis_backup.ipynb",
)


def _historical_scn_names() -> list[str]:
    names: set[str] = set()
    for path in (repo_root() / "notebooks/agarose_gel").glob("*.ipynb"):
        text = path.read_text(errors="ignore")
        for token in text.replace('"', "'").split("'"):
            if token.lower().endswith(".scn") and "*" not in token:
                names.add(Path(token).name)
    return sorted(names)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    baseline = Path(str(request.get("baseline", repo_root() / "data/tracked_path_baseline_v2.json")))
    registry = Path(str(request.get("artifact_registry", repo_root() / "data/artifact_registry_v2.json")))
    for path in (baseline, registry):
        if not path.is_file():
            raise FileNotFoundError(path)
    return {"status": "passed", "baseline": str(baseline), "registry": str(registry)}


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    inputs = preflight(context, request)
    baseline = json.loads(Path(inputs["baseline"]).read_text())
    rows: list[dict[str, Any]] = []
    for relative in baseline["paths"]:
        path = repo_root() / relative
        rows.append({"check": "historical_path", "target": relative, "status": "passed" if path.exists() else "failed", "detail": "" if path.exists() else "baseline path missing"})
    registry = ArtifactRegistry(inputs["registry"])
    for record in registry.list():
        try:
            path = registry.verify(record.artifact_id) if record.repo_path else None
            status, detail = ("passed", str(path or "release-backed_not_downloaded"))
        except Exception as exc:
            status, detail = "failed", f"{type(exc).__name__}: {exc}"
        rows.append({"check": "artifact_registry", "target": record.artifact_id, "status": status, "detail": detail})
    for manifest in sorted(context.workspace_root.glob("*/workspace_manifest.json")):
        try:
            payload = json.loads(manifest.read_text())
            status = "passed" if payload.get("schema_version") == "hurdler-notebook-workspace-v2" else "failed"
            detail = payload.get("backend_id", "")
        except Exception as exc:
            status, detail = "failed", f"{type(exc).__name__}: {exc}"
        rows.append({"check": "workspace_manifest", "target": str(manifest), "status": status, "detail": detail})
    for relative in DEFERRED_LONG:
        if (repo_root() / relative).is_file():
            rows.append({
                "check": "historical_workflow", "target": relative,
                "status": "deferred_long",
                "detail": "Historical monolithic workflow exceeds the two-hour policy; use the matching V2 sharded production bundle.",
                "rerun_command": f"python studies/hurdler_validation/scripts/execute_notebook.py {relative} EXECUTED.ipynb REPORT.html --cwd . --timeout 7200",
            })
    scn_names = _historical_scn_names()
    if scn_names and not list(repo_root().rglob("*.scn")):
        rows.append({
            "check": "historical_input", "target": ", ".join(scn_names),
            "status": "blocked_missing_input",
            "detail": "Historical SCN files are absent; no gel image was fabricated.",
            "rerun_command": "Upload the named SCN files in notebook 12 and rerun Full mode.",
        })
    for notebook in sorted((repo_root() / "notebooks/v2").glob("*.ipynb")):
        try:
            payload = json.loads(notebook.read_text())
            code_cells = [cell for cell in payload.get("cells", []) if cell.get("cell_type") == "code"]
            clean = all(not cell.get("outputs") and cell.get("execution_count") is None for cell in code_cells)
            status, detail = ("passed", "source notebook is output-free") if clean else ("failed", "source notebook contains execution output")
        except Exception as exc:
            status, detail = "failed", f"{type(exc).__name__}: {exc}"
        rows.append({"check": "v2_source_notebook", "target": str(notebook.relative_to(repo_root())), "status": status, "detail": detail})
    frame = pd.DataFrame(rows)
    table_paths = write_frame(frame, context.directory("tables") / "reproducibility_status")
    counts = {str(key): int(value) for key, value in frame.status.value_counts().items()}
    report = context.directory("reports") / "reproducibility_dashboard.html"
    report.write_text(
        "<!doctype html><meta charset='utf-8'><title>HURDLER reproducibility</title>"
        "<style>body{font:14px sans-serif;margin:2rem}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:5px}</style>"
        f"<h1>HURDLER reproducibility dashboard</h1><p>Passed: {counts.get('passed',0)}; "
        f"Blocked: {counts.get('blocked_missing_input',0)}; Deferred: {counts.get('deferred_long',0)}; "
        f"Failed: {counts.get('failed',0)}</p>"
        + frame.to_html(index=False, escape=True)
    )
    summary_path = context.directory("reports") / "reproducibility_summary.json"
    summary_path.write_text(json.dumps({
        "schema_version": "hurdler-reproducibility-v2",
        "status_counts": counts,
        "checks": len(frame),
        "missing_scn_files": scn_names,
        "failed_targets": frame.loc[frame.status.eq("failed"), "target"].astype(str).tolist(),
    }, indent=2, sort_keys=True) + "\n")
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[*table_paths, report, summary_path],
        metrics={
            "checks": len(frame), "passed": int(counts.get("passed", 0)),
            "blocked_missing_input": int(counts.get("blocked_missing_input", 0)),
            "deferred_long": int(counts.get("deferred_long", 0)),
            "failed": int(counts.get("failed", 0)),
        },
        status="passed" if not counts.get("failed") else "failed",
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
