"""Backend for V2 notebook 01: reference database builder."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from ..reference import REFERENCE_FILES, REFERENCE_METADATA_FILES, build_reference_manifest
from .common import BackendSpec, preflight_files, repo_root, result_from_paths, zip_paths


SPEC = BackendSpec(
    "01_reference_database",
    "Reference database builder",
    "Validate and package the frozen REBASE, methylation, buffer, fidelity, codon and plasmid references.",
    default_request={"reference_dir": "data/reference_output"},
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _reference_dir(request: Mapping[str, Any]) -> Path:
    value = Path(str(request.get("reference_dir", "data/reference_output")))
    return value if value.is_absolute() else repo_root() / value


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    root = _reference_dir(request)
    return preflight_files(root / name for name in (*REFERENCE_FILES, *REFERENCE_METADATA_FILES))


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    preflight(context, request)
    if progress_callback:
        progress_callback({"stage": "reference", "status": "started"})
    reference = _reference_dir(request)
    manifest_path = context.directory("reports") / "reference_manifest_v2.json"
    manifest = build_reference_manifest(reference, manifest_path)
    funnel_rows: list[dict[str, Any]] = []
    for entry in manifest["files"]:
        if entry["rows"] is not None:
            funnel_rows.append(
                {
                    "source_table": entry["name"],
                    "input_rows": int(entry["rows"]),
                    "retained_rows": int(entry["rows"]),
                    "rejected_rows": 0,
                    "rejection_reason": "snapshot_validation_only",
                    "sha256": entry["sha256"],
                }
            )
    funnel = pd.DataFrame(funnel_rows)
    funnel_path = context.directory("tables") / "reference_selection_funnel.csv"
    funnel.to_csv(funnel_path, index=False)
    archive = context.run_root / "reference_database_v2.zip"
    zip_paths(
        archive,
        [reference / name for name in (*REFERENCE_FILES, *REFERENCE_METADATA_FILES)]
        + [manifest_path, funnel_path],
        base=repo_root(),
    )
    result = result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[manifest_path, funnel_path, archive],
        metrics={
            "reference_file_count": len(manifest["files"]),
            "total_csv_rows": int(funnel.input_rows.sum()),
            "source_mode": context.source_mode,
        },
        next_notebooks=["02_lookup_plasmid"],
        limitations=(
            ["Restricted-source refresh requires user-supplied files; snapshots were not replaced."]
            if context.source_mode != "refresh"
            else []
        ),
    )
    if progress_callback:
        progress_callback({"stage": "reference", "status": "completed", "files": len(manifest["files"])})
    return result


def write_outputs(
    context: NotebookContext, result: NotebookResult
) -> dict[str, Any]:
    return result.to_dict()
