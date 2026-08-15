"""Backend for V2 notebook 05: repeat-protein designer orchestration."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from ..vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    design_construct_v2,
    design_query,
    write_design_outputs_v2,
)
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "05_repeat_designer",
    "Repeat-protein designer V3",
    "Select an annotation-safe HURDLER route and optionally run the existing GA/IDT design engine.",
    production_workflows=("module-stage2",),
    default_request={
        "query": {
            "schema_version": DESIGN_SCHEMA_VERSION_V2,
            "input_mode": "split",
            "sequence_id": "tutorial_DHR12",
            "n_cap": "",
            "repeat_module": "TDDEEIARIIAYAARQTT",
            "c_cap": "",
            "repeat_copies": 4,
            "max_restoration_length_bp": 100,
        },
        "run_ga": False,
    },
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _query(request: Mapping[str, Any]) -> CompatibilityQuery:
    payload = dict(SPEC.default_request["query"])
    payload.update(dict(request.get("query", {})))
    return CompatibilityQuery.from_dict(payload)


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    query = _query(request)
    index = Path(str(request.get("protein_index_dir", repo_root() / "data/artifacts/vector-aware-hurdler-v2")))
    plasmids = Path(str(request.get("plasmid_reference", repo_root() / "data/reference_output/plasmid_reference_v2.json")))
    for path in (index / "metadata.json", plasmids):
        if not path.is_file():
            raise FileNotFoundError(path)
    return {"status": "passed", "sequence_id": query.sequence_id, "run_ga": bool(request.get("run_ga", False))}


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    preflight(context, request)
    index = Path(str(request.get("protein_index_dir", repo_root() / "data/artifacts/vector-aware-hurdler-v2")))
    plasmids = Path(str(request.get("plasmid_reference", repo_root() / "data/reference_output/plasmid_reference_v2.json")))
    query = _query(request)
    paths: list[Path] = []
    if request.get("run_ga"):
        if "design_request" not in request:
            raise ValueError("run_ga=True requires a confirmed design_request with DesignSelection")
        design_request = DesignRequestV2.from_dict(dict(request["design_request"]))
        if design_request.validation_mode == "api":
            raise ValueError("The backend accepts an injected in-memory IDT scorer; credentials never belong in request data")
        result = design_construct_v2(
            design_request,
            protein_index_dir=index,
            plasmid_reference_path=plasmids,
            progress_callback=progress_callback,
        )
        written = write_design_outputs_v2(result, context.run_root / "design")
        paths.extend(Path(value) for value in written.values())
        status = result.status
        route_count = len(result.vector_routes)
        selected = bool(result.selected_route)
    else:
        result = design_query(query, protein_index_dir=index, plasmid_reference_path=plasmids)
        result_json = context.directory("reports") / "repeat_design_query.json"
        result_json.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True, default=str) + "\n")
        paths.append(result_json)
        route_paths = write_frame(pd.DataFrame(result.vector_routes), context.directory("tables") / "repeat_design_routes")
        paths.extend(route_paths)
        status = result.status
        route_count = len(result.vector_routes)
        selected = False
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=paths,
        metrics={"designer_status": status, "route_count": route_count, "selected_route": selected},
        next_notebooks=["07_production_builder", "09_module_result_analysis"],
        status="passed" if status not in {"failed", "error"} else "failed",
        limitations=[] if request.get("run_ga") else ["Tutorial mode performs route selection only; GA is opt-in."],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
