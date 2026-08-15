"""Backend for V2 notebook 06: arbitrary exact-DNA HURDLER design."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from ..exact_dna_design import (
    EXACT_DNA_SCHEMA_VERSION,
    ExactDNAQuery,
    ExactDNASelection,
    confirm_best_exact_dna_route,
    confirm_exact_dna_route,
    query_exact_dna,
    write_exact_dna_outputs,
)
from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "06_exact_dna_designer",
    "Exact-DNA designer V2",
    "Plan exact active/latent-RE seed-to-target assembly without changing any final target base.",
    production_workflows=("exact-dna-routes", "exact-dna-purchase"),
    default_request={
        "query": {
            "schema_version": EXACT_DNA_SCHEMA_VERSION,
            "input_mode": "array",
            "sequence_id": "tutorial_regulatory_array",
            "repeat_unit": "GCGTATACGCGTATACGCGTATAC",
            "spacer": "",
            "repeat_copies": 2,
            "max_restoration_length_bp": 100,
            "max_states": 250,
            "timeout_seconds": 5,
            "paths_per_state": 2,
            "max_complete_routes": 5,
        },
        "confirm": False,
    },
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _query(request: Mapping[str, Any]) -> ExactDNAQuery:
    payload = dict(SPEC.default_request["query"])
    payload.update(dict(request.get("query", {})))
    return ExactDNAQuery.from_dict(payload)


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    query = _query(request)
    reference = Path(str(request.get("reference_dir", repo_root() / "data/reference_output")))
    plasmids = Path(str(request.get("plasmid_reference", reference / "plasmid_reference_v2.json")))
    for path in (reference / "restriction_enzyme.csv", plasmids):
        if not path.is_file():
            raise FileNotFoundError(path)
    return {"status": "passed", "sequence_id": query.sequence_id, "target_length_bp": len(query.target_sequence)}


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    preflight(context, request)
    query = _query(request)
    reference = Path(str(request.get("reference_dir", repo_root() / "data/reference_output")))
    artifacts = Path(str(request.get("artifact_dir", repo_root() / "data/artifacts")))
    plasmids = Path(str(request.get("plasmid_reference", reference / "plasmid_reference_v2.json")))
    result = query_exact_dna(
        query,
        reference_dir=reference,
        artifact_dir=artifacts,
        plasmid_reference_path=plasmids,
        progress_callback=progress_callback,
    )
    if request.get("confirm"):
        if request.get("selection"):
            selection = ExactDNASelection(**dict(request["selection"]))
            result = confirm_exact_dna_route(result, selection, progress_callback=progress_callback)
        else:
            if not result.route_candidates:
                raise ValueError("No exact-DNA route is available to confirm")
            selection = ExactDNASelection(
                str(result.route_candidates[0]["route_id"]), "batch"
            )
            result = confirm_best_exact_dna_route(
                result, selection, progress_callback=progress_callback
            )
        written = write_exact_dna_outputs(result, context.run_root / "design")
        paths = [Path(value) for value in written.values()]
    else:
        report = context.directory("reports") / "exact_dna_query.json"
        report.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True, default=str) + "\n")
        paths = [report, *write_frame(pd.DataFrame(result.route_candidates), context.directory("tables") / "exact_dna_routes")]
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=paths,
        metrics={
            "designer_status": result.status,
            "target_length_bp": result.target_length_bp,
            "pair_candidate_count": len(result.pair_candidates),
            "route_candidate_count": len(result.route_candidates),
            "final_target_exact": bool(result.independent_verification.get("passed")),
        },
        next_notebooks=["07_production_builder", "10_exact_dna_result_analysis"],
        limitations=[] if request.get("confirm") else ["Tutorial search does not perform IDT scoring or ordering."],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
