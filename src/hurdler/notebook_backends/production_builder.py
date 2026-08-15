"""Backend for V2 notebook 07: production bundle builder/importer."""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import Any, Mapping

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from ..production_bundle import (
    WORKFLOWS,
    ProductionBundleRequest,
    build_production_bundle,
    validate_production_bundle,
)
from .common import BackendSpec, result_from_paths, zip_paths
from .common import repo_root


SPEC = BackendSpec(
    "07_production_builder",
    "Production bundle builder",
    "Generate auditable Digs task files, submission scripts, recovery and finalization without submitting from Colab.",
    production_workflows=tuple(WORKFLOWS),
    default_request={
        "production_request": {
            "workflow_id": "reports",
            "parameter_version": "v2",
            "repo_commit": "AUTO",
            "cluster_profile": {
                "repo_root": "AUTO",
                "scratch_root": "AUTO",
                "conda_prefix": "AUTO",
                "taskrunner": "/net/software/taskrunner/taskrunner",
                "partition": "cpu",
                "cpu_per_task": 1,
                "memory": "8G",
                "walltime": "02:00:00",
                "array_throttle": 1
            },
            "inputs": [{"artifact_id": "tutorial-report", "path": "AUTO", "sha256": "AUTO"}],
            "shard_count": 1,
            "scientific_parameters": {},
            "random_seed": 42,
            "idt_mode": "batch",
            "output_dir": "AUTO",
            "scratch_dir": "AUTO"
        }
    },
)


def get_spec() -> dict[str, Any]:
    payload = SPEC.to_dict()
    payload["workflows"] = {key: value.__dict__ for key, value in WORKFLOWS.items()}
    return payload


def _production_request(
    context: NotebookContext, request: Mapping[str, Any]
) -> ProductionBundleRequest:
    payload = json.loads(json.dumps(request["production_request"]))
    profile = payload["cluster_profile"]
    if profile.get("repo_root") == "AUTO":
        profile["repo_root"] = str(repo_root())
    if profile.get("scratch_root") == "AUTO":
        profile["scratch_root"] = str(context.run_root / "production_scratch")
    if profile.get("conda_prefix") == "AUTO":
        profile["conda_prefix"] = sys.prefix
    if payload.get("repo_commit") == "AUTO":
        payload["repo_commit"] = context.repo_commit
    if payload.get("output_dir") == "AUTO":
        payload["output_dir"] = str(context.run_root / "production_output")
    if payload.get("scratch_dir") == "AUTO":
        payload["scratch_dir"] = str(context.run_root / "production_scratch")
    for item in payload.get("inputs", []):
        if item.get("path") == "AUTO":
            item["path"] = str(repo_root() / "notebooks/v2/11_reproducibility.ipynb")
        if item.get("sha256") == "AUTO":
            from ..io import sha256_file

            item["sha256"] = sha256_file(item["path"])
    return ProductionBundleRequest.from_dict(payload)


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    production = _production_request(context, request)
    spec = production.validate()
    return {"status": "passed", "workflow": spec.workflow_id, "default_shards": spec.default_shards}


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    preflight(context, request)
    production = _production_request(context, request)
    bundle = context.run_root / "production_bundle"
    build_production_bundle(production, bundle)
    validation = validate_production_bundle(bundle)
    archive = zip_paths(context.run_root / f"{production.workflow_id}_production_bundle.zip", bundle.rglob("*"), base=bundle)
    preview = context.directory("reports") / "production_bundle_preview.json"
    tasks = (bundle / "tasks.txt").read_text().splitlines()
    preview.write_text(json.dumps({
        **validation,
        "first_tasks": tasks[:3],
        "last_tasks": tasks[-3:],
        "contains_credentials": False,
        "submitted": False,
    }, indent=2, sort_keys=True) + "\n")
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[archive, preview, bundle / "request.json", bundle / "task_index.csv", bundle / "submit_digs.sh"],
        metrics={**validation, "submitted": False},
        next_notebooks=["08_success_landscape_analysis", "09_module_result_analysis", "10_exact_dna_result_analysis", "11_reproducibility"],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
