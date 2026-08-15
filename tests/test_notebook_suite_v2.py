from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import nbformat
import pandas as pd
import pytest
from nbclient import NotebookClient

from hurdler.artifacts import ArtifactRegistry
from hurdler.notebook_backends import BACKEND_MODULES, load_backend
from hurdler.notebook_workspace import (
    NotebookContext,
    NotebookResult,
    export_workspace,
    import_workspace,
    write_run_manifest,
)
from hurdler.production_bundle import (
    ClusterProfile,
    ProductionBundleRequest,
    build_production_bundle,
    directory_manifest_sha256,
    validate_production_bundle,
)


REPO = Path(__file__).resolve().parents[1]


def test_tracked_path_baseline_is_intact():
    payload = json.loads((REPO / "data/tracked_path_baseline_v2.json").read_text())
    assert payload["tracked_path_count"] == 1119
    assert len(payload["paths"]) == len(set(payload["paths"]))
    missing = [path for path in payload["paths"] if not (REPO / path).exists()]
    assert missing == []
    text = "\n".join(payload["paths"]) + "\n"
    assert hashlib.sha256(text.encode()).hexdigest() == payload["sorted_paths_sha256"]


def test_workspace_round_trip_and_secret_rejection(tmp_path: Path):
    context = NotebookContext("workspace_test", workspace_root=tmp_path).prepare()
    table = context.directory("tables") / "small.csv"
    table.write_text("x\n1\n")
    result = NotebookResult("passed", metrics={"rows": 1})
    write_run_manifest(context, backend_id="test", request={"mode": "fixture"}, result=result)
    archive = export_workspace(context)
    imported = import_workspace(archive, tmp_path / "imported")
    assert imported.run_id == context.run_id
    with pytest.raises(ValueError, match="Secret-bearing"):
        write_run_manifest(
            context,
            backend_id="test",
            request={"access_token": "must-never-be-written"},
            result=NotebookResult("failed"),
        )


def test_artifact_registry_verifies_every_bundled_artifact():
    registry = ArtifactRegistry(REPO / "data/artifact_registry_v2.json")
    assert len(registry.list()) >= 12
    for record in registry.list():
        if record.repo_path:
            assert registry.verify(record.artifact_id).is_file()


def _report_request(tmp_path: Path) -> ProductionBundleRequest:
    source = REPO / "notebooks/v2/11_reproducibility.ipynb"
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    return ProductionBundleRequest(
        workflow_id="reports",
        parameter_version="v2",
        repo_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        cluster_profile=ClusterProfile(
            repo_root=str(REPO),
            scratch_root=str(tmp_path / "scratch-root"),
            conda_prefix=sys.prefix,
            taskrunner="/bin/true",
            partition="cpu",
            cpu_per_task=1,
            memory="8G",
            walltime="02:00:00",
            array_throttle=1,
        ),
        inputs=({"artifact_id": "notebook", "path": str(source), "sha256": digest},),
        shard_count=1,
        idt_mode="batch",
        output_dir=str(tmp_path / "results"),
        scratch_dir=str(tmp_path / "scratch"),
    )


def test_production_bundle_is_complete_safe_and_stable(tmp_path: Path):
    request = _report_request(tmp_path)
    first = build_production_bundle(request, tmp_path / "bundle-a")
    result = validate_production_bundle(first)
    assert result == {"status": "passed", "workflow": "reports", "task_count": 1}
    assert (first / "tasks.txt").read_text().startswith("/")
    assert "submit)" in (first / "submit_digs.sh").read_text()
    assert "Refusing to overwrite taskrunner state" in (first / "submit_digs.sh").read_text()
    assert json.loads((first / "request.json").read_text())["cluster_profile"]["idt_env_path"] == ""
    with pytest.raises(FileExistsError):
        build_production_bundle(request, first)


def test_all_nine_production_workflows_generate_expected_default_tasks(tmp_path: Path):
    table = tmp_path / "inventory.parquet"
    pd.DataFrame(
        {"module_id": [f"m{index}" for index in range(10)],
         "hurdler_compatible": [True] * 9 + [False]}
    ).to_parquet(table, index=False)
    raw = tmp_path / "raw-routes"
    raw.mkdir()
    (raw / "manifest.txt").write_text("fixture\n")
    notebook = REPO / "notebooks/v2/11_reproducibility.ipynb"
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    expected = {
        "success-landscape": 1,
        "repeatsdb-natural": 128,
        "designed-structure": 10,
        "missing-af3": 10,
        "module-stage1": 128,
        "module-stage2": 3,
        "exact-dna-routes": 512,
        "exact-dna-purchase": 1,
        "reports": 1,
    }
    for workflow, task_count in expected.items():
        source = notebook if workflow == "reports" else table
        inputs = ({
            "artifact_id": "fixture",
            "path": str(source),
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        },)
        if workflow == "success-landscape":
            inputs = ()
        elif workflow == "exact-dna-purchase":
            inputs = ({
                "artifact_id": "raw-route-directory",
                "path": str(raw),
                "sha256": directory_manifest_sha256(raw),
                "kind": "directory",
            },)
        spec_profile = {
            "success-landscape": (16, "32G", "02:00:00", 1),
            "repeatsdb-natural": (8, "16G", "02:00:00", 16),
            "designed-structure": (4, "8G", "02:00:00", 16),
            "missing-af3": (4, "32G", "08:00:00", 4),
            "module-stage1": (1, "8G", "02:00:00", 16),
            "module-stage2": (1, "8G", "02:00:00", 16),
            "exact-dna-routes": (1, "8G", "02:00:00", 16),
            "exact-dna-purchase": (1, "8G", "00:30:00", 1),
            "reports": (1, "8G", "02:00:00", 1),
        }[workflow]
        cpu, memory, walltime, throttle = spec_profile
        request = ProductionBundleRequest(
            workflow_id=workflow,
            parameter_version="v2",
            repo_commit=commit,
            cluster_profile=ClusterProfile(
                repo_root=str(REPO), scratch_root=str(tmp_path / "cluster-scratch"),
                conda_prefix=sys.prefix, taskrunner="/bin/true",
                af3_runner="/bin/true", gpu="a100" if workflow == "missing-af3" else "",
                cpu_per_task=cpu, memory=memory, walltime=walltime,
                array_throttle=throttle,
            ),
            inputs=inputs,
            idt_mode="batch",
            output_dir=str(tmp_path / workflow / "results"),
            scratch_dir=str(tmp_path / workflow / "scratch"),
        )
        bundle = build_production_bundle(request, tmp_path / workflow / "bundle")
        assert validate_production_bundle(bundle)["task_count"] == task_count
        assert "OMP_NUM_THREADS=1" in (bundle / "tasks.txt").read_text()


def test_all_v2_backend_interfaces_are_uniform():
    assert len(BACKEND_MODULES) == 14
    for notebook_id in BACKEND_MODULES:
        backend = load_backend(notebook_id)
        assert backend.get_spec()["notebook_id"] == notebook_id
        assert callable(backend.preflight)
        assert callable(backend.run)
        assert callable(backend.write_outputs)


def test_generated_notebooks_are_stable_output_free_and_portable(tmp_path: Path):
    command = [sys.executable, str(REPO / "scripts/generate_notebook_suite_v2.py"), "--output-dir", str(tmp_path / "first")]
    subprocess.run(command, cwd=REPO, check=True)
    command[-1] = str(tmp_path / "second")
    subprocess.run(command, cwd=REPO, check=True)
    first = sorted((tmp_path / "first").glob("*.ipynb"))
    second = sorted((tmp_path / "second").glob("*.ipynb"))
    assert len(first) == len(second) == 14
    for left, right in zip(first, second, strict=True):
        assert left.name == right.name
        assert left.read_bytes() == right.read_bytes()
        notebook = nbformat.read(left, as_version=4)
        assert len({cell.id for cell in notebook.cells}) == len(notebook.cells)
        assert all(cell.get("execution_count") is None for cell in notebook.cells if cell.cell_type == "code")
        assert all(not cell.get("outputs") for cell in notebook.cells if cell.cell_type == "code")
        text = left.read_text()
        assert "/home/wendai" not in text
        assert "/net/scratch/wendai" not in text
        assert "agent/vector-aware-designer-v2" not in text


def test_compact_scientific_artifacts_have_frozen_counts():
    success = pd.read_parquet(REPO / "data/results/success_landscape_compact_v2.parquet")
    assert len(success) == 60 * 8
    assert success.groupby("module_length").plasmid.nunique().eq(8).all()
    for length in range(1, 6):
        assert set(success.loc[success.module_length.eq(length), "tests"]) == {20**length}
    modules = pd.read_parquet(REPO / "data/results/module_analysis_compact_v2.parquet")
    assert modules.collection.value_counts().to_dict() == {"Natural": 25_913, "Designed": 182}
    three = pd.read_csv(REPO / "data/results/repeatsdb_designed_hurdler_3mer_results.csv")
    compatible = three.hurdler_compatible.astype(bool)
    assert len(three) == len(modules) == 26_095
    assert three.loc[compatible, ["selected_re_pair", "site_i_3mer_aa", "site_ii_3mer_aa"]].ne("").all().all()
    exact_targets = pd.read_parquet(REPO / "data/results/exact_dna_target_analysis_compact_v2.parquet")
    exact_elements = pd.read_parquet(REPO / "data/results/exact_dna_element_matrix_compact_v2.parquet")
    exact_routes = pd.read_parquet(REPO / "data/results/exact_dna_selected_routes_compact_v2.parquet")
    assert len(exact_targets) == 145_210
    assert len(exact_elements) == 29_042
    assert len(exact_routes) == int(exact_targets.complete_route_verified.sum()) == 15_535
    assert exact_targets.groupby(["source_database", "element_id"]).target_copy_count.agg(set).eq({2, 4, 8, 16, 32}).all()
    assert not (exact_targets.fragment_rescued_by_hurdler & ~exact_targets.complete_route_verified).any()


@pytest.mark.parametrize("notebook_id,filename", [
    (entry["id"], entry["file"])
    for entry in json.loads((REPO / "notebooks/v2/catalog.json").read_text())["notebooks"]
])
def test_v2_tutorial_notebook_executes_in_clean_kernel(
    notebook_id: str, filename: str, tmp_path: Path
):
    notebook = nbformat.read(REPO / "notebooks/v2" / filename, as_version=4)
    parameter = next(
        cell for cell in notebook.cells
        if cell.cell_type == "code" and "parameters" in cell.metadata.get("tags", [])
    )
    lines = []
    for line in parameter.source.splitlines():
        if line.startswith("RUN_ID ="):
            line = f"RUN_ID = {notebook_id!r}"
        elif line.startswith("WORKSPACE_ROOT ="):
            line = f"WORKSPACE_ROOT = {str(tmp_path / 'workspace')!r}"
        lines.append(line)
    parameter.source = "\n".join(lines)
    client = NotebookClient(
        notebook,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(REPO)}},
    )
    executed = client.execute()
    assert all(
        output.get("output_type") != "error"
        for cell in executed.cells if cell.cell_type == "code"
        for output in cell.get("outputs", [])
    )
