from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

from hurdler.cli import main
from hurdler.design import bundled_index_dir


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "workflows" / "01_interactive_hurdler_designer.ipynb"
COLAB_NOTEBOOK = ROOT / "notebooks" / "workflows" / "02_colab_hurdler_designer.ipynb"


def test_notebook_source_is_output_free_and_has_no_widget_state():
    payload = json.loads(NOTEBOOK.read_text())
    assert payload["nbformat"] == 4
    assert "widgets" not in payload.get("metadata", {})
    assert any("parameters" in cell.get("metadata", {}).get("tags", []) for cell in payload["cells"])
    for cell in payload["cells"]:
        if cell["cell_type"] == "code":
            assert cell.get("outputs", []) == []
            assert cell.get("execution_count") is None
    text = NOTEBOOK.read_text()
    assert "IDT_CLIENT_SECRET=" not in text
    assert "IDT_PASSWORD=" not in text
    assert ".config/hurdler/idt.env" not in text


def test_colab_cells_are_named_hidden_forms():
    payload = json.loads(COLAB_NOTEBOOK.read_text())
    ids = []
    for cell in payload["cells"]:
        cell_id = cell.get("metadata", {}).get("id")
        assert cell_id
        ids.append(cell_id)
        if cell["cell_type"] != "code":
            continue
        assert cell["metadata"]["cellView"] == "form"
        assert cell["metadata"]["colab"] == {}
        assert "".join(cell["source"]).startswith("#@title ")
        assert cell.get("outputs", []) == []
        assert cell.get("execution_count") is None
    assert len(ids) == len(set(ids))


def test_readme_notebook_links_resolve():
    assert NOTEBOOK.is_file()
    assert COLAB_NOTEBOOK.is_file()
    assert (ROOT / "apps" / "hurdler_designer.py").is_file()
    assert (ROOT / "scripts" / "start_hurdler_web.sh").is_file()
    assert (ROOT / "notebooks" / "README.md").is_file()
    readme = (ROOT / "README.md").read_text()
    assert "notebooks/workflows/01_interactive_hurdler_designer.ipynb" in readme
    assert "./scripts/start_hurdler_web.sh" in readme
    assert "codespaces" not in readme.lower()


def test_bundled_index_checksum_manifest():
    root = bundled_index_dir()
    assert (root / "pattern_index.npz").stat().st_size > 3_000_000
    entries = [line.split(None, 1) for line in (root / "SHA256SUMS").read_text().splitlines()]
    assert len(entries) == 13
    for expected, relative in entries:
        path = root / relative
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


def test_design_construct_cli_writes_topology_bundle(tmp_path, capsys):
    output = tmp_path / "cli"
    assert main(
        [
            "design-construct",
            "--request",
            str(ROOT / "data" / "example_design_request.json"),
            "--output-dir",
            str(output),
            "--legacy-v1",
        ]
    ) == 0
    response = json.loads(capsys.readouterr().out)
    assert response["status"] == "hurdler_incompatible"
    assert (output / "design_summary.json").is_file()
    assert not (output / "optimized_construct.fasta").exists()


def test_unversioned_v1_request_requires_explicit_compatibility_flag(tmp_path):
    with pytest.raises(SystemExit, match="--legacy-v1"):
        main(
            [
                "design-construct",
                "--request",
                str(ROOT / "data" / "example_design_request.json"),
                "--output-dir",
                str(tmp_path / "rejected"),
            ]
        )


@pytest.mark.parametrize("notebook_path", [NOTEBOOK, COLAB_NOTEBOOK])
def test_notebook_headless_mock_executes_clean_kernel(
    tmp_path, monkeypatch, notebook_path
):
    monkeypatch.setenv("HURDLER_NOTEBOOK_SMOKE", "1")
    monkeypatch.setenv("HURDLER_NOTEBOOK_SMOKE_OUTPUT", str(tmp_path / "design"))
    notebook = nbformat.read(notebook_path, as_version=4)
    executed = NotebookClient(
        notebook,
        timeout=180,
        kernel_name="python3",
        resources={"metadata": {"path": str(ROOT)}},
    ).execute(cwd=str(ROOT))
    assert all(
        output.get("output_type") != "error"
        for cell in executed.cells
        for output in cell.get("outputs", [])
    )
    assert (tmp_path / "design" / "optimized_construct.fasta").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.csv").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.tsv").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.fasta").stat().st_size > 0
    summary = json.loads((tmp_path / "design" / "design_summary.json").read_text())
    assert summary["status"] == "optimized_unvalidated_batch"
    assert summary["idt_audit"] == []
