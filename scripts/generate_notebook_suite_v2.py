#!/usr/bin/env python3
"""Generate all 14 stable, output-free V2 notebooks from the catalog."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import pprint
from pathlib import Path
from typing import Any

import nbformat


REPO = Path(__file__).resolve().parents[1]
CATALOG = REPO / "notebooks/v2/catalog.json"
NOTEBOOK_DIR = CATALOG.parent
GITHUB_REPO = "https://github.com/Wenzhao-protein/clone_repeat_protein"


def _cell_id(notebook_id: str, role: str) -> str:
    return hashlib.sha1(f"{notebook_id}:{role}".encode()).hexdigest()[:12]


def _markdown(notebook_id: str, role: str, source: str):
    cell = nbformat.v4.new_markdown_cell(source)
    cell.id = _cell_id(notebook_id, role)
    return cell


def _code(notebook_id: str, role: str, source: str, *, tags: list[str] | None = None):
    cell = nbformat.v4.new_code_cell(source)
    cell.id = _cell_id(notebook_id, role)
    cell.execution_count = None
    cell.outputs = []
    if tags:
        cell.metadata["tags"] = tags
    return cell


def _default_request(module_name: str) -> dict[str, Any]:
    module = __import__(f"hurdler.notebook_backends.{module_name}", fromlist=["get_spec"])
    return dict(module.get_spec().get("default_request") or {})


def _designer_notebook(entry: dict[str, str], *, git_ref: str):
    """Reuse the maintained interactive Colab generators for V2 designers."""
    if entry["id"] == "05_repeat_designer":
        path = REPO / "scripts/generate_vector_designer_notebooks.py"
    elif entry["id"] == "06_exact_dna_designer":
        path = REPO / "scripts/generate_exact_dna_colab.py"
    else:
        raise ValueError(entry["id"])
    module_spec = importlib.util.spec_from_file_location(
        f"hurdler_v2_generator_{entry['id']}", path
    )
    if module_spec is None or module_spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    book = module.notebook(colab=True) if entry["id"] == "05_repeat_designer" else module.notebook()
    notebook_id = entry["id"]
    backend = entry["backend"]
    colab_url = f"https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/{git_ref}/notebooks/v2/{entry['file']}"
    banner = _markdown(
        notebook_id,
        "v2-banner",
        f"# {notebook_id}: {entry['title']}\n\n"
        "This is the interactive application in the HURDLER Notebook Suite V2. "
        "Its widgets are generated from the maintained designer UI while route, GA, IDT, "
        "GenBank and verification logic remains in the tested Python package.\n\n"
        f"[Open the current `main` version in Colab]({colab_url})",
    )
    parameters = _code(
        notebook_id,
        "parameters",
        "import datetime, pathlib\n"
        "MODE = 'tutorial'\n"
        f"RUN_ID = '{notebook_id}_' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')\n"
        "WORKSPACE_ROOT = '/content/hurdler_workspace' if pathlib.Path('/content').exists() else '/tmp/hurdler_workspace'\n"
        "RANDOM_SEED = 42",
        tags=["parameters"],
    )
    backend_contract = _code(
        notebook_id,
        "v2-backend-contract",
        f"from hurdler.notebook_backends import {backend} as v2_backend\n"
        "V2_BACKEND_SPEC = v2_backend.get_spec()\n"
        "try:\n"
        "    from IPython.display import JSON, display\n"
        "    display(JSON(V2_BACKEND_SPEC, expanded=False))\n"
        "except Exception:\n"
        "    pass",
    )
    application_cells = list(book.cells)
    for index, cell in enumerate(application_cells):
        stable_id = _cell_id(notebook_id, f"application-{index:02d}")
        cell.id = stable_id
        cell.metadata["id"] = stable_id
    book.cells = [banner, parameters, *application_cells, backend_contract]
    for cell in book.cells:
        if cell.cell_type == "code":
            cell.execution_count = None
            cell.outputs = []
    book.metadata["colab"] = {"name": entry["file"], "provenance": [], "toc_visible": True}
    book.metadata["hurdler"] = {
        "notebook_id": notebook_id, "backend": backend,
        "suite_version": "2026.1", "source_output_free": True,
        "interactive_application": True,
    }
    return book


def build_notebook(entry: dict[str, str], *, git_ref: str) -> Any:
    if entry["id"] in {"05_repeat_designer", "06_exact_dna_designer"}:
        return _designer_notebook(entry, git_ref=git_ref)
    notebook_id = entry["id"]
    backend = entry["backend"]
    title = entry["title"]
    default_request = _default_request(backend)
    request_python = pprint.pformat(default_request, sort_dicts=True, width=88)
    colab_url = f"https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/{git_ref}/notebooks/v2/{entry['file']}"
    cells = [
        _markdown(
            notebook_id,
            "title",
            f"# {notebook_id}: {title}\n\n"
            "This tutorial is one authoritative HURDLER V2 entry point. It uses a thin notebook and a tested Python backend, records input hashes and limitations, and exports a credential-free workspace.\n\n"
            f"[Open the current `main` version in Colab]({colab_url})",
        ),
        _markdown(
            notebook_id,
            "modes",
            "## 1. Choose a mode and data policy\n\n"
            "`tutorial` uses committed fixtures; `colab_full` uses frozen snapshots or explicit refresh/upload; "
            "`production_bundle` writes cluster files but never submits; `analyze` consumes finalized compact results. "
            "Do not place IDT, Google, GitHub or cluster secrets in `REQUEST`.",
        ),
        _code(
            notebook_id,
            "bootstrap",
            "import datetime, importlib.util, os, pathlib, subprocess, sys\n"
            "if importlib.util.find_spec('hurdler') is None:\n"
            "    checkout = pathlib.Path('/content/clone_repeat_protein')\n"
            "    if not checkout.exists():\n"
            f"        subprocess.check_call(['git', 'clone', '--depth', '1', '--branch', '{git_ref}', '{GITHUB_REPO}.git', str(checkout)])\n"
            "    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', '-e', str(checkout) + '[notebooks]'])\n"
            "    os.chdir(checkout)\n"
            "from IPython.display import JSON, Markdown, display\n"
            "from hurdler.notebook_workspace import NotebookContext, export_workspace\n"
            f"from hurdler.notebook_backends import {backend} as backend\n"
            "display(Markdown('**HURDLER backend ready.**'))",
        ),
        _markdown(notebook_id, "inputs_header", "## 2. Inputs\n\nEdit only this parameter cell for a normal run. Paths may be uploaded Colab files or artifacts imported from a previous workspace."),
        _code(
            notebook_id,
            "parameters",
            "MODE = 'tutorial'  # tutorial | colab_full | production_bundle | analyze\n"
            "SOURCE_MODE = 'snapshot'  # snapshot | refresh | upload\n"
            f"RUN_ID = '{notebook_id}_' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')\n"
            "WORKSPACE_ROOT = '/content/hurdler_workspace' if pathlib.Path('/content').exists() else '/tmp/hurdler_workspace'\n"
            f"REQUEST = {request_python}\n"
            "RANDOM_SEED = 42",
            tags=["parameters"],
        ),
        _code(
            notebook_id,
            "context",
            "context = NotebookContext(\n"
            "    run_id=RUN_ID, mode=MODE, source_mode=SOURCE_MODE,\n"
            "    workspace_root=WORKSPACE_ROOT, random_seed=RANDOM_SEED,\n"
            ").prepare()\n"
            "display(JSON(backend.get_spec(), expanded=False))",
        ),
        _markdown(notebook_id, "preflight_header", "## 3. Preflight\n\nThis stops on a missing input, checksum/schema mismatch or unavailable required external tool. It never silently substitutes another artifact."),
        _code(
            notebook_id,
            "preflight",
            "preflight_result = backend.preflight(context, REQUEST)\n"
            "display(JSON(preflight_result, expanded=False))",
        ),
        _markdown(notebook_id, "run_header", "## 4. Run\n\nTutorial work runs in Colab. Heavy production work is exported through notebook 07 and executed on Digs."),
        _code(
            notebook_id,
            "run",
            "progress_events = []\n"
            "def on_progress(event):\n"
            "    progress_events.append(event.to_dict() if hasattr(event, 'to_dict') else dict(event))\n"
            "result = backend.run(context, REQUEST, progress_callback=on_progress)\n"
            "backend.write_outputs(context, result)\n"
            "display(JSON(result.to_dict(), expanded=False))",
        ),
        _markdown(notebook_id, "production", "## 5. Production export\n\nIf this workflow declares a production requirement, open notebook 07, select its workflow ID, preview every task, and download the bundle. Colab does not submit jobs or SSH to Digs."),
        _code(
            notebook_id,
            "progress",
            "display(JSON({'progress_event_count': len(progress_events), 'last_event': progress_events[-1] if progress_events else None}, expanded=False))",
        ),
        _markdown(notebook_id, "results_header", "## 6. Results and provenance\n\nAll tables and figures are generated by the backend. The manifest records the repository commit, mode, source policy, warnings, limitations and next notebook IDs."),
        _code(
            notebook_id,
            "download",
            "workspace_zip = export_workspace(context)\n"
            "display(Markdown(f'Workspace ready: `{workspace_zip}`'))\n"
            "try:\n"
            "    from google.colab import files\n"
            "except ImportError:\n"
            "    files = None\n"
            "# In Colab, uncomment only when you are ready to download:\n"
            "# if files is not None: files.download(str(workspace_zip))",
        ),
        _markdown(notebook_id, "next", "## 7. Continue\n\nUse `next_notebook_ids` in the result manifest. Keep the exported workspace when moving between notebooks or runtimes."),
    ]
    notebook = nbformat.v4.new_notebook(cells=cells)
    notebook.metadata = {
        "colab": {"name": entry["file"], "provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
        "hurdler": {"notebook_id": notebook_id, "backend": backend, "suite_version": "2026.1", "source_output_free": True},
    }
    return notebook


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--output-dir", type=Path, default=NOTEBOOK_DIR)
    args = parser.parse_args()
    payload = json.loads(args.catalog.read_text())
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    for entry in payload["notebooks"]:
        notebook = build_notebook(entry, git_ref=payload.get("default_git_ref", "main"))
        nbformat.write(notebook, output / entry["file"])
    print(json.dumps({"generated": len(payload["notebooks"]), "output_dir": str(output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
