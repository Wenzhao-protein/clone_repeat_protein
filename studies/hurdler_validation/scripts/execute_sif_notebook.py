#!/usr/bin/env python3
"""Execute a notebook with the shared Papermill-SIF helper and export HTML."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import nbformat
from nbconvert import HTMLExporter
from traitlets.config import Config


HELPER = Path(
    "/home/wendai/projects/skills/papermill-sif-notebook/scripts/run_papermill_sif.sh"
)
SIF = Path("/net/software/containers/universal.sif")


def parse_parameter_value(value: str):
    lowered = value.strip().lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered in {"none", "null"}:
        return None
    try:
        return ast.literal_eval(value)
    except (SyntaxError, ValueError):
        return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("executed", type=Path)
    parser.add_argument("html", type=Path)
    parser.add_argument("--cwd", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument(
        "--parameter",
        nargs=2,
        action="append",
        default=[],
        metavar=("NAME", "VALUE"),
        help="Papermill parameter override; may be repeated",
    )
    args = parser.parse_args()
    args.executed.parent.mkdir(parents=True, exist_ok=True)
    args.html.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(
        PAPERMILL_SIF=str(SIF),
        PAPERMILL_BIND=f"{args.cwd}:{args.cwd},/net/scratch:/net/scratch,/home:/home",
        PAPERMILL_IN_CONTAINER_CWD=str(args.cwd),
        PAPERMILL_TIMEOUT=str(args.timeout),
    )
    started = time.monotonic()
    execution_source = args.source
    temporary_source: Path | None = None
    if args.parameter:
        parameterized = nbformat.read(args.source, as_version=4)
        parameter_cell_index = next(
            (
                index
                for index, cell in enumerate(parameterized.cells)
                if "parameters" in cell.get("metadata", {}).get("tags", [])
            ),
            None,
        )
        if parameter_cell_index is None:
            raise RuntimeError("Notebook has no cell tagged 'parameters'")
        assignments = "\n".join(
            f"{name} = {parse_parameter_value(value)!r}"
            for name, value in args.parameter
        )
        injected = nbformat.v4.new_code_cell(
            assignments,
            metadata={"tags": ["injected-parameters"]},
        )
        parameterized.cells.insert(parameter_cell_index + 1, injected)
        temporary_source = args.executed.with_name(
            f".{args.executed.stem}.parameterized-source.ipynb"
        )
        nbformat.write(parameterized, temporary_source)
        execution_source = temporary_source
    try:
        subprocess.run(
            [str(HELPER), str(execution_source), str(args.executed)],
            check=True,
            cwd=args.cwd,
            env=environment,
            timeout=args.timeout + 300,
        )
    finally:
        if temporary_source is not None:
            temporary_source.unlink(missing_ok=True)
    notebook = nbformat.read(args.executed, as_version=4)
    errors = [
        output
        for cell in notebook.cells
        if cell.cell_type == "code"
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    if errors:
        raise RuntimeError(f"Executed notebook contains {len(errors)} errors")
    config = Config()
    config.HTMLExporter.exclude_input_prompt = True
    config.HTMLExporter.exclude_output_prompt = True
    body, _ = HTMLExporter(config=config).from_notebook_node(notebook)
    args.html.write_text(body)
    payload = {
        "source": str(args.source.resolve()),
        "source_sha256": sha256(args.source),
        "executed": str(args.executed.resolve()),
        "executed_sha256": sha256(args.executed),
        "html": str(args.html.resolve()),
        "html_sha256": sha256(args.html),
        "sif": str(SIF),
        "runtime_seconds": round(time.monotonic() - started, 3),
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "error_output_count": 0,
        "parameters": {name: value for name, value in args.parameter},
        "status": "passed",
    }
    args.executed.with_suffix(".manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
