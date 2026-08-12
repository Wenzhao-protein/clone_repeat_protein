#!/usr/bin/env python3
"""Execute one source notebook, render HTML, and validate both artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import nbformat
import papermill as pm
from nbconvert import HTMLExporter
from traitlets.config import Config


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
    parser.add_argument("--parameters-json", default="{}")
    parser.add_argument("--timeout", type=int, default=7200)
    args = parser.parse_args()
    parameters = json.loads(args.parameters_json)
    args.executed.parent.mkdir(parents=True, exist_ok=True)
    args.html.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    pm.execute_notebook(
        str(args.source),
        str(args.executed),
        parameters=parameters,
        kernel_name="hurdler",
        cwd=str(args.cwd),
        execution_timeout=args.timeout,
        progress_bar=False,
        log_output=True,
    )
    book = nbformat.read(args.executed, as_version=4)
    errors = [
        output
        for cell in book.cells
        if cell.cell_type == "code"
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    if errors:
        raise RuntimeError(f"Executed notebook contains {len(errors)} error outputs")
    config = Config()
    config.HTMLExporter.exclude_input_prompt = True
    config.HTMLExporter.exclude_output_prompt = True
    body, _resources = HTMLExporter(config=config).from_notebook_node(book)
    args.html.write_text(body)
    if args.executed.stat().st_size == 0 or args.html.stat().st_size == 0:
        raise RuntimeError("Notebook or HTML output is empty")
    payload = {
        "source": str(args.source.resolve()),
        "source_sha256": sha256(args.source),
        "executed": str(args.executed.resolve()),
        "executed_sha256": sha256(args.executed),
        "html": str(args.html.resolve()),
        "html_sha256": sha256(args.html),
        "parameters": parameters,
        "runtime_seconds": round(time.monotonic() - started, 3),
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "error_output_count": 0,
        "status": "passed",
    }
    manifest = args.executed.with_suffix(".manifest.json")
    manifest.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
