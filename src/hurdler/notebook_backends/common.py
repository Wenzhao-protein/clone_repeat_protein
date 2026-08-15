"""Shared implementation helpers; scientific algorithms stay in hurdler.*."""

from __future__ import annotations

import json
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

from ..io import sha256_file
from ..notebook_workspace import (
    NotebookContext,
    NotebookResult,
    artifact_record,
    write_run_manifest,
)
from ..paths import ProjectPaths


@dataclass(frozen=True)
class BackendSpec:
    notebook_id: str
    title: str
    purpose: str
    production_workflows: tuple[str, ...] = ()
    default_request: Mapping[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "notebook_id": self.notebook_id,
            "title": self.title,
            "purpose": self.purpose,
            "production_workflows": list(self.production_workflows),
            "default_request": dict(self.default_request or {}),
        }


def repo_root() -> Path:
    return ProjectPaths.discover().root


def resolve_path(value: str | Path | None, default: str | Path) -> Path:
    return Path(value or default).expanduser().absolute()


def copy_file(source: str | Path, destination: str | Path) -> Path:
    src, dst = Path(source), Path(destination)
    if not src.is_file():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def write_frame(frame: pd.DataFrame, stem: Path) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    parquet = stem.with_suffix(".parquet")
    csv = stem.with_suffix(".csv")
    frame.to_parquet(parquet, index=False)
    frame.to_csv(csv, index=False)
    return [parquet, csv]


def zip_paths(output: Path, files: Iterable[Path], *, base: Path | None = None) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    paths = sorted({Path(path).absolute() for path in files if Path(path).is_file()})
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            if base is not None:
                try:
                    name = path.relative_to(base.absolute())
                except ValueError:
                    name = Path(path.name)
            else:
                name = Path(path.name)
            archive.write(path, name)
    return output


def result_from_paths(
    context: NotebookContext,
    *,
    backend_id: str,
    request: Mapping[str, Any],
    paths: Iterable[Path],
    metrics: Mapping[str, Any],
    next_notebooks: Iterable[str] = (),
    status: str = "passed",
    warnings: Iterable[str] = (),
    limitations: Iterable[str] = (),
) -> NotebookResult:
    records = [
        artifact_record(path, artifact_id=f"{context.run_id}:{path.name}", role=backend_id)
        for path in paths
    ]
    result = NotebookResult(
        status=status,
        metrics=dict(metrics),
        artifacts=records,
        warnings=list(warnings),
        limitations=list(limitations),
        next_notebook_ids=list(next_notebooks),
    )
    write_run_manifest(context, backend_id=backend_id, request=dict(request), result=result)
    return result


def load_table(path: str | Path, **kwargs: Any) -> pd.DataFrame:
    source = Path(path)
    if source.suffix.lower() == ".parquet":
        return pd.read_parquet(source, **kwargs)
    if source.suffix.lower() in {".csv", ".tsv"}:
        return pd.read_csv(source, sep="\t" if source.suffix.lower() == ".tsv" else ",", **kwargs)
    raise ValueError(f"Unsupported table format: {source}")


def preflight_files(paths: Iterable[Path]) -> dict[str, Any]:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + ", ".join(missing))
    return {
        "status": "passed",
        "inputs": [
            {
                "path": str(path.absolute()),
                "sha256": sha256_file(path) if path.is_file() else None,
                "size_bytes": path.stat().st_size if path.is_file() else None,
            }
            for path in paths
        ],
    }


def json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, default=str)
