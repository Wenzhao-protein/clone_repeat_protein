"""Portable, credential-free workspaces shared by the V2 notebooks.

The notebooks intentionally contain presentation code only.  This module owns
the durable hand-off contract so a tutorial run, a Colab run and a returned
production bundle all look the same to downstream notebooks.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping

from .io import sha256_file, utc_now, write_json_atomic


WORKSPACE_SCHEMA_VERSION = "hurdler-notebook-workspace-v2"
VALID_MODES = {"tutorial", "colab_full", "production_bundle", "analyze"}
VALID_SOURCE_MODES = {"snapshot", "refresh", "upload"}
WORKSPACE_DIRECTORIES = ("tables", "figures", "fasta", "genbank", "reports")
_SECRET_KEY = re.compile(
    r"(password|secret|access[_-]?token|authorization|credential_contents)", re.I
)


def _repo_commit(start: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(start), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _assert_secret_free(value: Any, trail: tuple[str, ...] = ()) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if _SECRET_KEY.search(str(key)) and item not in (None, "", False):
                raise ValueError(
                    "Secret-bearing values may not be written to a notebook workspace: "
                    + ".".join((*trail, str(key)))
                )
            _assert_secret_free(item, (*trail, str(key)))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _assert_secret_free(item, (*trail, str(index)))


@dataclass(frozen=True)
class NotebookContext:
    """Execution context shared by every authoritative V2 backend."""

    run_id: str
    mode: str = "tutorial"
    source_mode: str = "snapshot"
    workspace_root: Path = Path("/content/hurdler_workspace")
    repo_commit: str = ""
    random_seed: int = 42
    artifact_registry: Path | None = None

    def __post_init__(self) -> None:
        if self.mode not in VALID_MODES:
            raise ValueError(f"Unsupported notebook mode: {self.mode}")
        if self.source_mode not in VALID_SOURCE_MODES:
            raise ValueError(f"Unsupported source mode: {self.source_mode}")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", self.run_id):
            raise ValueError("run_id may contain only letters, digits, dot, dash and underscore")
        root = Path(self.workspace_root).expanduser().absolute()
        object.__setattr__(self, "workspace_root", root)
        if not self.repo_commit:
            object.__setattr__(self, "repo_commit", _repo_commit(Path(__file__).parents[2]))
        if self.artifact_registry is not None:
            object.__setattr__(
                self, "artifact_registry", Path(self.artifact_registry).expanduser().absolute()
            )

    @property
    def run_root(self) -> Path:
        return self.workspace_root / self.run_id

    def prepare(self) -> "NotebookContext":
        self.run_root.mkdir(parents=True, exist_ok=True)
        for name in WORKSPACE_DIRECTORIES:
            (self.run_root / name).mkdir(parents=True, exist_ok=True)
        return self

    def directory(self, name: str) -> Path:
        if name not in WORKSPACE_DIRECTORIES:
            raise ValueError(f"Unknown workspace directory: {name}")
        return self.run_root / name


@dataclass
class NotebookResult:
    """Uniform result returned by all V2 notebook backends."""

    status: str
    metrics: dict[str, Any] = field(default_factory=dict)
    artifacts: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    next_notebook_ids: list[str] = field(default_factory=list)
    run_manifest_path: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        _assert_secret_free(payload)
        return payload


def artifact_record(path: str | Path, *, artifact_id: str, role: str) -> dict[str, Any]:
    source = Path(path).absolute()
    if not source.is_file():
        raise FileNotFoundError(source)
    return {
        "artifact_id": artifact_id,
        "role": role,
        "path": str(source),
        "sha256": sha256_file(source),
        "size_bytes": source.stat().st_size,
    }


def write_run_manifest(
    context: NotebookContext,
    *,
    backend_id: str,
    request: Mapping[str, Any],
    result: NotebookResult,
) -> Path:
    """Write the canonical per-notebook manifest after rejecting secrets."""
    context.prepare()
    request_payload = dict(request)
    _assert_secret_free(request_payload)
    payload = {
        "schema_version": WORKSPACE_SCHEMA_VERSION,
        "created_at": utc_now(),
        "backend_id": backend_id,
        "context": {
            "run_id": context.run_id,
            "mode": context.mode,
            "source_mode": context.source_mode,
            "repo_commit": context.repo_commit,
            "random_seed": context.random_seed,
        },
        "request": request_payload,
        "result": result.to_dict(),
    }
    destination = context.run_root / "workspace_manifest.json"
    write_json_atomic(payload, destination)
    result.run_manifest_path = str(destination)
    return destination


def export_workspace(context: NotebookContext, output: str | Path | None = None) -> Path:
    """Create a deterministic-path, credential-free ZIP of one run."""
    context.prepare()
    manifest = context.run_root / "workspace_manifest.json"
    if not manifest.exists():
        raise FileNotFoundError("workspace_manifest.json is required before export")
    _assert_secret_free(json.loads(manifest.read_text()))
    destination = Path(output or context.workspace_root / f"hurdler_workspace_{context.run_id}.zip")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".zip", delete=False, dir=destination.parent) as handle:
        temporary = Path(handle.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(context.run_root.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(context.run_root.parent))
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def import_workspace(archive: str | Path, destination: str | Path) -> NotebookContext:
    """Import a workspace ZIP after preventing traversal and checking its manifest."""
    source = Path(archive)
    root = Path(destination).expanduser().absolute()
    with zipfile.ZipFile(source) as bundle:
        members = bundle.namelist()
        for member in members:
            candidate = (root / member).resolve()
            if root not in candidate.parents and candidate != root:
                raise ValueError(f"Unsafe workspace member: {member}")
        bundle.extractall(root)
    manifests = sorted(root.glob("*/workspace_manifest.json"))
    if len(manifests) != 1:
        raise ValueError("Workspace ZIP must contain exactly one workspace_manifest.json")
    payload = json.loads(manifests[0].read_text())
    _assert_secret_free(payload)
    if payload.get("schema_version") != WORKSPACE_SCHEMA_VERSION:
        raise ValueError("Unsupported workspace schema")
    data = payload["context"]
    return NotebookContext(
        run_id=data["run_id"],
        mode=data["mode"],
        source_mode=data["source_mode"],
        workspace_root=root,
        repo_commit=data["repo_commit"],
        random_seed=int(data["random_seed"]),
    )


ProgressCallback = Callable[[dict[str, Any]], None]
