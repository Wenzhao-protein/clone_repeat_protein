"""Versioned artifact registry used by the V2 notebook suite."""

from __future__ import annotations

import json
import shutil
import urllib.parse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import requests

from .io import sha256_file
from .paths import ProjectPaths


ARTIFACT_REGISTRY_SCHEMA = "hurdler-artifact-registry-v2"
VALID_LEVELS = {"fixture", "snapshot", "compact_result", "production_raw"}


@dataclass(frozen=True)
class ArtifactRecord:
    artifact_id: str
    version: str
    role: str
    schema_version: str
    level: str
    sha256: str
    size_bytes: int
    media_type: str
    repo_path: str = ""
    release_url: str = ""
    generated_by: str = ""
    source: str = ""
    download_date: str = ""
    citation: str = ""
    license: str = ""
    row_count: int | None = None

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "ArtifactRecord":
        record = cls(**payload)
        if record.level not in VALID_LEVELS:
            raise ValueError(f"Invalid artifact level for {record.artifact_id}: {record.level}")
        if bool(record.repo_path) == bool(record.release_url):
            raise ValueError(
                f"Artifact {record.artifact_id} must define exactly one of repo_path/release_url"
            )
        if len(record.sha256) != 64:
            raise ValueError(f"Artifact {record.artifact_id} has an invalid SHA256")
        return record


class ArtifactRegistry:
    def __init__(self, path: str | Path | None = None, *, repo_root: str | Path | None = None):
        paths = ProjectPaths.discover(repo_root)
        self.repo_root = paths.root
        self.path = Path(path or paths.root / "data" / "artifact_registry_v2.json")
        payload = json.loads(self.path.read_text())
        if payload.get("schema_version") != ARTIFACT_REGISTRY_SCHEMA:
            raise ValueError(f"Unsupported artifact registry: {self.path}")
        records = [ArtifactRecord.from_dict(item) for item in payload.get("artifacts", [])]
        if len({record.artifact_id for record in records}) != len(records):
            raise ValueError("Artifact IDs must be unique")
        self._records = {record.artifact_id: record for record in records}

    def list(self, *, level: str | None = None) -> list[ArtifactRecord]:
        records = sorted(self._records.values(), key=lambda item: item.artifact_id)
        return [record for record in records if level is None or record.level == level]

    def get(self, artifact_id: str) -> ArtifactRecord:
        try:
            return self._records[artifact_id]
        except KeyError as exc:
            raise KeyError(f"Unknown artifact ID: {artifact_id}") from exc

    def resolve_local(self, artifact_id: str) -> Path:
        record = self.get(artifact_id)
        if not record.repo_path:
            raise FileNotFoundError(f"Artifact {artifact_id} is release-backed and not bundled")
        return (self.repo_root / record.repo_path).resolve()

    def verify(self, artifact_id: str, path: str | Path | None = None) -> Path:
        record = self.get(artifact_id)
        candidate = Path(path) if path is not None else self.resolve_local(artifact_id)
        if not candidate.is_file():
            raise FileNotFoundError(candidate)
        observed_size = candidate.stat().st_size
        if observed_size != record.size_bytes:
            raise ValueError(
                f"Artifact {artifact_id} size mismatch: {observed_size} != {record.size_bytes}"
            )
        observed = sha256_file(candidate)
        if observed != record.sha256:
            raise ValueError(f"Artifact {artifact_id} checksum mismatch")
        if record.row_count is not None:
            suffix = candidate.suffix.lower()
            if suffix == ".parquet":
                rows = len(pd.read_parquet(candidate))
            elif suffix == ".csv":
                rows = len(pd.read_csv(candidate))
            else:
                rows = record.row_count
            if rows != record.row_count:
                raise ValueError(f"Artifact {artifact_id} row-count mismatch")
        return candidate

    def fetch(
        self,
        artifact_id: str,
        destination: str | Path | None = None,
        *,
        allow_production_raw: bool = False,
        timeout: int = 120,
    ) -> Path:
        record = self.get(artifact_id)
        if record.level == "production_raw" and not allow_production_raw:
            raise PermissionError("Production-raw artifacts require explicit opt-in")
        if record.repo_path:
            source = self.verify(artifact_id)
            if destination is None:
                return source
            target = Path(destination)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            return self.verify(artifact_id, target)
        name = Path(urllib.parse.urlparse(record.release_url).path).name or artifact_id
        target = Path(destination or Path.home() / ".cache" / "hurdler" / "artifacts" / name)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            return self.verify(artifact_id, target)
        temporary = target.with_suffix(target.suffix + ".part")
        with requests.get(record.release_url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        temporary.replace(target)
        return self.verify(artifact_id, target)


def registry_rows(records: Iterable[ArtifactRecord]) -> list[dict[str, Any]]:
    return [record.__dict__.copy() for record in records]
