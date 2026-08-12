#!/usr/bin/env python3
"""Safely extract supplementary archives and record every file hash."""

from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path


def safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            target = (destination / member.filename).resolve()
            if destination.resolve() not in target.parents and target != destination.resolve():
                raise ValueError(f"Unsafe archive member: {member.filename}")
        handle.extractall(destination)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    safe_extract(args.archive, args.destination)
    for nested in sorted(args.destination.rglob("*.zip")):
        nested_destination = nested.with_suffix("")
        safe_extract(nested, nested_destination)
    rows = [
        {"path": str(path.relative_to(args.destination)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(args.destination.rglob("*"))
        if path.is_file()
    ]
    (args.destination / "source_manifest.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps({"files": len(rows), "destination": str(args.destination)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
