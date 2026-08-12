"""Reference-data validation and versioned snapshots."""

from __future__ import annotations

import csv
from pathlib import Path

from .io import sha256_file, utc_now, write_json_atomic


REFERENCE_FILES = (
    "codon_usage.csv",
    "methylation_check.csv",
    "neb_buffer_activity_cleaned.csv",
    "orthogonality.csv",
    "plasmid_digest_check.csv",
    "restriction_enzyme_seamless_insert.csv",
    "restriction_enzyme_slient_mutation.csv",
)
REFERENCE_METADATA_FILES = ("codon_usage.metadata.json",)


def _csv_shape(path: Path) -> tuple[int, int]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = sum(1 for _ in reader)
    return rows, len(header)


def build_reference_manifest(reference_dir: str | Path, output_path: str | Path) -> dict[str, object]:
    root = Path(reference_dir)
    entries: list[dict[str, object]] = []
    for name in REFERENCE_FILES:
        path = root / name
        if not path.exists():
            raise FileNotFoundError(path)
        rows, columns = _csv_shape(path)
        entries.append(
            {
                "name": name,
                "path": str(path.absolute()),
                "bytes": path.stat().st_size,
                "rows": rows,
                "columns": columns,
                "sha256": sha256_file(path),
            }
        )
    for name in REFERENCE_METADATA_FILES:
        path = root / name
        if not path.exists():
            raise FileNotFoundError(path)
        entries.append(
            {
                "name": name,
                "path": str(path.absolute()),
                "bytes": path.stat().st_size,
                "rows": None,
                "columns": None,
                "sha256": sha256_file(path),
            }
        )
    manifest = {
        "schema_version": 1,
        "created_at": utc_now(),
        "reference_dir": str(root.absolute()),
        "files": entries,
    }
    write_json_atomic(manifest, output_path)
    return manifest
