#!/usr/bin/env python3
"""Convert a PDF to layout-preserving text with an auditable manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    executable = Path(sys.executable).with_name("pdftotext")
    if not executable.is_file():
        raise FileNotFoundError(f"pdftotext is not installed beside the active Python: {executable}")
    command = [str(executable), "-layout", str(args.source), str(args.destination)]
    subprocess.run(command, check=True)
    payload = {
        "source": str(args.source.resolve()),
        "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "destination": str(args.destination.resolve()),
        "destination_sha256": hashlib.sha256(args.destination.read_bytes()).hexdigest(),
        "destination_bytes": args.destination.stat().st_size,
        "converted_at": datetime.now(timezone.utc).isoformat(),
        "command": command,
    }
    manifest = args.destination.with_suffix(args.destination.suffix + ".manifest.json")
    manifest.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
