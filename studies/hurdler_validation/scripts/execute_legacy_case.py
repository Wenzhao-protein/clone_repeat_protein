#!/usr/bin/env python3
"""Run one historical notebook in an isolated scratch compatibility overlay."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import nbformat


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def materialize_compatibility_copy(
    case_id: str, source: Path, destination: Path
) -> tuple[Path, list[str]]:
    """Copy a notebook and apply only audited, scratch-local compatibility fixes."""
    book = nbformat.read(source, as_version=4)
    patches: list[str] = []
    if case_id == "archive_get_re_site":
        # The archive notebook contains this exact function body in a
        # triple-quoted cell before its first call.  Activate the same body in
        # the scratch copy without changing the archived source.
        book.cells.insert(
            0,
            nbformat.v4.new_code_cell(
                """import itertools
from Bio.Seq import Seq

def re_to_aa(seq):
    nucleotides = ['A', 'T', 'G', 'C']
    three_bp_sequences = [''.join(item) for item in itertools.product(nucleotides, repeat=3)]
    seq_rc = str(Seq(seq).reverse_complement())
    seq_ex = [seq, seq_rc]
    seq_ex += [item[0] + seq + item[1] + item[2] for item in three_bp_sequences]
    seq_ex += [item[0] + item[1] + seq + item[2] for item in three_bp_sequences]
    seq_ex += [item[0] + seq_rc + item[1] + item[2] for item in three_bp_sequences]
    seq_ex += [item[0] + item[1] + seq_rc + item[2] for item in three_bp_sequences]
    return sorted({str(Seq(item).translate()) for item in seq_ex if '*' not in str(Seq(item).translate())})
"""
            ),
        )
        patches.append("activated_existing_triple_quoted_re_to_aa_in_scratch_copy")
    destination.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(book, destination)
    return destination, patches


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--cwd-kind", choices=["root", "utils", "sec", "codon", "archive_get_re"], default="root")
    parser.add_argument("--timeout", type=int, default=6600)
    args = parser.parse_args()
    started = time.monotonic()
    overlay = args.scratch_root / args.case_id / "overlay"
    prepare = args.repo / "studies" / "hurdler_validation" / "scripts" / "prepare_legacy_overlay.py"
    subprocess.run(
        [
            "/home/wendai/.conda/envs/hurdler/bin/python",
            str(prepare),
            "--repo",
            str(args.repo),
            "--overlay",
            str(overlay),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    cwd = overlay if args.cwd_kind == "root" else overlay / args.cwd_kind
    case_artifacts = args.artifact_root / args.case_id
    executed = case_artifacts / f"{args.source.stem}_executed.ipynb"
    html = case_artifacts / f"{args.source.stem}.html"
    log = overlay.parent / "execution.log"
    execution_source, compatibility_patches = materialize_compatibility_copy(
        args.case_id,
        args.source,
        overlay.parent / f"{args.source.stem}_compatibility.ipynb",
    )
    command = [
        "/home/wendai/.conda/envs/hurdler/bin/python",
        str(args.repo / "studies" / "hurdler_validation" / "scripts" / "execute_notebook.py"),
        str(execution_source),
        str(executed),
        str(html),
        "--cwd",
        str(cwd),
        "--timeout",
        str(args.timeout),
    ]
    status = "failed"
    error = ""
    return_code: int | None = None
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=args.timeout + 60,
            check=False,
        )
        return_code = completed.returncode
        log.write_text(completed.stdout + "\n--- STDERR ---\n" + completed.stderr)
        if completed.returncode == 0:
            status = "passed"
        else:
            error = (completed.stderr or completed.stdout)[-4000:]
    except subprocess.TimeoutExpired as exc:
        status = "deferred_long"
        error = f"Historical notebook exceeded the {args.timeout}-second per-case profiling limit"
        log.write_text((exc.stdout or "") + "\n--- STDERR ---\n" + (exc.stderr or ""))

    payload = {
        "workflow": args.case_id,
        "workflow_type": "legacy_notebook",
        "source": str(args.source.resolve()),
        "source_sha256": sha256(args.source),
        "executed_source": str(execution_source),
        "compatibility_patches": compatibility_patches,
        "status": status,
        "return_code": return_code,
        "runtime_seconds": round(time.monotonic() - started, 3),
        "output_path": str(case_artifacts.resolve()),
        "scratch_overlay": str(overlay),
        "log": str(log),
        "error": error,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "rerun_command": " ".join(command),
    }
    case_artifacts.mkdir(parents=True, exist_ok=True)
    (case_artifacts / "status.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    # Capturing a notebook failure is a successful audit task. The case-level
    # status remains `failed` and is never inferred from Taskrunner completion.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
