#!/usr/bin/env python3
"""Create scratch-only compatibility trees for historical relative paths."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


def link(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() != source.resolve():
            raise RuntimeError(f"Unexpected existing symlink: {destination}")
        return
    if destination.exists():
        raise RuntimeError(f"Refusing to replace existing path: {destination}")
    destination.symlink_to(source, target_is_directory=source.is_dir())


def copy_once(source: Path, destination: Path) -> None:
    if destination.exists():
        return
    shutil.copytree(source, destination)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--overlay", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    root = args.overlay.resolve()
    root.mkdir(parents=True, exist_ok=True)

    # Historical inputs are exposed read-only through symlinks; paths that old
    # notebooks overwrite are copied into scratch first.
    (root / "utils").mkdir(exist_ok=True)
    link(repo / "data" / "reference_input", root / "utils" / "input")
    copy_once(repo / "data" / "reference_output", root / "utils" / "output")
    copy_once(repo / "output", root / "output")
    (root / "hurdler_analysis").mkdir(exist_ok=True)
    link(repo / "data" / "hurdler_analysis_input", root / "hurdler_analysis" / "input")

    sec = root / "sec"
    sec.mkdir(exist_ok=True)
    (sec / "input").mkdir(exist_ok=True)
    sec_mat = repo / "SEC" / "input" / "20250530_sec_result.mat"
    link(sec_mat, sec / "input" / "20250530_sec_result.mat")
    # The canonical notebook retained the pre-restructure `_total` filename;
    # both names point to the same tracked chromatogram without copying data.
    link(sec_mat, sec / "input" / "20250530_sec_result_total.mat")
    link(repo / "SEC" / "src", sec / "sec_utils")
    (sec / "output").mkdir(exist_ok=True)
    (sec / "temp").mkdir(exist_ok=True)

    codon = root / "codon"
    codon.mkdir(exist_ok=True)
    link(
        repo / "codon_opt_benchmark_extended" / "results" / "well_color_iterations.csv",
        codon / "well_color_iterations.csv",
    )
    archive_get_re = root / "archive_get_re"
    archive_get_re.mkdir(exist_ok=True)
    (archive_get_re / "output").mkdir(exist_ok=True)

    manifest = {
        "repo": str(repo),
        "overlay": str(root),
        "read_only_links": [
            str(root / "utils" / "input"),
            str(root / "hurdler_analysis" / "input"),
            str(sec / "input" / "20250530_sec_result.mat"),
            str(sec / "input" / "20250530_sec_result_total.mat"),
            str(sec / "sec_utils"),
            str(codon / "well_color_iterations.csv"),
        ],
        "scratch_copies": [str(root / "utils" / "output"), str(root / "output")],
    }
    (root / "overlay_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
