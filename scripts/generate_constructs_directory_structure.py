#!/usr/bin/env python3
"""Regenerate the deposited SI file-tree listing."""

from __future__ import annotations

import argparse
from pathlib import Path

from hurdler.data_deposition import write_directory_structure


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "directory",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "constructs_and_sequencing_result",
    )
    args = parser.parse_args()
    print(write_directory_structure(args.directory))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
