"""Compatibility module for module-level HURDLER queries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .index import PatternIndex
from .matching import materialize_best_solution, query_all_plasmids
from .paths import ProjectPaths


def query_module(module: str, index_dir: str | Path) -> list[dict[str, object]]:
    index = PatternIndex.load(index_dir)
    return [materialize_best_solution(result, index) for result in query_all_plasmids(module, index)]


def main(argv: list[str] | None = None) -> int:
    root = ProjectPaths.discover()
    parser = argparse.ArgumentParser(description="Query a repeat module against the HURDLER index")
    parser.add_argument("--module", required=True)
    parser.add_argument("--index-dir", type=Path, default=root.output / "artifacts" / "legacy-optimized-v1")
    args = parser.parse_args(argv)
    print(json.dumps(query_module(args.module, args.index_dir), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
