#!/usr/bin/env python3
"""Run independent refine-ga shards concurrently inside one Digs task."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import subprocess
import time
from pathlib import Path


def run_command(command: list[str]) -> dict[str, object]:
    started = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    return {
        "command": command,
        "returncode": result.returncode,
        "runtime_seconds": round(time.monotonic() - started, 3),
        "stdout": result.stdout[-4000:],
        "stderr": result.stderr[-4000:],
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hurdler", type=Path, required=True)
    parser.add_argument("--constructs", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--codon-usage", type=Path, required=True)
    parser.add_argument("--restriction-sites", type=Path, required=True)
    parser.add_argument("--group-index", type=int, required=True)
    parser.add_argument("--group-count", type=int, required=True)
    parser.add_argument("--total-shards", type=int, required=True)
    parser.add_argument("--shards-per-group", type=int, default=4)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--population-size", type=int, default=64)
    parser.add_argument("--short-generations", type=int, default=10)
    parser.add_argument(
        "--generation-schedule", type=int, nargs="+", default=[10, 20, 40, 60, 80, 100]
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--use-idt", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.group_index < args.group_count:
        raise ValueError("group-index must satisfy 0 <= index < group-count")
    if args.max_workers < 1 or args.shards_per_group < 1:
        raise ValueError("max-workers and shards-per-group must be positive")
    shard_indices = list(
        range(args.group_index, args.total_shards, args.group_count)
    )[: args.shards_per_group]
    if not shard_indices:
        raise ValueError("group does not contain a shard")

    commands: list[list[str]] = []
    for shard_index in shard_indices:
        command = [
            str(args.hurdler),
            "refine-ga",
            "--constructs",
            str(args.constructs),
            "--output-dir",
            str(args.output_root / f"shard_{shard_index:03d}"),
            "--codon-usage",
            str(args.codon_usage),
            "--restriction-sites",
            str(args.restriction_sites),
            "--shard-index",
            str(shard_index),
            "--shard-count",
            str(args.total_shards),
            "--population-size",
            str(args.population_size),
            "--adaptive-copy-search",
            "--short-generations",
            str(args.short_generations),
            "--generation-schedule",
            *(str(value) for value in args.generation_schedule),
            "--seed",
            str(args.seed),
        ]
        if args.use_idt:
            command.append("--use-idt")
        commands.append(command)

    started = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        results = list(pool.map(run_command, commands))
    runtime = round(time.monotonic() - started, 3)
    failures = [result for result in results if result["returncode"] != 0]
    output_hashes = {}
    for shard_index in shard_indices:
        output = args.output_root / f"shard_{shard_index:03d}" / "optimized_constructs_ga.parquet"
        if output.is_file():
            output_hashes[str(shard_index)] = sha256(output)
    payload = {
        "group_index": args.group_index,
        "group_count": args.group_count,
        "total_shards": args.total_shards,
        "shard_indices": shard_indices,
        "max_workers": args.max_workers,
        "runtime_seconds": runtime,
        "failed_commands": len(failures),
        "output_hashes": output_hashes,
        "results": results,
        "passed": not failures and len(output_hashes) == len(shard_indices),
    }
    destination = args.output_root / "groups"
    destination.mkdir(parents=True, exist_ok=True)
    validation = destination / f"group_{args.group_index:03d}.json"
    validation.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: payload[key] for key in payload if key != "results"}, indent=2))
    if not payload["passed"]:
        raise RuntimeError(f"Parallel refine group failed; see {validation}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
