#!/usr/bin/env python3
"""Merge and validate independently recoverable module-optimization shards."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_frames(shard_dirs: list[Path], name: str) -> list[pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for directory in shard_dirs:
        path = directory / name
        if not path.is_file():
            raise FileNotFoundError(f"Missing shard output: {path}")
        frame = pd.read_parquet(path)
        if not frame.empty:
            frame["source_shard"] = directory.name
            frames.append(frame)
    return frames


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-root", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--table-dir", type=Path, required=True)
    parser.add_argument("--candidate-output", type=Path, required=True)
    parser.add_argument("--expected-shards", type=int, required=True)
    args = parser.parse_args()

    shard_dirs = [args.shard_root / f"shard_{index:03d}" for index in range(args.expected_shards)]
    catalog = pd.read_parquet(args.catalog)
    results = pd.concat(read_frames(shard_dirs, "module_hurdler_results.parquet"), ignore_index=True)
    constructs = pd.concat(read_frames(shard_dirs, "optimized_constructs.parquet"), ignore_index=True)
    candidate_frames = read_frames(shard_dirs, "module_hurdler_candidates.parquet")
    candidates = pd.concat(candidate_frames, ignore_index=True) if candidate_frames else pd.DataFrame()
    if "failure_reason" not in constructs:
        constructs["failure_reason"] = ""
    constructs["failure_reason"] = constructs["failure_reason"].fillna("")
    no_solution = constructs.optimization_status.eq("no_hurdler_solution")
    constructs.loc[no_solution, "failure_reason"] = (
        "No compatible legacy-optimized-v1 HURDLER solution across the eight plasmids"
    )

    expected_modules = set(catalog.module_id)
    result_modules = set(results.module_id)
    construct_modules = set(constructs.module_id)
    duplicate_result_rows = int(results.duplicated(["module_id", "plasmid"]).sum())
    duplicate_construct_rows = int(constructs.duplicated(["module_id", "fragment_limit_bp"]).sum())
    validation = {
        "expected_shards": args.expected_shards,
        "observed_shards": len(shard_dirs),
        "catalog_modules": len(catalog),
        "result_modules": len(result_modules),
        "construct_modules": len(construct_modules),
        "result_rows": len(results),
        "construct_rows": len(constructs),
        "candidate_rows": len(candidates),
        "missing_result_modules": sorted(expected_modules - result_modules),
        "missing_construct_modules": sorted(expected_modules - construct_modules),
        "unexpected_result_modules": sorted(result_modules - expected_modules),
        "unexpected_construct_modules": sorted(construct_modules - expected_modules),
        "duplicate_result_rows": duplicate_result_rows,
        "duplicate_construct_rows": duplicate_construct_rows,
    }
    validation["passed"] = bool(
        validation["observed_shards"] == validation["expected_shards"]
        and len(results) == len(catalog) * 8
        and len(constructs) == len(catalog) * 2
        and not validation["missing_result_modules"]
        and not validation["missing_construct_modules"]
        and not validation["unexpected_result_modules"]
        and not validation["unexpected_construct_modules"]
        and duplicate_result_rows == 0
        and duplicate_construct_rows == 0
    )
    if not validation["passed"]:
        raise RuntimeError(json.dumps(validation, indent=2))

    results = results.sort_values(["collection", "module_id", "plasmid"]).reset_index(drop=True)
    constructs = constructs.sort_values(["collection", "module_id", "fragment_limit_bp"]).reset_index(drop=True)
    if not candidates.empty:
        candidates = candidates.sort_values(["collection", "module_id", "plasmid", "candidate_pair_id"]).reset_index(drop=True)
        join_keys = [
            "module_id", "plasmid", "direction", "site_i_position", "site_ii_position",
            "site_i_enzyme", "site_ii_enzyme",
        ]
        selected_columns = [
            *join_keys,
            "site_i_9mer_bp", "site_ii_9mer_bp_mutated",
            "site_i_recognition_site", "site_ii_recognition_site", "site_iii_sites",
        ]
        solution_metadata = candidates[selected_columns].drop_duplicates(join_keys, keep="first")
        constructs = constructs.merge(solution_metadata, on=join_keys, how="left", validate="many_to_one")

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.candidate_output.parent.mkdir(parents=True, exist_ok=True)
    results.to_parquet(args.table_dir / "module_hurdler_results.parquet", index=False)
    results.to_csv(args.table_dir / "module_hurdler_results.csv", index=False)
    constructs.to_parquet(args.table_dir / "optimized_constructs.parquet", index=False)
    constructs.to_csv(args.table_dir / "optimized_constructs.csv", index=False)
    candidates.to_parquet(args.candidate_output, index=False)

    fasta_path = args.table_dir / "optimized_constructs.fasta"
    with fasta_path.open("w") as output:
        for row in constructs.itertuples(index=False):
            dna = getattr(row, "dna_sequence", None)
            copies = getattr(row, "verified_max_copies", 0)
            if not isinstance(dna, str) or not dna:
                continue
            output.write(f">{row.module_id}|cap={row.fragment_limit_bp}|copies={copies}\n")
            for start in range(0, len(dna), 80):
                output.write(dna[start : start + 80] + "\n")

    validation["catalog_sha256"] = sha256(args.catalog)
    validation["result_sha256"] = sha256(args.table_dir / "module_hurdler_results.parquet")
    validation["construct_sha256"] = sha256(args.table_dir / "optimized_constructs.parquet")
    validation["candidate_sha256"] = sha256(args.candidate_output)
    (args.table_dir / "optimization_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
