#!/usr/bin/env python3
"""Validate and compact the completed run101/run103 production artifacts.

The command is deliberately narrow: it only touches the active
expanded-middle-repeatsdb-foldseek-v1 tables and the two completed Digs raw
directories named below.  Run without ``--apply`` for a read-only audit.
"""

from __future__ import annotations

import argparse
import getpass
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import BinaryIO

import duckdb
import pyarrow.parquet as pq

from hurdler.optimization import translate_dna


CORPUS = "expanded-middle-repeatsdb-foldseek-v1"
EXPECTED = {
    "modules": 26_095,
    "compatible": 25_168,
    "module_capacity_results": 50_336,
    "trace_rows": 806_726,
    "accepted": 48_163,
    "idt_audits": 546_608,
    "no_two_copy_capacity": 2,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def path_inventory(path: Path) -> dict[str, int]:
    if path.is_file():
        return {"files": 1, "logical_bytes": path.stat().st_size}
    files = 0
    logical_bytes = 0
    for base, _directories, names in os.walk(path):
        for name in names:
            target = Path(base) / name
            if target.is_file():
                files += 1
                logical_bytes += target.stat().st_size
    return {"files": files, "logical_bytes": logical_bytes}


def display_path(path: Path, repo: Path, scratch_repo: Path) -> str:
    for root, prefix in ((repo, "$REPO"), (scratch_repo, "$SCRATCH")):
        try:
            return f"{prefix}/{path.relative_to(root)}"
        except ValueError:
            continue
    return str(path)


def scalar(connection: duckdb.DuckDBPyConnection, sql: str) -> int:
    return int(connection.execute(sql).fetchone()[0])


def validate_tables(table_dir: Path, public_csv: Path) -> dict[str, object]:
    compatibility = table_dir / "module_compatibility.parquet"
    results = table_dir / "maximum_copy_results.parquet"
    traces = table_dir / "adaptive_copy_search_trace.parquet"
    validation = table_dir / "adaptive_copy_validation.parquet"
    summary = table_dir / "module_final_summary.parquet"
    for required in (compatibility, results, traces, validation, summary, public_csv):
        if not required.is_file():
            raise FileNotFoundError(required)

    con = duckdb.connect()
    q = lambda path: str(path).replace("'", "''")
    compatibility_rows = scalar(
        con, f"SELECT count(*) FROM read_parquet('{q(compatibility)}')"
    )
    compatible_rows = scalar(
        con,
        f"SELECT count(*) FROM read_parquet('{q(compatibility)}') "
        "WHERE hurdler_compatible",
    )
    result_rows = scalar(con, f"SELECT count(*) FROM read_parquet('{q(results)}')")
    unique_result_rows = scalar(
        con,
        f"SELECT count(*) FROM (SELECT DISTINCT collection, module_id, "
        f"fragment_limit_bp FROM read_parquet('{q(results)}'))",
    )
    trace_rows = scalar(con, f"SELECT count(*) FROM read_parquet('{q(traces)}')")
    accepted_rows = scalar(
        con,
        f"SELECT count(*) FROM read_parquet('{q(validation)}') "
        "WHERE validation_passed",
    )
    summary_rows = scalar(con, f"SELECT count(*) FROM read_parquet('{q(summary)}')")
    public_rows = scalar(
        con,
        f"SELECT count(*) FROM read_csv_auto('{q(public_csv)}', header=true, "
        "all_varchar=true, sample_size=-1)",
    )
    no_capacity = scalar(
        con,
        f"SELECT count(*) FROM read_parquet('{q(results)}') "
        "WHERE fragment_limit_bp=1800 AND mathematical_max_copies=1 "
        "AND stage2_preparation_status='no_two_copy_capacity'",
    )
    invalid_accepted = scalar(
        con,
        f"SELECT count(*) FROM read_parquet('{q(results)}') "
        "WHERE validation_passed AND (verified_max_copies < 2 "
        "OR idt_complexity_score >= 10 OR idt_complexity_score IS NULL "
        "OR selected_pair_re_site_excess <> 0 "
        "OR dna_length > fragment_limit_bp - external_deduction_bp "
        "OR idt_score_policy <> 'idt-rule-score-sum-lt10-v1')",
    )
    invalid_capacity_proof = scalar(
        con,
        f"SELECT count(*) FROM read_parquet('{q(results)}') "
        "WHERE validation_passed "
        "AND adaptive_maximum_proof_status='capacity_limit_reached' "
        "AND verified_max_copies <> adaptive_search_upper_bound_copies",
    )
    missing_next_copy_proof = scalar(
        con,
        f"SELECT count(*) FROM read_parquet('{q(results)}') r "
        "WHERE r.validation_passed "
        "AND r.adaptive_maximum_proof_status='next_copy_failed_at_100' "
        "AND NOT EXISTS (SELECT 1 FROM "
        f"read_parquet('{q(traces)}') t WHERE t.collection=r.collection "
        "AND t.module_id=r.module_id "
        "AND t.fragment_limit_bp=r.fragment_limit_bp "
        "AND t.copies=r.verified_max_copies+1 "
        "AND t.generations=100 AND t.passed=false)",
    )
    observed = {
        "modules": compatibility_rows,
        "compatible": compatible_rows,
        "module_capacity_results": result_rows,
        "unique_module_capacity_results": unique_result_rows,
        "trace_rows": trace_rows,
        "accepted": accepted_rows,
        "summary_rows": summary_rows,
        "public_csv_rows": public_rows,
        "no_two_copy_capacity": no_capacity,
        "invalid_accepted_rows": invalid_accepted,
        "invalid_capacity_proofs": invalid_capacity_proof,
        "missing_next_copy_100_generation_proofs": missing_next_copy_proof,
    }
    required_values = {
        "modules": EXPECTED["modules"],
        "compatible": EXPECTED["compatible"],
        "module_capacity_results": EXPECTED["module_capacity_results"],
        "unique_module_capacity_results": EXPECTED["module_capacity_results"],
        "trace_rows": EXPECTED["trace_rows"],
        "accepted": EXPECTED["accepted"],
        "summary_rows": EXPECTED["modules"],
        "public_csv_rows": EXPECTED["modules"],
        "no_two_copy_capacity": EXPECTED["no_two_copy_capacity"],
        "invalid_accepted_rows": 0,
        "invalid_capacity_proofs": 0,
        "missing_next_copy_100_generation_proofs": 0,
    }
    mismatches = {
        key: {"expected": expected, "observed": observed[key]}
        for key, expected in required_values.items()
        if observed[key] != expected
    }
    if mismatches:
        raise RuntimeError(f"Digs result integrity gate failed: {mismatches}")

    translated = 0
    for batch in pq.ParquetFile(results).iter_batches(
        batch_size=512,
        columns=[
            "validation_passed",
            "unit_sequence",
            "verified_max_copies",
            "dna_sequence",
        ],
    ):
        for row in batch.to_pylist():
            if not row["validation_passed"]:
                continue
            expected_protein = row["unit_sequence"] * int(row["verified_max_copies"])
            if translate_dna(row["dna_sequence"]) != expected_protein:
                raise RuntimeError("Accepted Stage-2 DNA failed translation validation")
            translated += 1
    if translated != EXPECTED["accepted"]:
        raise RuntimeError(
            f"Translated accepted rows differ: {translated} != {EXPECTED['accepted']}"
        )
    observed["translated_accepted_rows"] = translated
    return observed


def scan_audit(
    source: Path,
    *,
    compressed_output: Path | None = None,
) -> dict[str, object]:
    digest = hashlib.sha256()
    hashes: set[str] = set()
    records = 0
    input_handle: BinaryIO
    is_gzip = source.name.endswith(".gz") or source.name.endswith(".gz.tmp")
    input_handle = gzip.open(source, "rb") if is_gzip else source.open("rb")
    raw_output = None
    compressed_handle = None
    try:
        if compressed_output is not None:
            raw_output = compressed_output.open("wb")
            compressed_handle = gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_output,
                compresslevel=6,
                mtime=0,
            )
        for line_number, line in enumerate(input_handle, start=1):
            if not line.strip():
                continue
            digest.update(line)
            if compressed_handle is not None:
                compressed_handle.write(line)
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RuntimeError(f"Invalid IDT audit at line {line_number}") from exc
            response_sha = str(
                record.get("response_sha256")
                or record.get("summary", {}).get("idt_response_sha256")
                or ""
            )
            if not re.fullmatch(r"[0-9a-f]{64}", response_sha):
                raise RuntimeError(f"IDT audit line {line_number} lacks a SHA256")
            hashes.add(response_sha)
            records += 1
    finally:
        input_handle.close()
        if compressed_handle is not None:
            compressed_handle.close()
        if raw_output is not None:
            raw_output.close()
    return {
        "records": records,
        "unique_response_hashes": len(hashes),
        "uncompressed_sha256": digest.hexdigest(),
        "response_hashes": hashes,
    }


def referenced_response_hashes(results: Path) -> set[str]:
    con = duckdb.connect()
    quoted = str(results).replace("'", "''")
    rows = con.execute(
        f"SELECT DISTINCT idt_response_sha256 FROM read_parquet('{quoted}') "
        "WHERE validation_passed"
    ).fetchall()
    return {str(row[0]) for row in rows if row[0]}


def compact_results(results: Path) -> dict[str, object]:
    before_bytes = results.stat().st_size
    temp = results.with_name(results.name + ".compact.tmp")
    temp.unlink(missing_ok=True)
    con = duckdb.connect()
    source = str(results).replace("'", "''")
    target = str(temp).replace("'", "''")
    columns = pq.read_schema(results).names
    if "adaptive_search_trace_json" in columns:
        con.execute(
            f"COPY (SELECT * EXCLUDE (adaptive_search_trace_json) "
            f"FROM read_parquet('{source}')) TO '{target}' "
            "(FORMAT PARQUET, COMPRESSION ZSTD)"
        )
        if "adaptive_search_trace_json" in pq.read_schema(temp).names:
            raise RuntimeError("Compact Stage-2 results still contain embedded trace")
        if pq.ParquetFile(temp).metadata.num_rows != EXPECTED["module_capacity_results"]:
            raise RuntimeError("Compacted Stage-2 result row count changed")
        os.replace(temp, results)
    return {
        "before_bytes": before_bytes,
        "after_bytes": results.stat().st_size,
        "sha256": sha256_file(results),
        "embedded_trace_removed": "adaptive_search_trace_json" not in pq.read_schema(results).names,
    }


def state_record(path: Path) -> dict[str, object]:
    status = subprocess.run(
        ["/net/software/taskrunner/taskrunner", f"--state-file={path}", "status"],
        text=True,
        capture_output=True,
        check=False,
        timeout=30,
    )
    match = re.search(r"Slurm JobID:\s+(\d+)", status.stdout)
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "job_id": int(match.group(1)) if match else None,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true", help="Perform validated compaction/deletion")
    parser.add_argument(
        "--repo",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--scratch-repo",
        type=Path,
        default=(
            Path("/net/scratch")
            / getpass.getuser()
            / "projects/hurdler/clone_repeat_protein"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo = args.repo.resolve()
    scratch_repo = args.scratch_repo.resolve()
    table_dir = (
        repo
        / "studies/hurdler_validation/step04_module_optimization/tables"
        / CORPUS
    )
    public_csv = repo / "data/results/natural_designed_repeat_protein_hurdler_idt.csv"
    results = table_dir / "maximum_copy_results.parquet"
    audit_plain = table_dir / "idt_audit_records.jsonl"
    audit_gzip = table_dir / "idt_audit_records.jsonl.gz"
    manifest_path = repo / "studies/hurdler_validation/DIGS_CLEANUP_MANIFEST.json"

    validation = validate_tables(table_dir, public_csv)
    source_audit = audit_plain if audit_plain.is_file() else audit_gzip
    if not source_audit.is_file():
        raise FileNotFoundError("Neither plain nor compressed IDT audit exists")
    referenced_hashes = referenced_response_hashes(results)

    if not args.apply:
        print(json.dumps({"status": "validated_dry_run", "validation": validation}, indent=2))
        return 0

    original_audit = scan_audit(
        source_audit,
        compressed_output=(audit_gzip.with_suffix(audit_gzip.suffix + ".tmp") if source_audit == audit_plain else None),
    )
    if original_audit["records"] != EXPECTED["idt_audits"]:
        raise RuntimeError(
            f"IDT audit rows differ: {original_audit['records']} != {EXPECTED['idt_audits']}"
        )
    missing_hashes = referenced_hashes - original_audit["response_hashes"]
    if missing_hashes:
        raise RuntimeError(f"Accepted results reference {len(missing_hashes)} missing IDT responses")

    if source_audit == audit_plain:
        audit_temp = audit_gzip.with_suffix(audit_gzip.suffix + ".tmp")
        compressed_check = scan_audit(audit_temp)
        if (
            compressed_check["records"] != original_audit["records"]
            or compressed_check["uncompressed_sha256"]
            != original_audit["uncompressed_sha256"]
            or compressed_check["response_hashes"] != original_audit["response_hashes"]
        ):
            raise RuntimeError("Compressed IDT audit failed round-trip validation")
        os.replace(audit_temp, audit_gzip)
    else:
        compressed_check = original_audit

    compact = compact_results(results)
    post_compaction = validate_tables(table_dir, public_csv)
    if post_compaction != validation:
        raise RuntimeError("Stage-2 table metrics changed during compaction")

    home_duplicates = [
        table_dir / "maximum_copy_results.csv",
        table_dir / "adaptive_copy_search_trace.csv",
        table_dir / "module_final_summary.csv",
        table_dir / "module_compatibility.csv",
    ]
    if audit_plain.is_file():
        home_duplicates.append(audit_plain)
    scratch_targets = [
        scratch_repo
        / "studies/hurdler_validation/step04_module_optimization/runs"
        / "run101_expanded_stage1_compatibility/raw",
        scratch_repo
        / "studies/hurdler_validation/step04_module_optimization/runs"
        / "run103_expanded_stage2_adaptive/raw",
    ]
    state_dirs = [
        repo
        / "studies/hurdler_validation/step04_module_optimization/runs"
        / run
        / "taskfiles"
        for run in (
            "run101_expanded_stage1_compatibility",
            "run103_expanded_stage2_adaptive",
        )
    ]
    state_paths = sorted(path for directory in state_dirs for path in directory.glob("*.state"))
    states = [state_record(path) for path in state_paths]
    for record in states:
        record["path"] = display_path(Path(str(record["path"])), repo, scratch_repo)

    deletion_targets = [path for path in [*home_duplicates, *scratch_targets, *state_paths] if path.exists()]
    deleted_records = []
    for target in deletion_targets:
        deleted_records.append(
            {
                "path": display_path(target, repo, scratch_repo),
                **path_inventory(target),
                "kind": "directory" if target.is_dir() else "file",
            }
        )

    manifest = {
        "schema_version": "digs-cleanup-manifest-v1",
        "status": "validated_pending_cleanup",
        "corpus_version": CORPUS,
        "created_unix_time": time.time(),
        "validation": validation,
        "compaction": compact,
        "idt_audit": {
            "path": display_path(audit_gzip, repo, scratch_repo),
            "records": original_audit["records"],
            "unique_response_hashes": original_audit["unique_response_hashes"],
            "referenced_accepted_response_hashes": len(referenced_hashes),
            "uncompressed_sha256": original_audit["uncompressed_sha256"],
            "compressed_sha256": sha256_file(audit_gzip),
            "compressed_bytes": audit_gzip.stat().st_size,
        },
        "retained_results": {
            display_path(path, repo, scratch_repo): {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in (
                public_csv,
                table_dir / "module_compatibility.parquet",
                table_dir / "module_compatibility_candidates.parquet",
                results,
                table_dir / "adaptive_copy_search_trace.parquet",
                table_dir / "adaptive_copy_validation.parquet",
                table_dir / "module_final_summary.parquet",
                table_dir / "optimized_constructs.fasta",
                table_dir / "optimized_constructs.protein.fasta",
            )
        },
        "taskrunner_states": states,
        "deleted": deleted_records,
        "recovery_commands": [
            "cd $REPO && python scripts/create_module_experiment_tasks.py stage1 --catalog $REPO/studies/hurdler_validation/step03_module_corpus/tables/expanded-middle-repeatsdb-foldseek-v1/module_catalog.parquet --run-dir $REPO/studies/hurdler_validation/step04_module_optimization/runs/run101_expanded_stage1_compatibility --scratch-run-dir $SCRATCH/studies/hurdler_validation/step04_module_optimization/runs/run101_expanded_stage1_compatibility --final-output-dir $REPO/studies/hurdler_validation/step04_module_optimization/tables/expanded-middle-repeatsdb-foldseek-v1 --shards 128",
            "cd $REPO && python scripts/create_module_experiment_tasks.py stage2 --compatibility $REPO/studies/hurdler_validation/step04_module_optimization/tables/expanded-middle-repeatsdb-foldseek-v1/module_compatibility.parquet --run-dir $REPO/studies/hurdler_validation/step04_module_optimization/runs/run103_expanded_stage2_adaptive --scratch-run-dir $SCRATCH/studies/hurdler_validation/step04_module_optimization/runs/run103_expanded_stage2_adaptive --final-output-dir $REPO/studies/hurdler_validation/step04_module_optimization/tables/expanded-middle-repeatsdb-foldseek-v1 --modules-per-task 4",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    for target in deletion_targets:
        if target.is_dir():
            shutil.rmtree(target)
        else:
            target.unlink()

    if audit_plain.exists() or any(target.exists() for target in deletion_targets):
        raise RuntimeError("One or more validated deletion targets remain")
    final_audit = scan_audit(audit_gzip)
    if final_audit["uncompressed_sha256"] != original_audit["uncompressed_sha256"]:
        raise RuntimeError("Retained compressed IDT audit changed after cleanup")
    final_validation = validate_tables(table_dir, public_csv)
    manifest["status"] = "passed"
    manifest["completed_unix_time"] = time.time()
    manifest["validation_after_cleanup"] = final_validation
    manifest["deleted_files"] = sum(item["files"] for item in deleted_records)
    manifest["deleted_logical_bytes"] = sum(item["logical_bytes"] for item in deleted_records)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "passed",
        "manifest": str(manifest_path),
        "deleted_files": manifest["deleted_files"],
        "deleted_logical_bytes": manifest["deleted_logical_bytes"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
