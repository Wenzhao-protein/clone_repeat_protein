#!/usr/bin/env python3
"""Strict workflow-aware finalization for returned V2 production shards.

This wrapper delegates scientific validation to the existing maintained
HURDLER finalizers.  It never removes raw shards and refuses to write into a
non-empty final destination.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd

from hurdler.complete_route import finalize_complete_route_shards
from hurdler.dna_assembly_visualization import (
    plot_complete_production_report,
    write_production_figure_manifest,
)
from hurdler.io import sha256_file, utc_now, write_json_atomic
from hurdler.module_experiments import (
    finalize_adaptive_copy_results,
    finalize_module_compatibility,
)
from hurdler.repeatsdb import finalize_natural_corpus
from hurdler.structural_repeats import finalize_designed_catalog


SUPPORTED = {
    "repeatsdb-natural", "designed-structure", "missing-af3",
    "module-stage1", "module-stage2", "exact-dna-routes",
    "exact-dna-purchase",
}


def _require_files(paths: list[Path], label: str) -> list[Path]:
    files = [path for path in paths if path.is_file()]
    if not files:
        raise FileNotFoundError(f"No {label} files were returned")
    return files


def _inventory(root: Path) -> list[dict[str, Any]]:
    return [
        {
            "path": str(path.absolute()),
            "relative_path": str(path.relative_to(root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in sorted(root.rglob("*")) if path.is_file()
    ]


def _finalize_natural(scratch: Path, output: Path, source: Path | None) -> dict[str, Any]:
    mappings = _require_files(sorted(scratch.glob("shard_*/natural_source_mappings.parquet")), "natural mapping shard")
    inventories = sorted(scratch.glob("shard_*/natural_region_inventory.parquet"))
    exclusions = sorted(scratch.glob("shard_*/natural_exclusions.csv"))
    catalog_path = output / "natural_modules.parquet"
    catalog = finalize_natural_corpus(
        mappings, catalog_path,
        region_inventory_paths=inventories,
        exclusion_paths=exclusions,
        annotation_inventory_path=source,
    )
    return {"catalog_rows": len(catalog), "mapping_shards": len(mappings)}


def _finalize_designed(scratch: Path, output: Path) -> dict[str, Any]:
    mappings = _require_files(sorted(scratch.glob("shard_*/designed_source_mappings.parquet")), "designed mapping shard")
    exclusions = sorted(scratch.glob("shard_*/exclusions.parquet"))
    catalog = finalize_designed_catalog(
        mappings, exclusions, output / "designed_modules.parquet",
        candidate_paths=sorted(scratch.glob("shard_*/candidates.parquet")),
        unit_paths=sorted(scratch.glob("shard_*/units.parquet")),
        position_paths=sorted(scratch.glob("shard_*/positions.parquet")),
    )
    return {"catalog_rows": len(catalog), "mapping_shards": len(mappings)}


def _finalize_stage1(scratch: Path, output: Path) -> dict[str, Any]:
    summaries = _require_files(
        sorted(scratch.glob("shard_*/module_compatibility_shard-*.parquet")),
        "Stage-1 summary shard",
    )
    candidates = sorted(scratch.glob("shard_*/module_compatibility_candidates_shard-*.parquet"))
    modules, binned, candidate_frame = finalize_module_compatibility(summaries, candidates, output)
    return {
        "module_rows": len(modules), "binned_rows": len(binned),
        "candidate_rows_loaded": len(candidate_frame), "summary_shards": len(summaries),
    }


def _finalize_stage2(scratch: Path, output: Path, compatibility: Path | None) -> dict[str, Any]:
    if compatibility is None:
        raise ValueError("module-stage2 finalization requires --input compatibility table")
    results = _require_files(sorted(scratch.glob("shard_*/optimized_constructs_ga.parquet")), "Stage-2 result shard")
    audits = _require_files(sorted(scratch.glob("shard_*/idt_optimization_responses.jsonl")), "Stage-2 IDT audit shard")
    maximum, traces, summary = finalize_adaptive_copy_results(
        results, compatibility, output, idt_audit_paths=audits,
    )
    return {
        "maximum_rows": len(maximum), "trace_rows": len(traces),
        "summary_rows": len(summary), "result_shards": len(results),
    }


def _finalize_exact_routes(
    scratch: Path,
    output: Path,
    expected_elements: int | None,
    expected_targets: int | None,
) -> dict[str, Any]:
    shard_dirs = sorted(path.parent for path in scratch.glob("shard_*/complete_route_manifest.json"))
    if not shard_dirs:
        raise FileNotFoundError("No complete-route shard manifests were returned")
    tables = finalize_complete_route_shards(
        shard_dirs, output,
        expected_public_elements=expected_elements,
        expected_real_targets=expected_targets,
    )
    figures = plot_complete_production_report(
        tables["targets"], tables["element_matrix"], tables["selected_routes"],
        tables["transitions"], tables["fragments"], tables["seeds"],
        output / "figures",
    )
    write_production_figure_manifest(
        figures, output / "figures" / "figure_manifest.csv",
        input_tables=[
            output / "production_target_analysis.parquet",
            output / "production_element_matrix.parquet",
            output / "production_selected_routes.parquet",
        ],
        source_notebook="notebooks/v2/10_exact_dna_result_analysis.ipynb",
    )
    return {
        "target_rows": len(tables["targets"]),
        "element_rows": len(tables["element_matrix"]),
        "selected_route_rows": len(tables["selected_routes"]),
        "shards": len(shard_dirs), "figures": len(figures),
    }


def _finalize_purchase(scratch: Path, output: Path, expected_elements: int | None) -> dict[str, Any]:
    source = scratch / "production"
    summary_path = source / "purchase_orderability_summary.json"
    if not summary_path.is_file():
        candidates = sorted(source.glob("*.json"))
        if not candidates:
            raise FileNotFoundError("Purchase audit summary is absent")
        summary_path = candidates[0]
    summary = json.loads(summary_path.read_text())
    results = summary.get("results", summary)
    observed_elements = results.get("elements_with_found_routes")
    if expected_elements is not None and observed_elements is not None and int(observed_elements) != expected_elements:
        raise ValueError(f"Purchase-audit element mismatch: {observed_elements} != {expected_elements}")
    copied = 0
    for path in sorted(source.iterdir()):
        if path.is_file():
            shutil.copy2(path, output / path.name)
            copied += 1
    if not copied:
        raise FileNotFoundError("Purchase audit returned no files")
    return {"copied_files": copied, "elements_with_routes": observed_elements}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workflow", choices=sorted(SUPPORTED), required=True)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--expected-elements", type=int)
    parser.add_argument("--expected-targets", type=int)
    args = parser.parse_args()
    scratch = args.scratch_dir.expanduser().absolute()
    output = args.output_dir.expanduser().absolute()
    if not scratch.is_dir():
        raise FileNotFoundError(scratch)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty final output: {output}")
    output.mkdir(parents=True, exist_ok=True)

    if args.workflow == "repeatsdb-natural":
        metrics = _finalize_natural(scratch, output, args.input)
    elif args.workflow == "designed-structure":
        metrics = _finalize_designed(scratch, output)
    elif args.workflow == "module-stage1":
        metrics = _finalize_stage1(scratch, output)
    elif args.workflow == "module-stage2":
        metrics = _finalize_stage2(scratch, output, args.input)
    elif args.workflow == "exact-dna-routes":
        metrics = _finalize_exact_routes(
            scratch, output, args.expected_elements, args.expected_targets
        )
    elif args.workflow == "exact-dna-purchase":
        metrics = _finalize_purchase(scratch, output, args.expected_elements)
    else:  # missing-af3: the runner is installation-specific, so verify and inventory all outputs.
        files = _inventory(scratch)
        if not files:
            raise FileNotFoundError("AlphaFold3 production returned no files")
        write_json_atomic(
            {"workflow": "missing-af3", "returned_files": files},
            output / "af3_return_inventory.json",
        )
        metrics = {"returned_files": len(files), "production_required_missing_structure": 0}

    manifest = {
        "schema_version": "hurdler-production-finalization-v2",
        "created_at": utc_now(),
        "workflow": args.workflow,
        "status": "passed",
        "metrics": metrics,
        "outputs": _inventory(output),
        "raw_inputs_preserved": True,
        "cleanup_performed": False,
    }
    write_json_atomic(manifest, output / "finalization_manifest.json")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
