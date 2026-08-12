#!/usr/bin/env python3
"""Curate DHR and THR exact repeat units from primary supplements."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from hurdler.modules import (
    merge_module_catalogs,
    parse_dhr_supplement,
    parse_fasta_modules,
    parse_pdb_exact_modules,
    refine_module_boundaries,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dhr-text", type=Path, required=True)
    parser.add_argument("--thr-fasta", type=Path, required=True)
    parser.add_argument("--thr-pdb-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--boundary-workers", type=int, default=1)
    args = parser.parse_args()
    source_url = "https://www.nature.com/articles/s41586-024-07188-4"
    exclusions = args.output.with_name("thr_fasta_exclusions.csv")
    frames = [
        parse_dhr_supplement(args.dhr_text),
        parse_pdb_exact_modules(
            sorted(args.thr_pdb_dir.glob("THR*_design.pdb")),
            family="THR",
            source_url=source_url,
            experimental_accessions=[f"THR{number}" for number in range(1, 13)],
            structure_accessions=["THR1", "THR2", "THR3", "THR5", "THR6"],
        ),
        parse_fasta_modules(
            [args.thr_fasta],
            family="THR-associated constructs",
            evidence_tier="B",
            source_url=source_url,
            exclusions_path=exclusions,
        ),
    ]
    refined = refine_module_boundaries(
        pd.concat(frames, ignore_index=True),
        audit_path=args.output.with_name("designed_boundary_audit.parquet"),
        candidates_path=args.output.with_name("designed_period_candidates.parquet"),
        unit_alignment_path=args.output.with_name("designed_unit_alignment.parquet"),
        position_variability_path=args.output.with_name("designed_position_variability.parquet"),
        workers=args.boundary_workers,
    )
    catalog = merge_module_catalogs([refined], args.output)
    dhr_count = int(catalog["family"].eq("DHR").sum())
    designed_count = int(catalog["collection"].eq("designed_all").sum())
    primary_count = int(catalog["in_designed_primary100"].sum())
    validation = {
        "dhr_count": dhr_count,
        "designed_unique_count": designed_count,
        "designed_primary100_count": primary_count,
        "dhr1_83_complete": dhr_count == 83,
        "designed_at_least_100": designed_count >= 100,
        "primary100_complete": primary_count == 100,
        "output": str(args.output.resolve()),
    }
    validation["passed"] = all(
        validation[key]
        for key in ("dhr1_83_complete", "designed_at_least_100", "primary100_complete")
    )
    args.output.with_name("designed_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    if not validation["passed"]:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
