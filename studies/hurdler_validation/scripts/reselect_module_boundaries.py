#!/usr/bin/env python3
"""Apply the current boundary selector to frozen period-candidate scores."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from hurdler.modules import merge_module_catalogs, reselect_module_boundaries
from hurdler.periodicity import BOUNDARY_METHOD_VERSION


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = pd.read_parquet(args.catalog)
    candidates = pd.read_parquet(args.candidates)
    refined = reselect_module_boundaries(
        source,
        candidates,
        output_path=args.output.with_name(args.output.stem + "_boundary_audit.parquet"),
        unit_alignment_path=args.output.with_name(args.output.stem + "_unit_alignment.parquet"),
        position_variability_path=args.output.with_name(
            args.output.stem + "_position_variability.parquet"
        ),
    )
    catalog = merge_module_catalogs([refined], args.output)
    validation = {
        "boundary_method_version": BOUNDARY_METHOD_VERSION,
        "source_rows": len(source),
        "catalog_rows": len(catalog),
        "all_reselected": bool(refined.boundary_reselected_from_candidates.fillna(False).all()),
        "all_boundaries_resolved": bool(
            refined.boundary_refinement_status.isin(["refined", "source_prior_fallback"]).all()
        ),
        "shortened_source_units": int(
            (refined.primitive_period < refined.prior_unit_length).sum()
        ),
        "output": str(args.output),
    }
    validation["passed"] = bool(
        validation["source_rows"] == validation["catalog_rows"]
        and validation["all_reselected"]
        and validation["all_boundaries_resolved"]
    )
    args.output.with_name(args.output.stem + "_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0 if validation["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
