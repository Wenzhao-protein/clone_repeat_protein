#!/usr/bin/env python3
"""Merge independently recoverable natural and designed curation runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from hurdler.modules import merge_module_catalogs


def load(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--natural", type=Path, required=True)
    parser.add_argument("--designed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    catalog = merge_module_catalogs([load(args.natural), load(args.designed)], args.output)
    natural_count = int(catalog["collection"].eq("natural100").sum())
    designed_count = int(catalog["collection"].eq("designed_all").sum())
    primary_count = int(catalog["in_designed_primary100"].sum())
    standard_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    validation = {
        "natural_count": natural_count,
        "designed_unique_count": designed_count,
        "designed_primary100_count": primary_count,
        "all_standard_amino_acids": bool(
            catalog["unit_sequence"].map(lambda value: not (set(value) - standard_amino_acids)).all()
        ),
        "unique_within_collection": not bool(
            catalog.duplicated(["collection", "unit_sequence"]).any()
        ),
        "download_dates_complete": bool(catalog["download_date"].astype(str).str.len().gt(0).all()),
        "output": str(args.output.resolve()),
    }
    validation["passed"] = (
        natural_count == 100
        and designed_count >= 100
        and primary_count == 100
        and validation["all_standard_amino_acids"]
        and validation["unique_within_collection"]
        and validation["download_dates_complete"]
    )
    args.output.with_name("module_catalog_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )
    print(json.dumps(validation, indent=2))
    if not validation["passed"]:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
