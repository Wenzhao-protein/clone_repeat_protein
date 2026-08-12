#!/usr/bin/env python3
"""Validate primitive repeat boundaries and build the periodic-v3 catalog."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from hurdler.modules import merge_module_catalogs
from hurdler.periodicity import BOUNDARY_METHOD_VERSION, MODULE_SELECTION_POLICY


PALETTE = {"Natural": "#4B2E83", "Designed": "#E57200"}


def load(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)


def sequence_coordinate_valid(row: pd.Series) -> bool:
    required = (row.full_sequence_origin, row.unit_start, row.unit_end, row.full_sequence, row.unit_sequence)
    if any(value is None or (not isinstance(value, str) and pd.isna(value)) for value in required):
        return False
    mapping_json = row.get("full_sequence_auth_mapping_json")
    if isinstance(mapping_json, str) and mapping_json.startswith("["):
        mapping = [
            int(match.group()) if (match := re.match(r"^-?\d+", str(value))) else None
            for value in json.loads(mapping_json)
        ]
        indices = [
            index
            for index, coordinate in enumerate(mapping)
            if coordinate is not None and int(row.unit_start) <= coordinate <= int(row.unit_end)
        ]
        observed = "".join(row.full_sequence[index] for index in indices)
        return bool(indices) and observed == row.unit_sequence
    origin = int(row.full_sequence_origin)
    local_start = int(row.unit_start) - origin
    local_end = int(row.unit_end) - origin + 1
    return (
        0 <= local_start < local_end <= len(row.full_sequence)
        and row.full_sequence[local_start:local_end] == row.unit_sequence
    )


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_figure(figure: plt.Figure, base: Path) -> list[Path]:
    outputs = [base.with_suffix(".png"), base.with_suffix(".pdf")]
    for output in outputs:
        figure.savefig(output, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return outputs


def combine_optional_tables(paths: list[Path | None], output: Path) -> int:
    existing = [path for path in paths if path is not None and path.exists()]
    if not existing:
        return 0
    frame = pd.concat([load(path) for path in existing], ignore_index=True)
    frame.to_parquet(output, index=False)
    frame.to_csv(output.with_suffix(".csv"), index=False)
    return len(frame)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--natural", type=Path, required=True)
    parser.add_argument("--designed", type=Path, required=True)
    parser.add_argument("--natural-candidates", type=Path)
    parser.add_argument("--designed-candidates", type=Path)
    parser.add_argument("--natural-units", type=Path)
    parser.add_argument("--designed-units", type=Path)
    parser.add_argument("--natural-positions", type=Path)
    parser.add_argument("--designed-positions", type=Path)
    parser.add_argument("--natural-ss-candidates", type=Path)
    parser.add_argument("--designed-ss-candidates", type=Path)
    parser.add_argument("--natural-ss-residues", type=Path)
    parser.add_argument("--designed-ss-residues", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    natural = load(args.natural)
    designed = load(args.designed)
    source = pd.concat([natural, designed], ignore_index=True)
    source["module_type"] = source.collection.map(
        {"natural100": "Natural", "designed_all": "Designed"}
    )
    source["length_ratio"] = pd.to_numeric(source.prior_unit_length, errors="coerce") / pd.to_numeric(
        source.primitive_period, errors="coerce"
    )
    source["coordinate_valid"] = source.apply(sequence_coordinate_valid, axis=1)
    source["prior_overlaps_region"] = (
        (source.prior_unit_start <= source.repeat_region_end)
        & (source.prior_unit_end >= source.repeat_region_start)
    ).fillna(False)
    source["integer_harmonic"] = np.isclose(
        source.length_ratio, source.length_ratio.round(), atol=0.08
    )
    source["qa_flags"] = [
        ";".join(
            flag
            for flag, active in (
                (
                    "source_prior_fallback",
                    row.boundary_refinement_status == "source_prior_fallback",
                ),
                (
                    "source_annotation_middle_unit",
                    row.boundary_refinement_status
                    == "source_annotation_middle_unit",
                ),
                (
                    "unresolved_boundary",
                    row.boundary_refinement_status
                    not in {
                        "refined",
                        "source_prior_fallback",
                        "source_annotation_middle_unit",
                    },
                ),
                ("low_confidence", row.periodicity_confidence == "low"),
                (
                    "primitive_lt6",
                    pd.notna(row.primitive_period) and int(row.primitive_period) < 6,
                ),
                ("coordinate_mismatch", not bool(row.coordinate_valid)),
                ("prior_outside_region", not bool(row.prior_overlaps_region)),
                (
                    "shorter_noninteger_period",
                    pd.notna(row.primitive_period)
                    and row.primitive_period < row.prior_unit_length
                    and not bool(row.integer_harmonic),
                ),
                ("two_copy_only", pd.notna(row.repeat_count) and int(row.repeat_count) == 2),
                (
                    "secondary_structure_unavailable",
                    getattr(row, "secondary_structure_status", None) != "passed",
                ),
            )
            if active
        )
        for row in source.itertuples(index=False)
    ]
    source.to_parquet(args.output_dir / "module_boundary_audit.parquet", index=False)
    source.to_csv(args.output_dir / "module_boundary_audit.csv", index=False)
    source.loc[source.qa_flags.ne("")].to_csv(
        args.output_dir / "module_boundary_manual_review.csv", index=False
    )

    catalog_path = args.output_dir / "module_catalog_periodic_v4_middle.parquet"
    catalog = merge_module_catalogs([natural, designed], catalog_path)
    source_counts = source.groupby("module_type").size().to_dict()
    catalog_counts = catalog.groupby("collection").size().to_dict()

    raw_dir = args.natural.parent
    combined_counts = {
        "period_candidates": combine_optional_tables(
            [
                args.natural_candidates
                or raw_dir / "natural100_period_candidates.parquet",
                args.designed_candidates
                or args.designed.parent / "designed_period_candidates.parquet",
            ],
            args.output_dir / "module_period_candidates.parquet",
        ),
        "unit_alignment": combine_optional_tables(
            [
                args.natural_units or raw_dir / "natural100_unit_alignment.parquet",
                args.designed_units
                or args.designed.parent / "designed_unit_alignment.parquet",
            ],
            args.output_dir / "module_unit_alignment.parquet",
        ),
        "position_variability": combine_optional_tables(
            [
                args.natural_positions
                or raw_dir / "natural100_position_variability.parquet",
                args.designed_positions
                or args.designed.parent / "designed_position_variability.parquet",
            ],
            args.output_dir / "module_position_variability.parquet",
        ),
        "secondary_structure_candidates": combine_optional_tables(
            [
                args.natural_ss_candidates,
                args.designed_ss_candidates,
            ],
            args.output_dir / "module_secondary_structure_candidates.parquet",
        ),
        "secondary_structure_residues": combine_optional_tables(
            [
                args.natural_ss_residues,
                args.designed_ss_residues,
            ],
            args.output_dir / "module_secondary_structure_residues.parquet",
        ),
    }

    sns.set_theme(style="whitegrid")
    figure_outputs: list[Path] = []
    figure, axis = plt.subplots(figsize=(7.2, 6.2))
    for module_type, part in source.groupby("module_type"):
        axis.scatter(
            part.prior_unit_length,
            part.primitive_period,
            s=34,
            alpha=0.72,
            color=PALETTE[module_type],
            label=f"{module_type} (n={len(part)})",
        )
    maximum = int(max(source.prior_unit_length.max(), source.primitive_period.max()))
    axis.plot([0, maximum], [0, maximum], linestyle="--", color="#444444", linewidth=1)
    axis.set(xlabel="Source candidate unit length (AA)", ylabel="Inferred primitive period (AA)")
    axis.legend(frameon=True)
    axis.set_title("Source candidate units versus full-protein primitive periods")
    figure_outputs += save_figure(figure, args.figure_dir / "source_vs_primitive_module_length")

    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    sns.histplot(
        data=source,
        x="length_ratio",
        hue="module_type",
        palette=PALETTE,
        bins=np.arange(0.75, max(4.25, source.length_ratio.max() + 0.25), 0.25),
        multiple="layer",
        element="step",
        stat="count",
        common_norm=False,
        ax=axis,
    )
    axis.axvline(1, linestyle="--", color="#444444", linewidth=1)
    axis.set(xlabel="Source-unit length / primitive period", ylabel="Proteins")
    axis.set_title("Harmonics revealed by complete-sequence analysis")
    figure_outputs += save_figure(figure, args.figure_dir / "module_harmonic_ratio")

    source["fixed_fraction"] = source.fixed_mask.fillna("").str.count("F") / source.primitive_period
    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    sns.scatterplot(
        data=source,
        x="primitive_period",
        y="fixed_fraction",
        hue="module_type",
        style="periodicity_confidence",
        palette=PALETTE,
        s=56,
        alpha=0.75,
        ax=axis,
    )
    axis.set(
        xlabel="Primitive module length (AA)",
        ylabel="Fixed-position fraction (conservation ≥ 0.8)",
        ylim=(-0.02, 1.02),
    )
    axis.set_title("Fixed and variable ranges across inferred repeat copies")
    figure_outputs += save_figure(figure, args.figure_dir / "module_fixed_fraction")

    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    sns.scatterplot(
        data=source,
        x="primitive_period",
        y="secondary_structure_known_fraction",
        hue="module_type",
        style="secondary_structure_selected_support",
        palette=PALETTE,
        s=58,
        alpha=0.78,
        ax=axis,
    )
    axis.axhline(0.70, linestyle="--", color="#444444", linewidth=1)
    axis.set(
        xlabel="Selected module length (AA)",
        ylabel="Fraction of full sequence with residue-level SS annotation",
        ylim=(-0.02, 1.02),
    )
    axis.set_title("Secondary-structure coverage used in boundary assessment")
    figure_outputs += save_figure(figure, args.figure_dir / "secondary_structure_coverage")

    ss_candidate_path = args.output_dir / "module_secondary_structure_candidates.parquet"
    if ss_candidate_path.exists():
        ss_candidates = load(ss_candidate_path)
        module_types = source.set_index("module_id")["module_type"]
        ss_candidates["module_type"] = ss_candidates.module_id.map(module_types)
        figure, axis = plt.subplots(figsize=(6.4, 5.5))
        sns.scatterplot(
            data=ss_candidates,
            x="sequence_score",
            y="score",
            hue="module_type",
            palette=PALETTE,
            alpha=0.45,
            s=24,
            ax=axis,
        )
        axis.axvline(0.55, linestyle="--", color="#666666", linewidth=1)
        axis.axhline(0.62, linestyle="--", color="#666666", linewidth=1)
        axis.set(
            xlabel="Amino-acid periodicity score",
            ylabel="Secondary-structure periodicity score",
            xlim=(-0.02, 1.02),
            ylim=(-0.02, 1.02),
        )
        axis.set_title("Independent sequence and structure support per candidate")
        figure_outputs += save_figure(
            figure, args.figure_dir / "sequence_vs_secondary_structure_evidence"
        )

    examples = source.loc[source.boundary_refinement_status.eq("refined")].sort_values(
        ["length_ratio", "periodicity_score", "module_id"], ascending=[False, False, True]
    ).head(16)
    figure, axis = plt.subplots(figsize=(10.5, 7.2))
    for y, row in enumerate(examples.itertuples(index=False)):
        origin = int(row.full_sequence_origin)
        full_end = origin + len(row.full_sequence) - 1
        axis.plot([origin, full_end], [y, y], color="#BBBBBB", linewidth=3, solid_capstyle="butt")
        axis.plot(
            [row.repeat_region_start, row.repeat_region_end],
            [y, y],
            color=PALETTE[row.module_type],
            linewidth=7,
            solid_capstyle="butt",
        )
        axis.plot(
            [row.prior_unit_start, row.prior_unit_end],
            [y + 0.18, y + 0.18],
            color="#111111",
            linewidth=2,
            solid_capstyle="butt",
        )
        axis.plot(
            [row.selected_module_start, row.selected_module_end],
            [y - 0.18, y - 0.18],
            color="#B7A57A",
            linewidth=3.2,
            solid_capstyle="butt",
        )
        for boundary in range(
            int(row.repeat_region_start), int(row.repeat_region_end) + 1, int(row.primitive_period)
        ):
            axis.plot([boundary, boundary], [y - 0.25, y + 0.25], color="white", linewidth=0.8)
    axis.set_yticks(range(len(examples)), examples.module_id)
    axis.invert_yaxis()
    axis.set_xlabel("Full-protein coordinate (AA)")
    axis.set_title(
        "Full sequence (gray), repeat region (color), source unit (black), selected middle unit (gold)"
    )
    figure_outputs += save_figure(figure, args.figure_dir / "module_boundary_examples")

    summary = (
        source.groupby("module_type")
        .agg(
            proteins=("module_id", "nunique"),
            median_source_length=("prior_unit_length", "median"),
            median_primitive_length=("primitive_period", "median"),
            shorter_than_source=("length_ratio", lambda values: int((values > 1.08).sum())),
            high_confidence=("periodicity_confidence", lambda values: int((values == "high").sum())),
            medium_confidence=("periodicity_confidence", lambda values: int((values == "medium").sum())),
            low_confidence=("periodicity_confidence", lambda values: int((values == "low").sum())),
            secondary_structure_passed=(
                "secondary_structure_status", lambda values: int((values == "passed").sum())
            ),
            jointly_selected=(
                "secondary_structure_selected_support", lambda values: int(values.fillna(False).sum())
            ),
        )
        .reset_index()
    )
    summary.to_csv(args.output_dir / "module_boundary_summary.csv", index=False)
    validation = {
        "boundary_method_version": BOUNDARY_METHOD_VERSION,
        "source_rows": len(source),
        "source_counts": source_counts,
        "catalog_rows": len(catalog),
        "catalog_counts": catalog_counts,
        "designed_primary100_count": int(catalog.in_designed_primary100.sum()),
        "all_full_sequences_present": bool(source.full_sequence.astype(str).str.len().gt(0).all()),
        "all_boundaries_resolved": bool(
            source.boundary_refinement_status.isin(
                [
                    "refined",
                    "source_prior_fallback",
                    "source_annotation_middle_unit",
                ]
            ).all()
        ),
        "source_prior_fallback_count": int(
            source.boundary_refinement_status.eq("source_prior_fallback").sum()
        ),
        "source_annotation_fallback_count": int(
            source.boundary_refinement_status.eq("source_annotation_middle_unit").sum()
        ),
        "all_coordinates_valid": bool(source.coordinate_valid.all()),
        "all_prior_units_overlap_repeat_region": bool(source.prior_overlaps_region.all()),
        "module_selection_policy": MODULE_SELECTION_POLICY,
        "all_selected_modules_are_middle_policy": bool(
            source.selected_module_policy.eq(MODULE_SELECTION_POLICY).all()
        ),
        "all_unit_sequences_are_selected_module": bool(
            source.unit_sequence.eq(source.selected_module_sequence).all()
            and source.unit_start.eq(source.selected_module_start).all()
            and source.unit_end.eq(source.selected_module_end).all()
        ),
        "low_confidence_count": int(source.periodicity_confidence.eq("low").sum()),
        "secondary_structure_passed_count": int(
            source.secondary_structure_status.eq("passed").sum()
        ),
        "secondary_structure_unavailable_count": int(
            source.secondary_structure_status.ne("passed").sum()
        ),
        "jointly_selected_count": int(
            source.secondary_structure_selected_support.fillna(False).sum()
        ),
        "manual_review_count": int(source.qa_flags.ne("").sum()),
        "combined_table_rows": combined_counts,
    }
    validation["passed"] = bool(
        validation["source_counts"].get("Natural", 0) == 100
        and validation["source_counts"].get("Designed", 0) >= 100
        and validation["catalog_counts"].get("natural100", 0) == 100
        and validation["catalog_counts"].get("designed_all", 0) >= 100
        and validation["designed_primary100_count"] == 100
        and validation["all_full_sequences_present"]
        and validation["all_boundaries_resolved"]
        and validation["all_coordinates_valid"]
        and validation["all_prior_units_overlap_repeat_region"]
        and validation["all_selected_modules_are_middle_policy"]
        and validation["all_unit_sequences_are_selected_module"]
        and validation["secondary_structure_unavailable_count"] == 0
    )
    (args.output_dir / "module_boundary_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )
    manifest_rows = [
        {
            "path": str(path),
            "sha256": file_sha256(path),
            "bytes": path.stat().st_size,
            "status": "passed" if path.stat().st_size > 0 else "failed_empty",
        }
        for path in figure_outputs
    ]
    pd.DataFrame(manifest_rows).to_csv(args.output_dir / "module_boundary_figure_manifest.csv", index=False)
    print(json.dumps({**validation, "summary": summary.to_dict(orient="records")}, indent=2))
    if not validation["passed"]:
        raise RuntimeError("Boundary validation did not pass; canonical catalog remains untouched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
