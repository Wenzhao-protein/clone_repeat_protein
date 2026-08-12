"""Versioned two-stage repeat-module experiment helpers.

Stage 1 is deliberately independent from codon optimization: it evaluates the
exact middle module against the frozen lookup, preserves every solution, and
materializes both per-module and binned tables.  Those tables are the only
inputs required by the reporting notebook.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .constants import PLASMIDS, RULE_PROFILE_NAME, validate_protein_sequence
from .index import PatternIndex
from .matching import enumerate_module_solutions, expand_short_module
from .optimization import (
    _candidate_rank,
    _construct_metrics,
    load_codon_weights,
    recognition_site_count,
    translate_dna,
)

CORPUS_VERSION = "expanded-middle-repeatsdb-foldseek-v1"
COMPATIBILITY_METHOD_VERSION = "middle-module-all-plasmids-legacy-optimized-v1"

_UW_PURPLE = "#4B2E83"
_UW_GOLD = "#B7A57A"
_NATURAL = "#0072B2"
_DESIGNED = "#D55E00"
_CANDIDATE_PANDAS_ROW_LIMIT = 250_000

# Stage-1 is an experiment table, not a second copy of the complete structure
# catalog.  Keep the identifiers and boundary evidence needed for analysis;
# large full-chain sequences, DSSP strings, alignments and source maps remain
# normalized in the authoritative corpus artifacts.
_STAGE1_PROVENANCE_FIELDS = (
    "module_id",
    "module_type",
    "family",
    "unit_sequence",
    "unit_length",
    "evidence_tier",
    "source_name",
    "source_url",
    "source_accession",
    "source_chain",
    "source_annotation_id",
    "annotation_uuid",
    "structure_accession",
    "structure_chain",
    "structure_source",
    "protein_key",
    "unit_start",
    "unit_end",
    "selected_module_start",
    "selected_module_end",
    "selected_module_index",
    "selected_module_count",
    "repeat_region_start",
    "repeat_region_end",
    "region_start",
    "region_end",
    "boundary_method",
    "boundary_method_version",
    "boundary_refinement_status",
    "module_selection_policy",
    "sequence_coordinate_method",
    "strict_dual_evidence_passed",
    "dssp_state_agreement",
    "dssp_transition_agreement",
    "dssp_median_transitions_per_unit",
    "foldseek_3di_identity",
    "foldseek_median_min_tm",
    "foldseek_median_lddt",
    "foldseek_median_coverage",
    "structure_source_type",
    "structure_sequence_sha256",
    "sequence_sha256",
    "source_sha256",
    "retrieved_date",
    "download_date",
    "citation",
    "license_name",
)

STAGE2_INPUT_COLUMNS = (
    *_STAGE1_PROVENANCE_FIELDS,
    "collection",
    "hurdler_compatible",
    "selected_solution_json",
    "selected_plasmid",
    "selected_site_i_enzyme",
    "selected_site_ii_enzyme",
    "selected_site_i_recognition_site",
    "selected_site_ii_recognition_site",
    "selected_direction",
    "selected_site_i_position",
    "selected_site_ii_position",
    "selected_candidate_pair_id",
)


def _read_table(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    if source.suffix == ".parquet":
        return pd.read_parquet(source)
    if source.suffix == ".csv":
        return pd.read_csv(source)
    raise ValueError(f"Expected a Parquet or CSV table: {source}")


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _length_bin(length: int) -> tuple[int, int, str]:
    lower = ((int(length) - 1) // 10) * 10 + 1
    upper = lower + 9
    return lower, upper, f"{lower}–{upper}"


def _collection_label(row: dict[str, Any]) -> str:
    module_type = str(row.get("module_type", "")).strip().lower()
    collection = str(row.get("collection", "")).strip().lower()
    if module_type == "natural" or collection.startswith("natural"):
        return "Natural"
    if module_type == "designed" or collection.startswith("designed"):
        return "Designed"
    raise ValueError(
        f"Module {row.get('module_id')} is neither Natural nor Designed"
    )


def _require_stage1_catalog(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"module_id", "unit_sequence", "unit_length", "collection"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Module catalog is missing required columns: {missing}")
    normalized = frame.copy()
    normalized["unit_sequence"] = normalized.unit_sequence.map(
        validate_protein_sequence
    )
    if not normalized.unit_length.astype(int).eq(
        normalized.unit_sequence.str.len()
    ).all():
        raise ValueError("unit_length does not match unit_sequence")
    normalized["experiment_collection"] = [
        _collection_label(row) for row in normalized.to_dict(orient="records")
    ]
    duplicated = normalized.duplicated(
        ["experiment_collection", "unit_sequence"], keep=False
    )
    if duplicated.any():
        ids = normalized.loc[duplicated, "module_id"].astype(str).head(5).tolist()
        raise ValueError(
            "Stage 1 requires one exact middle-module sequence per collection; "
            f"duplicates include {ids}"
        )
    return normalized.sort_values(
        ["experiment_collection", "module_id"], kind="mergesort"
    ).reset_index(drop=True)


def run_module_compatibility(
    catalog_path: str | Path,
    index_dir: str | Path,
    output_dir: str | Path,
    *,
    shard_index: int = 0,
    shard_count: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate one deterministic catalog shard and preserve all candidates."""
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must satisfy 0 <= shard_index < shard_count")
    catalog = _require_stage1_catalog(_read_table(catalog_path))
    catalog = catalog.iloc[shard_index::shard_count].copy()
    index = PatternIndex.load(index_dir)
    summaries: list[dict[str, Any]] = []
    candidates: list[dict[str, Any]] = []
    for module in catalog.to_dict(orient="records"):
        sequence = str(module["unit_sequence"])
        effective, expansion_copies = (
            expand_short_module(sequence) if len(sequence) < 6 else (sequence, 1)
        )
        solutions = enumerate_module_solutions(sequence, index)
        ordered = sorted(solutions, key=_candidate_rank)
        selected = ordered[0] if ordered else {}
        for rank, solution in enumerate(ordered, start=1):
            candidate = {
                "module_id": module["module_id"],
                "collection": module["experiment_collection"],
                "unit_sequence": sequence,
                "unit_length": len(sequence),
                "candidate_rank": rank,
                "selected_candidate": rank == 1,
                "corpus_version": str(
                    module.get("corpus_version") or CORPUS_VERSION
                ),
                "rules_version": RULE_PROFILE_NAME,
                **solution,
            }
            candidates.append(candidate)
        lower, upper, label = _length_bin(len(sequence))
        compatible_plasmids = sorted(
            {str(value["plasmid"]) for value in ordered},
            key=PLASMIDS.index,
        )
        provenance = {
            field: module[field]
            for field in _STAGE1_PROVENANCE_FIELDS
            if field in module
        }
        summary = {
            **provenance,
            "collection": module["experiment_collection"],
            "corpus_version": str(module.get("corpus_version") or CORPUS_VERSION),
            "compatibility_method_version": COMPATIBILITY_METHOD_VERSION,
            "rules_version": RULE_PROFILE_NAME,
            "effective_module_sequence": effective,
            "effective_module_length": len(effective),
            "short_module_expansion_copies": expansion_copies,
            "hurdler_compatible": bool(ordered),
            "compatible_plasmid_count": len(compatible_plasmids),
            "compatible_plasmids_json": json.dumps(compatible_plasmids),
            "candidate_solution_count": len(ordered),
            "selected_solution_json": json.dumps(
                selected, sort_keys=True, default=str
            ),
            "length_bin_lower": lower,
            "length_bin_upper": upper,
            "length_bin": label,
            "selected_plasmid": selected.get("plasmid"),
            "selected_site_i_enzyme": selected.get("site_i_enzyme"),
            "selected_site_ii_enzyme": selected.get("site_ii_enzyme"),
            "selected_site_i_recognition_site": selected.get(
                "site_i_recognition_site"
            ),
            "selected_site_ii_recognition_site": selected.get(
                "site_ii_recognition_site"
            ),
            "selected_direction": selected.get("direction"),
            "selected_site_i_position": selected.get("site_i_position"),
            "selected_site_ii_position": selected.get("site_ii_position"),
            "selected_candidate_pair_id": selected.get(
                "candidate_pair_id", selected.get("best_pair_id")
            ),
        }
        summaries.append(summary)

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    suffix = f"shard-{shard_index:05d}-of-{shard_count:05d}"
    summary_frame = pd.DataFrame(summaries)
    candidate_frame = pd.DataFrame(candidates)
    summary_path = destination / f"module_compatibility_{suffix}.parquet"
    candidate_path = destination / f"module_compatibility_candidates_{suffix}.parquet"
    summary_frame.to_parquet(summary_path, index=False)
    candidate_frame.to_parquet(candidate_path, index=False)
    manifest = {
        "corpus_version": CORPUS_VERSION,
        "method_version": COMPATIBILITY_METHOD_VERSION,
        "rules_version": RULE_PROFILE_NAME,
        "catalog": str(Path(catalog_path).resolve()),
        "catalog_sha256": _sha256(catalog_path),
        "index_dir": str(Path(index_dir).resolve()),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "input_rows": len(catalog),
        "compatible_rows": int(
            summary_frame.get("hurdler_compatible", pd.Series(dtype=bool)).sum()
        ),
        "candidate_rows": len(candidate_frame),
    }
    (destination / f"module_compatibility_{suffix}.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return summary_frame, candidate_frame


def summarize_compatibility(per_module: pd.DataFrame) -> pd.DataFrame:
    """Build complete shared 10-AA bins for both collections."""
    if per_module.empty:
        raise ValueError("Cannot summarize an empty compatibility table")
    maximum = int(per_module.unit_length.max())
    maximum_upper = _length_bin(maximum)[1]
    rows: list[dict[str, Any]] = []
    for collection in ("Natural", "Designed"):
        subset = per_module.loc[per_module.collection.eq(collection)]
        for lower in range(1, maximum_upper + 1, 10):
            upper = lower + 9
            values = subset.loc[subset.unit_length.between(lower, upper)]
            compatible = int(values.hurdler_compatible.astype(bool).sum())
            total = len(values)
            incompatible = total - compatible
            rows.append(
                {
                    "collection": collection,
                    "length_bin_lower": lower,
                    "length_bin_upper": upper,
                    "length_bin": f"{lower}–{upper}",
                    "compatible_count": compatible,
                    "incompatible_count": incompatible,
                    "total_count": total,
                    "compatible_fraction": compatible / total if total else 0.0,
                    "incompatible_fraction": incompatible / total if total else 0.0,
                }
            )
    result = pd.DataFrame(rows)
    if not (
        result.compatible_count + result.incompatible_count
    ).eq(result.total_count).all():
        raise AssertionError("Compatibility counts do not sum to bin totals")
    return result


def plot_compatibility(summary: pd.DataFrame, output_stem: str | Path) -> None:
    """Create the requested absolute and 100%-stacked four-panel figure."""
    collections = ("Natural", "Designed")
    bins = (
        summary.sort_values("length_bin_lower").length_bin.drop_duplicates().tolist()
    )
    x = np.arange(len(bins))
    figure, axes = plt.subplots(
        2, 2, figsize=(max(12, 0.55 * len(bins)), 10), sharex=True
    )
    for column, collection in enumerate(collections):
        values = (
            summary.loc[summary.collection.eq(collection)]
            .set_index("length_bin")
            .reindex(bins)
            .fillna(0)
        )
        incompatible = values.incompatible_count.to_numpy(dtype=float)
        compatible = values.compatible_count.to_numpy(dtype=float)
        total = values.total_count.to_numpy(dtype=float)
        top = axes[0, column]
        top.bar(
            x, incompatible, color=_UW_GOLD, edgecolor="white", label="Incompatible"
        )
        top.bar(
            x,
            compatible,
            bottom=incompatible,
            color=_UW_PURPLE,
            edgecolor="white",
            hatch="//",
            label="HURDLER-compatible",
        )
        for index, (no, yes, count) in enumerate(
            zip(incompatible, compatible, total, strict=True)
        ):
            if no:
                top.text(index, no / 2, f"{int(no)}", ha="center", va="center", fontsize=7)
            if yes:
                top.text(index, no + yes / 2, f"{int(yes)}", ha="center", va="center", fontsize=7, color="white")
            top.text(index, count, f"n={int(count)}", ha="center", va="bottom", fontsize=7)
        top.set_title(f"{collection}: counts")
        top.set_ylabel("Unique middle modules")
        top.grid(axis="y", alpha=0.2)

        bottom = axes[1, column]
        no_fraction = np.divide(
            incompatible, total, out=np.zeros_like(incompatible), where=total > 0
        )
        yes_fraction = np.divide(
            compatible, total, out=np.zeros_like(compatible), where=total > 0
        )
        bottom.bar(x, no_fraction, color=_UW_GOLD, edgecolor="white")
        bottom.bar(
            x,
            yes_fraction,
            bottom=no_fraction,
            color=_UW_PURPLE,
            edgecolor="white",
            hatch="//",
        )
        for index, (no, yes, count) in enumerate(
            zip(no_fraction, yes_fraction, total, strict=True)
        ):
            if count:
                if no >= 0.08:
                    bottom.text(index, no / 2, f"{no:.0%}\n({int(incompatible[index])})", ha="center", va="center", fontsize=7)
                if yes >= 0.08:
                    bottom.text(index, no + yes / 2, f"{yes:.0%}\n({int(compatible[index])})", ha="center", va="center", fontsize=7, color="white")
            bottom.text(
                index,
                1.01,
                f"n={int(count)}",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=5.5,
            )
        bottom.set_title(f"{collection}: proportions")
        bottom.set_ylim(0, 1.17)
        bottom.set_ylabel("Fraction within length bin")
        bottom.set_xlabel("Middle-module length (AA)")
        bottom.set_xticks(x, bins, rotation=45, ha="right")
        bottom.grid(axis="y", alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.968),
        ncol=2,
        frameon=False,
    )
    figure.suptitle(
        "HURDLER compatibility of RepeatsDB-direct and DSSP/Foldseek middle modules",
        y=0.997,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def finalize_module_compatibility(
    summary_paths: Iterable[str | Path],
    candidate_paths: Iterable[str | Path],
    output_dir: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Merge stable shards, validate cardinality, and generate Stage-1 outputs."""
    summary_files = [Path(path) for path in summary_paths]
    candidate_files = [Path(path) for path in candidate_paths]
    if not summary_files:
        raise ValueError("At least one summary shard is required")
    per_module = pd.concat(
        [pd.read_parquet(path) for path in summary_files], ignore_index=True
    )
    duplicate = per_module.duplicated(["collection", "unit_sequence"], keep=False)
    if duplicate.any():
        raise ValueError("Duplicate collection/unit_sequence rows across Stage-1 shards")
    per_module = per_module.sort_values(
        ["collection", "unit_length", "module_id"], kind="mergesort"
    ).reset_index(drop=True)
    binned = summarize_compatibility(per_module)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    for name, frame in (
        ("module_compatibility", per_module),
        ("module_compatibility_binned", binned),
    ):
        frame.to_parquet(destination / f"{name}.parquet", index=False)
        frame.to_csv(destination / f"{name}.csv", index=False)

    candidate_row_count = sum(
        int(pq.ParquetFile(path).metadata.num_rows) for path in candidate_files
    )
    candidate_output = destination / "module_compatibility_candidates.parquet"
    candidate_csv_written = candidate_row_count <= _CANDIDATE_PANDAS_ROW_LIMIT
    if candidate_csv_written:
        candidates = (
            pd.concat(
                [pd.read_parquet(path) for path in candidate_files],
                ignore_index=True,
            )
            if candidate_files
            else pd.DataFrame()
        )
        candidates.to_parquet(candidate_output, index=False)
        candidates.to_csv(candidate_output.with_suffix(".csv"), index=False)
    else:
        # Candidate solutions are a normalized large table.  Stream batches
        # into one Parquet artifact rather than constructing a many-gigabyte
        # pandas DataFrame or duplicating it as CSV.
        writer: pq.ParquetWriter | None = None
        output_schema: pa.Schema | None = None
        try:
            for path in candidate_files:
                parquet = pq.ParquetFile(path)
                if parquet.metadata.num_rows == 0 or parquet.metadata.num_columns == 0:
                    continue
                if output_schema is None:
                    output_schema = parquet.schema_arrow
                    writer = pq.ParquetWriter(
                        candidate_output,
                        output_schema,
                        compression="zstd",
                    )
                elif parquet.schema_arrow != output_schema:
                    raise ValueError(f"Candidate shard schema mismatch: {path}")
                assert writer is not None
                for batch in parquet.iter_batches(batch_size=100_000):
                    writer.write_table(pa.Table.from_batches([batch], output_schema))
        finally:
            if writer is not None:
                writer.close()
        if writer is None:
            pd.DataFrame().to_parquet(candidate_output, index=False)
        candidates = pd.DataFrame(columns=list(output_schema.names) if output_schema else [])
        candidates.attrs["candidate_row_count"] = candidate_row_count
    plot_compatibility(
        binned, destination / "module_compatibility_by_length"
    )
    manifest = {
        "corpus_version": CORPUS_VERSION,
        "method_version": COMPATIBILITY_METHOD_VERSION,
        "summary_shards": [str(path.resolve()) for path in summary_files],
        "candidate_shards": [str(path.resolve()) for path in candidate_files],
        "module_rows": len(per_module),
        "candidate_rows": candidate_row_count,
        "candidate_table_format": "normalized_parquet",
        "candidate_csv_written": candidate_csv_written,
        "compatible_rows": int(per_module.hurdler_compatible.sum()),
        "collection_counts": per_module.collection.value_counts().to_dict(),
    }
    (destination / "module_compatibility.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return per_module, candidates, binned


def prepare_adaptive_copy_frame(
    source: pd.DataFrame,
    *,
    codon_weights: dict[str, float],
    fragment_limits: tuple[int, ...] = (1800, 3000),
    external_deduction_bp: int = 0,
) -> pd.DataFrame:
    """Prepare compact Stage-2 rows using the frozen Stage-1 selected pair."""
    if external_deduction_bp < 0:
        raise ValueError("external_deduction_bp must be non-negative")
    required = {
        "module_id",
        "collection",
        "unit_sequence",
        "unit_length",
        "hurdler_compatible",
        "selected_solution_json",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise ValueError(f"Compatibility table is missing columns: {missing}")
    rows: list[dict[str, Any]] = []
    for module in source.loc[source.hurdler_compatible.astype(bool)].to_dict(
        orient="records"
    ):
        solution = json.loads(str(module["selected_solution_json"]))
        unit = validate_protein_sequence(str(module["unit_sequence"]))
        provenance = {
            field: module[field]
            for field in STAGE2_INPUT_COLUMNS
            if field in module
        }
        for cap in fragment_limits:
            available_bp = max(0, int(cap) - external_deduction_bp)
            mathematical_max = available_bp // (3 * len(unit))
            minimum = max(
                2,
                int(
                    np.ceil(
                        (
                            max(
                                int(solution["site_i_position"]),
                                int(solution["site_ii_position"]),
                            )
                            + 3
                        )
                        / len(unit)
                    )
                ),
            )
            row = {
                **provenance,
                **solution,
                "collection": _collection_label(module),
                "unit_sequence": unit,
                "unit_length": len(unit),
                "fragment_limit_bp": int(cap),
                "external_deduction_bp": external_deduction_bp,
                "available_coding_bp": available_bp,
                "mathematical_max_copies": mathematical_max,
                "verified_max_copies": None,
                "stage2_minimum_copies": minimum,
                "stage2_selected_pair_frozen": True,
                "stage2_preparation_status": "pending",
            }
            if mathematical_max < minimum:
                row.update(
                    dna_sequence=None,
                    dna_length=0,
                    stage2_preparation_status="no_two_copy_capacity",
                    failure_reason="Fragment cap cannot contain two complete modules and both locked windows",
                )
            else:
                try:
                    metrics = _construct_metrics(
                        unit,
                        minimum,
                        solution,
                        codon_weights,
                        validate_hard_constraints=False,
                    )
                    row.update(metrics)
                    row["stage2_preparation_status"] = "prepared"
                    row["failure_reason"] = ""
                except Exception as exc:
                    row.update(
                        dna_sequence=None,
                        dna_length=0,
                        stage2_preparation_status="preparation_failed",
                        failure_reason=f"{type(exc).__name__}: {exc}",
                    )
            rows.append(row)
    return pd.DataFrame(rows)


def prepare_adaptive_copy_inputs(
    compatibility_path: str | Path,
    output_path: str | Path,
    *,
    codon_usage_path: str | Path,
    fragment_limits: tuple[int, ...] = (1800, 3000),
    external_deduction_bp: int = 0,
) -> pd.DataFrame:
    """Prepare and materialize compact Stage-2 rows from a compatibility table."""
    source = _read_table(compatibility_path)
    frame = prepare_adaptive_copy_frame(
        source,
        codon_weights=load_codon_weights(codon_usage_path),
        fragment_limits=fragment_limits,
        external_deduction_bp=external_deduction_bp,
    )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(destination, index=False)
    frame.to_csv(destination.with_suffix(".csv"), index=False)
    return frame


def plot_maximum_copy_scatter(
    results: pd.DataFrame, output_stem: str | Path
) -> None:
    """Plot only independently verified maxima of at least two copies."""
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharex=True)
    for axis, cap in zip(axes, (1800, 3000), strict=True):
        cap_rows = results.loc[results.fragment_limit_bp.astype(int).eq(cap)].copy()
        verified = pd.to_numeric(cap_rows.verified_max_copies, errors="coerce")
        plotted = cap_rows.loc[verified.ge(2)].copy()
        plotted["verified_max_copies"] = pd.to_numeric(
            plotted.verified_max_copies
        )
        for collection, color, marker in (
            ("Natural", _NATURAL, "o"),
            ("Designed", _DESIGNED, "^"),
        ):
            values = plotted.loc[plotted.collection.eq(collection)]
            axis.scatter(
                values.unit_length,
                values.verified_max_copies,
                color=color,
                marker=marker,
                alpha=0.72,
                s=30,
                edgecolors="white",
                linewidths=0.4,
                label=collection,
            )
        failures = int(verified.lt(2).sum() + verified.isna().sum())
        capacity_limited = int(
            cap_rows.get(
                "adaptive_maximum_proof_status",
                pd.Series("", index=cap_rows.index),
            )
            .astype(str)
            .eq("capacity_limit_reached")
            .sum()
        )
        axis.text(
            0.02,
            0.98,
            f"compatible inputs={len(cap_rows)}\nplotted={len(plotted)}\nIDT/GA failures={failures}\ncapacity-limited={capacity_limited}",
            transform=axis.transAxes,
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
        )
        axis.set_title(f"{cap:,}-bp fragment limit")
        axis.set_xlabel("Middle-module length (AA)")
        axis.set_ylabel("Maximum verified repeat copies")
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    figure.tight_layout()
    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def _open_audit_text(source: Path):
    if source.suffix == ".gz":
        return gzip.open(source, "rt", encoding="utf-8")
    return source.open("r", encoding="utf-8")


def _load_idt_audits(paths: Iterable[str | Path]) -> tuple[list[dict[str, Any]], set[str]]:
    records: dict[str, dict[str, Any]] = {}
    for source in (Path(path) for path in paths):
        if not source.is_file():
            raise FileNotFoundError(f"IDT audit file does not exist: {source}")
        with _open_audit_text(source) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(
                        f"Invalid IDT audit JSON at {source}:{line_number}"
                    ) from exc
                response_sha = str(
                    record.get("response_sha256")
                    or record.get("summary", {}).get("idt_response_sha256")
                    or ""
                )
                if len(response_sha) != 64:
                    raise ValueError(
                        f"IDT audit record lacks a SHA256 at {source}:{line_number}"
                    )
                records.setdefault(response_sha, record)
    return list(records.values()), set(records)


def _validate_adaptive_result(
    row: dict[str, Any], idt_response_hashes: set[str]
) -> tuple[bool, list[str]]:
    """Machine-check one Stage-2 row and return all failure reasons."""
    reasons: list[str] = []
    maximum = pd.to_numeric(
        pd.Series([row.get("verified_max_copies")]), errors="coerce"
    ).iloc[0]
    final_passed = row.get("final_passed")
    accepted = bool(
        pd.notna(maximum)
        and int(maximum) >= 2
        and pd.notna(final_passed)
        and bool(final_passed)
    )
    if not accepted:
        return False, reasons

    copies = int(maximum)
    unit = validate_protein_sequence(str(row.get("unit_sequence", "")))
    dna = str(row.get("dna_sequence") or "")
    if not dna or translate_dna(dna) != unit * copies:
        reasons.append("translation_mismatch")
    if len(dna) != 3 * len(unit) * copies:
        reasons.append("coding_length_mismatch")
    available = int(
        row.get(
            "available_coding_bp",
            int(row.get("fragment_limit_bp", 0))
            - int(row.get("external_deduction_bp", 0)),
        )
    )
    if len(dna) > available:
        reasons.append("fragment_capacity_exceeded")

    site_i = str(row.get("site_i_recognition_site") or "")
    site_ii = str(row.get("site_ii_recognition_site") or "")
    if not site_i or not site_ii or site_i == site_ii:
        reasons.append("invalid_selected_pair")
    else:
        selected_excess = max(0, recognition_site_count(dna, site_i) - 1)
        selected_excess += max(0, recognition_site_count(dna, site_ii))
        if selected_excess != 0:
            reasons.append("selected_pair_excess_site")
        reported_excess = pd.to_numeric(
            pd.Series([row.get("selected_pair_re_site_excess")]),
            errors="coerce",
        ).iloc[0]
        if pd.isna(reported_excess) or int(reported_excess) != selected_excess:
            reasons.append("selected_pair_excess_metric_mismatch")

    score = pd.to_numeric(
        pd.Series([row.get("idt_complexity_score")]), errors="coerce"
    ).iloc[0]
    if pd.isna(score) or float(score) >= 10.0:
        reasons.append("idt_score_not_below_10")
    if row.get("idt_score_policy") != "idt-rule-score-sum-lt10-v1":
        reasons.append("idt_policy_mismatch")
    sequence_sha = hashlib.sha256(dna.encode()).hexdigest() if dna else ""
    if row.get("idt_scored_sequence_sha256") != sequence_sha:
        reasons.append("idt_sequence_hash_mismatch")
    response_sha = str(row.get("idt_response_sha256") or "")
    if len(response_sha) != 64 or response_sha not in idt_response_hashes:
        reasons.append("idt_response_audit_missing")
    proof_status = row.get("adaptive_maximum_proof_status")
    if proof_status not in {
        "capacity_limit_reached",
        "next_copy_failed_at_100",
    }:
        reasons.append("maximum_proof_missing")
    else:
        try:
            trace = json.loads(str(row.get("adaptive_search_trace_json") or "[]"))
        except json.JSONDecodeError:
            trace = None
            reasons.append("maximum_proof_trace_invalid_json")
        if not isinstance(trace, list):
            trace = []
            if "maximum_proof_trace_invalid_json" not in reasons:
                reasons.append("maximum_proof_trace_not_list")
        if proof_status == "capacity_limit_reached":
            upper = pd.to_numeric(
                pd.Series(
                    [
                        row.get(
                            "adaptive_search_upper_bound_copies",
                            row.get("mathematical_max_copies"),
                        )
                    ]
                ),
                errors="coerce",
            ).iloc[0]
            if pd.isna(upper) or copies != int(upper):
                reasons.append("capacity_limit_proof_mismatch")
        elif proof_status == "next_copy_failed_at_100":
            next_copy_proof = any(
                isinstance(item, dict)
                and int(item.get("copies", -1)) == copies + 1
                and int(item.get("generations", -1)) == 100
                and item.get("passed") is False
                for item in trace
            )
            if not next_copy_proof:
                reasons.append("next_copy_100_generation_failure_missing")
    return True, reasons


def finalize_adaptive_copy_results(
    result_paths: Iterable[str | Path],
    compatibility_path: str | Path,
    output_dir: str | Path,
    *,
    idt_audit_paths: Iterable[str | Path],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Finalize Stage 2, expand traces, prove maxima, and write the summary."""
    result_files = [Path(path) for path in result_paths]
    if not result_files:
        raise ValueError("At least one adaptive result shard is required")
    results = pd.concat(
        [pd.read_parquet(path) for path in result_files], ignore_index=True
    )
    required = {
        "module_id",
        "collection",
        "unit_sequence",
        "unit_length",
        "fragment_limit_bp",
        "verified_max_copies",
        "adaptive_search_trace_json",
    }
    missing = sorted(required - set(results.columns))
    if missing:
        raise ValueError(f"Adaptive results are missing columns: {missing}")
    duplicated = results.duplicated(
        ["collection", "module_id", "fragment_limit_bp"], keep=False
    )
    if duplicated.any():
        raise ValueError("Duplicate module/capacity rows across adaptive shards")

    compatibility = _read_table(compatibility_path)
    base_key = ["collection", "module_id"]
    compatibility_required = {*base_key, "hurdler_compatible", "unit_sequence"}
    compatibility_missing = sorted(
        compatibility_required - set(compatibility.columns)
    )
    if compatibility_missing:
        raise ValueError(
            "Compatibility table is missing columns: "
            f"{compatibility_missing}"
        )
    if compatibility.duplicated(base_key).any():
        raise ValueError("Compatibility table contains duplicate module rows")

    fragment_limits = tuple(
        sorted(results.fragment_limit_bp.astype(int).unique().tolist())
    )
    if fragment_limits != (1800, 3000):
        raise ValueError(
            "Stage 2 requires exactly the 1800-bp and 3000-bp capacities; "
            f"observed {fragment_limits}"
        )
    compatible = compatibility.loc[
        compatibility.hurdler_compatible.astype(bool),
        [*base_key, "unit_sequence"],
    ].copy()
    expected_keys = {
        (str(row.collection), str(row.module_id), int(capacity))
        for row in compatible.itertuples(index=False)
        for capacity in fragment_limits
    }
    observed_keys = {
        (str(row.collection), str(row.module_id), int(row.fragment_limit_bp))
        for row in results.itertuples(index=False)
    }
    missing_keys = sorted(expected_keys - observed_keys)
    extra_keys = sorted(observed_keys - expected_keys)
    if missing_keys or extra_keys:
        raise ValueError(
            "Adaptive Stage-2 results are incomplete or out of scope: "
            f"expected={len(expected_keys)}, observed={len(observed_keys)}, "
            f"missing={len(missing_keys)} examples={missing_keys[:5]}, "
            f"extra={len(extra_keys)} examples={extra_keys[:5]}"
        )
    compatibility_sequences = compatible.rename(
        columns={"unit_sequence": "compatibility_unit_sequence"}
    )
    sequence_check = results.merge(
        compatibility_sequences,
        on=base_key,
        how="left",
        validate="many_to_one",
    )
    if not sequence_check.unit_sequence.astype(str).eq(
        sequence_check.compatibility_unit_sequence.astype(str)
    ).all():
        raise ValueError(
            "Stage-2 unit_sequence differs from the Stage-1 compatibility table"
        )

    audit_records, audit_hashes = _load_idt_audits(idt_audit_paths)
    validation_rows: list[dict[str, Any]] = []
    for row in results.to_dict(orient="records"):
        accepted, reasons = _validate_adaptive_result(row, audit_hashes)
        validation_rows.append(
            {
                "module_id": row["module_id"],
                "collection": row["collection"],
                "fragment_limit_bp": int(row["fragment_limit_bp"]),
                "accepted_repeat_construct": accepted,
                "validation_passed": accepted and not reasons,
                "validation_reasons_json": json.dumps(reasons),
            }
        )
    validation = pd.DataFrame(validation_rows)
    bad_accepted = validation.loc[
        validation.accepted_repeat_construct & ~validation.validation_passed
    ]
    if not bad_accepted.empty:
        raise ValueError(
            "Accepted adaptive constructs failed validation: "
            + bad_accepted.head(5).to_json(orient="records")
        )
    results = results.merge(
        validation,
        on=["module_id", "collection", "fragment_limit_bp"],
        how="left",
        validate="one_to_one",
    )

    trace_rows: list[dict[str, Any]] = []
    for result in results.to_dict(orient="records"):
        try:
            trace = json.loads(str(result["adaptive_search_trace_json"]))
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Invalid adaptive trace for {result['module_id']}"
            ) from exc
        if not isinstance(trace, list):
            raise ValueError(f"Adaptive trace is not a list for {result['module_id']}")
        for evaluation_index, item in enumerate(trace):
            if not isinstance(item, dict):
                raise ValueError(
                    f"Adaptive trace item is not an object for {result['module_id']}"
                )
            trace_rows.append(
                {
                    "module_id": result["module_id"],
                    "collection": result["collection"],
                    "unit_length": int(result["unit_length"]),
                    "fragment_limit_bp": int(result["fragment_limit_bp"]),
                    "evaluation_index": evaluation_index,
                    **item,
                }
            )
    traces = pd.DataFrame(trace_rows)

    value_columns = [
        "verified_max_copies",
        "mathematical_max_copies",
        "adaptive_maximum_proof_status",
        "adaptive_stop_reason",
        "optimization_status",
        "failure_reason",
        "final_ga_weights_json",
        "idt_status",
        "idt_complexity_score",
        "idt_score_policy",
        "idt_positive_score_names_json",
        "idt_rule_details_json",
        "idt_response_sha256",
        "selected_pair_re_site_excess",
        "dna_sequence",
        "dna_length",
        "validation_passed",
    ]
    present_values = [column for column in value_columns if column in results]
    wide_parts: list[pd.DataFrame] = []
    for cap in sorted(results.fragment_limit_bp.astype(int).unique()):
        cap_rows = results.loc[
            results.fragment_limit_bp.astype(int).eq(cap),
            [*base_key, *present_values],
        ].copy()
        cap_rows = cap_rows.rename(
            columns={column: f"cap_{cap}_{column}" for column in present_values}
        )
        wide_parts.append(cap_rows)
    summary = compatibility.copy()
    for part in wide_parts:
        summary = summary.merge(part, on=base_key, how="left", validate="one_to_one")

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    compact_results = results.drop(
        columns=["adaptive_search_trace_json"], errors="ignore"
    )
    for name, frame in (
        ("maximum_copy_results", compact_results),
        ("adaptive_copy_search_trace", traces),
        ("adaptive_copy_validation", validation),
        ("module_final_summary", summary),
    ):
        frame.to_parquet(destination / f"{name}.parquet", index=False)
    # The validation CSV is compact enough for manual review.  Large result,
    # trace, and summary CSV mirrors are deliberately not generated.
    validation.to_csv(destination / "adaptive_copy_validation.csv", index=False)

    audit_output = destination / "idt_audit_records.jsonl.gz"
    with audit_output.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            compresslevel=6,
            mtime=0,
        ) as compressed_handle:
            with io.TextIOWrapper(compressed_handle, encoding="utf-8") as handle:
                for record in sorted(
                    audit_records,
                    key=lambda item: str(item.get("response_sha256", "")),
                ):
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
    with (destination / "optimized_constructs.fasta").open("w") as dna_handle, (
        destination / "optimized_constructs.protein.fasta"
    ).open("w") as protein_handle:
        accepted = results.loc[results.validation_passed.astype(bool)]
        for row in accepted.itertuples(index=False):
            header = (
                f"{row.collection}|{row.module_id}|cap={int(row.fragment_limit_bp)}|"
                f"copies={int(row.verified_max_copies)}"
            )
            dna = str(row.dna_sequence)
            protein = str(row.unit_sequence) * int(row.verified_max_copies)
            dna_handle.write(f">{header}\n")
            protein_handle.write(f">{header}\n")
            for start in range(0, len(dna), 80):
                dna_handle.write(dna[start : start + 80] + "\n")
            for start in range(0, len(protein), 80):
                protein_handle.write(protein[start : start + 80] + "\n")

    plot_maximum_copy_scatter(
        results, destination / "maximum_verified_repeat_copies"
    )
    manifest = {
        "corpus_version": CORPUS_VERSION,
        "result_shards": [str(path.resolve()) for path in result_files],
        "compatibility": str(Path(compatibility_path).resolve()),
        "result_rows": len(results),
        "trace_rows": len(traces),
        "compatibility_rows": len(compatibility),
        "accepted_rows": int(validation.validation_passed.sum()),
        "idt_audit_records": len(audit_records),
        "artifacts": {
            "maximum_copy_results": "maximum_copy_results.parquet",
            "adaptive_copy_search_trace": "adaptive_copy_search_trace.parquet",
            "adaptive_copy_validation": "adaptive_copy_validation.parquet",
            "module_final_summary": "module_final_summary.parquet",
            "idt_audit": "idt_audit_records.jsonl.gz",
        },
        "status_counts": results.get(
            "optimization_status", pd.Series(dtype=str)
        )
        .fillna("missing")
        .astype(str)
        .value_counts(dropna=False)
        .to_dict(),
    }
    (destination / "adaptive_copy_search.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n"
    )
    return results, traces, summary
