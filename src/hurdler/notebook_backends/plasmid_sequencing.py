"""Backend for V2 notebook 14: plasmid sequencing QC."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any, Mapping

import matplotlib
import numpy as np
import pandas as pd
from Bio import SeqIO

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "14_plasmid_sequencing",
    "Plasmid sequencing QC",
    "Summarize retained per-position sequencing statistics and verify reference sequence annotations.",
    production_workflows=("reports",),
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _folder(request: Mapping[str, Any]) -> Path:
    return Path(str(request.get("construct_dir", repo_root() / "plasmid_sequencing_result/Na4M13A")))


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    folder = _folder(request)
    stats = Path(str(request.get("stats_csv", next(iter(sorted(folder.glob("*_stats.csv"))), ""))))
    genbank = Path(str(request.get("reference_genbank", next(iter(sorted(folder.glob("*.gbk"))), ""))))
    if not stats.is_file() or not genbank.is_file():
        raise FileNotFoundError(f"Construct folder requires *_stats.csv and *.gbk: {folder}")
    return {"status": "passed", "construct_dir": str(folder), "stats_csv": str(stats), "reference_genbank": str(genbank)}


def _read_all_positions(path: Path) -> pd.DataFrame:
    text = path.read_text(errors="replace")
    marker = "All positions:\n"
    if marker not in text:
        return pd.read_csv(path)
    return pd.read_csv(io.StringIO(text.split(marker, 1)[1]))


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    info = preflight(context, request)
    positions = _read_all_positions(Path(info["stats_csv"]))
    reference = SeqIO.read(info["reference_genbank"], "genbank")
    required = {"pos", "reads_all", "matches", "mismatches", "deletions", "insertions"}
    missing = sorted(required - set(positions.columns))
    if missing:
        raise ValueError(f"Sequencing stats are missing columns: {missing}")
    positions["coverage"] = positions.reads_all.astype(float)
    positions["mismatch_fraction"] = positions.mismatches / positions.reads_all.replace(0, np.nan)
    positions["deletion_fraction"] = positions.deletions / positions.reads_all.replace(0, np.nan)
    positions["insertion_fraction"] = positions.insertions / positions.reads_all.replace(0, np.nan)
    table_paths = write_frame(positions, context.directory("tables") / "plasmid_sequencing_positions")
    summary = pd.DataFrame([{
        "construct": Path(info["construct_dir"]).name,
        "reference_length_bp": len(reference.seq),
        "covered_positions": int(positions.coverage.gt(0).sum()),
        "mean_coverage": float(positions.coverage.mean()),
        "median_coverage": float(positions.coverage.median()),
        "positions_mismatch_gt_5pct": int(positions.mismatch_fraction.gt(0.05).sum()),
        "cds_feature_count": sum(feature.type == "CDS" for feature in reference.features),
        "restriction_site_feature_count": sum("restriction" in feature.type.lower() for feature in reference.features),
    }])
    summary_paths = write_frame(summary, context.directory("tables") / "plasmid_sequencing_summary")
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.5), sharex=True)
    axes[0].plot(positions.pos, positions.coverage, color="#4B2E83", linewidth=0.8)
    axes[0].set_ylabel("Coverage")
    axes[1].plot(positions.pos, positions.mismatch_fraction * 100, color="#D55E00", linewidth=0.7)
    axes[1].set_ylabel("Mismatch (%)")
    axes[1].set_xlabel("Reference position (bp)")
    fig.suptitle(f"Plasmid sequencing QC: {summary.iloc[0].construct}")
    fig.tight_layout()
    figure_paths = []
    for suffix in ("png", "pdf", "svg"):
        path = context.directory("figures") / f"plasmid_sequencing_qc.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white")
        figure_paths.append(path)
    plt.close(fig)
    copied_gb = context.directory("genbank") / "annotated_reference.gbk"
    copied_gb.write_bytes(Path(info["reference_genbank"]).read_bytes())
    report = context.directory("reports") / "plasmid_sequencing_qc.html"
    report.write_text("<!doctype html><meta charset='utf-8'><h1>Plasmid sequencing QC</h1>" + summary.to_html(index=False))
    return result_from_paths(context, backend_id=SPEC.notebook_id, request=request,
                             paths=[*table_paths, *summary_paths, *figure_paths, copied_gb, report],
                             metrics=summary.iloc[0].to_dict(), limitations=["Large raw-read realignment should be exported as a production report workflow."])


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
