"""Backend for V2 notebook 08: the validated 1--60 AA success landscape."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ..constants import PLASMIDS
from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from ..short_screen import SHORT_MOTIF_COUNTS
from .common import BackendSpec, load_table, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "08_success_landscape_analysis",
    "Success landscape analysis",
    "Validate exhaustive 1-5 AA and seeded 6-60 AA results and regenerate the historical square 1-50 AA plot.",
    production_workflows=("success-landscape",),
    default_request={"max_plot_length": 50},
)
UW_PURPLE = "#4B2E83"
UW_GOLD = "#B7A57A"
PLASMID_COLORS = ["#4B2E83", "#85754D", "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9"]


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _source(request: Mapping[str, Any]) -> Path:
    return Path(str(request.get("summary", repo_root() / "data/results/success_landscape_compact_v2.parquet")))


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    source = _source(request)
    if not source.is_file():
        raise FileNotFoundError(
            f"{source}; generate the success-landscape production bundle with notebook 07"
        )
    frame = load_table(source)
    required = {"module_length", "plasmid", "successes", "tests", "success_rate", "method"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Success summary is missing columns: {missing}")
    return {"status": "passed", "rows": len(frame)}


def _validate(frame: pd.DataFrame) -> dict[str, Any]:
    for length, expected in SHORT_MOTIF_COUNTS.items():
        rows = frame.loc[frame.module_length.eq(length)]
        if len(rows) != len(PLASMIDS):
            raise ValueError(f"Expected eight plasmid rows for {length}AA")
        if set(rows.tests.astype(int)) != {expected}:
            raise ValueError(f"{length}AA denominator is not exactly 20**{length}")
        if not rows.method.astype(str).str.contains("exhaustive").all():
            raise ValueError(f"{length}AA must use exhaustive enumeration")
    random = frame.loc[frame.module_length.between(6, 60)]
    if random.empty or random.groupby("module_length").plasmid.nunique().min() != len(PLASMIDS):
        raise ValueError("6-60AA random results are incomplete")
    if not random.loc[random.module_length.ge(7), "seed"].fillna(42).astype(int).eq(42).all():
        raise ValueError("7-60AA rows do not preserve seed 42")
    return {
        "short_total_motifs": sum(SHORT_MOTIF_COUNTS.values()),
        "length_min": int(frame.module_length.min()),
        "length_max": int(frame.module_length.max()),
        "plasmid_count": int(frame.plasmid.nunique()),
    }


def _plot(frame: pd.DataFrame, output: Path, max_length: int) -> list[Path]:
    view = frame.loc[frame.module_length.between(1, max_length)].copy()
    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    for color, plasmid in zip(PLASMID_COLORS, PLASMIDS, strict=True):
        rows = view.loc[view.plasmid.eq(plasmid)].sort_values("module_length")
        ax.plot(rows.module_length, rows.success_rate * 100, marker="o", markersize=2.4,
                linewidth=1.35, label=plasmid, color=color)
    ax.set_title("3-mer Probability vs Sequence Length", fontsize=13)
    ax.set_xlabel("Module length (AA)")
    ax.set_ylabel("HURDLER success rate (%)")
    ax.set_xlim(1, max_length)
    ax.set_ylim(0, 100)
    ax.grid(axis="y", color="#D1D5DB", linewidth=0.5, alpha=0.65)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="lower right")
    fig.tight_layout()
    output.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for suffix in ("png", "pdf", "svg"):
        path = output / f"success_landscape_1_{max_length}.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    preflight(context, request)
    frame = load_table(_source(request))
    metrics = _validate(frame)
    table_paths = write_frame(frame, context.directory("tables") / "success_landscape_1_60")
    figure_paths = _plot(frame, context.directory("figures"), int(request.get("max_plot_length", 50)))
    comparison_paths: list[Path] = []
    if request.get("three_copy_summary"):
        three = load_table(str(request["three_copy_summary"]))
        keys = ["module_length", "plasmid"]
        comparison = frame.merge(three, on=keys, suffixes=("_2x", "_3x"))
        comparison["improvement"] = comparison.success_rate_3x - comparison.success_rate_2x
        comparison_paths = write_frame(comparison, context.directory("tables") / "success_landscape_2x_vs_3x")
        metrics["improved_module_lengths"] = sorted(
            comparison.loc[comparison.improvement.gt(0), "module_length"].astype(int).unique().tolist()
        )
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[*table_paths, *comparison_paths, *figure_paths],
        metrics=metrics,
        next_notebooks=["11_reproducibility"],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
