"""Backend for V2 notebook 09: module, maximum-copy, GA/codon and IDT analysis."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "09_module_result_analysis",
    "Repeat-module, GA, codon and IDT analysis",
    "Rebuild Natural/Designed compatibility, 3-mer and maximum-copy figures from compact production tables.",
    production_workflows=("module-stage1", "module-stage2"),
)
COLORS = {"Natural": "#4B2E83", "Designed": "#E69F00"}


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _module_path(request: Mapping[str, Any]) -> Path:
    return Path(str(request.get("module_results", repo_root() / "data/results/module_analysis_compact_v2.parquet")))


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    source = _module_path(request)
    three = Path(str(request.get("three_mer_results", repo_root() / "data/results/repeatsdb_designed_hurdler_3mer_results.csv")))
    codon = Path(str(request.get(
        "codon_benchmark",
        repo_root() / "codon_opt_benchmark_extended/results/well_color_iterations.csv",
    )))
    for path in (source, three, codon):
        if not path.is_file():
            raise FileNotFoundError(path)
    return {
        "status": "passed", "module_results": str(source),
        "three_mer_results": str(three), "codon_benchmark": str(codon),
    }


def _read_modules(path: Path) -> pd.DataFrame:
    columns = [
        "module_id", "collection", "middle_module_sequence_aa", "middle_module_length_aa",
        "hurdler_compatible", "selected_plasmid", "selected_site_i_enzyme", "selected_site_ii_enzyme",
        "cap1800_maximum_verified_copies", "cap1800_result_status", "cap1800_final_ga_weights",
        "cap1800_mathematical_max_copies", "cap1800_stop_reason", "cap1800_failure_reason",
        "cap1800_idt_status", "cap1800_idt_score_sum", "cap1800_idt_rule_reasons", "cap3000_maximum_verified_copies",
        "cap3000_result_status", "cap3000_final_ga_weights", "cap3000_idt_score_sum",
        "cap3000_mathematical_max_copies", "cap3000_stop_reason", "cap3000_failure_reason",
        "cap3000_idt_status", "cap3000_idt_rule_reasons", "middle_module_number_one_based", "boundary_method",
    ]
    if path.suffix.lower() == ".parquet":
        frame = pd.read_parquet(path)
        return frame[[value for value in columns if value in frame.columns]]
    return pd.read_csv(path, usecols=lambda value: value in columns, low_memory=False)


def _compatibility_summary(frame: pd.DataFrame) -> pd.DataFrame:
    maximum = int(np.ceil(frame.middle_module_length_aa.max() / 10.0) * 10)
    bins = list(range(0, maximum + 10, 10))
    labels = [f"{start + 1}-{start + 10}" for start in bins[:-1]]
    work = frame.copy()
    work["length_bin"] = pd.cut(work.middle_module_length_aa, bins=bins, labels=labels, include_lowest=True)
    rows = (
        work.groupby(["collection", "length_bin", "hurdler_compatible"], observed=False)
        .size().rename("count").reset_index()
    )
    totals = rows.groupby(["collection", "length_bin"], observed=False)["count"].transform("sum")
    rows["fraction"] = np.where(totals > 0, rows["count"] / totals, 0.0)
    rows["total_n"] = totals
    return rows


def _plot_compatibility(summary: pd.DataFrame, output: Path) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.8), sharex="col")
    for column, collection in enumerate(("Natural", "Designed")):
        rows = summary.loc[summary.collection.eq(collection)].copy()
        labels = [str(value) for value in rows.length_bin.drop_duplicates()]
        x = np.arange(len(labels))
        incompatible = rows.loc[~rows.hurdler_compatible.astype(bool)].set_index("length_bin")
        compatible = rows.loc[rows.hurdler_compatible.astype(bool)].set_index("length_bin")
        no = np.array([incompatible["count"].get(label, 0) for label in labels])
        yes = np.array([compatible["count"].get(label, 0) for label in labels])
        total = no + yes
        axes[0, column].bar(x, no, color="#9CA3AF", label="Incompatible")
        axes[0, column].bar(x, yes, bottom=no, color=COLORS[collection], label="Compatible")
        no_fraction = np.divide(no, total, out=np.zeros_like(no, dtype=float), where=total > 0)
        yes_fraction = np.divide(yes, total, out=np.zeros_like(yes, dtype=float), where=total > 0)
        axes[1, column].bar(x, no_fraction, color="#9CA3AF")
        axes[1, column].bar(x, yes_fraction, bottom=no_fraction, color=COLORS[collection])
        axes[0, column].set_title(collection)
        axes[1, column].set_xticks(x, labels, rotation=90, fontsize=6)
        for index, value in enumerate(total):
            if value:
                axes[0, column].text(index, value, f"n={value}", ha="center", va="bottom", fontsize=6)
                if no[index]:
                    axes[1, column].text(
                        index, no_fraction[index] / 2,
                        f"{no[index]}\n{no_fraction[index] * 100:.0f}%",
                        ha="center", va="center", fontsize=4.8,
                    )
                if yes[index]:
                    axes[1, column].text(
                        index, no_fraction[index] + yes_fraction[index] / 2,
                        f"{yes[index]}\n{yes_fraction[index] * 100:.0f}%",
                        ha="center", va="center", fontsize=4.8, color="white",
                    )
    axes[0, 0].set_ylabel("Modules")
    axes[1, 0].set_ylabel("Proportion")
    axes[1, 0].set_ylim(0, 1)
    axes[1, 1].set_ylim(0, 1)
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    return _save(fig, output, "module_compatibility_by_length")


def _save(fig: Any, output: Path, stem: str) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix in ("png", "pdf", "svg"):
        path = output / f"{stem}.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def _plot_scatter(frame: pd.DataFrame, output: Path) -> list[Path]:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.1), sharex=True)
    for ax, capacity in zip(axes, (1800, 3000), strict=True):
        column = f"cap{capacity}_maximum_verified_copies"
        for collection in ("Natural", "Designed"):
            rows = frame.loc[
                frame.collection.eq(collection)
                & frame.hurdler_compatible.astype(bool)
                & pd.to_numeric(frame[column], errors="coerce").ge(2)
            ]
            ax.scatter(rows.middle_module_length_aa, rows[column], s=9, alpha=0.55,
                       color=COLORS[collection], label=collection)
        ax.set_title(f"{capacity:,} bp capacity")
        ax.set_xlabel("Middle-module length (AA)")
        ax.set_ylabel("Maximum verified module copies")
        ax.grid(color="#D1D5DB", linewidth=0.4, alpha=0.6)
        compatible = int(frame.hurdler_compatible.astype(bool).sum())
        plotted = int(pd.to_numeric(frame[column], errors="coerce").ge(2).sum())
        failed = compatible - plotted
        mathematical = pd.to_numeric(
            frame[f"cap{capacity}_mathematical_max_copies"], errors="coerce"
        )
        verified = pd.to_numeric(frame[column], errors="coerce")
        capacity_limited = int((verified.ge(2) & verified.eq(mathematical)).sum())
        ax.text(
            0.98, 0.98,
            f"compatible={compatible:,}\nplotted={plotted:,}\nGA/IDT failed={failed:,}\ncapacity-limited={capacity_limited:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5,
            bbox={"facecolor": "white", "edgecolor": "#D1D5DB", "alpha": 0.9},
        )
    axes[0].legend(frameon=False)
    fig.tight_layout()
    return _save(fig, output, "module_length_vs_maximum_verified_copies")


def _parse_weights(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for capacity in (1800, 3000):
        column = f"cap{capacity}_final_ga_weights"
        status_column = f"cap{capacity}_result_status"
        for record in frame[["module_id", "collection", column, status_column]].itertuples(index=False, name=None):
            module_id, collection, encoded, status = record
            if not isinstance(encoded, str) or not encoded or encoded == "None":
                continue
            for field in encoded.split(";"):
                name, separator, value = field.partition("=")
                if not separator:
                    continue
                try:
                    numeric = float(value)
                except ValueError:
                    continue
                rows.append({
                    "module_id": module_id, "collection": collection,
                    "capacity_bp": capacity, "result_status": status,
                    "weight_name": name, "weight_value": numeric,
                })
    return pd.DataFrame(rows)


def _parse_idt_reasons(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for capacity in (1800, 3000):
        column = f"cap{capacity}_idt_rule_reasons"
        status_column = f"cap{capacity}_idt_status"
        for reasons, status in frame[[column, status_column]].itertuples(index=False, name=None):
            if not isinstance(reasons, str) or not reasons or reasons == "None":
                continue
            for encoded in reasons.split("|"):
                rule_name = encoded.split(":score=", 1)[0].strip()
                if rule_name:
                    rows.append({"capacity_bp": capacity, "idt_status": status, "rule_name": rule_name})
    if not rows:
        return pd.DataFrame(columns=["capacity_bp", "idt_status", "rule_name", "count"])
    return (
        pd.DataFrame(rows).groupby(["capacity_bp", "idt_status", "rule_name"])
        .size().rename("count").reset_index()
    )


def _ga_weight_summary(weights: pd.DataFrame) -> pd.DataFrame:
    if weights.empty:
        return pd.DataFrame()
    return (
        weights.groupby(["capacity_bp", "result_status", "weight_name"])
        .weight_value.agg(n="size", minimum="min", q25=lambda x: x.quantile(.25), median="median", q75=lambda x: x.quantile(.75), maximum="max")
        .reset_index()
    )


def _status_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for capacity in (1800, 3000):
        grouped = frame.groupby(["collection", f"cap{capacity}_result_status"]).size()
        for (collection, status), count in grouped.items():
            rows.append({"capacity_bp": capacity, "collection": collection, "result_status": status, "count": int(count)})
    return pd.DataFrame(rows)


def _plot_optimization_qc(
    frame: pd.DataFrame,
    reasons: pd.DataFrame,
    weight_summary: pd.DataFrame,
    output: Path,
) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.8))
    status = _status_summary(frame)
    statuses = list(status.result_status.drop_duplicates())
    x = np.arange(4)
    labels = ["Natural\n1,800", "Designed\n1,800", "Natural\n3,000", "Designed\n3,000"]
    bottoms = np.zeros(4)
    status_colors = ["#4B2E83", "#B7A57A", "#9CA3AF", "#D55E00"]
    for status_name, color in zip(statuses, status_colors, strict=False):
        values = []
        for capacity, collection in ((1800, "Natural"), (1800, "Designed"), (3000, "Natural"), (3000, "Designed")):
            matched = status.loc[
                status.capacity_bp.eq(capacity) & status.collection.eq(collection) & status.result_status.eq(status_name), "count"
            ]
            values.append(int(matched.iloc[0]) if len(matched) else 0)
        axes[0, 0].bar(x, values, bottom=bottoms, color=color, label=str(status_name))
        bottoms += np.asarray(values)
    axes[0, 0].set_xticks(x, labels)
    axes[0, 0].set_ylabel("Modules")
    axes[0, 0].set_title("A  Optimization outcomes", loc="left", weight="bold")
    axes[0, 0].legend(frameon=False, fontsize=6)

    for capacity, color in ((1800, "#4B2E83"), (3000, "#E69F00")):
        values = pd.to_numeric(frame[f"cap{capacity}_idt_score_sum"], errors="coerce").dropna()
        axes[0, 1].hist(values, bins=np.linspace(0, 10, 26), histtype="step", linewidth=1.5, color=color, label=f"{capacity:,} bp")
    axes[0, 1].axvline(10, color="#D55E00", linewidth=1)
    axes[0, 1].set_xlabel("Final IDT rule-score sum")
    axes[0, 1].set_ylabel("Accepted modules")
    axes[0, 1].set_title("B  Final IDT scores (<10)", loc="left", weight="bold")
    axes[0, 1].legend(frameon=False)

    top_reasons = reasons.groupby("rule_name")["count"].sum().nlargest(10).index if not reasons.empty else []
    reason_values = reasons.loc[reasons.rule_name.isin(top_reasons)].groupby("rule_name")["count"].sum().sort_values()
    axes[1, 0].barh(reason_values.index, reason_values.values, color="#0072B2")
    axes[1, 0].set_xlabel("Recorded positive-rule occurrences")
    axes[1, 0].set_title("C  Most common IDT rule feedback", loc="left", weight="bold")
    axes[1, 0].tick_params(axis="y", labelsize=6)

    selected = weight_summary.loc[
        weight_summary.weight_name.isin(["repeated_re_site_excess", "repeated_8mer", "repeated_13mer", "repeated_14mer", "gc_window_soft_violation", "hairpin_10mer_proxy"])
        & weight_summary.result_status.eq("idt_accepted")
    ]
    for capacity, color in ((1800, "#4B2E83"), (3000, "#E69F00")):
        local = selected.loc[selected.capacity_bp.eq(capacity)].sort_values("weight_name")
        axes[1, 1].plot(local.weight_name, local["median"], marker="o", color=color, label=f"{capacity:,} bp")
    axes[1, 1].set_yscale("symlog", linthresh=1)
    axes[1, 1].tick_params(axis="x", rotation=35, labelsize=6)
    axes[1, 1].set_ylabel("Median final GA weight")
    axes[1, 1].set_title("D  Adaptive GA weights", loc="left", weight="bold")
    axes[1, 1].legend(frameon=False)
    fig.tight_layout()
    return _save(fig, output, "module_ga_idt_qc")


def _plot_three_mer(three: pd.DataFrame, output: Path) -> list[Path]:
    compatible = three.loc[three.hurdler_compatible.astype(bool)].copy()
    if compatible[["selected_re_pair", "site_i_3mer_aa", "site_ii_3mer_aa"]].isna().any().any():
        raise ValueError("A compatible module lacks its selected pair or Site-I/II 3-mer")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.2))
    pairs = compatible.groupby(["collection", "selected_re_pair"]).size().rename("count").reset_index()
    top_pairs = pairs.groupby("selected_re_pair")["count"].sum().nlargest(12).index
    pivot = pairs.loc[pairs.selected_re_pair.isin(top_pairs)].pivot(index="selected_re_pair", columns="collection", values="count").fillna(0)
    pivot = pivot.reindex(pivot.sum(axis=1).sort_values().index)
    pivot.plot.barh(ax=axes[0], color=[COLORS.get(value, "#777777") for value in pivot.columns])
    axes[0].set_xlabel("Compatible modules")
    axes[0].set_ylabel("Selected Site-I/Site-II pair")
    axes[0].set_title("A  Most-used RE pairs", loc="left", weight="bold")
    axes[0].legend(frameon=False)
    triplets = compatible.groupby(["site_i_3mer_aa", "site_ii_3mer_aa"]).size().unstack(fill_value=0)
    common_i = triplets.sum(axis=1).nlargest(15).index
    common_ii = triplets.sum(axis=0).nlargest(15).index
    image = axes[1].imshow(np.log1p(triplets.loc[common_i, common_ii]), cmap="Purples", aspect="auto")
    axes[1].set_xticks(range(len(common_ii)), common_ii, rotation=90, fontsize=6)
    axes[1].set_yticks(range(len(common_i)), common_i, fontsize=6)
    axes[1].set_xlabel("Site-II 3-mer AA")
    axes[1].set_ylabel("Site-I 3-mer AA")
    axes[1].set_title("B  Selected 3-mer combinations", loc="left", weight="bold")
    fig.colorbar(image, ax=axes[1], label="log(1 + modules)")
    fig.tight_layout()
    return _save(fig, output, "module_selected_re_pair_and_3mer")


def _plot_codon_benchmark(codon: pd.DataFrame, output: Path) -> list[Path]:
    color_map = {"red": "#D55E00", "yellow": "#E69F00", "green": "#009E73"}
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    ax.scatter(codon.iteration, np.ones(len(codon)), c=codon.color.map(color_map).fillna("#9CA3AF"), s=18)
    ax.set_yticks([])
    ax.set_xlabel("Historical codon-optimization iteration")
    ax.set_title("Committed codon benchmark well classification", loc="left", weight="bold")
    ax.set_xlim(codon.iteration.min() - 1, codon.iteration.max() + 1)
    fig.tight_layout()
    return _save(fig, output, "historical_codon_benchmark_well_classes")


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    inputs = preflight(context, request)
    modules = _read_modules(Path(inputs["module_results"]))
    three = pd.read_csv(inputs["three_mer_results"])
    codon = pd.read_csv(inputs["codon_benchmark"])
    if len(modules) != 26_095:
        raise ValueError(f"Expected 26,095 active module rows, observed {len(modules)}")
    counts = modules.collection.value_counts().to_dict()
    if counts != {"Natural": 25_913, "Designed": 182}:
        raise ValueError(f"Unexpected active collection counts: {counts}")
    if len(three) != len(modules):
        raise ValueError("3-mer and module result row counts differ")
    if three.sequence_id.duplicated().any():
        raise ValueError("3-mer result contains duplicate sequence IDs")
    summary = _compatibility_summary(modules)
    summary_paths = write_frame(summary, context.directory("tables") / "module_compatibility_bins")
    compact_columns = [column for column in modules.columns if "accepted_dna" not in column]
    compact_paths = write_frame(modules[compact_columns], context.directory("tables") / "module_analysis_compact")
    weights = _parse_weights(modules)
    weight_summary = _ga_weight_summary(weights)
    reasons = _parse_idt_reasons(modules)
    status = _status_summary(modules)
    analysis_paths = [
        *write_frame(weight_summary, context.directory("tables") / "module_ga_weight_summary"),
        *write_frame(reasons, context.directory("tables") / "module_idt_rule_summary"),
        *write_frame(status, context.directory("tables") / "module_optimization_status_summary"),
        *write_frame(codon, context.directory("tables") / "historical_codon_benchmark_well_classes"),
    ]
    figures = [
        *_plot_compatibility(summary, context.directory("figures")),
        *_plot_scatter(modules, context.directory("figures")),
        *_plot_optimization_qc(modules, reasons, weight_summary, context.directory("figures")),
        *_plot_three_mer(three, context.directory("figures")),
        *_plot_codon_benchmark(codon, context.directory("figures")),
    ]
    compatible = int(modules.hurdler_compatible.astype(bool).sum())
    plotted_1800 = int(pd.to_numeric(modules.cap1800_maximum_verified_copies, errors="coerce").ge(2).sum())
    plotted_3000 = int(pd.to_numeric(modules.cap3000_maximum_verified_copies, errors="coerce").ge(2).sum())
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[*summary_paths, *compact_paths, *analysis_paths, *figures],
        metrics={
            "module_rows": len(modules), "natural_rows": counts["Natural"], "designed_rows": counts["Designed"],
            "hurdler_compatible": compatible, "hurdler_compatible_fraction": compatible / len(modules),
            "scatter_1800_rows": plotted_1800, "scatter_3000_rows": plotted_3000,
            "ga_weight_records": len(weights), "idt_rule_groups": len(reasons),
            "codon_benchmark_iterations": len(codon),
        },
        next_notebooks=["11_reproducibility"],
        limitations=[
            "Compact production tables retain final GA weights and final IDT feedback; per-attempt trajectories require the optional raw production trace.",
            "The historical codon benchmark repository snapshot retains well classifications and published PDFs, but its gitignored raw per-run workbooks are unavailable for re-aggregation.",
        ],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
