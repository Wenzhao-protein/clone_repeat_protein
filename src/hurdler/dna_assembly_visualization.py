"""Publication figures for complete-route regulatory-array production runs."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .complete_route import TARGET_COPY_COUNTS
from .dna_assembly import DNA_COMPLETE_ROUTE_VERSION
from .io import sha256_file, utc_now


UW_PURPLE = "#4B2E83"
UW_GOLD = "#B7A57A"
BLUE = "#2D7DD2"
ORANGE = "#F45D01"
GREEN = "#2A9D8F"
RED = "#C44E52"
GREY = "#8C8C8C"
SOURCE_ORDER = ["CRISPRCasdb", "Rfam", "Ribocentre_Aptamer"]
SOURCE_LABELS = {
    "CRISPRCasdb": "CRISPRCasdb",
    "Rfam": "Rfam",
    "Ribocentre_Aptamer": "Ribocentre",
}
SOURCE_COLORS = {
    "CRISPRCasdb": UW_PURPLE,
    "Rfam": ORANGE,
    "Ribocentre_Aptamer": BLUE,
}
LENGTH_BINS = [0, 20, 30, 40, 60, 100, 200, np.inf]
LENGTH_LABELS = ["<20", "20–29", "30–39", "40–59", "60–99", "100–199", "≥200"]


def _prepare_style() -> None:
    import matplotlib as mpl
    import seaborn as sns

    sns.set_theme(style="white", context="paper")
    mpl.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8,
            "figure.titlesize": 13,
            "axes.linewidth": 0.8,
            "savefig.facecolor": "white",
        }
    )


def _save(figure, output_dir: str | Path, stem: str) -> list[Path]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for suffix in ("png", "pdf", "svg"):
        path = destination / f"{stem}.{suffix}"
        figure.savefig(path, dpi=300, facecolor="white", bbox_inches="tight")
        outputs.append(path)
    return outputs


def _source_values(frame: pd.DataFrame) -> list[str]:
    observed = set(frame.source_database.dropna().astype(str))
    return [source for source in SOURCE_ORDER if source in observed]


def plot_public_element_outcomes(
    targets: pd.DataFrame,
    elements: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Figure 1: corpus size and strict five-target element outcomes."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    _prepare_style()
    sources = _source_values(elements)
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.8), facecolor="white")

    counts = elements.groupby("source_database").size().reindex(sources)
    axes[0, 0].barh(
        [SOURCE_LABELS.get(value, value) for value in sources],
        counts.values,
        color=[SOURCE_COLORS[value] for value in sources],
    )
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_xlabel("Unique elements (log scale)")
    axes[0, 0].set_title("A  Public source coverage", loc="left", weight="bold")
    for y_value, count in enumerate(counts):
        axes[0, 0].text(count * 1.08, y_value, f"n={count:,}", va="center", fontsize=8)

    distribution = (
        elements.groupby(["source_database", "successful_target_count"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=sources, columns=range(6), fill_value=0)
    )
    palette = ["#D9D9D9", "#C7B9D9", "#AB91C8", "#8D68B4", "#6C4598", UW_PURPLE]
    bottoms = np.zeros(len(sources))
    totals = distribution.sum(axis=1).to_numpy()
    for success_count in range(6):
        values = distribution[success_count].to_numpy()
        fractions = np.divide(values, totals, out=np.zeros_like(values, dtype=float), where=totals > 0)
        axes[0, 1].bar(
            range(len(sources)), fractions, bottom=bottoms,
            color=palette[success_count], label=str(success_count), width=0.72,
        )
        for x_value, (bottom, fraction, count) in enumerate(zip(bottoms, fractions, values)):
            if count:
                axes[0, 1].text(
                    x_value, bottom + fraction / 2,
                    f"{count:,}\n{fraction*100:.1f}%",
                    ha="center", va="center", fontsize=5.4,
                    color="white" if success_count >= 4 else "black",
                )
        bottoms += fractions
    axes[0, 1].set_xticks(range(len(sources)), [SOURCE_LABELS.get(v, v) for v in sources], rotation=18)
    axes[0, 1].set_ylim(0, 1)
    axes[0, 1].set_ylabel("Fraction of unique elements")
    axes[0, 1].set_title("B  Number of exact target lengths completed", loc="left", weight="bold")
    axes[0, 1].legend(title="Passed copy levels", ncol=3, frameon=False, loc="upper center")
    for x_value, total in enumerate(totals):
        axes[0, 1].text(x_value, 1.02, f"n={total:,}", ha="center", va="bottom", fontsize=7)

    rates = (
        targets.groupby(["source_database", "target_copy_count"])
        .complete_route_verified.agg(["sum", "count"])
        .reset_index()
    )
    rates["percent"] = 100 * rates["sum"] / rates["count"]
    for source in sources:
        subset = rates.loc[rates.source_database.eq(source)].sort_values("target_copy_count")
        axes[1, 0].plot(
            subset.target_copy_count, subset.percent,
            marker="o", linewidth=1.7, markersize=4.5,
            color=SOURCE_COLORS[source], label=SOURCE_LABELS.get(source, source),
        )
        for row in subset.itertuples(index=False):
            axes[1, 0].annotate(
                f"{row.percent:.1f}%",
                (row.target_copy_count, row.percent),
                xytext=(0, 5), textcoords="offset points", ha="center", fontsize=6.2,
            )
    axes[1, 0].set_xticks(TARGET_COPY_COUNTS)
    axes[1, 0].set_ylim(0, 103)
    axes[1, 0].set_xlabel("Exact target copy number")
    axes[1, 0].set_ylabel("Complete-route success (%)")
    axes[1, 0].set_title("C  Exact-target success", loc="left", weight="bold")
    axes[1, 0].legend(frameon=False)

    length_frame = elements.copy()
    length_frame["unit_length_bin"] = pd.cut(
        length_frame.unit_length_bp,
        bins=LENGTH_BINS,
        labels=LENGTH_LABELS,
        right=False,
    )
    length_summary = (
        length_frame.groupby("unit_length_bin", observed=False)
        .all_five_complete.agg(["sum", "count"])
        .reindex(LENGTH_LABELS, fill_value=0)
    )
    passed_values = length_summary["sum"].to_numpy()
    failed_values = (length_summary["count"] - length_summary["sum"]).to_numpy()
    axes[1, 1].bar(range(len(LENGTH_LABELS)), passed_values, color=UW_PURPLE, label="All five complete")
    axes[1, 1].bar(
        range(len(LENGTH_LABELS)), failed_values,
        bottom=passed_values, color=UW_GOLD, label="Not all five",
    )
    for index, total in enumerate(length_summary["count"]):
        axes[1, 1].text(index, total, f"n={int(total):,}", ha="center", va="bottom", fontsize=6, rotation=90)
        if total:
            axes[1, 1].text(
                index,
                passed_values[index] / 2,
                f"{100 * passed_values[index] / total:.1f}%",
                ha="center",
                va="center",
                fontsize=5.8,
                color="white" if passed_values[index] else "black",
            )
    axes[1, 1].set_xticks(range(len(LENGTH_LABELS)), LENGTH_LABELS, rotation=30, ha="right")
    axes[1, 1].set_xlabel("Unit length (bp)")
    axes[1, 1].set_ylabel("Unique elements")
    axes[1, 1].set_title("D  Five-level completion by unit length", loc="left", weight="bold")
    axes[1, 1].legend(frameon=False)
    sns.despine(fig=figure)
    figure.tight_layout()
    outputs = _save(figure, output_dir, "production_public_element_outcomes")
    plt.close(figure)
    return outputs


def _heatmap_table(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    local = frame.copy()
    local["unit_length_bin"] = pd.cut(
        local.unit_length_bp, bins=LENGTH_BINS, labels=LENGTH_LABELS, right=False
    )
    grouped = (
        local.groupby(["unit_length_bin", "target_copy_count"], observed=False)
        .complete_route_verified.agg(["sum", "count"])
        .reset_index()
    )
    grouped["percent"] = np.where(grouped["count"] > 0, 100 * grouped["sum"] / grouped["count"], np.nan)
    values = grouped.pivot(index="unit_length_bin", columns="target_copy_count", values="percent").reindex(index=LENGTH_LABELS, columns=TARGET_COPY_COUNTS)
    annotation = grouped.assign(
        label=grouped.apply(
            lambda row: "n=0" if row["count"] == 0 else f"{row['percent']:.1f}%\n{int(row['sum'])}/{int(row['count'])}",
            axis=1,
        )
    ).pivot(index="unit_length_bin", columns="target_copy_count", values="label").reindex(index=LENGTH_LABELS, columns=TARGET_COPY_COUNTS)
    return values, annotation


def plot_exact_target_landscape(targets: pd.DataFrame, output_dir: str | Path) -> list[Path]:
    """Figure 2: source-specific unit-length by copy-number heatmaps."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    _prepare_style()
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.8), facecolor="white")
    panels = [*SOURCE_ORDER, "All public sources"]
    for axis, source in zip(axes.ravel(), panels):
        subset = targets if source == "All public sources" else targets.loc[targets.source_database.eq(source)]
        values, annotations = _heatmap_table(subset)
        sns.heatmap(
            values, annot=annotations, fmt="", cmap="Purples", vmin=0, vmax=100,
            linewidths=0.5, linecolor="white", cbar=source == "All public sources",
            cbar_kws={"label": "Complete-route success (%)"}, ax=axis,
            annot_kws={"fontsize": 6.2},
        )
        axis.set_title(SOURCE_LABELS.get(source, source), weight="bold")
        axis.set_xlabel("Target copy number")
        axis.set_ylabel("Unit length (bp)")
        axis.tick_params(axis="y", rotation=0)
    figure.suptitle("Exact-target scalability landscape")
    figure.tight_layout()
    outputs = _save(figure, output_dir, "production_exact_target_landscape")
    plt.close(figure)
    return outputs


def plot_failure_and_idt_evidence(targets: pd.DataFrame, output_dir: str | Path) -> list[Path]:
    """Figure 3: mutually exclusive failure classes and synthesis evidence."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    _prepare_style()
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.8), facecolor="white")
    local = targets.copy()
    local["outcome"] = np.where(
        local.complete_route_verified,
        "Complete exact route",
        local.failure_reason.replace(
            {
                "seed_unavailable": "Seed unavailable",
                "no_active_latent_pair": "No active/latent pair",
                "no_exact_repeat_gain_pair": "No exact repeat-gain pair",
                "vector_or_digest_failure": "Vector/digest failure",
                "purchase_or_idt_failure": "Purchase/IDT failure",
                "no_reachable_precursor_state": "No reachable precursor",
            }
        ),
    )
    outcome_counts = local.outcome.value_counts().sort_values()
    axes[0, 0].barh(outcome_counts.index, outcome_counts.values, color=[UW_PURPLE if value == "Complete exact route" else UW_GOLD for value in outcome_counts.index])
    for y_value, count in enumerate(outcome_counts.values):
        axes[0, 0].text(count, y_value, f" {count:,}", va="center", fontsize=7)
    axes[0, 0].set_xlabel("Exact targets")
    axes[0, 0].set_title("A  Mutually exclusive outcomes", loc="left", weight="bold")

    idt_order = [
        "passed", "failed", "not_applicable_primer_pair_under_90bp",
        "api_unclassified", "api_failure", "not_required", "not_run",
    ]
    idt_counts = local.whole_target_idt_status.value_counts().reindex(idt_order).dropna()
    axes[0, 1].barh(idt_counts.index, idt_counts.values, color=BLUE)
    for y_value, count in enumerate(idt_counts.values):
        axes[0, 1].text(count, y_value, f" {int(count):,}", va="center", fontsize=7)
    axes[0, 1].set_xlabel("Exact targets")
    axes[0, 1].set_title("B  Whole-target IDT status", loc="left", weight="bold")

    compatible = local.loc[local.complete_route_verified]
    evidence_counts = compatible.idt_evidence_tier.value_counts().sort_values()
    axes[1, 0].barh(evidence_counts.index, evidence_counts.values, color=[GREEN, BLUE, UW_PURPLE][:len(evidence_counts)])
    for y_value, count in enumerate(evidence_counts.values):
        axes[1, 0].text(count, y_value, f" {count:,}", va="center", fontsize=7)
    axes[1, 0].set_xlabel("Complete exact targets")
    axes[1, 0].set_title("C  Purchase-evidence tier", loc="left", weight="bold")

    rescued = (
        local.groupby(["source_database", "target_copy_count"])
        .fragment_rescued_by_hurdler.sum()
        .reset_index(name="rescued")
    )
    for source in _source_values(local):
        subset = rescued.loc[rescued.source_database.eq(source)]
        axes[1, 1].plot(
            subset.target_copy_count, subset.rescued,
            marker="o", color=SOURCE_COLORS[source],
            label=SOURCE_LABELS.get(source, source),
        )
    axes[1, 1].set_xticks(TARGET_COPY_COUNTS)
    axes[1, 1].set_xlabel("Target copy number")
    axes[1, 1].set_ylabel("Valid whole-target rescues")
    axes[1, 1].set_title("D  Valid HURDLER rescue", loc="left", weight="bold")
    axes[1, 1].legend(frameon=False)
    sns.despine(fig=figure)
    figure.tight_layout()
    outputs = _save(figure, output_dir, "production_failure_and_idt_evidence")
    plt.close(figure)
    return outputs


def plot_route_complexity(
    targets: pd.DataFrame,
    elements: pd.DataFrame,
    routes: pd.DataFrame,
    transitions: pd.DataFrame,
    fragments: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Figure 4: experimental complexity and selected molecular choices."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    _prepare_style()
    figure, axes = plt.subplots(2, 3, figsize=(9.6, 6.8), facecolor="white")
    passed = targets.loc[targets.complete_route_verified].copy()
    if not passed.empty:
        sns.boxplot(
            data=passed, x="target_copy_count", y="hurdler_step_count",
            color=UW_PURPLE, fliersize=1.5, linewidth=0.8, ax=axes[0, 0],
        )
    axes[0, 0].set_title("A  HURDLER cycles", loc="left", weight="bold")
    axes[0, 0].set_xlabel("Target copy number")
    axes[0, 0].set_ylabel("Digest–ligation steps")

    max_counts = elements.maximum_verified_copy_count.value_counts().sort_index()
    axes[0, 1].bar(max_counts.index.astype(str), max_counts.values, color=UW_GOLD)
    axes[0, 1].set_title("B  Maximum verified copies", loc="left", weight="bold")
    axes[0, 1].set_xlabel("Maximum exact target")
    axes[0, 1].set_ylabel("Unique elements")

    plasmids = (
        routes.plasmid.value_counts().sort_values()
        if not routes.empty else pd.Series(dtype=int)
    )
    axes[0, 2].barh(plasmids.index, plasmids.values, color=BLUE)
    axes[0, 2].set_title("C  Fixed plasmid", loc="left", weight="bold")
    axes[0, 2].set_xlabel("Selected complete routes")

    pair_changes = (
        routes.pair_change_count.value_counts().sort_index()
        if not routes.empty else pd.Series(dtype=int)
    )
    axes[1, 0].bar(pair_changes.index.astype(str), pair_changes.values, color=ORANGE)
    axes[1, 0].set_title("D  RE-pair changes", loc="left", weight="bold")
    axes[1, 0].set_xlabel("Pair changes per route")
    axes[1, 0].set_ylabel("Selected complete routes")

    if not transitions.empty:
        pair_table = pd.crosstab(transitions.site_i_enzyme, transitions.site_ii_enzyme)
        common_i = pair_table.sum(axis=1).nlargest(8).index
        common_ii = pair_table.sum(axis=0).nlargest(8).index
        sns.heatmap(
            pair_table.loc[common_i, common_ii], cmap="Purples", ax=axes[1, 1],
            cbar_kws={"label": "Transitions"}, linewidths=0.3,
        )
    axes[1, 1].set_title("E  Most-used RE pairs", loc="left", weight="bold")
    axes[1, 1].set_xlabel("Site II")
    axes[1, 1].set_ylabel("Site I")

    product_counts = fragments.product_type.value_counts().sort_values() if not fragments.empty else pd.Series(dtype=int)
    axes[1, 2].barh(product_counts.index, product_counts.values, color=GREEN)
    if not fragments.empty:
        product_lengths = fragments.groupby("product_type").purchase_length_bp.agg(
            ["median", "min", "max"]
        ).reindex(product_counts.index)
        for y_value, (product, count) in enumerate(product_counts.items()):
            length = product_lengths.loc[product]
            axes[1, 2].text(
                count,
                y_value,
                f"  n={int(count):,}; median {length['median']:.0f} bp "
                f"({length['min']:.0f}–{length['max']:.0f})",
                va="center",
                fontsize=6.2,
            )
    axes[1, 2].set_title("F  Purchase products", loc="left", weight="bold")
    axes[1, 2].set_xlabel("Selected fragment occurrences")
    sns.despine(fig=figure)
    figure.tight_layout()
    outputs = _save(figure, output_dir, "production_route_complexity")
    plt.close(figure)
    return outputs


def plot_worked_regulatory_example(
    targets: pd.DataFrame,
    routes: pd.DataFrame,
    transitions: pd.DataFrame,
    seeds: pd.DataFrame,
    output_dir: str | Path,
    *,
    element_id_contains: str = "3cc020d8d3025df7",
) -> list[Path]:
    """Figure 5: exact complete route and IDT evidence for RF00050."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    _prepare_style()
    matched_targets = targets.loc[
        targets.element_id.astype(str).str.contains(element_id_contains, regex=False)
    ]
    if matched_targets.empty:
        raise ValueError("RF00050 worked-example element is absent from production")
    target4 = matched_targets.loc[matched_targets.target_copy_count.eq(4)].iloc[0]
    figure, axes = plt.subplots(2, 1, figsize=(7.2, 6.8), facecolor="white")
    axis = axes[0]
    axis.set_xlim(0, 12)
    axis.set_ylim(0, 4.5)
    axis.axis("off")
    axis.set_title("A  RF00050 complete assembly route", loc="left", weight="bold")
    if bool(target4.complete_route_verified):
        route = routes.loc[routes.complete_route_id.eq(target4.complete_route_id)].iloc[0]
        route_transitions = transitions.loc[
            transitions.complete_route_id.eq(target4.complete_route_id)
        ].sort_values("transition_index")
        x_positions = np.linspace(0.6, 10.6, len(route_transitions) + 1)
        seed_copy = int(route.seed_copy_count)
        copy_values = [seed_copy, *route_transitions.result_copy_count.astype(int).tolist()]
        for index, (x_value, copies) in enumerate(zip(x_positions, copy_values)):
            width = min(2.2, 0.45 + 0.22 * copies)
            axis.add_patch(
                FancyBboxPatch(
                    (x_value - width / 2, 1.55), width, 0.8,
                    boxstyle="round,pad=0.05", facecolor="#ECE9F2",
                    edgecolor=UW_PURPLE, linewidth=1.2,
                )
            )
            axis.text(x_value, 1.95, f"{copies} copies", ha="center", va="center", fontsize=8)
            if index < len(route_transitions):
                edge = route_transitions.iloc[index]
                next_x = x_positions[index + 1]
                axis.add_patch(
                    FancyArrowPatch(
                        (x_value + width / 2 + 0.08, 1.95),
                        (next_x - 0.55, 1.95), arrowstyle="-|>",
                        mutation_scale=12, color="#333333",
                    )
                )
                axis.text(
                    (x_value + next_x) / 2, 2.55,
                    f"+{int(edge.donor_copy_count)} donor\n{edge.site_i_enzyme}/{edge.site_ii_enzyme}",
                    ha="center", va="center", fontsize=7,
                )
        axis.text(
            0.5, 0.75,
            f"Plasmid: {route.plasmid}   Final SHA exact: yes   "
            f"HURDLER steps: {int(route.hurdler_step_count)}",
            fontsize=8.5, color=UW_PURPLE,
        )
    else:
        axis.text(
            0.6, 2.0,
            f"RF00050 failed the complete-route criterion:\n{target4.failure_reason}",
            color=RED, fontsize=10,
        )

    axis = axes[1]
    axis.set_title("B  IDT evidence and exact-copy reachability", loc="left", weight="bold")
    scores = matched_targets.set_index("target_copy_count").whole_target_idt_score.reindex(TARGET_COPY_COUNTS)
    plot_scores = scores.fillna(0)
    axis.bar(np.arange(len(TARGET_COPY_COUNTS)) - 0.16, plot_scores, width=0.32, color=UW_GOLD, label="Whole target score")
    max_scores = matched_targets.set_index("target_copy_count").maximum_idt_score.reindex(TARGET_COPY_COUNTS).fillna(0)
    axis.bar(np.arange(len(TARGET_COPY_COUNTS)) + 0.16, max_scores, width=0.32, color=UW_PURPLE, label="Max purchase score")
    axis.axhline(10, color=RED, linewidth=1, label="IDT threshold (<10)")
    axis.set_xticks(range(len(TARGET_COPY_COUNTS)), [str(value) for value in TARGET_COPY_COUNTS])
    axis.set_xlabel("Exact target copy number")
    axis.set_ylabel("IDT rule-score sum")
    axis.legend(frameon=False)
    figure.tight_layout()
    outputs = _save(figure, output_dir, "production_rf00050_complete_route")
    plt.close(figure)
    return outputs


def plot_production_qc(
    metrics: pd.DataFrame,
    benchmark: dict,
    output_dir: str | Path,
) -> list[Path]:
    """Supplement: shard completeness and measured worker scaling."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    _prepare_style()
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.8), facecolor="white")
    shard_count = int(metrics.shard_count.iloc[0])
    side = int(math.ceil(math.sqrt(shard_count))) if shard_count else 1
    grid = np.zeros(side * side)
    for index in metrics.shard_index.astype(int):
        grid[index] = 1
    sns.heatmap(
        grid.reshape(side, side), cmap=["#EFEFEF", GREEN], vmin=0, vmax=1,
        cbar=False, square=True, xticklabels=False, yticklabels=False, ax=axes[0, 0],
    )
    axes[0, 0].set_title(f"A  Completed shards ({len(metrics)}/{shard_count})", loc="left", weight="bold")
    axes[0, 1].hist(metrics.target_rows, bins=min(20, max(2, len(metrics) // 4)), color=UW_PURPLE)
    axes[0, 1].set_title("B  Target rows per shard", loc="left", weight="bold")
    axes[0, 1].set_xlabel("Target rows")
    axes[0, 1].set_ylabel("Shards")
    runs = pd.DataFrame(benchmark.get("runs", []))
    if not runs.empty:
        rate_column = (
            "elements_per_second"
            if "elements_per_second" in runs.columns else
            "targets_per_second"
        )
        axes[1, 0].plot(
            runs.workers, runs[rate_column], marker="o", color=BLUE
        )
        axes[1, 0].set_xticks(runs.workers)
    axes[1, 0].set_title("C  Worker benchmark", loc="left", weight="bold")
    axes[1, 0].set_xlabel("Workers")
    axes[1, 0].set_ylabel("Elements per second")
    axes[1, 1].axis("off")
    axes[1, 1].text(
        0.03, 0.95,
        "Production checks\n"
        f"Version: {DNA_COMPLETE_ROUTE_VERSION}\n"
        f"Elements: {int(metrics.element_rows.sum()):,}\n"
        f"Targets: {int(metrics.target_rows.sum()):,}\n"
        f"Missing shards: {max(0, shard_count-len(metrics))}\n"
        f"Scientific equivalence: {benchmark.get('scientific_equivalence', 'not recorded')}",
        va="top", fontsize=9, linespacing=1.5,
    )
    figure.tight_layout()
    outputs = _save(figure, output_dir, "supplement_production_qc")
    plt.close(figure)
    return outputs


def plot_synthetic_factorial_landscape(
    synthetic_targets: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Supplement: architecture × unit length × GC × copy-number success."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    required = {
        "architecture", "synthetic_unit_length_bp", "gc_target",
        "target_copy_count", "complete_route_verified",
    }
    missing = sorted(required - set(synthetic_targets.columns))
    if missing:
        raise ValueError(f"Synthetic analysis lacks columns: {missing}")
    architectures = [
        "exact_tandem", "fixed_spacer", "alternating_ab",
        "nonrepetitive_control",
    ]
    _prepare_style()
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.8), facecolor="white")
    for axis, architecture in zip(axes.ravel(), architectures):
        local = synthetic_targets.loc[
            synthetic_targets.architecture.eq(architecture)
        ].copy()
        local["length_gc"] = local.apply(
            lambda row: (
                f"{int(row.synthetic_unit_length_bp)} bp / "
                f"{100 * float(row.gc_target):.0f}% GC"
            ),
            axis=1,
        )
        grouped = local.groupby(
            ["length_gc", "target_copy_count"], sort=False
        ).complete_route_verified.agg(["sum", "count"]).reset_index()
        grouped["percent"] = 100 * grouped["sum"] / grouped["count"]
        values = grouped.pivot(
            index="length_gc", columns="target_copy_count", values="percent"
        ).reindex(columns=TARGET_COPY_COUNTS)
        labels = grouped.assign(
            label=grouped.apply(
                lambda row: (
                    f"{row['percent']:.0f}%\n"
                    f"{int(row['sum'])}/{int(row['count'])}"
                ),
                axis=1,
            )
        ).pivot(
            index="length_gc", columns="target_copy_count", values="label"
        ).reindex(index=values.index, columns=TARGET_COPY_COUNTS)
        sns.heatmap(
            values, annot=labels, fmt="", vmin=0, vmax=100, cmap="Purples",
            linewidths=0.4, linecolor="white", cbar=architecture == architectures[-1],
            cbar_kws={"label": "Complete-route success (%)"},
            annot_kws={"fontsize": 5.6}, ax=axis,
        )
        axis.set_title(architecture.replace("_", " ").title(), weight="bold")
        axis.set_xlabel("Target copy number")
        axis.set_ylabel("Unit length / GC target")
        axis.tick_params(axis="y", rotation=0)
    figure.suptitle("Synthetic factorial complete-route landscape")
    figure.tight_layout()
    outputs = _save(figure, output_dir, "supplement_synthetic_factorial")
    plt.close(figure)
    return outputs


def write_production_figure_manifest(
    figures: Iterable[str | Path],
    destination: str | Path,
    *,
    input_tables: Iterable[str | Path],
    source_notebook: str = "notebooks/tasks/08_long_repetitive_dna_assembly.ipynb",
) -> pd.DataFrame:
    input_payload = {
        str(Path(path).absolute()): sha256_file(path)
        for path in input_tables
        if Path(path).is_file()
    }
    input_hash = hashlib.sha256(
        json.dumps(input_payload, sort_keys=True).encode()
    ).hexdigest()
    rows = []
    for value in figures:
        path = Path(value)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty production figure: {path}")
        rows.append(
            {
                "version": DNA_COMPLETE_ROUTE_VERSION,
                "figure": str(path.absolute()),
                "filename": path.name,
                "format": path.suffix.lstrip("."),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "input_table_bundle_sha256": input_hash,
                "input_tables_json": json.dumps(input_payload, sort_keys=True),
                "source_notebook": source_notebook,
                "generated_at": utc_now(),
                "status": "passed",
            }
        )
    frame = pd.DataFrame(rows)
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    return frame


def plot_complete_production_report(
    targets: pd.DataFrame,
    elements: pd.DataFrame,
    routes: pd.DataFrame,
    transitions: pd.DataFrame,
    fragments: pd.DataFrame,
    seeds: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Generate the production-first main figures from finalized tables."""
    outputs: list[Path] = []
    outputs.extend(plot_public_element_outcomes(targets, elements, output_dir))
    outputs.extend(plot_exact_target_landscape(targets, output_dir))
    outputs.extend(plot_failure_and_idt_evidence(targets, output_dir))
    outputs.extend(
        plot_route_complexity(
            targets, elements, routes, transitions, fragments, output_dir
        )
    )
    outputs.extend(
        plot_worked_regulatory_example(
            targets, routes, transitions, seeds, output_dir
        )
    )
    return outputs
