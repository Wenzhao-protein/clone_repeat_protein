"""Credential-free IDT score histories and trajectory figures."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


IDT_SCORE_THRESHOLD = 10.0
UW_PURPLE = "#4B2E83"
UW_GOLD = "#B7A57A"
PASS_GREEN = "#2D6A4F"
FAIL_RED = "#B31B1B"
UNCLASSIFIED_GREY = "#6B7280"
RULE_COLORS = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#D55E00",
    "#F0E442",
    "#332288",
    "#44AA99",
    "#AA4499",
)


HISTORY_COLUMNS = (
    "evaluation_index",
    "fragment_id",
    "fragment_kind",
    "repeat_copies",
    "ga_generations",
    "feedback_round",
    "request_length_bp",
    "idt_total_score",
    "idt_classification",
    "idt_cache_hit",
    "idt_response_sha256",
    "positive_rule_names_json",
    "positive_rule_reasons_json",
    "rule_scores_json",
)


def _finite_number(value: Any) -> float | None:
    if (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    ):
        return float(value)
    return None


def _json_value(value: Any, fallback: Any) -> Any:
    if isinstance(value, type(fallback)):
        return value
    try:
        parsed = json.loads(str(value))
    except (TypeError, ValueError, json.JSONDecodeError):
        return fallback
    return parsed if isinstance(parsed, type(fallback)) else fallback


def _audit_reasons(audit: Mapping[str, Any]) -> dict[str, str]:
    reasons: dict[str, str] = {}
    details = _json_value(audit.get("idt_rule_details_json", "[]"), [])
    for detail in details:
        if not isinstance(detail, dict):
            continue
        score = _finite_number(detail.get("score"))
        if not (bool(detail.get("is_violated")) or (score is not None and score > 0)):
            continue
        name = str(detail.get("name") or "unnamed_rule")
        reasons[name] = str(detail.get("display_text") or "").strip() or "positive IDT rule score"
    return reasons


def idt_score_history_rows(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Normalize IDT audit rows or ``idt_fragment_scored`` progress events."""
    rows: list[dict[str, Any]] = []
    for ordinal, record in enumerate(records, start=1):
        is_event = (
            str(record.get("event_name") or "") == "idt_fragment_scored"
            or str(record.get("status") or "") == "fragment_scored"
        )
        raw_scores = (
            record.get("idt_rule_scores", {})
            if is_event
            else _json_value(record.get("idt_rule_scores_json", "{}"), {})
        )
        rule_scores = {
            str(name): float(score)
            for name, score in dict(raw_scores).items()
            if _finite_number(score) is not None
        }
        score = _finite_number(
            record.get("idt_score") if is_event else record.get("idt_complexity_score")
        )
        if is_event:
            classification = str(record.get("idt_classification") or "unclassified")
            reasons = {
                str(name): str(reason)
                for name, reason in dict(record.get("idt_rule_reasons") or {}).items()
            }
            positive_names = [
                str(value) for value in record.get("idt_positive_rules", ())
            ]
            response_sha = str(record.get("idt_response_sha256") or "")
            details = dict(record.get("details") or {})
        else:
            if record.get("idt_explicit_pass") is True and score is not None and score < IDT_SCORE_THRESHOLD:
                classification = "passed"
            elif bool(record.get("idt_score_complete")) and score is not None:
                classification = "rejected"
            else:
                classification = "unclassified"
            reasons = _audit_reasons(record)
            positive_names = sorted(name for name, value in rule_scores.items() if value > 0)
            response_sha = str(record.get("response_sha256") or "")
            details = record
        rows.append(
            {
                "evaluation_index": int(record.get("idt_evaluation_index") or ordinal),
                "fragment_id": str(
                    record.get("idt_fragment_name") if is_event else record.get("fragment_id")
                    or "fragment"
                ),
                "fragment_kind": str(record.get("fragment_kind") or ""),
                "repeat_copies": record.get("copies") if is_event else record.get("repeat_copies"),
                "ga_generations": record.get("generations") if is_event else record.get("ga_generations"),
                "feedback_round": record.get("feedback_round"),
                "request_length_bp": details.get("request_length_bp"),
                "idt_total_score": score,
                "idt_classification": classification,
                "idt_cache_hit": bool(record.get("idt_cache_hit", False)),
                "idt_response_sha256": response_sha,
                "positive_rule_names_json": json.dumps(positive_names, sort_keys=True),
                "positive_rule_reasons_json": json.dumps(reasons, sort_keys=True),
                "rule_scores_json": json.dumps(rule_scores, sort_keys=True),
            }
        )
    return rows


def plot_idt_score_trajectory(
    rows: Sequence[Mapping[str, Any]],
    *,
    title: str = "IDT complexity score trajectory",
):
    """Create the shared live/static total-and-rule score figure."""
    normalized = (
        [dict(row) for row in rows]
        if rows and all("idt_total_score" in row for row in rows)
        else idt_score_history_rows(rows)
    )
    if not normalized:
        raise ValueError("At least one IDT score record is required")
    width = min(24.0, max(8.0, 0.48 * len(normalized)))
    fig, (total_ax, rule_ax) = plt.subplots(
        2,
        1,
        figsize=(width, 6.6),
        sharex=True,
        gridspec_kw={"height_ratios": (1.0, 1.25)},
        constrained_layout=True,
    )
    x = np.array([int(row["evaluation_index"]) for row in normalized])
    totals = np.array(
        [
            np.nan if _finite_number(row.get("idt_total_score")) is None else float(row["idt_total_score"])
            for row in normalized
        ]
    )
    total_ax.plot(x, totals, color=UW_PURPLE, marker="o", linewidth=2.2, label="Total score")
    total_ax.axhline(
        IDT_SCORE_THRESHOLD,
        color=UW_GOLD,
        linewidth=1.8,
        label="Acceptance threshold (<10)",
    )
    marker_colors = {
        "passed": PASS_GREEN,
        "rejected": FAIL_RED,
        "unclassified": UNCLASSIFIED_GREY,
    }
    for row, x_value, score in zip(normalized, x, totals, strict=True):
        classification = str(row["idt_classification"])
        color = marker_colors.get(classification, UNCLASSIFIED_GREY)
        if math.isfinite(float(score)):
            total_ax.scatter([x_value], [score], color=color, edgecolor="white", linewidth=0.8, s=55, zorder=4)
            total_ax.annotate(
                f"{score:g}",
                (x_value, score),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=7,
                color=color,
            )
        else:
            total_ax.scatter([x_value], [0], color=color, marker="x", s=55, zorder=4)
            total_ax.annotate("unclassified", (x_value, 0), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=7, color=color)
        if classification == "rejected":
            reasons = _json_value(row.get("positive_rule_reasons_json", "{}"), {})
            names = list(reasons) or _json_value(row.get("positive_rule_names_json", "[]"), [])
            if names and math.isfinite(float(score)):
                total_ax.annotate(
                    ", ".join(str(value) for value in names[:3]),
                    (x_value, score),
                    xytext=(5, -16),
                    textcoords="offset points",
                    fontsize=6.5,
                    color=FAIL_RED,
                )
    total_ax.set_title(title, color=UW_PURPLE, fontweight="bold")
    total_ax.set_ylabel("IDT score sum")
    total_ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    total_ax.legend(loc="upper left", frameon=False)

    parsed_scores = [
        _json_value(row.get("rule_scores_json", "{}"), {}) for row in normalized
    ]
    rule_names = sorted({str(name) for scores in parsed_scores for name in scores})
    for index, name in enumerate(rule_names):
        values = [
            float(scores[name]) if name in scores and _finite_number(scores[name]) is not None else np.nan
            for scores in parsed_scores
        ]
        rule_ax.plot(
            x,
            values,
            marker="o",
            markersize=3.5,
            linewidth=1.4,
            color=RULE_COLORS[index % len(RULE_COLORS)],
            label=name,
        )
    rule_ax.set_ylabel("Rule score")
    rule_ax.set_xlabel("Chronological IDT evaluation")
    rule_ax.set_xticks(x)
    rule_ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    if rule_names:
        rule_ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=7)
    labels = [
        f"#{row['evaluation_index']} {row['fragment_kind'] or 'fragment'}\n"
        f"c={row['repeat_copies'] if row['repeat_copies'] is not None else '—'} "
        f"fb={row['feedback_round'] if row['feedback_round'] is not None else '—'}"
        for row in normalized
    ]
    rule_ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    for axis in (total_ax, rule_ax):
        axis.spines[["top", "right"]].set_visible(False)
    return fig


def write_idt_score_trajectory(
    audits: Sequence[Mapping[str, Any]], output_dir: str | Path
) -> dict[str, str]:
    """Write the canonical score table and figures when API scores exist."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    rows = idt_score_history_rows(audits)
    csv_path = destination / "idt_score_history.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HISTORY_COLUMNS)
        writer.writeheader()
        writer.writerows({key: row.get(key) for key in HISTORY_COLUMNS} for row in rows)
    paths = {"idt_score_history_csv": str(csv_path)}
    if not rows:
        return paths
    figure = plot_idt_score_trajectory(rows)
    try:
        for suffix in ("png", "pdf", "svg"):
            path = destination / f"idt_score_trajectory.{suffix}"
            figure.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight", facecolor="white")
            paths[f"idt_score_trajectory_{suffix}"] = str(path)
    finally:
        plt.close(figure)
    return paths
