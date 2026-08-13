from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

from hurdler.idt_trajectory import (
    idt_score_history_rows,
    plot_idt_score_trajectory,
    write_idt_score_trajectory,
)


def _audit(
    fragment: str,
    score: float | None,
    *,
    complete: bool,
    passed: bool,
    rule_scores: dict[str, float],
    cache_hit: bool = False,
):
    return {
        "fragment_id": fragment,
        "fragment_kind": "secondary",
        "repeat_copies": 12,
        "ga_generations": 10,
        "feedback_round": 1,
        "request_length_bp": 500,
        "idt_complexity_score": score,
        "idt_score_complete": complete,
        "idt_explicit_pass": passed,
        "idt_cache_hit": cache_hit,
        "response_sha256": fragment * 8,
        "idt_rule_scores_json": json.dumps(rule_scores),
        "idt_rule_details_json": json.dumps(
            [
                {
                    "name": name,
                    "score": value,
                    "is_violated": value > 0,
                    "display_text": f"{name} reason",
                }
                for name, value in rule_scores.items()
            ]
        ),
    }


def test_idt_history_preserves_missing_components_and_classification():
    rows = idt_score_history_rows(
        [
            _audit("a", 3.0, complete=True, passed=True, rule_scores={"Repeat": 3.0, "GC": 0.0}),
            _audit("b", 10.0, complete=True, passed=False, rule_scores={"Repeat": 10.0}),
            _audit("c", None, complete=False, passed=False, rule_scores={}),
        ]
    )
    assert [row["evaluation_index"] for row in rows] == [1, 2, 3]
    assert [row["idt_classification"] for row in rows] == [
        "passed", "rejected", "unclassified"
    ]
    assert json.loads(rows[1]["rule_scores_json"]) == {"Repeat": 10.0}
    assert "GC" not in json.loads(rows[1]["rule_scores_json"])
    assert "Repeat reason" in rows[1]["positive_rule_reasons_json"]


def test_idt_progress_events_keep_fragment_order_cache_and_missing_rules():
    events = [
        {
            "event_name": "idt_fragment_scored",
            "status": "fragment_scored",
            "idt_evaluation_index": 7,
            "idt_fragment_name": "primary_7",
            "fragment_kind": "primary",
            "copies": 7,
            "feedback_round": 2,
            "idt_score": 4.0,
            "idt_classification": "passed",
            "idt_rule_scores": {"Repeat": 4.0, "GC": 0.0},
            "idt_rule_reasons": {"Repeat": "repeat reason"},
            "idt_cache_hit": False,
        },
        {
            "event_name": "idt_fragment_scored",
            "status": "fragment_scored",
            "idt_evaluation_index": 8,
            "idt_fragment_name": "secondary_7",
            "fragment_kind": "secondary",
            "copies": 7,
            "feedback_round": 2,
            "idt_score": None,
            "idt_classification": "unclassified",
            "idt_rule_scores": {"Repeat": 2.0},
            "idt_rule_reasons": {},
            "idt_cache_hit": True,
        },
    ]
    rows = idt_score_history_rows(events)
    assert [row["evaluation_index"] for row in rows] == [7, 8]
    assert [row["fragment_id"] for row in rows] == ["primary_7", "secondary_7"]
    assert rows[1]["idt_cache_hit"] is True
    assert rows[1]["idt_classification"] == "unclassified"
    assert "GC" not in json.loads(rows[1]["rule_scores_json"])
    figure = plot_idt_score_trajectory(rows)
    try:
        assert len(figure.axes) == 2
    finally:
        plt.close(figure)


def test_idt_trajectory_has_only_the_required_threshold_and_exports(tmp_path):
    audits = [
        _audit("a", 2.0, complete=True, passed=True, rule_scores={"Repeat": 2.0}),
        _audit("b", 12.0, complete=True, passed=False, rule_scores={"Repeat": 9.0, "GC": 3.0}, cache_hit=True),
    ]
    figure = plot_idt_score_trajectory(audits)
    try:
        assert len(figure.axes) == 2
        assert any(line.get_label() == "Acceptance threshold (<10)" for line in figure.axes[0].lines)
        assert not any(line.get_linestyle() == "--" for axis in figure.axes for line in axis.lines)
    finally:
        plt.close(figure)

    paths = write_idt_score_trajectory(audits, tmp_path)
    assert set(paths) == {
        "idt_score_history_csv",
        "idt_score_trajectory_png",
        "idt_score_trajectory_pdf",
        "idt_score_trajectory_svg",
    }
    for value in paths.values():
        assert Path(value).stat().st_size > 0
    with Path(paths["idt_score_history_csv"]).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[1]["idt_cache_hit"] == "True"
    assert "Sequence" not in rows[0]


def test_empty_idt_history_writes_table_only(tmp_path):
    paths = write_idt_score_trajectory([], tmp_path)
    assert set(paths) == {"idt_score_history_csv"}
    assert Path(paths["idt_score_history_csv"]).read_text().startswith("evaluation_index,")
