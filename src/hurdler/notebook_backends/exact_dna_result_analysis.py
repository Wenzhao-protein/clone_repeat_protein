"""Backend for V2 notebook 10: exact-DNA production analysis.

The plotting code remains in :mod:`hurdler.dna_assembly_visualization`.  This
backend only resolves the versioned compact tables, validates their shared
denominators, and registers regenerated tables, figures, and reviewer text in
the notebook workspace.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from ..dna_assembly_visualization import (
    plot_complete_production_report,
    plot_production_qc,
    write_production_figure_manifest,
)
from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "10_exact_dna_result_analysis",
    "Exact-DNA production analysis",
    "Analyze complete routes, purchase evidence, HURDLER rescue and reviewer-response counts.",
    production_workflows=("exact-dna-routes", "exact-dna-purchase"),
)

TABLE_FILENAMES = {
    "targets": "exact_dna_target_analysis_compact_v2.parquet",
    "elements": "exact_dna_element_matrix_compact_v2.parquet",
    "routes": "exact_dna_selected_routes_compact_v2.parquet",
    "transitions": "exact_dna_transitions_compact_v2.parquet",
    "fragments": "exact_dna_fragments_compact_v2.parquet",
    "seeds": "exact_dna_seeds_compact_v2.parquet",
    "metrics": "exact_dna_run_metrics_compact_v2.parquet",
    "purchases": "exact_dna_purchase_compact_v2.parquet",
}


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _paths(request: Mapping[str, Any]) -> dict[str, Path]:
    bundled = repo_root() / "data/results"
    production_dir_value = request.get("production_dir")
    production_dir = Path(str(production_dir_value)).expanduser() if production_dir_value else None
    production_names = {
        "targets": "production_target_analysis.parquet",
        "elements": "production_element_matrix.parquet",
        "routes": "production_selected_routes.parquet",
        "transitions": "production_transitions.parquet",
        "fragments": "production_fragments.parquet",
        "seeds": "production_seeds.parquet",
        "metrics": "production_run_metrics.parquet",
    }
    resolved: dict[str, Path] = {}
    for key, filename in TABLE_FILENAMES.items():
        explicit = request.get(f"{key}_path")
        if explicit:
            resolved[key] = Path(str(explicit)).expanduser()
        elif production_dir is not None and key in production_names:
            resolved[key] = production_dir / production_names[key]
        else:
            resolved[key] = bundled / filename
    resolved["headline"] = Path(
        str(request.get("headline_summary", bundled / "exact_dna_production_headline_v2.json"))
    ).expanduser()
    return resolved


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    paths = _paths(request)
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing exact-DNA analysis inputs: " + ", ".join(missing))
    return {"status": "passed", **{key: str(value) for key, value in paths.items()}}


def _validate(
    targets: pd.DataFrame,
    elements: pd.DataFrame,
    routes: pd.DataFrame,
    metrics: pd.DataFrame,
) -> None:
    if len(targets) != 145_210:
        raise ValueError(f"Expected 145,210 public targets, observed {len(targets):,}")
    if len(elements) != 29_042:
        raise ValueError(f"Expected 29,042 public elements, observed {len(elements):,}")
    expected_copies = {2, 4, 8, 16, 32}
    observed = targets.groupby(["source_database", "element_id"]).target_copy_count.agg(
        lambda values: set(map(int, values))
    )
    if not observed.map(lambda values: values == expected_copies).all():
        raise ValueError("Every public element must have exact 2/4/8/16/32-copy rows")
    passed = targets.loc[targets.complete_route_verified.astype(bool)]
    if not passed.final_target_exact.astype(bool).all():
        raise ValueError("A complete route lacks exact final-target verification")
    if set(passed.complete_route_id.astype(str)) != set(routes.complete_route_id.astype(str)):
        raise ValueError("Passing targets and selected routes are not one-to-one")
    invalid_rescue = targets.fragment_rescued_by_hurdler.astype(bool) & ~targets.complete_route_verified.astype(bool)
    if invalid_rescue.any():
        raise ValueError("An incomplete route is incorrectly marked as HURDLER-rescued")
    if len(metrics) != 512 or set(metrics.shard_index.astype(int)) != set(range(512)):
        raise ValueError("The complete-route production shard set is incomplete")


def _analysis_tables(targets: pd.DataFrame, elements: pd.DataFrame) -> dict[str, pd.DataFrame]:
    source_copy = (
        targets.groupby(["source_database", "target_copy_count"])
        .complete_route_verified.agg(successful_targets="sum", total_targets="count")
        .reset_index()
    )
    source_copy["success_rate"] = source_copy.successful_targets / source_copy.total_targets
    outcomes = (
        targets.assign(
            outcome=targets.failure_reason.where(
                ~targets.complete_route_verified.astype(bool), "complete_exact_route"
            )
        )
        .groupby(["source_database", "target_copy_count", "outcome"])
        .size().rename("count").reset_index()
    )
    element_counts = (
        elements.groupby(["source_database", "successful_target_count", "all_five_complete"])
        .size().rename("element_count").reset_index()
    )
    return {
        "exact_target_success_by_source_copy": source_copy,
        "exact_target_outcomes": outcomes,
        "exact_element_completion": element_counts,
    }


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    inputs = preflight(context, request)
    frames = {
        key: pd.read_parquet(inputs[key])
        for key in ("targets", "elements", "routes", "transitions", "fragments", "seeds", "metrics", "purchases")
    }
    headline = json.loads(Path(inputs["headline"]).read_text())
    _validate(frames["targets"], frames["elements"], frames["routes"], frames["metrics"])
    if progress_callback:
        progress_callback({"stage": "validated", "targets": len(frames["targets"])})

    artifacts: list[Path] = []
    for stem, frame in _analysis_tables(frames["targets"], frames["elements"]).items():
        artifacts.extend(write_frame(frame, context.directory("tables") / stem))
    artifacts.extend(
        write_frame(frames["purchases"], context.directory("tables") / "exact_dna_purchase_analysis")
    )
    figures = plot_complete_production_report(
        frames["targets"], frames["elements"], frames["routes"],
        frames["transitions"], frames["fragments"], frames["seeds"],
        context.directory("figures"),
    )
    figures.extend(plot_production_qc(frames["metrics"], {}, context.directory("figures")))
    artifacts.extend(figures)
    manifest_path = context.directory("figures") / "exact_dna_figure_manifest.csv"
    write_production_figure_manifest(
        figures,
        manifest_path,
        input_tables=[Path(inputs[key]) for key in ("targets", "elements", "routes", "transitions", "fragments", "seeds")],
        source_notebook="notebooks/v2/10_exact_dna_result_analysis.ipynb",
    )
    artifacts.append(manifest_path)

    reviewer_text = str(headline.get("reviewer_response_text", "")).strip()
    report = context.directory("reports") / "exact_dna_reviewer_numbers.md"
    report.write_text(
        "# Exact-DNA HURDLER production evidence\n\n"
        f"{reviewer_text}\n\n"
        f"- Public elements: **{len(frames['elements']):,}**\n"
        f"- Exact targets: **{len(frames['targets']):,}**\n"
        f"- Complete exact routes: **{int(frames['targets'].complete_route_verified.sum()):,}**\n"
        f"- Elements complete at all five copy counts: **{int(frames['elements'].all_five_complete.sum()):,}**\n"
        f"- Valid whole-target rescues: **{int(frames['targets'].fragment_rescued_by_hurdler.sum()):,}**\n"
        "- `<90 bp` primer pairs are not assigned a gBlock complexity score.\n"
        "- IDT score `<10` is a complexity-screen result, not a formal ordering guarantee.\n"
        "- The historical 53.67% result is a final-step-only baseline and is excluded.\n"
    )
    artifacts.append(report)
    complete = int(frames["targets"].complete_route_verified.sum())
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=artifacts,
        metrics={
            "unique_public_elements": len(frames["elements"]),
            "real_exact_targets": len(frames["targets"]),
            "real_exact_targets_complete": complete,
            "complete_route_fraction": complete / len(frames["targets"]),
            "elements_all_five_complete": int(frames["elements"].all_five_complete.sum()),
            "valid_hurdler_rescues": int(frames["targets"].fragment_rescued_by_hurdler.sum()),
            "selected_routes": len(frames["routes"]),
            "production_shards": len(frames["metrics"]),
        },
        next_notebooks=["11_reproducibility"],
        limitations=[
            "IDT complexity scores are not quotes or wet-lab ordering guarantees.",
            "Synthetic factorial cases are intentionally excluded from reviewer N/X.",
        ],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
