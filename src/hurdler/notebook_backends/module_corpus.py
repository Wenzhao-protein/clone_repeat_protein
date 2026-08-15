"""Backend for V2 notebook 03: natural/direct and designed/structural modules."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from ..repeatsdb import annotation_repeat_regions, select_longest_region_per_protein
from ..structural_repeats import lag_evidence, scan_dual_evidence_periods
from .common import BackendSpec, load_table, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "03_module_corpus",
    "Repeat-module corpus",
    "Use RepeatsDB unit coordinates for natural proteins and strict DSSP/Foldseek evidence for designed proteins.",
    production_workflows=("repeatsdb-natural", "designed-structure", "missing-af3"),
    default_request={"tutorial_rows_per_collection": 12},
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    source = Path(str(request.get("catalog", repo_root() / "data/results/module_corpus_tutorial_fixture_v2.parquet")))
    if not source.is_file():
        raise FileNotFoundError(source)
    return {"status": "passed", "catalog": str(source.absolute())}


def _tutorial_catalog(source: Path, rows: int) -> pd.DataFrame:
    wanted = [
        "module_id", "display_name", "collection", "family", "middle_module_sequence_aa",
        "middle_module_length_aa", "full_protein_sequence_aa", "source_accession", "source_url",
        "middle_module_start", "middle_module_end", "repeat_region_start", "repeat_region_end",
        "repeat_count", "boundary_method", "boundary_status", "strict_dual_evidence_passed",
        "dssp_state_agreement", "foldseek_3di_identity", "structure_source_type", "corpus_version",
    ]
    if source.suffix.lower() == ".parquet":
        frame = pd.read_parquet(source)
        frame = frame[[name for name in wanted if name in frame.columns]]
    else:
        frame = pd.read_csv(source, usecols=lambda name: name in wanted, low_memory=False)
    selected = (
        frame.sort_values(["collection", "module_id"], kind="mergesort")
        .groupby("collection", sort=True, group_keys=False)
        .head(rows)
        .reset_index(drop=True)
    )
    selected["middle_unit_selection_qc"] = "earlier_middle_from_source_or_strict_block"
    selected["natural_boundaries_modified"] = False
    return selected


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    preflight(context, request)
    source = Path(str(request.get("catalog", repo_root() / "data/results/module_corpus_tutorial_fixture_v2.parquet")))
    if context.mode == "tutorial":
        frame = _tutorial_catalog(source, int(request.get("tutorial_rows_per_collection", 12)))
    else:
        frame = load_table(source)
    natural = frame.loc[frame.collection.astype(str).str.casefold() == "natural"]
    if not natural.empty and "boundary_method" in natural:
        invalid = ~natural.boundary_method.astype(str).str.contains("repeatsdb", case=False, na=False)
        if invalid.any():
            raise ValueError("Natural V2 rows must retain direct RepeatsDB boundaries")
    designed = frame.loc[frame.collection.astype(str).str.casefold() == "designed"]
    if not designed.empty and "strict_dual_evidence_passed" in designed:
        accepted = designed.boundary_status.astype(str).str.contains("accepted|passed", case=False, na=False)
        dual_evidence = designed.strict_dual_evidence_passed.eq(True)
        if (accepted & ~dual_evidence).any():
            raise ValueError("Accepted designed rows require strict DSSP/Foldseek dual evidence")
    paths = write_frame(frame, context.directory("tables") / "expanded_middle_repeatsdb_foldseek_v2")
    exclusions = pd.DataFrame(columns=["module_id", "collection", "exclusion_reason", "rerun_workflow"])
    exclusion_paths = write_frame(exclusions, context.directory("tables") / "module_exclusions_v2")
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[*paths, *exclusion_paths],
        metrics={
            "catalog_rows": len(frame),
            "natural_rows": len(natural),
            "designed_rows": len(designed),
            "unique_middle_modules": int(frame.middle_module_sequence_aa.nunique()),
            "boundary_policy": "repeatsdb-direct-natural_strict-dssp-foldseek-designed",
        },
        next_notebooks=["07_production_builder", "09_module_result_analysis"],
        limitations=["Full structure processing and missing AF3 predictions require a production bundle."],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
