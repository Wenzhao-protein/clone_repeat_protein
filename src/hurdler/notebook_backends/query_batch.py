"""Backend for V2 notebook 04: single and batch HURDLER screens."""

from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
from Bio import SeqIO

from ..constants import validate_protein_sequence
from ..index import PatternIndex
from ..matching import materialize_best_solution, query_all_plasmids
from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "04_query_batch",
    "Query and batch screen",
    "Screen one amino-acid module, FASTA records or a CSV table against all eight plasmid profiles.",
    production_workflows=("module-stage1",),
    default_request={"module": "TDDEEIARIIAYAARQTT"},
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _sequences(request: Mapping[str, Any]) -> list[tuple[str, str]]:
    if request.get("module"):
        return [(str(request.get("sequence_id", "module_001")), validate_protein_sequence(str(request["module"])))]
    source = Path(str(request["input_path"]))
    if source.suffix.lower() in {".fa", ".faa", ".fasta"}:
        return [(str(record.id), validate_protein_sequence(str(record.seq))) for record in SeqIO.parse(source, "fasta")]
    frame = pd.read_csv(source)
    id_column = str(request.get("id_column", "sequence_id"))
    sequence_column = str(request.get("sequence_column", "sequence"))
    if id_column not in frame or sequence_column not in frame:
        raise ValueError(f"CSV requires {id_column!r} and {sequence_column!r}")
    return [(str(row[id_column]), validate_protein_sequence(str(row[sequence_column]))) for _, row in frame.iterrows()]


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    index = Path(str(request.get("index_dir", repo_root() / "data/artifacts/legacy-optimized-v1")))
    if not (index / "metadata.json").is_file():
        raise FileNotFoundError(index / "metadata.json")
    sequences = _sequences(request)
    if not sequences:
        raise ValueError("No protein sequences were supplied")
    return {"status": "passed", "sequence_count": len(sequences), "index_dir": str(index)}


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    info = preflight(context, request)
    index_dir = Path(info["index_dir"])
    index = PatternIndex.load(index_dir)
    summaries: list[dict[str, Any]] = []
    candidates: list[dict[str, Any]] = []
    fasta_path = context.directory("fasta") / "query_modules.fasta"
    fasta_lines: list[str] = []
    for ordinal, (sequence_id, module) in enumerate(_sequences(request), 1):
        fasta_lines.extend([f">{sequence_id}", module])
        matches = query_all_plasmids(module, index)
        materialized = [materialize_best_solution(item, index) for item in matches]
        successes = [item for item in materialized if bool(item.get("success"))]
        selected = successes[0] if successes else {}
        summaries.append(
            {
                "sequence_id": sequence_id,
                "selected_module_sequence_aa": module,
                "module_length_aa": len(module),
                "hurdler_compatible": bool(successes),
                "compatible_plasmid_count": len(successes),
                "selected_plasmid": selected.get("plasmid", ""),
                "selected_re_pair": (
                    f"{selected.get('site_i_enzyme')}/{selected.get('site_ii_enzyme')}"
                    if successes else ""
                ),
                "site_i_enzyme": selected.get("site_i_enzyme", ""),
                "site_i_3mer_aa": selected.get("site_i_3mer", ""),
                "site_ii_enzyme": selected.get("site_ii_enzyme", ""),
                "site_ii_3mer_aa": selected.get("site_ii_3mer", ""),
                "direction": selected.get("direction", ""),
                "site_i_start_aa_0based": selected.get("site_i_position"),
                "site_ii_start_aa_0based": selected.get("site_ii_position"),
            }
        )
        for item in materialized:
            candidates.append({"sequence_id": sequence_id, **item})
        if progress_callback:
            progress_callback({"stage": "query", "status": "progress", "completed": ordinal, "total": info["sequence_count"]})
    fasta_path.write_text("\n".join(fasta_lines) + "\n")
    summary_paths = write_frame(pd.DataFrame(summaries), context.directory("tables") / "hurdler_query_summary")
    candidate_paths = write_frame(pd.DataFrame(candidates), context.directory("tables") / "hurdler_query_candidates")
    compatible = sum(bool(row["hurdler_compatible"]) for row in summaries)
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[*summary_paths, *candidate_paths, fasta_path],
        metrics={"sequence_count": len(summaries), "compatible_count": compatible, "compatible_fraction": compatible / len(summaries)},
        next_notebooks=["05_repeat_designer", "07_production_builder"],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
