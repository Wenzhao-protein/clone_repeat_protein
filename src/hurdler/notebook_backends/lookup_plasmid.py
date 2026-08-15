"""Backend for V2 notebook 02: sparse lookups and annotated plasmids."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from ..index import PatternIndex, build_pattern_index
from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from ..plasmid_reference import (
    build_plasmid_reference,
    bundled_plasmid_reference_path,
    load_plasmid_reference,
    validate_plasmid_reference,
)
from ..protein_index import ProteinPatternIndex, build_protein_pattern_index
from ..qc import legacy_qc
from .common import BackendSpec, preflight_files, repo_root, result_from_paths, zip_paths


SPEC = BackendSpec(
    "02_lookup_plasmid",
    "Lookup and plasmid database",
    "Build or validate the frozen sparse lookup, vector-aware protein index and seven-vector/eight-profile database.",
    default_request={"rebuild": False, "rules": "legacy-optimized-v1"},
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _paths(request: Mapping[str, Any]) -> tuple[Path, Path, Path, Path]:
    root = repo_root()
    source = Path(str(request.get("source_dir", root / "output"))).expanduser()
    reference = Path(str(request.get("reference_dir", root / "data/reference_output"))).expanduser()
    legacy = Path(str(request.get("legacy_index_dir", root / "data/artifacts/legacy-optimized-v1"))).expanduser()
    protein = Path(str(request.get("protein_index_dir", root / "data/artifacts/vector-aware-hurdler-v2"))).expanduser()
    return source.absolute(), reference.absolute(), legacy.absolute(), protein.absolute()


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    source, reference, legacy, protein = _paths(request)
    required = [reference / "orthogonality.csv", bundled_plasmid_reference_path()]
    if request.get("rebuild"):
        required.extend(
            source / name
            for name in (
                "hurdler_site_i_dataframe.csv", "hurdler_site_ii_dataframe.csv",
                "selected_site_iii_enzymes.csv", "site_i_site_ii_pairing_matrix.csv",
            )
        )
    else:
        required.extend([legacy / "metadata.json", protein / "metadata.json"])
    return preflight_files(required)


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    preflight(context, request)
    source, reference, legacy, protein = _paths(request)
    if request.get("rebuild"):
        legacy = context.run_root / "artifacts" / "legacy-optimized-v1"
        protein = context.run_root / "artifacts" / "vector-aware-hurdler-v2"
        build_pattern_index(source, legacy, orthogonality_path=reference / "orthogonality.csv")
        build_protein_pattern_index(source, protein, orthogonality_path=reference / "orthogonality.csv")
    legacy_index = PatternIndex.load(legacy)
    protein_index = ProteinPatternIndex.load(protein)
    plasmids = load_plasmid_reference()
    plasmid_qc = validate_plasmid_reference(plasmids)
    qc_path = context.directory("reports") / "lookup_qc_v2.json"
    qc = legacy_qc(source, legacy, qc_path)
    qc["protein_pattern_count"] = int(len(protein_index.keys))
    qc["plasmid_reference"] = plasmid_qc
    qc_path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
    archive = context.run_root / "lookup_plasmid_database_v2.zip"
    files = list(legacy.glob("*")) + list(protein.glob("*")) + [bundled_plasmid_reference_path(), qc_path]
    zip_paths(archive, [path for path in files if path.is_file()], base=repo_root())
    return result_from_paths(
        context,
        backend_id=SPEC.notebook_id,
        request=request,
        paths=[qc_path, archive],
        metrics={
            "legacy_pattern_count": int(len(legacy_index.keys)),
            "protein_pattern_count": int(len(protein_index.keys)),
            "enzyme_pair_count": int(len(legacy_index.pair_table)),
            **plasmid_qc,
        },
        next_notebooks=["03_module_corpus", "04_query_batch", "05_repeat_designer", "06_exact_dna_designer"],
    )


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
