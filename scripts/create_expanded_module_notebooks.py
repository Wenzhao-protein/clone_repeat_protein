#!/usr/bin/env python3
"""Create the three authoritative expanded-module notebooks."""

from __future__ import annotations

import hashlib
from pathlib import Path

import nbformat as nbf


REPO = Path(__file__).resolve().parents[1]
STUDY = REPO / "studies" / "hurdler_validation"
VERSION = "expanded-middle-repeatsdb-foldseek-v1"
TABLES3 = STUDY / "step03_module_corpus" / "tables" / VERSION
TABLES4 = STUDY / "step04_module_optimization" / "tables" / VERSION
FIGURES4 = STUDY / "step04_module_optimization" / "figures" / VERSION
SCRATCH = Path(
    "/net/scratch/wendai/projects/hurdler/clone_repeat_protein/"
    "studies/hurdler_validation"
)


def code(value: str, *, parameters: bool = False) -> nbf.NotebookNode:
    metadata = {"tags": ["parameters"]} if parameters else {}
    return nbf.v4.new_code_cell(value, metadata=metadata)


def markdown(value: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(value)


def book(title: str, introduction: str, parameters: str, cells: list[object]):
    result = nbf.v4.new_notebook()
    result.metadata.kernelspec = {
        "display_name": "HURDLER",
        "language": "python",
        "name": "hurdler",
    }
    result.metadata.language_info = {"name": "python", "version": "3.11"}
    result.cells = [
        markdown(f"# {title}\n\n{introduction}"),
        code(parameters, parameters=True),
        code(
            "from pathlib import Path\n"
            "import hashlib, json\n"
            "import pandas as pd\n"
            "VERSION = 'expanded-middle-repeatsdb-foldseek-v1'\n"
            "def sha256(path):\n"
            "    path = Path(path)\n"
            "    if not path.is_file(): return None\n"
            "    h = hashlib.sha256()\n"
            "    with path.open('rb') as handle:\n"
            "        for chunk in iter(lambda: handle.read(1024*1024), b''): h.update(chunk)\n"
            "    return h.hexdigest()\n"
            "def read_optional(path):\n"
            "    path = Path(path)\n"
            "    if not path.is_file(): return None\n"
            "    return pd.read_parquet(path) if path.suffix == '.parquet' else pd.read_csv(path)\n"
            "run_context = {'corpus_version': VERSION, 'rules_version': 'legacy-optimized-v1', 'inputs': {}, 'row_counts': {}, 'filter_flow': [], 'limitations': [], 'status': 'passed'}"
        ),
        *cells,
        code("run_context"),
    ]
    return result


def write(path: Path, notebook: nbf.NotebookNode) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    relative = path.relative_to(REPO).as_posix()
    for index, cell in enumerate(notebook.cells):
        identity = f"{relative}\0{index}\0{cell.cell_type}\0{cell.source}"
        cell["id"] = hashlib.sha256(identity.encode()).hexdigest()[:16]
    nbf.write(notebook, path)


def main() -> int:
    reference_parameters = (
        f"REPO = {str(REPO)!r}\n"
        f"ANNOTATION_INVENTORY = {str(SCRATCH / 'step03_module_corpus/runs/run93_repeatsdb_inventory/raw/repeatsdb_annotations.parquet')!r}\n"
        f"NATURAL_CATALOG = {str(TABLES3 / 'natural_module_catalog.parquet')!r}\n"
        f"NATURAL_MAPPINGS = {str(TABLES3 / 'natural_module_catalog_source_mappings.parquet')!r}\n"
        f"NATURAL_REGIONS = {str(TABLES3 / 'natural_module_catalog_all_region_source_mappings.parquet')!r}\n"
        f"DESIGNED_INVENTORY = {str(TABLES3 / 'designed_structure_inventory_expanded.parquet')!r}\n"
        f"AF3_VALIDATION = {str(TABLES3 / 'designed_af3_validation.parquet')!r}\n"
        f"DESIGNED_CATALOG = {str(TABLES3 / 'designed_module_catalog.parquet')!r}\n"
        f"DESIGNED_EXCLUSIONS = {str(TABLES3 / 'designed_module_catalog_exclusions.csv')!r}\n"
        f"DESIGNED_CANDIDATES = {str(TABLES3 / 'designed_module_catalog_boundary_candidates.parquet')!r}\n"
        f"DESIGNED_UNITS = {str(TABLES3 / 'designed_module_catalog_unit_alignment.parquet')!r}\n"
        f"DESIGNED_POSITIONS = {str(TABLES3 / 'designed_module_catalog_position_variability.parquet')!r}"
    )
    write(
        REPO / "notebooks/reference/03_expanded_middle_module_acquisition.ipynb",
        book(
            "RepeatsDB-direct natural and DSSP/Foldseek designed modules",
            "This is the authoritative acquisition and boundary notebook. Natural unit coordinates are copied directly from RepeatsDB and never replaced by periodic inference. Designed boundaries require independent Biotite DSSP and Foldseek 3Di/TM-align evidence. Source notebooks contain no execution outputs; missing production shards are reported explicitly.",
            reference_parameters,
            [
                markdown(
                    "## Boundary policy\n\n"
                    "Natural: group by canonical UniProt (else full-sequence SHA256), choose the longest annotated RepeatsDB region, sort its annotated units, and select index `(unit_count-1)//2`. DSSP/Foldseek are QC only.\n\n"
                    "Designed: scan lags 6..half-chain with eight-state DSSP and Foldseek 3Di, identify the dominant contiguous recurrent scale, globally validate adjacent fragments by Foldseek/TM-align, choose the smallest passing primitive period, then select the earlier middle copy. MAFFT defines fixed positions at ≥80% conservation."
                ),
                code(
                    "paths = [ANNOTATION_INVENTORY, NATURAL_CATALOG, NATURAL_MAPPINGS, NATURAL_REGIONS, DESIGNED_INVENTORY, AF3_VALIDATION, DESIGNED_CATALOG, DESIGNED_EXCLUSIONS, DESIGNED_CANDIDATES, DESIGNED_UNITS, DESIGNED_POSITIONS]\n"
                    "availability = pd.DataFrame({'artifact': [Path(p).name for p in paths], 'path': paths, 'exists': [Path(p).is_file() for p in paths], 'sha256': [sha256(p) for p in paths]})\n"
                    "run_context['inputs'] = dict(zip(availability.artifact, availability.sha256))\n"
                    "if not availability.exists.all():\n"
                    "    run_context['status'] = 'production_pending'\n"
                    "    run_context['limitations'].append('One or more Digs production/finalization artifacts are not complete; no values are imputed.')\n"
                    "availability"
                ),
                code(
                    "annotation_inventory = read_optional(ANNOTATION_INVENTORY)\n"
                    "natural = read_optional(NATURAL_CATALOG)\n"
                    "natural_mappings = read_optional(NATURAL_MAPPINGS)\n"
                    "natural_regions = read_optional(NATURAL_REGIONS)\n"
                    "designed_inventory = read_optional(DESIGNED_INVENTORY)\n"
                    "af3_validation = read_optional(AF3_VALIDATION)\n"
                    "designed = read_optional(DESIGNED_CATALOG)\n"
                    "exclusions = read_optional(DESIGNED_EXCLUSIONS)\n"
                    "candidates = read_optional(DESIGNED_CANDIDATES)\n"
                    "units = read_optional(DESIGNED_UNITS)\n"
                    "positions = read_optional(DESIGNED_POSITIONS)\n"
                    "run_context['row_counts'] = {name: (None if value is None else len(value)) for name, value in {'annotations':annotation_inventory,'natural_unique_units':natural,'natural_proteins':natural_mappings,'natural_annotated_regions':natural_regions,'designed_inventory':designed_inventory,'af3_validation':af3_validation,'designed_strict_units':designed,'designed_exclusions':exclusions,'period_candidates':candidates,'aligned_units':units,'position_rows':positions}.items()}\n"
                    "pd.DataFrame([run_context['row_counts']])"
                ),
                code(
                    "if natural is not None:\n"
                    "    assert natural.boundary_refinement_status.eq('source_annotation_middle_unit').all()\n"
                    "    assert natural.selected_module_index.eq((natural.selected_module_count.astype(int)-1)//2 + 1).all()\n"
                    "    assert natural.unit_sequence.str.len().eq(natural.unit_length.astype(int)).all()\n"
                    "    natural_qc = natural.groupby(['annotation_schema','structure_source']).agg(unique_middle_units=('unit_sequence','nunique'), proteins=('protein_key','nunique'), median_length=('unit_length','median')).reset_index()\n"
                    "else:\n"
                    "    natural_qc = pd.DataFrame({'status':['production_pending']})\n"
                    "natural_qc"
                ),
                code(
                    "if designed_inventory is not None:\n"
                    "    structure_flow = designed_inventory.groupby(['family','structure_inventory_status']).size().rename('proteins').reset_index()\n"
                    "else:\n"
                    "    structure_flow = pd.DataFrame({'status':['inventory_missing']})\n"
                    "structure_flow"
                ),
                code(
                    "if designed is not None:\n"
                    "    assert designed.strict_dual_evidence_passed.fillna(False).all()\n"
                    "    assert designed.selected_module_index.eq((designed.repeat_count.astype(int)-1)//2 + 1).all()\n"
                    "    designed_qc = designed[['module_id','family','period','repeat_count','selected_module_index','dssp_state_agreement','dssp_transition_agreement','foldseek_3di_identity','foldseek_median_min_tm','foldseek_median_lddt','foldseek_median_coverage']].sort_values(['family','module_id'])\n"
                    "else:\n"
                    "    designed_qc = pd.DataFrame({'status':['strict_DSSP_Foldseek_shards_pending']})\n"
                    "designed_qc.head(20)"
                ),
                code(
                    "if candidates is not None and 'module_id' in candidates and candidates.module_id.eq('designed_THR29').any():\n"
                    "    thr29 = candidates.loc[candidates.module_id.eq('designed_THR29') & candidates.period.isin([23,45,68,136,204]), ['period','dssp_state_agreement','dssp_transition_agreement','foldseek_3di_identity','repeat_block_recurrence_composite','dominant_recurrence_shortlist','foldseek_global_thresholds_passed']]\n"
                    "else:\n"
                    "    thr29 = pd.DataFrame({'golden_expectation':['THR29 period 68; DSSP≈0.994; Foldseek 3Di≈0.912']})\n"
                    "thr29"
                ),
                code(
                    "run_context['filter_flow'] = ['enumerate every RepeatsDB PDB and AlphaFoldDB annotation without class caps', 'one longest annotated region per biological protein', 'select exact earlier-middle RepeatsDB unit without inferred boundary replacement', 'match author/PDB structure then AlphaFoldDB then missing-only AF3', 'Biotite DsspApp eight-state annotation', 'Foldseek structureto3didescriptor lag self-alignment', 'dominant recurrence-scale filtering to reject local structural texture', 'global adjacent-copy Foldseek/TM-align validation', 'MAFFT fixed/variable position table', 'deduplicate exact middle AA sequence within Natural and Designed separately']\n"
                    "run_context['limitations'].extend(['designed rows without strict dual evidence remain in the exclusion inventory', 'AF3 is QC/boundary evidence for designed proteins only and cannot modify natural RepeatsDB coordinates'])"
                ),
            ],
        ),
    )

    stage1_parameters = (
        f"PER_MODULE = {str(TABLES4 / 'module_compatibility.parquet')!r}\n"
        f"BINNED = {str(TABLES4 / 'module_compatibility_binned.parquet')!r}\n"
        f"CANDIDATES = {str(TABLES4 / 'module_compatibility_candidates.parquet')!r}\n"
        f"FIGURE_STEM = {str(FIGURES4 / 'module_compatibility_by_length')!r}"
    )
    write(
        REPO / "notebooks/tasks/06_module_compatibility_by_length.ipynb",
        book(
            "Stage 1: HURDLER compatibility by middle-module length",
            "Every accepted unique middle module is evaluated against all eight plasmids with `legacy-optimized-v1`. The tables, not this notebook, contain all candidate solutions and the deterministic selected plasmid/RE pair.",
            stage1_parameters,
            [
                code(
                    "from hurdler.module_experiments import plot_compatibility\n"
                    "import pyarrow.parquet as pq\n"
                    "per_module = read_optional(PER_MODULE); binned = read_optional(BINNED)\n"
                    "candidate_rows = pq.ParquetFile(CANDIDATES).metadata.num_rows if Path(CANDIDATES).is_file() else None\n"
                    "for path in (PER_MODULE,BINNED,CANDIDATES): run_context['inputs'][Path(path).name] = sha256(path)\n"
                    "if per_module is None or binned is None or candidate_rows is None:\n"
                    "    run_context['status']='production_pending'; run_context['limitations'].append('Compatibility shards have not yet been finalized; no provisional bar heights are shown.')\n"
                    "    overview=pd.DataFrame({'status':['production_pending']})\n"
                    "else:\n"
                    "    assert ~per_module.duplicated(['collection','unit_sequence']).any()\n"
                    "    assert (binned.compatible_count+binned.incompatible_count).eq(binned.total_count).all()\n"
                    "    assert per_module.hurdler_compatible.sum() == binned.compatible_count.sum()\n"
                    "    overview=per_module.groupby('collection').agg(modules=('module_id','size'),compatible=('hurdler_compatible','sum'),median_length_aa=('unit_length','median')).reset_index()\n"
                    "    overview['compatible_fraction']=overview.compatible/overview.modules\n"
                    "    plot_compatibility(binned, FIGURE_STEM)\n"
                    "run_context['row_counts']={'modules':None if per_module is None else len(per_module),'candidate_solutions':candidate_rows,'length_bins':None if binned is None else len(binned)}\n"
                    "overview"
                ),
                code(
                    "binned if binned is not None else pd.DataFrame({'status':['production_pending']})"
                ),
                code(
                    "from IPython.display import Image, display\n"
                    "figure=Path(str(FIGURE_STEM)+'.png')\n"
                    "if figure.is_file(): display(Image(filename=str(figure)))\n"
                    "else: print('Figure pending finalized Stage-1 tables:', figure)"
                ),
                code(
                    "run_context['filter_flow']=['<6AA: repeat motif to shortest effective module >=6AA', 'scan effective module twice with frozen signed-overhang rules', 'evaluate all eight maintained plasmids', 'compatible iff at least one solution exists', 'retain all candidates and select by frozen deterministic ranking', 'shared 10-AA bins with empty intervening bins visible']\n"
                    "run_context['limitations'].append('Natural and Designed panels share bin order and compatibility encoding; corpus completeness depends on the reference notebook status.')"
                ),
            ],
        ),
    )

    stage2_parameters = (
        f"RESULTS = {str(TABLES4 / 'maximum_copy_results.parquet')!r}\n"
        f"TRACE = {str(TABLES4 / 'adaptive_copy_search_trace.parquet')!r}\n"
        f"FINAL_SUMMARY = {str(TABLES4 / 'module_final_summary.parquet')!r}\n"
        f"FIGURE_STEM = {str(FIGURES4 / 'module_length_vs_maximum_copies')!r}"
    )
    write(
        REPO / "notebooks/tasks/07_adaptive_selected_pair_capacity.ipynb",
        book(
            "Stage 2: adaptive selected-pair-clean maximum repeat count",
            "Only Stage-1-compatible modules enter this experiment. Each exact GA DNA is scored by the live IDT API; all finite rule Scores are summed and only totals strictly below 10 pass. Positive-score rules reweight the next GA attempt at the same copy count.",
            stage2_parameters,
            [
                code(
                    "from hurdler.module_experiments import plot_maximum_copy_scatter\n"
                    "import pyarrow.parquet as pq\n"
                    "results=read_optional(RESULTS); final_summary=read_optional(FINAL_SUMMARY)\n"
                    "trace_rows=pq.ParquetFile(TRACE).metadata.num_rows if Path(TRACE).is_file() else None\n"
                    "for path in (RESULTS,TRACE,FINAL_SUMMARY): run_context['inputs'][Path(path).name]=sha256(path)\n"
                    "if results is None or trace_rows is None or final_summary is None:\n"
                    "    run_context['status']='production_pending'; run_context['limitations'].append('Adaptive GA/IDT production shards are incomplete; no maximum is imputed as zero or one.')\n"
                    "    outcomes=pd.DataFrame({'status':['production_pending']})\n"
                    "else:\n"
                    "    accepted=results.final_passed.fillna(False)\n"
                    "    assert pd.to_numeric(results.loc[accepted,'verified_max_copies']).ge(2).all()\n"
                    "    assert results.loc[accepted,'selected_pair_re_site_excess'].eq(0).all()\n"
                    "    assert pd.to_numeric(results.loc[accepted,'idt_complexity_score']).lt(10).all()\n"
                    "    assert results.loc[accepted,'idt_response_sha256'].astype(str).str.len().eq(64).all()\n"
                    "    assert results.loc[accepted,'validation_passed'].fillna(False).all()\n"
                    "    proof=results.loc[accepted,'adaptive_maximum_proof_status'].isin(['capacity_limit_reached','next_copy_failed_at_100'])\n"
                    "    assert proof.all()\n"
                    "    outcomes=results.groupby(['collection','fragment_limit_bp','final_status']).size().rename('modules').reset_index()\n"
                    "    plot_maximum_copy_scatter(results,FIGURE_STEM)\n"
                    "run_context['row_counts']={'compatible_capacity_rows':None if results is None else len(results),'GA_IDT_attempts':trace_rows,'final_summary_rows':None if final_summary is None else len(final_summary),'verified_maxima':None if results is None else int(accepted.sum())}\n"
                    "outcomes"
                ),
                code(
                    "if trace_rows is not None:\n"
                    "    import duckdb\n"
                    "    escaped_trace=str(TRACE).replace(\"'\",\"''\")\n"
                    "    feedback=duckdb.connect().execute(f\"SELECT copies, generations, count(*) AS rejected_attempts FROM read_parquet('{escaped_trace}') WHERE try_cast(idt_complexity_score AS DOUBLE) >= 10 GROUP BY copies, generations ORDER BY copies, generations\").fetchdf()\n"
                    "else: feedback=pd.DataFrame({'status':['production_pending']})\n"
                    "feedback.head(30)"
                ),
                code(
                    "from IPython.display import Image, display\n"
                    "figure=Path(str(FIGURE_STEM)+'.png')\n"
                    "if figure.is_file(): display(Image(filename=str(figure)))\n"
                    "else: print('Scatter pending finalized Stage-2 tables:', figure)"
                ),
                code(
                    "run_context['filter_flow']=['freeze each Stage-1 selected plasmid and Site-I/Site-II pair', 'require exact translation of middle_unit × n', 'hard constraint: zero excess selected-pair sites', 'soft scores: Site III, nonselected RE sites, GC, repeats, hairpins, CAI', '10-generation binary lower-bound search', 'advance one copy at a time with 10/20/40/60/80/100 generations', 'score every completed GA DNA through live IDT API', 'sum every finite rule Score and accept only sum <10', 'use positive rule names/reasons to raise mapped GA weights before same-length retry', 'report >=2 only after cap or next-copy 100-generation proof']\n"
                    "run_context['limitations'].append('IDT API failures are unclassified and cannot pass; credentials and raw submitted DNA are never stored in Git.')"
                ),
            ],
        ),
    )
    print("created 3 expanded-middle authoritative notebooks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
