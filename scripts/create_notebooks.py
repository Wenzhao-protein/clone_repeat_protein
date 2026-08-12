#!/usr/bin/env python3
"""Generate thin, parameterized source notebooks for maintained workflows."""

from __future__ import annotations

import hashlib
from pathlib import Path

import nbformat as nbf


REPO = Path(__file__).resolve().parents[1]
STUDY = REPO / "studies" / "hurdler_validation"
SCRATCH_STUDY = Path("/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation")
INDEX = Path("/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step01_reference_lookup/runs/run01_production/raw/legacy-optimized-v1")


def notebook(
    title: str,
    description: str,
    parameters: str,
    cells: list[object],
    *,
    rules: str = "legacy-optimized-v1",
) -> nbf.NotebookNode:
    book = nbf.v4.new_notebook()
    book.metadata.kernelspec = {"display_name": "HURDLER", "language": "python", "name": "hurdler"}
    book.metadata.language_info = {"name": "python", "version": "3.11"}
    book.cells = [
        nbf.v4.new_markdown_cell(
            f"# {title}\n\n{description}\n\n"
            f"**Rules:** `{rules}`; **seed:** 42 unless explicitly noted."
        ),
        nbf.v4.new_code_cell(parameters, metadata={"tags": ["parameters"]}),
        nbf.v4.new_code_cell(
            "from pathlib import Path\n"
            "import hashlib, json\n"
            "import pandas as pd\n"
            "\n"
            "def sha256(path):\n"
            "    path = Path(path)\n"
            "    if not path.is_file(): return None\n"
            "    h = hashlib.sha256()\n"
            "    with path.open('rb') as handle:\n"
            "        for chunk in iter(lambda: handle.read(1024 * 1024), b''): h.update(chunk)\n"
            "    return h.hexdigest()\n"
            "\n"
            "run_context = {'rule_profile': RULE_PROFILE, 'input_hashes': {}, 'row_counts': {}, 'filter_flow': [], 'limitations': []}"
        ),
        *cells,
    ]
    return book


def code(source: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(source)


def markdown(source: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(source)


def write(path: Path, book: nbf.NotebookNode) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    relative = path.relative_to(REPO).as_posix()
    for index, cell in enumerate(book.cells):
        identity = f"{relative}\0{index}\0{cell.cell_type}\0{cell.source}"
        cell["id"] = hashlib.sha256(identity.encode()).hexdigest()[:16]
    nbf.write(book, path)


def main() -> int:
    common = f"REPO = {str(REPO)!r}\nRULE_PROFILE = 'legacy-optimized-v1'"
    write(
        REPO / "notebooks/reference/01_reference_manifest.ipynb",
        notebook(
            "Reference manifest",
            "Hashes the maintained reference database. Database generation is kept separate from HURDLER task calls.",
            common + f"\nREFERENCE_DIR = {str(REPO / 'data/reference_output')!r}\nOUTPUT = {str(STUDY / 'step01_reference_lookup/tables/reference_manifest.json')!r}",
            [
                code("from hurdler.reference import build_reference_manifest\nmanifest = build_reference_manifest(REFERENCE_DIR, OUTPUT)\nrun_context['input_hashes'] = {row['name']: row['sha256'] for row in manifest['files']}\nrun_context['row_counts']['reference_files'] = len(manifest['files'])\nrun_context['filter_flow'] = ['hash every maintained reference file; no implicit download or replacement']\npd.DataFrame(manifest['files'])"),
                markdown("The manifest records exact content hashes and row counts. It does not download or silently replace reference files."),
            ],
        ),
    )
    write(
        REPO / "notebooks/reference/02_lookup_qc.ipynb",
        notebook(
            "Sparse lookup and legacy QC",
            "Summarizes the versioned sparse index and quantifies known historical Site-II/Site-III overhang mismatches.",
            common + f"\nINDEX_DIR = {str(INDEX)!r}\nLEGACY_OUTPUT_DIR = {str(REPO / 'output')!r}\nQC_OUTPUT = {str(STUDY / 'step01_reference_lookup/tables/legacy_qc.json')!r}",
            [
                code("from hurdler.index import PatternIndex\nfrom hurdler.qc import legacy_qc\nindex = PatternIndex.load(INDEX_DIR)\nqc = legacy_qc(LEGACY_OUTPUT_DIR, INDEX_DIR, QC_OUTPUT)\nrun_context['input_hashes']['pattern_index.npz'] = sha256(Path(INDEX_DIR) / 'pattern_index.npz')\nrun_context['row_counts'] = {'patterns': len(index.keys), 'enzyme_pairs': len(index.pair_table)}\nrun_context['limitations'] = [qc['note']]\npd.DataFrame([qc])"),
                code("index.pair_table.groupby(['site_i_ovhg', 'site_ii_ovhg']).size().rename('enzyme_pairs').reset_index()"),
            ],
        ),
    )
    write(
        REPO / "notebooks/tasks/01_hurdler_query.ipynb",
        notebook(
            "HURDLER task query",
            "Runs one repeat module across the eight frozen plasmids using the same artifact as validation and success-rate analyses.",
            common + f"\nINDEX_DIR = {str(INDEX)!r}\nMODULE = 'VLA'",
            [
                code("from hurdler.index import PatternIndex\nfrom hurdler.matching import materialize_best_solution, query_all_plasmids\nindex = PatternIndex.load(INDEX_DIR)\nrows = [materialize_best_solution(result, index) for result in query_all_plasmids(MODULE, index)]\nrun_context['input_hashes']['pattern_index.npz'] = sha256(Path(INDEX_DIR) / 'pattern_index.npz')\nresults = pd.DataFrame(rows)\nrun_context['row_counts'] = {'plasmids_tested': len(results), 'successful_plasmids': int(results.success.sum())}\nrun_context['filter_flow'] = ['expand motifs shorter than 6AA', 'scan doubled module', 'return first frozen-order match per plasmid']\nresults"),
                code("results.groupby('success').size().rename('plasmids').to_frame()"),
            ],
        ),
    )
    success_run_2x = SCRATCH_STUDY / "step02_success_landscape/runs/run03_single_file_16core/raw"
    success_run_3x = SCRATCH_STUDY / "step02_success_landscape/runs/run06_three_copy_16core/raw"
    write(
        REPO / "notebooks/tasks/02_success_rate_1_60.ipynb",
        notebook(
            "HURDLER success rate, 1–60AA (figure through 50AA)",
            "Uses three repeated copies to complete every cyclic 3-mer window while retaining the original per-plasmid pattern population, distance criterion, random inputs, and get_re_sites.ipynb typography, lines, text, and plasmid order on a near-square 6×5-inch canvas. The historical two-copy result remains the exact comparison baseline; the authoritative plot displays 1–50AA.",
            f"REPO = {str(REPO)!r}\nRULE_PROFILE = 'historical-notebook-success-v1-three-copy-scan'"
            + f"\nTWO_SHORT_RESULTS = {str(success_run_2x / 'short_motifs_1_5.parquet')!r}"
            + f"\nTWO_RANDOM_RESULTS = {str(success_run_2x / 'random_modules_6_60.parquet')!r}"
            + f"\nSHORT_RESULTS = {str(success_run_3x / 'short_motifs_1_5.parquet')!r}"
            + f"\nRANDOM_RESULTS = {str(success_run_3x / 'random_modules_6_60.parquet')!r}"
            + f"\nCOMPARISON = {str(STUDY / 'step02_success_landscape/tables/scan_copy_comparison_2x_vs_3x.csv')!r}"
            + f"\nBY_LENGTH = {str(STUDY / 'step02_success_landscape/tables/scan_copy_improvement_by_length.csv')!r}"
            + f"\nLANDSCAPE_SCRIPT = {str(REPO / 'scripts/run_success_landscape_single_files.py')!r}"
            + f"\nFIGURE_DIR = {str(STUDY / 'step02_success_landscape/figures/scan_3x')!r}",
            [
                code("import importlib.util\nimport pyarrow.parquet as pq\nspec = importlib.util.spec_from_file_location('single_file_success_landscape', LANDSCAPE_SCRIPT)\nlandscape = importlib.util.module_from_spec(spec)\nspec.loader.exec_module(landscape)\npaths = {'two_short': TWO_SHORT_RESULTS, 'two_random': TWO_RANDOM_RESULTS, 'three_short': SHORT_RESULTS, 'three_random': RANDOM_RESULTS, 'comparison': COMPARISON, 'by_length': BY_LENGTH}\nrun_context['input_hashes'] = {name: sha256(path) for name, path in paths.items()}\nrates = landscape.success_summary(Path(SHORT_RESULTS), Path(RANDOM_RESULTS))\ntwo_copy_rates = landscape.success_summary(Path(TWO_SHORT_RESULTS), Path(TWO_RANDOM_RESULTS))\ncomparison = pd.read_csv(COMPARISON)\nby_length = pd.read_csv(BY_LENGTH)\nshort_rows = pq.ParquetFile(SHORT_RESULTS).metadata.num_rows\nrandom_rows = pq.ParquetFile(RANDOM_RESULTS).metadata.num_rows\nassert short_rows == 3_368_420\nassert random_rows == 440_000\nassert len(rates) == len(two_copy_rates) == len(comparison) == 60 * 8\nassert not {'ci95_low', 'ci95_high'}.intersection(rates.columns)\nassert comparison.success_delta.ge(0).all()\nassert rates.loc[rates.module_length.le(2), 'successes'].sum() == 0\nimproved_lengths = by_length.loc[by_length.any_improvement, 'module_length'].astype(int).tolist()\nrun_context['row_counts'] = {'rate_rows': len(rates), 'short_motifs': short_rows, 'sampled_sequences': random_rows, 'improved_length_plasmid_rows': int(comparison.improved.sum()), 'improved_module_lengths': len(improved_lengths)}\nrun_context['filter_flow'] = ['rebuild the original notebook plasmid-specific pattern population', 'keep identical exhaustive and random module inputs for 2x and 3x', 'repeat short motif to shortest >=6AA module', 'scan module+module+module with 5 <= d < effective module length', 'compare every length/plasmid to the historical two-copy baseline']\nrun_context['limitations'] = ['6-60AA values are seeded Monte Carlo point estimates; uncertainty intervals are intentionally not drawn']\nrates.head()"),
                code("by_length[['module_length','improved_plasmids','added_successes_all_plasmids','mean_success_rate_2x','mean_success_rate_3x','maximum_plasmid_rate_delta','improved_plasmid_names']]"),
                code("from IPython.display import Image, display\nfigures = landscape.plot_success_curve(rates, Path(FIGURE_DIR), file_stem='success_rate_1_60_scan_3x')\ndisplay(Image(filename=str(Path(FIGURE_DIR) / 'success_rate_1_60_scan_3x.png')))\nrates.groupby(['method','plasmid']).agg(lengths=('module_length','nunique'), tests=('tests','sum'), successes=('successes','sum')).reset_index()"),
                markdown("Every sequence is now scanned as `module + module + module`; the distance upper bound remains the single effective module length. The third copy completes the cyclic window whose Site-I 3-mer begins at the final residue of the first copy and whose Site-II start is `L-1` residues later. It adds no new biological motif or pattern rule. The 1AA and 2AA exhaustive rates remain zero for all eight plasmids. The plot uses a near-square 6×5-inch canvas with the original notebook's default font and line settings, plasmid order, exact title `3-mer Probability vs Sequence Length`, and 50AA upper limit. It shows one point estimate per length and plasmid, with no confidence band or method-divider annotation."),
            ],
            rules="historical-notebook-success-v1-three-copy-scan",
        ),
    )
    module_tables = STUDY / "step04_module_optimization/tables/periodic_v4"
    module_figures = STUDY / "step04_module_optimization/figures/periodic_v4"
    write(
        REPO / "notebooks/tasks/03_repeat_module_benchmark.ipynb",
        notebook(
            "Natural and designed repeat modules",
            "Compares HURDLER feasibility and IDT-orderable adaptive maximum fragment capacity for curated natural100, designed_all, and designed_primary100 collections. All downstream work uses the real middle copy of each inferred repeat region rather than its first copy. Every locally acceptable candidate is scored by the IDT API; rejection reasons reweight the GA before the same length is retried, and module count advances only after an explicit zero-violation result.",
            common
            + f"\nCATALOG = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_catalog_periodic_v4_middle.parquet')!r}"
            + f"\nRESULTS = {str(module_tables / 'module_hurdler_results.parquet')!r}"
            + f"\nCONSTRUCTS = {str(module_tables / 'optimized_constructs.parquet')!r}"
            + f"\nADAPTIVE_TRACE = {str(module_tables / 'adaptive_copy_search_trace.parquet')!r}"
            + f"\nFINAL_MODULE_SUMMARY = {str(module_tables / 'module_final_summary.parquet')!r}"
            + f"\nHURDLER_FRACTIONS = {str(module_tables / 'module_hurdler_usable_fraction.csv')!r}"
            + f"\nSCATTER_DATA = {str(module_tables / 'module_length_copy_scatter_data.csv')!r}"
            + f"\nSUMMARY_DIR = {str(module_tables)!r}"
            + f"\nFIGURE_DIR = {str(module_figures)!r}",
            [
                code("catalog = pd.read_parquet(CATALOG)\nresults = pd.read_parquet(RESULTS)\nconstructs = pd.read_parquet(CONSTRUCTS)\nassert catalog.selected_module_policy.eq('repeat-region-middle-unit-tie-earlier-v1').all()\nassert catalog.unit_sequence.eq(catalog.selected_module_sequence).all()\nrun_context['input_hashes'] = {Path(path).name: sha256(path) for path in (CATALOG, RESULTS, CONSTRUCTS)}\nresults_compare = pd.concat([results, results.loc[results.in_designed_primary100].assign(collection='designed_primary100')], ignore_index=True)\nconstructs_compare = pd.concat([constructs, constructs.loc[constructs.in_designed_primary100].assign(collection='designed_primary100')], ignore_index=True)\nga_applicable = constructs_compare.loc[constructs_compare.ga_status.eq('passed')].copy()\nga_applicable['ga_score_improvement'] = ga_applicable.ga_initial_score - ga_applicable.ga_score\nsummary = (results_compare.groupby(['collection','plasmid']).success.agg(successes='sum', modules='count').reset_index())\nsummary['success_rate'] = summary.successes/summary.modules\ncapacity = constructs_compare.groupby(['collection','fragment_limit_bp']).agg(modules=('module_id','nunique'), median_mathematical_max=('mathematical_max_copies','median'), median_verified_max=('verified_max_copies','median'), local_ga_passes=('ga_local_constraints_passed','sum'), idt_api_scored=('idt_api_called','sum'), idt_rule_violations=('idt_violation_count','sum'), final_idt_passes=('final_passed','sum')).reset_index()\nga_summary = ga_applicable.groupby(['collection','fragment_limit_bp']).agg(constructs=('module_id','count'), ga_improved=('ga_improved','sum'), initial_repeated_re_sites=('ga_initial_repeated_re_site_excess','sum'), final_repeated_re_sites=('repeated_re_site_excess','sum'), repeated_re_sites_removed=('ga_repeated_re_site_excess_removed','sum'), median_score_improvement=('ga_score_improvement','median')).reset_index()\nfailure_modes = constructs_compare.groupby(['collection','fragment_limit_bp','optimization_status']).size().rename('modules').reset_index()\nidt_summary = constructs_compare.groupby(['collection','idt_status']).agg(constructs=('module_id','size'), idt_rule_violations=('idt_violation_count','sum')).reset_index()\nrun_context['row_counts'] = {'catalog_modules': len(catalog), 'module_plasmid_rows': len(results), 'construct_cap_rows': len(constructs), 'ga_applicable': len(ga_applicable), 'natural100': int((catalog.collection=='natural100').sum()), 'designed_all': int((catalog.collection=='designed_all').sum()), 'designed_primary100': int(catalog.in_designed_primary100.sum())}\nrun_context['filter_flow'] = ['select the real middle repeat copy and preserve its exact variable residues', 'require middle-copy policy and selected-sequence equality', 'deduplicate exact units within collection', 'query eight plasmids', 'rank candidates by optimizability then frozen stable order', 'genetic refinement with repeated-RE-site fitness term', 'score every locally acceptable DNA with the IDT API before changing copy count', 'map IDT rejection reasons to GA weights and retry the same copy count', 'advance only after an explicit zero-violation/orderable result']\nrun_context['limitations'] = ['external vector/adapter deduction is configurable and is zero in this run', 'modules longer than 60AA remain in this report but not the 1-60AA curve', 'IDT scores require OAuth account credentials and are not imputed when unavailable', 'IDT orderability is a hard adaptive-search gate and can limit the achievable copy count even when all local constraints pass']\nPath(SUMMARY_DIR).mkdir(parents=True, exist_ok=True)\nsummary.to_csv(Path(SUMMARY_DIR)/'module_success_summary.csv', index=False)\ncapacity.to_csv(Path(SUMMARY_DIR)/'module_capacity_summary.csv', index=False)\nga_summary.to_csv(Path(SUMMARY_DIR)/'module_ga_summary.csv', index=False)\nfailure_modes.to_csv(Path(SUMMARY_DIR)/'module_failure_modes.csv', index=False)\nidt_summary.to_csv(Path(SUMMARY_DIR)/'module_idt_summary.csv', index=False)\nsummary"),
                code("import matplotlib.pyplot as plt\nimport seaborn as sns\nPath(FIGURE_DIR).mkdir(parents=True, exist_ok=True)\nsns.set_theme(style='whitegrid')\norder = ['natural100','designed_all','designed_primary100']\nfig, axes = plt.subplots(2, 2, figsize=(15, 11), facecolor='white')\naxes = axes.ravel()\nsns.barplot(data=summary, x='collection', y='success_rate', hue='plasmid', order=order, ax=axes[0], palette='colorblind')\nsns.boxplot(data=constructs_compare, x='collection', y='verified_max_copies', hue='fragment_limit_bp', order=order, ax=axes[1], palette=['#4B2E83','#B7A57A'])\nsns.boxplot(data=ga_applicable, x='collection', y='ga_repeated_re_site_excess_removed', hue='fragment_limit_bp', order=order, ax=axes[2], palette=['#2D7DD2','#F45D01'])\nsns.countplot(data=constructs_compare, x='collection', hue='optimization_status', order=order, ax=axes[3], palette=['#4B2E83','#B7A57A','#85754D','#999999'])\nfor ax in axes: ax.tick_params(axis='x', rotation=20)\naxes[0].set_title('HURDLER success by plasmid'); axes[1].set_title('Maximum verified full copies'); axes[2].set_title('Repeated RE sites removed by GA'); axes[3].set_title('Optimization outcomes (both caps)')\nsns.despine(); fig.tight_layout()\nfor suffix in ('png','pdf'): fig.savefig(Path(FIGURE_DIR) / f'module_benchmark.{suffix}', dpi=300, facecolor='white')\nfig"),
                code("capacity.merge(failure_modes.groupby(['collection','fragment_limit_bp']).modules.sum().rename('outcome_rows').reset_index(), on=['collection','fragment_limit_bp'])"),
                code("final_module_summary = pd.read_parquet(FINAL_MODULE_SUMMARY)\nhurdler_fractions = pd.read_csv(HURDLER_FRACTIONS)\nrun_context['input_hashes'][Path(FINAL_MODULE_SUMMARY).name] = sha256(FINAL_MODULE_SUMMARY)\nrun_context['input_hashes'][Path(HURDLER_FRACTIONS).name] = sha256(HURDLER_FRACTIONS)\nrun_context['row_counts']['final_module_summary_rows'] = len(final_module_summary)\nhurdler_fractions"),
                code("from IPython.display import Image, display\nscatter_data = pd.read_csv(SCATTER_DATA)\nscatter_figure = Path(FIGURE_DIR) / 'module_length_vs_max_orderable_copies.png'\nassert set(scatter_data.module_type) == {'Natural', 'Designed'}\nassert scatter_data.max_orderable_module_copies.ge(2).all()\nassert scatter_figure.stat().st_size > 0\nrun_context['input_hashes'][Path(SCATTER_DATA).name] = sha256(SCATTER_DATA)\nrun_context['row_counts']['idt_orderable_repeat_scatter_modules'] = scatter_data.module_id.nunique()\ndisplay(Image(filename=str(scatter_figure)))\nscatter_data.groupby('module_type').agg(modules=('module_id','nunique'), median_module_length_aa=('unit_length_aa','median'), median_orderable_copies=('max_orderable_module_copies','median')).reset_index()"),
                code("idt_violation_rows = []\nfor row in constructs_compare.itertuples(index=False):\n    value = getattr(row, 'idt_violation_names_json', None)\n    if not isinstance(value, str): continue\n    for rule_name in json.loads(value):\n        idt_violation_rows.append({'collection': row.collection, 'fragment_limit_bp': row.fragment_limit_bp, 'rule_name': rule_name})\nidt_violations = pd.DataFrame(idt_violation_rows, columns=['collection','fragment_limit_bp','rule_name']).groupby(['collection','fragment_limit_bp','rule_name']).size().rename('violations').reset_index()\nidt_violations.to_csv(Path(SUMMARY_DIR)/'module_idt_violations.csv', index=False)\nidt_violations.sort_values(['collection','fragment_limit_bp','violations'], ascending=[True,True,False])"),
                code("idt_summary"),
                code("""trace = pd.read_parquet(ADAPTIVE_TRACE)
run_context['input_hashes'][Path(ADAPTIVE_TRACE).name] = sha256(ADAPTIVE_TRACE)
trace_compare = pd.concat([trace, trace.loc[trace.in_designed_primary100].assign(collection='designed_primary100')], ignore_index=True)
constructs_compare['copy_change_vs_legacy'] = constructs_compare.verified_max_copies - constructs_compare.pre_adaptive_verified_max_copies
constructs_compare['reached_mathematical_bound'] = constructs_compare.verified_max_copies.eq(constructs_compare.adaptive_search_upper_bound_copies) & constructs_compare.final_passed.fillna(False)
constructs_compare['orderable_maximum'] = constructs_compare.final_passed.fillna(False) & constructs_compare.adaptive_orderable_passed.fillna(False) & constructs_compare.adaptive_boundary_proven.fillna(False)
adaptive_capacity = constructs_compare.groupby(['collection','fragment_limit_bp']).agg(modules=('module_id','nunique'), modules_with_maximum=('orderable_maximum','sum'), boundary_proven=('adaptive_boundary_proven','sum'), median_legacy_max=('pre_adaptive_verified_max_copies','median'), median_adaptive_max=('verified_max_copies','median'), maximum_adaptive_copies=('verified_max_copies','max'), recovered_copies=('copy_change_vs_legacy','sum'), reached_mathematical_bound=('reached_mathematical_bound','sum'), median_search_evaluations=('adaptive_search_evaluations','median'), median_winning_generations=('ga_generations','median')).reset_index()
adaptive_stops = constructs_compare.groupby(['collection','fragment_limit_bp','adaptive_stop_reason']).size().rename('constructs').reset_index()
trace_compare['idt_rejected'] = trace_compare.idt_explicit_pass.eq(False)
trace_compare['feedback_applied'] = trace_compare.idt_feedback_adjustments_json.fillna('[]').ne('[]')
generation_usage = trace_compare.groupby(['collection','fragment_limit_bp','phase','generations']).agg(evaluations=('module_id','size'), orderable_passes=('passed','sum'), idt_scored=('idt_api_called','sum'), idt_rejections=('idt_rejected','sum'), feedback_updates=('feedback_applied','sum')).reset_index()
rejection_rows = trace_compare.loc[trace_compare.idt_rejected, ['collection','fragment_limit_bp','copies','generations','idt_violation_names_json','idt_rule_scores_json','idt_feedback_adjustments_json']].copy()
rejection_rows['idt_reason'] = rejection_rows.idt_violation_names_json.map(json.loads)
idt_rejection_reasons = rejection_rows.explode('idt_reason').groupby(['collection','fragment_limit_bp','idt_reason']).size().rename('rejections').reset_index()
maximum_columns = ['module_id','collection','family','in_designed_primary100','fragment_limit_bp','unit_sequence','unit_length','mathematical_max_copies','pre_adaptive_verified_max_copies','verified_max_copies','adaptive_search_evaluations','adaptive_idt_scored_evaluations','ga_generations','adaptive_stop_reason','adaptive_boundary_proven','adaptive_boundary_evidence','adaptive_orderable_passed','plasmid','direction','site_i_position','site_ii_position','site_i_enzyme','site_ii_enzyme','site_iii_enzymes','dna_sequence','idt_status','idt_violation_count','idt_scored_sequence_sha256']
maximum_constructs = constructs.loc[constructs.final_passed.fillna(False) & constructs.adaptive_boundary_proven.fillna(False), maximum_columns].copy()
adaptive_capacity.to_csv(Path(SUMMARY_DIR)/'adaptive_maximum_summary.csv', index=False)
adaptive_stops.to_csv(Path(SUMMARY_DIR)/'adaptive_stop_reasons.csv', index=False)
generation_usage.to_csv(Path(SUMMARY_DIR)/'adaptive_generation_usage.csv', index=False)
idt_rejection_reasons.to_csv(Path(SUMMARY_DIR)/'adaptive_idt_rejection_reasons.csv', index=False)
rejection_rows.to_csv(Path(SUMMARY_DIR)/'adaptive_idt_feedback_trace.csv', index=False)
maximum_constructs.to_parquet(Path(SUMMARY_DIR)/'maximum_passed_constructs_notebook.parquet', index=False)
maximum_constructs.to_csv(Path(SUMMARY_DIR)/'maximum_passed_constructs_notebook.csv', index=False)
run_context['row_counts']['adaptive_trace_evaluations'] = len(trace)
run_context['row_counts']['idt_rejected_evaluations'] = int(trace_compare.idt_rejected.sum())
run_context['row_counts']['idt_feedback_updates'] = int(trace_compare.feedback_applied.sum())
run_context['row_counts']['maximum_constructs_with_dna'] = len(maximum_constructs)
run_context['row_counts']['boundary_proven_construct_caps'] = int(constructs.adaptive_boundary_proven.fillna(False).sum())
run_context['filter_flow'].extend(['set upper bound to floor((fragment cap - deduction)/(3*unit length))', 'binary search using the local plus IDT gate at 10 generations', 'after IDT rejection map its reasons to GA weights and retry the same copy count', 'from one copy above the short-search maximum, add exactly one module at a time', 'for each new copy count escalate generations through 10,20,40,60,80,100', 'accept a maximum only when it is IDT-orderable and reaches the mathematical cap, or the next copy remains non-orderable at 100 generations', 'return the explicitly orderable maximum count and exact DNA'])
run_context['limitations'].append('binary search treats short-generation orderability as monotone; the audited one-by-one phase and explicit 100-generation local plus IDT boundary proof determine the reported maximum')
adaptive_capacity"""),
                code("fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), facecolor='white')\nlong_capacity = constructs_compare.melt(id_vars=['collection','fragment_limit_bp'], value_vars=['pre_adaptive_verified_max_copies','verified_max_copies'], var_name='search_stage', value_name='copies')\nsns.boxplot(data=long_capacity, x='collection', y='copies', hue='search_stage', order=order, ax=axes[0], palette=['#B7A57A','#4B2E83'])\nsns.boxplot(data=constructs_compare.loc[constructs_compare.orderable_maximum], x='collection', y='ga_generations', hue='fragment_limit_bp', order=order, ax=axes[1], palette=['#4B2E83','#B7A57A'])\nsns.barplot(data=adaptive_capacity, x='collection', y='reached_mathematical_bound', hue='fragment_limit_bp', order=order, ax=axes[2], palette=['#2D7DD2','#F45D01'])\nfor ax in axes: ax.tick_params(axis='x', rotation=20)\naxes[0].set_title('Legacy maximum vs adaptive maximum'); axes[1].set_title('Generations used by winning DNA'); axes[2].set_title('Modules reaching fragment-length ceiling')\nsns.despine(); fig.tight_layout()\nfor suffix in ('png','pdf'): fig.savefig(Path(FIGURE_DIR) / f'adaptive_maximum_search.{suffix}', dpi=300, facecolor='white')\nfig"),
                code("run_context"),
            ],
        ),
    )
    write(
        REPO / "notebooks/tasks/05_module_boundary_inference.ipynb",
        notebook(
            "Full-protein sequence and secondary-structure repeat boundaries",
            "Audits source units against complete protein sequences and residue-level DSSP/author secondary structure. The real middle copy of the jointly supported primitive repeat is passed to HURDLER and codon optimization; structure-derived RepeatsDB units remain authoritative when a shorter natural harmonic lacks 3D confirmation.",
            common
            + f"\nBOUNDARY_TABLE = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_boundary_audit.parquet')!r}"
            + f"\nCANDIDATE_TABLE = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_period_candidates.parquet')!r}"
            + f"\nUNIT_TABLE = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_unit_alignment.parquet')!r}"
            + f"\nPOSITION_TABLE = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_position_variability.parquet')!r}"
            + f"\nSS_CANDIDATE_TABLE = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_secondary_structure_candidates.parquet')!r}"
            + f"\nSS_RESIDUE_TABLE = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_secondary_structure_residues.parquet')!r}"
            + f"\nVALIDATION = {str(STUDY / 'step03_module_corpus/tables/periodic_v4/module_boundary_validation.json')!r}"
            + f"\nFIGURE_DIR = {str(STUDY / 'step03_module_corpus/figures/periodic_v4')!r}",
            [
                code("boundary = pd.read_parquet(BOUNDARY_TABLE)\ncandidates = pd.read_parquet(CANDIDATE_TABLE)\nunits = pd.read_parquet(UNIT_TABLE)\npositions = pd.read_parquet(POSITION_TABLE)\nss_candidates = pd.read_parquet(SS_CANDIDATE_TABLE)\nss_residues = pd.read_parquet(SS_RESIDUE_TABLE)\nvalidation = json.loads(Path(VALIDATION).read_text())\nrun_context['input_hashes'] = {Path(path).name: sha256(path) for path in (BOUNDARY_TABLE, CANDIDATE_TABLE, UNIT_TABLE, POSITION_TABLE, SS_CANDIDATE_TABLE, SS_RESIDUE_TABLE, VALIDATION)}\nrun_context['row_counts'] = {'proteins': len(boundary), 'period_candidates': len(candidates), 'secondary_structure_candidate_scores': len(ss_candidates), 'secondary_structure_residues': len(ss_residues), 'aligned_units': len(units), 'module_positions': len(positions), 'natural': int(boundary.collection.eq('natural100').sum()), 'designed': int(boundary.collection.eq('designed_all').sum())}\nrun_context['filter_flow'] = ['read the complete protein/construct sequence', 'map RepeatsDB author chains exactly to RCSB label chains', 'run DSSP on natural and THR coordinates or use the DHR author H-loop-H-loop residue template', 'score every candidate independently by amino-acid Fourier/self-similarity and H/E/C state/transition periodicity', 'select the smallest harmonic passing both evidence gates', 'choose the real middle repeat copy, with an exact central tie resolved toward the earlier copy', 'pass that middle AA sequence to HURDLER and codon optimization', 'call positions fixed when conservation across inferred copies is at least 0.8']\nrun_context['limitations'] = ['a shorter natural harmonic is reported for review but does not replace a structure-derived RepeatsDB unit without 3D superposition', 'DHR secondary structure comes from the author residue-count template; THR and natural entries use DSSP', 'missing or mismatched structure chains are explicit failures, never inferred annotations']\nvalidation"),
                code("summary = boundary.groupby('module_type').agg(proteins=('module_id','nunique'), median_source_length=('prior_unit_length','median'), median_primitive_length=('primitive_period','median'), split_source_units=('length_ratio', lambda values: int((values > 1.08).sum())), secondary_structure_passed=('secondary_structure_status', lambda values: int((values == 'passed').sum())), jointly_selected=('secondary_structure_selected_support', lambda values: int(values.fillna(False).sum())), manual_review=('qa_flags', lambda values: int((values != '').sum()))).reset_index()\nsummary"),
                code("boundary[['module_id','module_type','prior_unit_start','prior_unit_end','repeat_region_start','repeat_region_end','first_module_start','first_module_end','selected_module_index','selected_module_start','selected_module_end','unit_sequence','prior_unit_length','primitive_period','repeat_count','periodicity_score','secondary_structure_known_fraction','secondary_structure_selected_support','selection_reason','qa_flags']].sort_values(['module_type','module_id'])"),
                code("from IPython.display import Image, display\nfor name in ['source_vs_primitive_module_length.png','module_harmonic_ratio.png','module_fixed_fraction.png','secondary_structure_coverage.png','sequence_vs_secondary_structure_evidence.png','module_boundary_examples.png']:\n    path = Path(FIGURE_DIR) / name\n    assert path.stat().st_size > 0\n    display(Image(filename=str(path)))"),
                markdown("The colored interval is the inferred repeat region within the full protein. White ticks are module boundaries; the black interval is the source unit. The AA sequence used downstream is the real middle copy of that region, with an exact central tie resolved toward the earlier copy. Designed harmonics are accepted only when both amino-acid and residue-level secondary-structure periodicity pass. RepeatsDB's structure-derived natural boundary is retained unless later 3D superposition establishes a smaller complete unit."),
                code("run_context"),
            ],
        ),
    )
    write(
        REPO / "notebooks/tasks/04_reproducibility_status.ipynb",
        notebook(
            "Repository reproducibility matrix",
            "Reports maintained, GUI, archive, blocked, deferred, and failed workflows without fabricating missing historical inputs.",
            common + f"\nSTATUS_TABLE = {str(STUDY / 'step05_reproducibility/tables/execution_status.csv')!r}",
            [
                code("status = pd.read_csv(STATUS_TABLE)\nrun_context['input_hashes']['execution_status.csv'] = sha256(STATUS_TABLE)\nrun_context['row_counts']['workflow_rows'] = len(status)\nrun_context['filter_flow'] = ['inventory canonical, maintained, adjacent, and archive notebook/script entry points', 'retain exact pass/block/defer/fail evidence']\nrun_context['limitations'] = ['blocked_missing_input results are never imputed']\nstatus.groupby(['status','workflow_type']).size().rename('count').reset_index()"),
                code("status[['workflow','status','job_id','runtime_seconds','output_path','rerun_command']].sort_values(['status','workflow'])"),
            ],
        ),
    )
    print("created 2 reference notebooks and 5 task notebooks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
