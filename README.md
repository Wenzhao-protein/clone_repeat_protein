# HURDLER repeat-protein cloning

**[Open the current preview in Colab](https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/agent/vector-aware-designer-v2/notebooks/workflows/02_colab_hurdler_designer.ipynb) · [Run the local web designer](#local-web-designer)**

[![Open current preview in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/agent/vector-aware-designer-v2/notebooks/workflows/02_colab_hurdler_designer.ipynb)

## Test in a browser

- **Colab:** use the badge above to open the review branch immediately. The
  notebook clones that exact branch, installs the locked project extras, and
  exposes the protein defaults as native Colab Forms. Choose **Runtime → Run
  all** once to install the package and reveal separate individual-enzyme,
  plasmid, RE-solution, and cut-scheme selectors. Review the route table and
  explicitly select Site I, Site II, Site III, plasmid, and cut scheme. GA/IDT
  controls appear only after that route is confirmed. Full-protein inputs without confirmed boundaries
  intentionally stop after proposing boundary candidates. Code remains hidden
  unless **Show code** is selected.
  After this branch is merged, the
  permanent [main-branch Colab link](https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/main/notebooks/workflows/02_colab_hurdler_designer.ipynb)
  provides the same entry point.
- **Local browser:** follow the commands below. The final script starts the
  Marimo app on `127.0.0.1:2718` and opens it in your browser.

Both interfaces use the same `hurdler.vector_design` implementation.
Credentials and IDT scoring stay inside the active Colab or local process; the
software produces design files only and never submits an order.

Live IDT scoring is the Colab default. A local runtime automatically tries
`~/.config/hurdler/idt.env`; hosted Colab instead displays a temporary env-file
upload. The file contains either `IDT_ACCESS_TOKEN`, or all of
`IDT_CLIENT_ID`, `IDT_CLIENT_SECRET`, `IDT_USERNAME`, and `IDT_PASSWORD`.
Store it outside the repository with mode 600. Colab Secrets and an invisible
runtime prompt remain alternatives; credentials are never notebook form parameters.

## Local web designer

```bash
git clone https://github.com/Wenzhao-protein/clone_repeat_protein.git
cd clone_repeat_protein
conda env create -f envs/hurdler.yml
conda activate hurdler
python -m pip install -e ".[notebooks,optimization]"
./scripts/start_hurdler_web.sh
```

The launcher checks that the active Python environment contains HURDLER and
Marimo before starting. Advanced users can select another port or suppress the
browser launch with `--port PORT` or `--no-browser`.

The public [Natural/Designed repeat-protein result catalog](data/results/README.md)
contains one manually searchable row for every active middle module. It includes
source and boundary evidence, Stage-1 HURDLER compatibility, the selected
plasmid/enzyme route, and the independently validated 1,800-bp and 3,000-bp
maximum constructs. An `*_idt_accepted_dna` cell is populated only when the
exact DNA translated correctly, had no selected-pair excess sites, carried a
maximum-copy proof, and received a live IDT rule-score sum strictly below 10.

The v2 interactive designer is the fastest user entry point. It first queries
all 776 protein-level Site-I/Site-II pairs, then evaluates the retained long
backbone of seven annotated physical vectors (eight selectable profiles) under
four MCS cut schemes. This avoids the old whole-plasmid prefilter. Default input
is `N-cap + module × n + C-cap`; complete-protein mode preserves every repeat
variant. Each route records restoration segments and the strict annotation-
aware cutter-silencing decision. For split inputs, the Colab workflow treats
the requested copy count as exact. It independently optimizes a cap-containing
primary and one reusable secondary, then proves
`primary + rounds × secondary = target` at every RDL intermediate. Live
progress reports fragment, copy count, generation, GA score, IDT result, and
elapsed time. Every successful bundle also contains IDT Bulk Input
CSV/TSV/FASTA. It never uses IDT codon optimization and never submits an order.
Run the local page with `hurdler web`.

This repository implements and validates the HURDLER three-site cloning
strategy for repeat proteins. The maintained path is an installable Python
package, thin parameterized notebooks, and recoverable Digs task manifests.
Historical notebooks and adjacent experimental workflows remain available for
traceability but do not define the current scientific result.

## Frozen result contract

The primary rule profile is `legacy-optimized-v1`:

- Site-I/Site-II orthogonality threshold: 1; missing fidelity is compatible.
- Site-II/Site-III use the same signed overhang.
- Plasmid compatibility requires both Site I and Site II.
- The active success landscape scans `module + module + module` with
  `5 <= distance < module_length`; the historical two-copy scan is retained as
  an immutable comparison baseline.
- 7--60AA Monte Carlo results preserve seed 42 and the historical loop order.
- 1--5AA results exhaust all `20**k` motifs without equivalence reduction.

The active module corpus is
`expanded-middle-repeatsdb-foldseek-v1`. Natural and designed boundaries have
different evidence contracts. Natural units are sliced directly from official
RepeatsDB PDB or AlphaFoldDB annotation coordinates: one longest annotated
region is selected per biological protein and unit index `(count - 1) // 2`
provides the earlier middle unit. DSSP/Foldseek are natural-QC fields only and
can never change those coordinates. Designed units have no annotation
fallback: eight-state Biotite/mkdssp periodicity and Foldseek 3Di plus
fragment-level TM/LDDT validation must independently support the same period
and repeat block. MAFFT then defines the fixed and variable positions, and the
earlier middle copy is passed to HURDLER.

The exact short-motif counts are 20, 400, 8,000, 160,000, and 3,200,000
(3,368,420 total). Each motif is repeated to the shortest module of at least
6AA before matching.

The maintained success-landscape entrypoints are
`scripts/run_success_landscape_single_files.py` and
`scripts/run_success_landscape_16core.sh`. The Digs wrapper requests 16 cores,
benchmarks serial and multiprocessing results for exact equality, and writes
the raw observations to exactly two files: `short_motifs_1_5.parquet`
(3,368,420 rows) and `random_modules_6_60.parquet` (440,000 rows). The 6AA
block uses seed 420006; 7--60AA preserve seed 42 and the historical
length→plasmid→test order. Notebook `02_success_rate_1_60.ipynb` derives the
eight-plasmid three-copy curve directly from those two files and compares it
row-for-row with the historical two-copy result. The pattern population,
distance rule and first-match behavior are unchanged. The third copy only
completes cyclic 3-mer windows that cannot be represented at the right edge of
the doubled string. The figure follows the original continuous-line settings
on a near-square 6 × 5 inch canvas, displays 1–50AA, uses the exact title
`3-mer Probability vs Sequence Length`, and contains no confidence band,
vertical method divider, or extra method annotation. After the replacement
files passed independent validation, the
superseded 404-file motif and hit shard directories were removed; the legacy
combined Parquet and two-copy single-file run remain available for regression
audit.

## Maintained layout

```text
src/hurdler/                 importable algorithms and unified schemas
notebooks/reference/         reference manifests and sparse-index QC
notebooks/tasks/             HURDLER queries and scientific reports
notebooks/workflows/         interactive end-user construct design
data/artifacts/              committed complete legacy-optimized-v1 lookup
data/reference_{input,output}/ versioned small reference data
envs/hurdler.yml             portable historical-style conda specification
envs/hurdler-linux-64.lock   explicit environment lock
studies/hurdler_validation/  Digs manifests, summaries, HTML, figures, status
tests/                       unit, property-style, and golden regressions
archive/                     read-only historical implementations
```

Large indexes, exhaustive Parquet tables, and downloads use the exact
`/net/scratch/wendai/projects/hurdler/clone_repeat_protein/...` mirror.
Reviewable code, task files, summary tables, and reports remain under `/home`.
Module HURDLER and adaptive-copy shards use the shared `/net/scratch` mirror;
only compact finalized tables, figures, notebooks, manifests, and reports are
promoted back to `/home`.

## Install and test

```bash
/net/software/conda/bin/conda env create \
  --prefix /home/wendai/.conda/envs/hurdler \
  --file envs/hurdler.yml
/home/wendai/.conda/envs/hurdler/bin/pip install -e .
/home/wendai/.conda/envs/hurdler/bin/python -m pytest -q
```

The executable exposes one artifact contract across all operations:

```bash
hurdler reference build --help
hurdler lookup build --rules legacy-optimized-v1 --help
hurdler lookup protein-build --help
hurdler plasmid-reference build --help
hurdler plasmid-reference validate --help
hurdler query --module VLA --help
hurdler screen-short --help
hurdler success-rate --help
hurdler curate-modules --help
hurdler infer-boundaries --help
hurdler designed-inventory --help
hurdler infer-designed-boundaries --help
hurdler module-compatibility --help
hurdler adaptive-copy-search --help
hurdler optimize-modules --help
hurdler refine-ga --help
hurdler design-construct --request request.json --output-dir output/my_design --help
hurdler design-query --request examples/vector_aware_query.json --help
hurdler web --help
hurdler dna-assembly build-corpus --help
hurdler dna-assembly plan --help
hurdler dna-assembly finalize --help
hurdler validate-run --help
```

The exact-DNA assembly workflow treats a functional donor shorter than 90 bp
as an annealed oligo product: it emits two complementary 5′→3′ primer
sequences with the required sticky ends exposed and does not call the IDT
gBlocks complexity endpoint for that fragment. Donors of 90 bp or longer
remain synthesis fragments and retain live IDT scoring.

A full-protein compatibility example is ready at
[`data/example_design_request.json`](data/example_design_request.json). Set
`optimize` to `true` only after confirming the coordinates and configuring the
live IDT scorer. It uses the 26-heptad *S. cerevisiae* Rpb1 CTD plus its
C-terminal tip. The confirmed 1--182 repeat region selects the earlier middle
heptad `YSPTSPS`; the frozen index reports that heptad incompatible across all
eight maintained plasmids, so no orderable DNA is produced. This negative
golden result is intentional rather than silently switching to a different
module or DNA-derived site.

`hurdler design-construct` and the interactive notebook share the strict
`DesignRequest`/`DesignResult` interface. A request without confirmed 1-based
inclusive repeat coordinates stops at `needs_boundary_confirmation`. A
confirmed request with optimization disabled produces only a topology draft
marked `not_orderable_not_for_purchase`; it never writes orderable DNA.

`hurdler refine-ga` includes repeated restriction-site occurrences directly in
the genetic fitness score. Reverse-complement aliases are canonicalized so one
physical site is not double-penalized. Translation and locked HURDLER codons are
invariant; selected-enzyme site counts are the local hard check. GC remains a
high-weight GA objective and a live-IDT rule rather than a second local veto.
Non-selected repeated RE sites are a soft optimization objective and never
prevent an otherwise locally valid candidate from being scored by IDT. This
contract is versioned as
`nonselected-re-sites-soft-score-selected-sites-hard-v2`. With
`--use-idt`, every locally acceptable candidate DNA is submitted to IDT's
gBlocks complexity screener before a copy-count decision is made. No
IDT-generated replacement sequence is requested or adopted. A position-matched
empty rule list has score zero. Otherwise every finite numeric rule `Score` is
summed and the exact GA DNA passes only when the total is strictly below 10
(`idt-rule-score-sum-lt10-v1`). `IsViolated` is diagnostic and never gates the
sequence by itself. Positive-score rule names, actual values and thresholds are
retained and mapped back to GA weights for GC, hairpins/palindromes,
homopolymers, terminal repeats, restriction sites, and 8/13/14-mer repeats.
Missing or non-numeric scores are unclassified and cannot pass. Credentials are
loaded only from a user-selected, repo-external mode-600 env file; values and
the resolved private path are never copied to code,
notebooks, manifests, logs, or Git.

Interactive purchase fragments enforce IDT's current
[125--3,000 bp gBlocks range](https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments).
Shorter 20--124 bp fragments are labeled as duplexed Ultramer candidates rather
than gBlocks; unsupported lengths cannot pass the purchase design.

For maximum-repeat searches, `hurdler refine-ga --adaptive-copy-search` uses
the fragment-length mathematical maximum as its upper bound. It first performs
a 10-generation binary search, then starts one copy above the short-search
maximum and increases the construct by exactly one module. Each new count is
attempted at 10, 20, 40, 60, 80, and 100 generations. A locally acceptable but
IDT-rejected candidate stays at the same copy count; its rejection reasons
raise the corresponding score weights before the next generation budget. Only
an explicit score-sum-below-10 IDT result permits the next module. The first count
that is still not orderable at 100 stops the search; the preceding orderable
count and its exact DNA are written to the versioned maximum-copy tables and
FASTA. The same route is applied to
every HURDLER-compatible natural and designed module.

The merge rejects any HURDLER-compatible row without an explicit boundary
proof: either the orderable count equals the mathematical fragment cap, or the
next count has a recorded failed local+IDT evaluation at exactly 100 generations. A
10/20/40/60/80-generation failure can never define the reported maximum.

The historical 249-row `periodic_v4` run and its parallel benchmark remain an
immutable baseline. New production shard counts are derived from the exhaustive
RepeatsDB and strict-dual-evidence designed inventories; no 100/149/249 row cap
is carried into the active corpus.

For Digs, copy `config/idt.env.example` to a repo-external private location,
fill it, and set mode `600`. Runtime taskfiles may receive that path through an
environment-specific launch configuration; the wrapper refuses group/world-readable
credential files and never prints their contents. The interactive designer
also accepts hidden manual OAuth fields, a hidden prompt for any repo-external
mode-600 env file, or a temporary in-memory env upload. Controls are cleared
immediately; credential values, upload contents, and the actual path are not
written to notebook state, design outputs, manifests, logs, or Git.

## Canonical reports

The completed historical handoff is
[`studies/hurdler_validation/FINAL_REPORT_PERIODIC_V4_MIDDLE.md`](studies/hurdler_validation/FINAL_REPORT_PERIODIC_V4_MIDDLE.md).
It records the authoritative artifacts, middle-module audit, HURDLER and IDT
outcomes, maximum-copy constructs, Digs job lineage, and the exact remaining
blocked/deferred workflows. The per-module machine-readable deliverable is
`step04_module_optimization/tables/periodic_v4/module_final_summary.parquet`;
it contains all 249 legacy tested AA sequences. It is retained for audit and is
not input to `expanded-middle-repeatsdb-foldseek-v1`.

- `notebooks/reference/01_reference_manifest.ipynb`
- `notebooks/reference/02_lookup_qc.ipynb`
- `notebooks/tasks/01_hurdler_query.ipynb`
- `notebooks/tasks/02_success_rate_1_60.ipynb`
- `notebooks/tasks/03_repeat_module_benchmark.ipynb`
- `notebooks/tasks/04_reproducibility_status.ipynb`
- `notebooks/tasks/05_module_boundary_inference.ipynb`
- `notebooks/reference/03_expanded_middle_module_acquisition.ipynb`
- `notebooks/tasks/06_module_compatibility_by_length.ipynb`
- `notebooks/tasks/07_adaptive_selected_pair_capacity.ipynb`
- `notebooks/tasks/08_long_repetitive_dna_assembly.ipynb`
- `notebooks/workflows/01_interactive_hurdler_designer.ipynb`

Notebook 08 uses `arbitrary-dna-complete-route-v2`: each public regulatory
element retains independent 2/4/8/16/32-copy outcomes, and a pass requires a
verified path from a purchasable exact seed. The earlier 53.67% one-step
baseline remains immutable QC and is excluded from reviewer conclusions.

Papermill-executed notebooks and HTML live under
`studies/hurdler_validation/step05_reproducibility/`. The complete run status,
including missing historical `.scn` inputs and deferred long workflows, is in
`step05_reproducibility/tables/execution_status.csv`.

The active method/corpus version is
`expanded-middle-repeatsdb-foldseek-v1`. `periodic_v4_middle_unit`,
`periodic_v3`, first-unit, and pre-run48 files are immutable legacy artifacts
and must not be combined with the new tables.

Adjacent SEC, agarose-gel, codon-optimization, and plasmid-sequencing projects
retain their existing top-level directories. See `docs/architecture.md` and
`studies/hurdler_validation/README.md` for data flow and run details.

## License

See [LICENSE](LICENSE).
