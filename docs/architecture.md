# Architecture

This repository supports the design and validation of the **HURDLER**
cloning strategy for repeat proteins, plus several adjacent experimental
analyses (codon optimisation, agarose-gel quantification, size-exclusion
chromatography, plasmid sequencing).

## Top-level layout

```
clone_repeat_protein/
├── README.md                  # High-level overview, quick start
├── LICENSE
├── envs/                      # Conda environment files
├── src/                       # Maintained Python source (hurdler package)
│   └── hurdler/
│       ├── pipeline.py        # Full data-generation pipeline (df1, df2, lookup)
│       ├── query.py           # CLI/library to query valid combinations
│       ├── validate.py        # Data-quality validation and statistics
│       ├── success_rate.py    # Success-rate analysis over random sequences
│       └── utils.py           # Shared enzyme/plasmid compatibility helpers
├── notebooks/                 # Topic-grouped Jupyter notebooks (analyses)
│   ├── hurdler/
│   ├── enzyme_selection/
│   ├── codon_optimization/
│   ├── agarose_gel/
│   ├── sec/
│   └── utils/                 # Reference-data generation notebooks
├── tests/                     # Smoke and regression tests
├── scripts/                   # CLI entry points / helpers (optional)
├── docs/                      # All maintained documentation
│   ├── architecture.md
│   ├── glossary.md
│   ├── contributing.md
│   ├── workflows/             # Per-workflow guides
│   └── reports/               # Historical writeups (read-only)
├── data/                      # Small, committed reference inputs/outputs
│   ├── reference_input/       # Source databases (REBASE, NEB, plasmids, …)
│   ├── reference_output/      # Curated reference CSVs derived from inputs
│   ├── hurdler_analysis_input/# Inputs for the HURDLER pipeline
│   └── example_batch_query.csv
├── output/                    # GIT-IGNORED generated artifacts (CSVs, plots)
├── codon_opt_benchmark_extended/  # Subproject: codon-optimisation benchmark
├── agarose_gel_analysis/      # Subproject: agarose-gel quantification
├── SEC/                       # Subproject: SEC chromatogram processing
├── plasmid_sequencing_result/ # Reference sequencing artifacts (Genbank, HTML)
└── archive/                   # Read-only historical/duplicate files
```

## Data flow (HURDLER)

```
data/reference_input/         (databases: REBASE, plasmids, …)
        │
        ▼
notebooks/utils/*.ipynb       (one-off curation; ouputs to data/reference_output/)
        │
        ▼
data/reference_output/        (curated enzyme/plasmid CSVs)
        │
        ▼
src/hurdler/pipeline.py       (build df1, df2, lookup)
        │                     │
        ▼                     ▼
output/*.csv, *.pkl, *.png    src/hurdler/query.py, validate.py, success_rate.py
                              (consume generated artifacts)
```

- **Reference inputs** in `data/reference_input/` are stable; only update
  them when the source databases change.
- **Reference outputs** in `data/reference_output/` are derived but small
  enough to commit; they back the maintained workflows.
- **Generated outputs** (large CSVs, pickles, PDFs) live in `output/` and
  are ignored by Git. They are reproducible by re-running the pipeline.
- **Archived outputs** in `archive/artifacts/` are historical snapshots kept
  only for traceability.

## Code/notebook responsibility split

- **Python modules under `src/hurdler/`** hold the stable, importable logic.
- **Notebooks under `notebooks/`** orchestrate, visualise, and document the
  analyses. They should call into `src/hurdler/` rather than re-implementing
  logic.
- **Tests under `tests/`** smoke-test the pipeline and query layer.

## Subprojects

The following subprojects are self-contained and have their own inputs,
outputs, and notebooks. They share only reference data with the HURDLER
core:

- `codon_opt_benchmark_extended/` — genetic-algorithm codon optimisation.
- `agarose_gel_analysis/` — agarose-gel band quantification (GUI + notebooks).
- `SEC/` — size-exclusion chromatography processing.
- `plasmid_sequencing_result/` — Genbank + annotated HTML for sequenced
  constructs.
