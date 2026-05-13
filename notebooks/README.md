# Notebooks

Topic-grouped Jupyter notebooks. Each folder mirrors a workflow that is
documented in [`../docs/workflows/`](../docs/workflows/) and (where
applicable) backed by a Python module under [`../src/hurdler/`](../src/hurdler/).

| Folder | Purpose | Backing code / docs |
|--------|---------|---------------------|
| [`hurdler/`](hurdler/) | HURDLER three-site combinations and success-rate analysis. | `src/hurdler/pipeline.py`, `src/hurdler/success_rate.py`, `docs/workflows/hurdler_*.md` |
| [`enzyme_selection/`](enzyme_selection/) | Selection and pairing analyses for restriction enzymes (3-mer coverage, plasmid compatibility, Sankey summaries). | `docs/workflows/enzyme_selection.md` |
| [`codon_optimization/`](codon_optimization/) | Codon-optimisation experiments and benchmark notebooks. | `codon_opt_benchmark_extended/README.md` |
| [`agarose_gel/`](agarose_gel/) | Agarose-gel image quantification. | `agarose_gel_analysis/` |
| [`sec/`](sec/) | Size-exclusion chromatography processing. | `SEC/` |
| [`utils/`](utils/) | One-off notebooks that curate `data/reference_output/` from raw databases (REBASE, NEB, codon usage, methylation, plasmids). | `data/reference_input/` → `data/reference_output/` |

## Conventions

- Each notebook should start with a markdown cell stating: purpose,
  required inputs (paths under `data/`), expected outputs (paths under
  `output/`), and the matching module/script (if any).
- Notebooks should call into `src/hurdler/` rather than re-implementing
  stable logic.
- Backup, executed, and exploratory copies belong in
  [`../archive/notebooks/`](../archive/notebooks/), not here.
