# clone_repeat_protein

Tools and analyses for cloning **repeat proteins** using the **HURDLER**
three-restriction-site strategy, plus several adjacent experimental
workflows used in the same project (codon optimisation, agarose-gel
quantification, size-exclusion chromatography, plasmid sequencing).

The repository is organised so that:

- a **new contributor** can understand the project from this README,
- a **maintainer** can locate canonical code, notebooks, and docs without
  guessing, and
- **experimental history** is preserved but clearly separated from the
  maintained path under [`archive/`](archive/).

---

## What lives where

```
README.md                  ← you are here
LICENSE
envs/                      Conda environment files (codon_opt, visualization)
src/hurdler/               Maintained Python package (pipeline, query, validate, …)
notebooks/                 Guided analyses, grouped by topic
tests/                     Smoke / regression tests
scripts/                   CLI helpers
docs/                      All maintained documentation (start at docs/README.md)
data/                      Small reference inputs and curated reference outputs
output/                    Generated artifacts (GIT-IGNORED)
codon_opt_benchmark_extended/  Subproject: codon-optimisation benchmark
agarose_gel_analysis/      Subproject: agarose-gel quantification
SEC/                       Subproject: size-exclusion chromatography
plasmid_sequencing_result/ Reference sequencing data (Genbank / annotated HTML)
archive/                   Historical, duplicate, and exploratory files (read-only)
```

For a deeper view see [`docs/architecture.md`](docs/architecture.md).
For terminology see [`docs/glossary.md`](docs/glossary.md).

---

## Subprojects

| Subproject | Status | Where to start |
|------------|--------|----------------|
| **HURDLER three-site combinations** | Stable | [`docs/workflows/hurdler_site_combinations.md`](docs/workflows/hurdler_site_combinations.md) → [`notebooks/hurdler/hurdler_site_combination_analysis.ipynb`](notebooks/hurdler/hurdler_site_combination_analysis.ipynb) |
| **HURDLER success-rate analysis** | Stable | [`docs/workflows/hurdler_success_rate.md`](docs/workflows/hurdler_success_rate.md) → [`notebooks/hurdler/hurdler_success_rate_analysis.ipynb`](notebooks/hurdler/hurdler_success_rate_analysis.ipynb) |
| **Restriction-enzyme selection** | Stable | [`docs/workflows/enzyme_selection.md`](docs/workflows/enzyme_selection.md) → [`notebooks/enzyme_selection/`](notebooks/enzyme_selection/) |
| **Reference-data curation** | Stable | [`notebooks/utils/`](notebooks/utils/) (writes into `data/reference_output/`) |
| **Codon-optimisation benchmark** | Exploratory | [`codon_opt_benchmark_extended/README.md`](codon_opt_benchmark_extended/README.md), [`notebooks/codon_optimization/`](notebooks/codon_optimization/) |
| **Agarose-gel analysis** | Exploratory | [`agarose_gel_analysis/`](agarose_gel_analysis/), [`notebooks/agarose_gel/`](notebooks/agarose_gel/) |
| **SEC analysis** | Exploratory | [`SEC/`](SEC/), [`notebooks/sec/`](notebooks/sec/) |
| **Plasmid sequencing results** | Reference data | [`plasmid_sequencing_result/`](plasmid_sequencing_result/) |

---

## Quick start

### 1. Create an environment

```bash
conda env create -f envs/codon_opt.yml          # core analysis + biopython
conda activate codon_opt
# Optional, for some visualisation notebooks:
conda env create -f envs/visualization.yml
```

Required Python packages (already declared in the env files):
`pandas`, `numpy`, `biopython`, `matplotlib`, `seaborn`, `tqdm`,
`scipy`.

### 2. Inspect curated reference data

Reference enzyme/plasmid CSVs derived from REBASE and NEB live in
[`data/reference_output/`](data/reference_output/). They back every
HURDLER workflow and are committed to the repository.

### 3. Run the HURDLER pipeline (generates `output/`)

```bash
# Full pipeline (df1, df2, lookup, success-rate plots)
python -m hurdler.pipeline                # from repo root, with src/ on PYTHONPATH
# or
PYTHONPATH=src python src/hurdler/pipeline.py
```

This populates `output/` with the generated CSVs, pickles, and plots
described in [`docs/workflows/hurdler_site_combinations.md`](docs/workflows/hurdler_site_combinations.md).

### 4. Query valid three-site combinations

```bash
PYTHONPATH=src python -m hurdler.query \
  --site-i-aa  "NEQ" \
  --site-ii-aa "IQA" \
  --plasmid    "pET-28a(+)"
```

Batch queries:

```bash
PYTHONPATH=src python -m hurdler.query \
  --batch  data/example_batch_query.csv \
  --output output/my_results.csv
```

### 5. Validate generated data

```bash
PYTHONPATH=src python -m hurdler.validate
```

---

## Recommended paths

### New users

1. Read this README.
2. Skim [`docs/architecture.md`](docs/architecture.md) and
   [`docs/glossary.md`](docs/glossary.md).
3. Open the canonical HURDLER notebook
   [`notebooks/hurdler/hurdler_site_combination_analysis.ipynb`](notebooks/hurdler/hurdler_site_combination_analysis.ipynb)
   and run it top to bottom.

### Maintainers

1. Read [`docs/contributing.md`](docs/contributing.md) for file-placement
   and naming rules.
2. Make changes to logic under [`src/hurdler/`](src/hurdler/); update the
   corresponding notebook only after the module change is stable.
3. Keep generated artifacts in `output/` (gitignored) and reference data
   in `data/`.
4. Anything you retire should move into [`archive/`](archive/) — never
   delete history.

---

## Status

- The maintained code path is **`src/hurdler/`** plus the notebooks under
  **`notebooks/`** that load reference data from **`data/`**.
- Documentation in **`docs/workflows/`** is the source of truth.
  Older summaries under **`docs/reports/`** are preserved for context
  only.
- Everything under **`archive/`** is read-only history. See
  [`archive/README.md`](archive/README.md) for what was moved and why.

## License

See [`LICENSE`](LICENSE).
