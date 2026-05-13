# Documentation Index

This folder is the entry point for all detailed documentation. The root
[`README.md`](../README.md) is the high-level overview; everything deeper
lives here.

## Start here

- [`architecture.md`](architecture.md) — How the code, notebooks, reference
  data, and generated outputs relate.
- [`glossary.md`](glossary.md) — Domain terminology (HURDLER, Site I/II/III,
  3-mer AA, overhang, methylation compatibility, Type IIS, …).
- [`contributing.md`](contributing.md) — File-placement and naming rules for
  new contributions.

## Workflow guides

Canonical, maintained workflows live in [`workflows/`](workflows/):

| Workflow | Document | Entry points |
|----------|----------|--------------|
| HURDLER three-site combinations (Site I + Site II + Site III) | [`workflows/hurdler_site_combinations.md`](workflows/hurdler_site_combinations.md) | `notebooks/hurdler/hurdler_site_combination_analysis.ipynb`, `src/hurdler/pipeline.py`, `src/hurdler/query.py`, `src/hurdler/validate.py` |
| HURDLER success-rate analysis | [`workflows/hurdler_success_rate.md`](workflows/hurdler_success_rate.md) | `notebooks/hurdler/hurdler_success_rate_analysis.ipynb`, `src/hurdler/success_rate.py` |
| Restriction-enzyme selection and pairing | [`workflows/enzyme_selection.md`](workflows/enzyme_selection.md) | `notebooks/enzyme_selection/*.ipynb` |

## Reports archive

Historical writeups, fix summaries, and one-off analysis reports are
preserved in [`reports/`](reports/). They are kept for traceability but
should **not** be used as onboarding material — start with the workflow
guides above instead.

## Subprojects with their own docs

- `codon_opt_benchmark_extended/README.md` — Codon-optimisation benchmark.
- `agarose_gel_analysis/` — Agarose-gel image analysis (notebook-driven).
- `SEC/` — Size-exclusion chromatography analysis (notebook-driven).
- `plasmid_sequencing_result/` — Reference sequencing data.
