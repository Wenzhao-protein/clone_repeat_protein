# Contributing

These conventions keep the repository navigable after the May 2026
reconstruction.

## File placement

| You are adding… | Put it under… |
|-----------------|---------------|
| Stable, importable Python logic | `src/hurdler/` (or a new `src/<package>/`). |
| A guided analysis / demo | `notebooks/<topic>/`. |
| A regression or smoke test | `tests/`. |
| A CLI entry point or shell helper | `scripts/`. |
| A small reference CSV / FASTA / database file | `data/reference_input/` or `data/reference_output/`. |
| A document for end users or maintainers | `docs/` (workflow guides in `docs/workflows/`). |
| A historical writeup or one-off report | `docs/reports/`. |
| A deprecated script, backup notebook, or executed copy | `archive/` (with a note in `archive/README.md`). |

Generated outputs (large CSVs, pickles, plots, executed notebooks) belong
in `output/` and **must not** be committed. See [`../.gitignore`](../.gitignore).
The sole maintained exception is the public, projected
`data/results/natural_designed_repeat_protein_hurdler_idt.csv`; its exporter
enforces the GitHub blob limit and excludes raw traces and private paths.

## Naming conventions

- Use **lowercase `snake_case`** for new files and folders.
- Avoid spaces and trailing copy markers (`-1`, `-copy`, `-backup`,
  `_v2`). If you need an alternative version, give it a meaningful
  suffix or branch.
- Filenames should match the workflow they belong to (e.g.
  `hurdler_success_rate_*`). Cross-reference filename, docs, and code so
  searches resolve to one place.

## Documentation conventions

- English is the primary maintenance language for docs and module
  docstrings. Chinese is welcome in companion docs (clearly labelled),
  but the root README, `docs/architecture.md`, `docs/glossary.md`, and
  workflow guides stay in English.
- Each maintained Python module should start with a docstring that
  states its purpose, inputs, outputs, and the workflows that depend on
  it.
- Each retained notebook should start with a markdown cell that
  describes purpose, required inputs, expected outputs, and the script
  or module it mirrors.
- Update `docs/README.md`, `docs/architecture.md`, and the relevant
  workflow guide whenever you add or rename a file in `src/`,
  `notebooks/`, or `data/`.

## Annotations / comments

- Prefer a clear module docstring + per-function docstrings over
  scattered inline comments.
- Avoid comments that refer to filenames or workflows that no longer
  exist; update or remove them.
- For biological context, explain *what* an enzyme/sequence constraint
  means — not only *how* the code expresses it.
