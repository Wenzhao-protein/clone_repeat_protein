# Archive

This folder preserves historical, duplicate, exploratory, and one-off files
that were removed from the main working paths during the May 2026 repository
reconstruction. Nothing here is part of the maintained workflows, but every
file is kept for traceability.

## Layout

| Path | What it holds | Why archived |
|------|---------------|--------------|
| `archive/scripts/` | Legacy generation, lookup, debug, demo and monitor scripts (`generate_df2_optimized_v2.py`, `create_hurdler_lookup*.py`, `debug_lookup*.py`, `demo_algorithm_difference.py`, `monitor_lookup*.sh`, …). | Superseded by `src/hurdler/pipeline.py`. |
| `archive/tests/` | Ad-hoc and one-off test scripts (`quick_test_*.py`, `simple_test.py`, `test_new_logic.py`, `test_plot_fix.py`, …). | Not exercised by the maintained workflows; useful only as historical references. |
| `archive/notebooks/` | Backup / standalone / executed notebook variants (`hurdler_standalone*.ipynb`, `*_backup.ipynb`, `inspect_site_candidates_executed.ipynb`, …). | Duplicates of the notebooks now living under `notebooks/`. |
| `archive/agarose_gel_analysis/` | Backups, GUI tests, and diagnostic helpers extracted from `agarose_gel_analysis/input/`. | Preserved alongside the maintained agarose-gel workflow. |
| `archive/get_re_dict/` | The early `get_re_ dict/` directory (note the space in the original name) containing first-iteration enzyme-reference notebooks and outputs. | Replaced by `notebooks/utils/` + `data/reference_*`. |
| `archive/artifacts/` | Large or one-off binary outputs that were committed at the repository root (`newplot.png`, `pae.*.npy`, `pred.*.cif`, an old `hurdler_three_site_combinations_df1.csv`). | Did not belong in the working tree. |
| `archive/docs/` | Historical reports and indices (`FILE_INDEX.md`, `OPTIMIZATION_NOTES.md`, GUI fix reports). | Superseded by `docs/README.md` and `docs/workflows/`. |

## Status

Files in `archive/` are **read-only references**. They may import paths or
filenames that no longer exist. Do not run them; consult the canonical
implementations under `src/hurdler/`, `notebooks/`, and `docs/` instead.
