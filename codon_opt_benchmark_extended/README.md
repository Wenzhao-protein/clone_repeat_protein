# codon_opt_benchmark_extended

Extended genetic-algorithm codon-optimisation benchmark. Iterates many
mutation seeds for each target protein and tracks the best score per
iteration.

## Layout

| Folder | Contents |
|--------|----------|
| [`src/`](src/) | Python optimiser (`codon_optimizer.py`) and SLURM-style batch drivers (`batch_optimization.sh`, `run_m_series.sh`). |
| [`tasks/`](tasks/) | Task CSVs that describe per-experiment input sequences and parameters (`batch_optimization_tasks.csv`, `m12_tasks.csv`, `m13_tasks.csv`). |
| [`results/`](results/) | Committed result figures and aggregated CSVs (`complexity_distribution*.pdf`, `scatter_plot*.pdf`, `well_color_iterations.csv`). Per-run outputs (`codon_opt_results_*/`) are gitignored. |

## Notebooks

Visualisation lives at
[`../notebooks/codon_optimization/mutation_grid_best_per_iteration.ipynb`](../notebooks/codon_optimization/mutation_grid_best_per_iteration.ipynb).

## Status

Exploratory subproject. Historical notes preserved at
[`../archive/docs/CODON_OPT_BENCHMARK_NOTES.md`](../archive/docs/CODON_OPT_BENCHMARK_NOTES.md).
