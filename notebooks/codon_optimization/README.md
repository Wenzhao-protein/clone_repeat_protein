# notebooks/codon_optimization

Notebooks for codon-optimisation experiments and the genetic-algorithm
benchmark.

| Notebook | Purpose |
|----------|---------|
| `reverse_translate.ipynb` | Reverse-translate protein sequences using a host codon-usage table. |
| `codon_opt_benchmark.ipynb` | Baseline benchmark notebook for the codon optimiser. |
| `mutation_grid_best_per_iteration.ipynb` | Visualises the GA mutation grid and best-per-iteration trajectories from the extended benchmark. |

The genetic-algorithm implementation itself lives at
[`../../codon_opt_benchmark_extended/src/codon_optimizer.py`](../../codon_opt_benchmark_extended/src/codon_optimizer.py).
