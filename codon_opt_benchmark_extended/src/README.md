# codon_opt_benchmark_extended/src

| File | Purpose |
|------|---------|
| `codon_optimizer.py` | Genetic-algorithm codon optimiser. Takes a task CSV (see `../tasks/`) and writes per-iteration best sequences. |
| `batch_optimization.sh` | Wrapper that runs `codon_optimizer.py` over an entire task CSV. |
| `run_m_series.sh` | Runs the M12 / M13 series benchmarks end-to-end. |

Per-run outputs land in `codon_opt_results_*/` folders that are
gitignored (see [`../../.gitignore`](../../.gitignore)).
