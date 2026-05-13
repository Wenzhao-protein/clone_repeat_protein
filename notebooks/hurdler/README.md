# notebooks/hurdler

Canonical notebooks for the HURDLER three-site cloning workflow.

| Notebook | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `hurdler_site_combination_analysis.ipynb` | End-to-end demonstration of generating df1, df2, and the lookup table. Mirrors [`../../src/hurdler/pipeline.py`](../../src/hurdler/pipeline.py). | `data/reference_output/*.csv` | `output/hurdler_three_site_combinations_df{1,2}*.csv`, lookup pickle |
| `hurdler_success_rate_analysis.ipynb` | Random-sequence success-rate analysis. Mirrors [`../../src/hurdler/success_rate.py`](../../src/hurdler/success_rate.py). | `output/hurdler_three_site_combinations_df2.csv` | success-rate CSVs + plots |
| `hurdler_success_rate_optimized.ipynb` | Pattern-based, faster variant of the success-rate analysis. | `output/hurdler_fast_match_package.pkl` | success-rate CSVs + plots |

Backups and earlier variants live in [`../../archive/notebooks/`](../../archive/notebooks/).
