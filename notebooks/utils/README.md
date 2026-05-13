# notebooks/utils

One-off curation notebooks that turn raw databases under
[`../../data/reference_input/`](../../data/reference_input/) into the
maintained CSVs under
[`../../data/reference_output/`](../../data/reference_output/).

| Notebook | Output(s) |
|----------|-----------|
| `get_codon_usage.ipynb` | `data/reference_output/codon_usage.csv` |
| `get_re_sites.ipynb` | `data/reference_output/restriction_enzyme*.csv` |
| `methylation_check.ipynb` | `data/reference_output/methylation_check.csv` |
| `neb_buffer_activity_check.ipynb` | `data/reference_output/neb_buffer_activity_cleaned.csv` |
| `plasmid_check.ipynb` | `data/reference_output/plasmid_check.csv`, `plasmid_digest_check.csv` |
| `re_pair_fidelity.ipynb` | Pairwise ligation-fidelity tables (Golden Gate). |

Re-run only when the underlying database (REBASE, NEB, …) is updated.
