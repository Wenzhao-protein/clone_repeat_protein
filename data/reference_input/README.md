# data/reference_input

Source databases consumed by the curation notebooks under
[`../../notebooks/utils/`](../../notebooks/utils/). These files are
checked in because they are small and define the reference universe for
all maintained workflows.

| Path | Source | Used by |
|------|--------|---------|
| `link_parsrefs.txt`, `link_withref.txt`, `rebase/link_withref.txt`, `rebase_parsref.txt` | REBASE restriction-enzyme database. | `get_re_sites.ipynb` |
| `neb_buffer_activity.csv` | NEB buffer-activity table. | `neb_buffer_activity_check.ipynb` |
| `golden_gate_fidelity/FileS*.{xlsx,csv}` | Published Golden Gate ligation-fidelity supplementary tables (T4 / T7 / HF / LF / DP / FP ligases at 25 °C / 37 °C / 18 h / 1 h). | `re_pair_fidelity.ipynb` |
| `plasmids/*.fa` | Cloning-host plasmid sequences (pCold I, pET-21a(+), pET-28a(+), pGEX-4T-1, pMAL-c5X, pQE-3, pUC18). | `plasmid_check.ipynb`, `src/hurdler/utils.load_plasmid_sequences` |

Refresh only when the underlying database publishes a new version.
