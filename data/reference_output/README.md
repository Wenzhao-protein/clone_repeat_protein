# data/reference_output

Curated reference CSVs derived from `../reference_input/` by the
notebooks in [`../../notebooks/utils/`](../../notebooks/utils/). These
files back the HURDLER pipeline and the enzyme-selection notebooks.

| File | Produced by | Consumed by |
|------|-------------|-------------|
| `restriction_enzyme.csv`, `available_restriction_enzyme.csv`, `selected_restriction_enzyme.csv` | `get_re_sites.ipynb` | `enzyme_selection_analysis.ipynb`, `src/hurdler/pipeline.py` |
| `restriction_enzyme_seamless_insert.csv` | `get_re_sites.ipynb` | Site I filter |
| `restriction_enzyme_slient_mutation.csv` | `get_re_sites.ipynb` | Site II filter |
| `seamless_insert.csv`, `slient_mutation.csv` | `get_re_sites.ipynb` | `src/hurdler/pipeline.py` |
| `orthogonality.csv` | `get_re_sites.ipynb` | Type-IIS pairing |
| `methylation_check.csv` | `methylation_check.ipynb` | All HURDLER filters |
| `neb_buffer_activity_cleaned.csv` | `neb_buffer_activity_check.ipynb` | All HURDLER filters |
| `plasmid_check.csv`, `plasmid_digest_check.csv` | `plasmid_check.ipynb` | Plasmid-compatibility filter |
| `codon_usage.csv` | `get_codon_usage.ipynb` | `src/hurdler/pipeline.py` (df2 generation) |

These CSVs are regeneratable. If you change a curation notebook,
re-run it to refresh the relevant file here.
