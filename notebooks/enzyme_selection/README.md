# notebooks/enzyme_selection

Notebooks that curate the enzyme universe used by the HURDLER workflow.

| Notebook | Purpose |
|----------|---------|
| `enzyme_selection_analysis.ipynb` | Applies per-site filters (methylation, NEB quality, plasmid compatibility, Type IIS, overhang length) and exports `selected_site_{i,ii,iii}_enzymes.csv`. |
| `re_plasmid_compatibility.ipynb` | Builds the enzyme × plasmid compatibility heatmap. |
| `re_3mer_analysis.ipynb` | Estimates sequence coverage of each enzyme via 3-mer AA windows. |
| `inspect_site_candidates.ipynb` | Manual inspection of candidate enzymes for each site. |

See [`../../docs/workflows/enzyme_selection.md`](../../docs/workflows/enzyme_selection.md)
for context.
