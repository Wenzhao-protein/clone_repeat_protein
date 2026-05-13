# Workflow: Restriction-Enzyme Selection

This workflow curates the set of restriction enzymes that are eligible
for HURDLER Site I, Site II, and Site III, and explores their pairwise
compatibility across the eight supported plasmids.

It is **upstream** of the HURDLER three-site combination workflow: it
produces the filtered enzyme lists that
[`hurdler_site_combinations.md`](hurdler_site_combinations.md) then
combines.

## Goal

Starting from the full REBASE enzyme catalogue plus NEB activity data,
produce:

- per-site eligible enzyme tables for Site I, Site II, and Site III, and
- per-plasmid compatibility data (no backbone cuts),

so that downstream HURDLER analyses can iterate over a clean, curated
universe rather than the raw enzyme database.

## Inputs

All committed under `data/reference_input/` and curated into
`data/reference_output/`:

| File | Role |
|------|------|
| `reference_output/restriction_enzyme.csv` | Master enzyme table from REBASE. |
| `reference_output/methylation_check.csv` | DH5α methylation compatibility. |
| `reference_output/neb_buffer_activity_cleaned.csv` | NEB buffer / ligation / star-activity. |
| `reference_output/plasmid_digest_check.csv` | Does each enzyme cut each plasmid backbone? |
| `reference_output/orthogonality.csv` | Sticky-end orthogonality between enzyme pairs. |
| `reference_input/plasmids/*.fa` | Plasmid sequences used for cut-site detection. |

## Notebooks

Located in [`../../notebooks/enzyme_selection/`](../../notebooks/enzyme_selection/):

| Notebook | Purpose |
|----------|---------|
| `enzyme_selection_analysis.ipynb` | Applies the site-by-site filters (methylation, NEB quality, plasmid compatibility, Type IIS, overhang length) and exports `selected_site_{i,ii,iii}_enzymes.csv`. |
| `re_plasmid_compatibility.ipynb` | Builds the enzyme × plasmid compatibility heatmap. |
| `re_3mer_analysis.ipynb` | Estimates the protein-sequence coverage of each enzyme as a function of how many 3-mer AA windows it can cut. |
| `inspect_site_candidates.ipynb` | Manual inspection of candidate enzyme pairs for each site. |

## Outputs

Written into the gitignored [`../../output/`](../../output/) folder:

- `selected_site_i_enzymes.csv`
- `selected_site_ii_enzymes.csv`
- `selected_site_iii_enzymes.csv`
- `selected_enzymes_summary.json`
- `re_plasmid_compatibility*.{csv,pdf,png}`
- `re_3mer_*.{csv,svg}`
- Funnel diagrams and Sankey HTML for the selection process.

## Related code

Shared helpers used by these notebooks live in
[`../../src/hurdler/utils.py`](../../src/hurdler/utils.py):
`check_neb_quality`, `build_enzyme_pairing_matrix`,
`group_enzymes_by_overhang`, `load_plasmid_sequences`,
`build_enzyme_plasmid_matrix`, `plot_enzyme_pairing_heatmap`.

## History

The original write-up of the export changes is preserved at
[`../reports/ENZYME_SELECTION_PIPELINE_CHANGES.md`](../reports/ENZYME_SELECTION_PIPELINE_CHANGES.md).
