# Citation guide for the 182 designed repeat-protein modules

The active `expanded-middle-repeatsdb-foldseek-v1` result table contains 182
unique designed middle-module sequences. They trace to **19 source records**:
**17 peer-reviewed papers with a DOI and two PDB-only records whose primary
citation is still “to be published.”** The row-level mapping is in
[`designed_repeat_protein_citations.csv`](../data/results/designed_repeat_protein_citations.csv),
and every individual module can be joined through `source_accession` in
[`natural_designed_repeat_protein_hurdler_idt.csv`](../data/results/natural_designed_repeat_protein_hurdler_idt.csv).

## Citations required for the whole designed collection

For a manuscript using all 182 rows, cite all 17 DOI papers listed in the CSV.
The largest blocks are:

- Brunette et al. (2015), which accounts for 82 DHR modules.
- Huddy et al. (2024), which accounts for 40 THR modules and 25 associated
  constructs (65 rows total).
- Nine designed-armadillo/DARPin papers, which collectively account for 25
  rows.
- Three TPR papers (four rows), two ankyrin-containing sources (five rows),
  the designed LRR paper (one row), and the iTHR ice-recrystallization paper
  (one row).

PDB 4HQD (OR265) and 4PQ8 (OR465) do not currently name a published paper in
their primary PDB citation. Cite the PDB records themselves and describe them
as structure records, not peer-reviewed article-derived designs.

## How the counts were derived

The counts are observation counts after strict DSSP/Foldseek boundary
acceptance and exact middle-module sequence deduplication. They are not counts
of papers, full constructs, or every sequence reported by those papers. The
65 rows associated with Huddy et al. are intentionally consolidated under one
paper even though 25 legacy catalog rows had the placeholder citation “See
source_url and source_accession”; their source URL and supplement are from the
same Nature article.
