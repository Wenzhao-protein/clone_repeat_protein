# Natural/Designed repeat-protein results

For the paper-level provenance of all 182 designed modules, see the
[designed-module citation guide](../../docs/designed_repeat_protein_citations.md)
and its machine-readable
[`designed_repeat_protein_citations.csv`](designed_repeat_protein_citations.csv).

[`natural_designed_repeat_protein_hurdler_idt.csv`](natural_designed_repeat_protein_hurdler_idt.csv)
is the public, spreadsheet-oriented export of the active
`expanded-middle-repeatsdb-foldseek-v1` corpus. It deliberately excludes the
249-row `periodic_v4` legacy corpus.

The table contains one row per unique middle-module AA sequence within the
Natural or Designed collection. Open it in pandas, LibreOffice, or Excel and
filter `collection`, `family`, `record_status`, `hurdler_compatible`, or either
`cap*_idt_passed` column. `display_name` is a short accession-oriented label;
`search_terms` combines the useful accession, family, structure, annotation,
module-sequence, plasmid, and enzyme identifiers for a simple text search.

## Result levels

- `hurdler_compatible` is a Stage-1 geometry result. It is true when scanning
  `middle_module × 2` finds at least one legal frozen
  `legacy-optimized-v1` pattern on any of the eight maintained plasmids.
- `cap1800_*` and `cap3000_*` are Stage-2 construct results. A populated
  `*_idt_accepted_dna` is the maximum verified construct for that capacity,
  not an arbitrary intermediate GA candidate.
- Stage 2 requires exact AA translation, zero selected Site-I/Site-II excess
  sites, fragment-length compliance, a live IDT response hash, IDT finite-rule
  score sum `<10`, and either a capacity-limit proof or failure of the next
  copy at 100 GA generations.
- An incompatible or unsuccessful row retains an explicit status/reason and
  leaves accepted-DNA and maximum-verified-copy fields empty. Zero or one is
  never presented as a verified repeat maximum.

The overall Stage-1 percentages should not be treated as the random-sequence
success landscape. This catalog tests curated biological modules, takes the
union across all eight plasmids, and uses the frozen two-copy classifier. The
success-landscape notebook estimates a separate per-length, per-plasmid
probability using exhaustive short motifs or Monte Carlo random sequences; its
current comparison scan uses three copies.

## Important column groups

| Group | Representative columns |
|---|---|
| Search and identity | `display_name`, `search_terms`, `collection`, `module_id`, `family` |
| Protein and module | `full_protein_sequence_aa`, `full_protein_length_aa`, `middle_module_sequence_aa`, `middle_module_length_aa`, repeat-region and middle-unit coordinates |
| Sources | UniProt/structure/annotation accessions, source URLs, citations, evidence tier, source-mapping count |
| Boundary evidence | boundary method/status, DSSP agreement, Foldseek 3Di identity, TM-score, LDDT, coverage |
| HURDLER | compatible plasmids, candidate count, selected plasmid, Site-I/Site-II/Site-III enzymes, recognition sites, positions, direction |
| Capacity result | mathematical and verified maxima, proof/stop reason, GA weights, IDT score/rules/hash, accepted DNA and DNA hash |
| Reproducibility | corpus/rule/IDT versions, input-bundle hash and generation timestamp; the exact repo-relative source paths are in the regeneration command below |

The CSV is UTF-8/RFC-4180 with single-line sequences. External text fields that
could be interpreted as spreadsheet formulas are prefixed safely. It contains
no credential, raw IDT response, private environment path, scratch path, or
complete adaptive-search trace.

## Regeneration

After the complete Stage-2 finalizer has produced `maximum_copy_results.parquet`:

```bash
hurdler export-module-results \
  --catalog studies/hurdler_validation/step03_module_corpus/tables/expanded-middle-repeatsdb-foldseek-v1/module_catalog.parquet \
  --source-mappings studies/hurdler_validation/step03_module_corpus/tables/expanded-middle-repeatsdb-foldseek-v1/natural_module_catalog_source_mappings.parquet studies/hurdler_validation/step03_module_corpus/tables/expanded-middle-repeatsdb-foldseek-v1/designed_module_catalog_source_mappings.parquet \
  --compatibility studies/hurdler_validation/step04_module_optimization/tables/expanded-middle-repeatsdb-foldseek-v1/module_compatibility.parquet \
  --maximum-results studies/hurdler_validation/step04_module_optimization/tables/expanded-middle-repeatsdb-foldseek-v1/maximum_copy_results.parquet
```
