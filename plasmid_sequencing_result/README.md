# plasmid_sequencing_result

Reference plasmid-sequencing artifacts. Each subfolder corresponds to a
sequenced construct and contains the raw read set plus annotation
exports (Genbank, annotated HTML, FASTA, statistics).

## Constructs

| Subfolder | Construct |
|-----------|-----------|
| `Na4M13A/` | Na4 + M13A repeat module. |
| `Na4M25A/` | Na4 + M25A repeat module. |
| `Na4M37A/` | Na4 + M37A repeat module. |
| `Na4M49A/` | Na4 + M49A repeat module. |
| `Na4M6Q6_7_MA/`, `Na4_M6Q6_2MA/`, `Na4_M6Q6_4MA/`, `Na4_M6Q6_5MA/` | Na4 + M6/Q6 repeat-module series. |
| `dArmRP_triangle_step6/` | Triangular dArmRP assembly, step 6. |

## File conventions

Within each construct folder you will find:

- `*.fasta` / `*.gbk` — sequence with annotations.
- `*_pLann.gbk` / `*_pLann.html` — pLannotate exports.
- `*.png` — quick-look maps.
- `*_stats.csv` — coverage statistics.
- `*.clc`, `*.fastq`, `*.distribution.tsv`, etc. — raw reads and
  intermediate analysis files (large, may be gitignored).

This folder is **reference data only** — no scripts here. Treat it as
read-only.
