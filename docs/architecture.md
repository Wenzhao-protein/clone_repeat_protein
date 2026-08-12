# Architecture

HURDLER is organized around a single versioned sparse index rather than
interchangeable legacy CSV/pickle bundles.

```text
data/reference_output
        |
        v
hurdler reference build ----> reference manifest + source hashes
        |
        v
hurdler lookup build --------> pattern_index.npz
                               enzyme/variant Parquet tables
                               partitioned solution catalog
                                      |
                 +--------------------+--------------------+
                 |                    |                    |
                 v                    v                    v
          query/validate       short + 6--60AA       module optimization
                                      |                    |
                                      v                    v
                               success report        catalog/results/FASTA
```

## Responsibilities

- `src/hurdler/reference.py`: reference parsing, hashes, and provenance.
- `rules.py`, `constants.py`: frozen scientific conventions.
- `index.py`: enzyme pairing and compact/normalized pattern artifacts.
- `matching.py`: shared query semantics for every downstream workflow.
- `short_screen.py`, `rate.py`: exhaustive 1--5AA and sampled 6--60AA runs.
- `periodicity.py`: full-protein sequence-period candidates. It combines
  Fourier power concentrated at repeat harmonics, lagged identity and
  BLOSUM62-positive similarity, and phase-wise residue conservation.
- `secondary_structure.py`: residue-level DSSP/author H-E-C annotations,
  state/transition periodicity, sequence-to-structure alignment, and the joint
  selector. A shorter designed module is accepted only when sequence and
  secondary structure independently support the same boundary.
- `repeatsdb.py`: exhaustive PDB/AlphaFoldDB annotation materialization,
  top-level and mapped-feature loci parsing, one-region-per-protein selection,
  exact earlier-middle annotated-unit slicing, and source-map-preserving
  sequence deduplication. It never infers a replacement natural boundary.
- `structural_repeats.py`: designed-only strict boundary inference from
  Biotite `DsspApp` eight-state DSSP and Foldseek 3Di lag agreement, followed
  by fragment Foldseek/TM-align validation and MAFFT fixed/variable positions.
- `module_experiments.py`: Stage-1 all-plasmid compatibility, complete shared
  length bins, Stage-2 input freezing, validated shard merges, expanded search
  traces, final summary/FASTA, and the count/proportion and maximum-copy plots.
- `module_results.py`: projected one-row-per-middle-module public CSV under
  `data/results/`; it exposes only validated maximum IDT-accepted DNA and
  removes raw responses, traces, credentials, and machine-specific paths.
- `modules.py`: natural/designed full-sequence curation, source-boundary
  provenance, primitive-unit selection, repeat-region/flank coordinates, and
  fixed/variable position tables.
- `optimization.py`: HURDLER-scheme selection and translation-preserving codon
  diversification under restriction-site, GC, repeat, CAI, and hairpin checks.
- `ga_optimization.py`: repeated-RE-site genetic scoring and adaptive capacity
  search. The search bounds at the mathematical fragment capacity, probes with
  a short binary stage, then audits the boundary one module at a time while
  escalating to 100 generations before declaring the next copy count failed.
  Every locally acceptable DNA is screened by IDT before the length advances;
  explicit rejection reasons reweight the implicated GA score components for
  the next attempt at that same copy count. Capacity searches obey a monotonic
  contract: a sequence already proved orderable at 1,800 bp is a valid lower
  bound for the same module at 3,000 bp. The larger-cap search imports that
  exact translated/locked DNA, obtains a fresh IDT score, and continues upward;
  independent stochastic searches therefore cannot report a smaller maximum
  for the larger capacity. Repeated non-selected RE sites and GC-window
  deviations remain soft score terms; only selected HURDLER Site-I/Site-II
  counts are local hard gates before every eligible candidate is sent to IDT.
- `idt.py`: credential-safe live gBlocks scoring and exact-DNA response cache.
  Every finite rule score is summed; only a complete total strictly below 10
  passes. `IsViolated` is retained as diagnostic evidence rather than a gate.
- `schemas.py`, `io.py`, `paths.py`: shared contracts, manifests, and paths.
- `cli.py`: the public command surface.

Notebooks contain parameters, hashes, summaries, and plots only. Reference
generation lives under `notebooks/reference`; HURDLER calls live under
`notebooks/tasks`. Source notebooks are never treated as completed outputs:
Papermill creates an `_executed.ipynb`, HTML, and an execution manifest.

For natural proteins, `unit_sequence` is exactly the official RepeatsDB unit
at sorted index `(unit_count - 1) // 2` inside the selected longest annotated
region. Top-level `content.loci` and mapped `RepeatsDB-*` feature loci are both
first-class annotation schemas. Insertions stay in the region audit but are
not units. DSSP, Foldseek, or AF3 can flag QC problems but cannot modify these
coordinates.

For designed proteins, `unit_sequence` is accepted only after independent
DSSP and Foldseek evidence agrees on the period and boundaries and adjacent
fragments meet the frozen TM-score, LDDT, coverage, DSSP-transition, and 3Di
thresholds. The earlier middle copy from the maximal passing block is selected.
Fixed positions have at least 0.8 conservation across the MAFFT alignment;
all other positions and contiguous variable ranges are reported explicitly.

Natural PDB chains are never case-normalized. RepeatsDB author chain IDs are
mapped through RCSB GraphQL to the exact `label_asym_id`; DSSP-to-full-sequence
identity below 0.75 is rejected as a probable cross-chain mapping error.

## Compute and recovery

Every scalable run has a local task file with one absolute command per line,
a task index, Taskrunner state, dry-run script, and scratch output directory.
Array-level concurrency is preferred. A failed shard is requeued independently;
large Cartesian products are Parquet partitions rather than monolithic CSVs.

Historical notebooks execute in isolated scratch overlays. Relative input
paths are compatibility symlinks, while writable legacy outputs are scratch
copies. Archive sources are not rewritten.

The active output namespace is `expanded-middle-repeatsdb-foldseek-v1`.
`FINAL_REPORT_PERIODIC_V4_MIDDLE.md` and run directories that mention
`periodic_v4`, `periodic_v3`, or a first-repeat selection are intentionally
kept for auditability but are not part of the active data path.
