# HURDLER validation and extension study

This multi-stage study rebuilds the versioned HURDLER lookup, reproduces the
legacy-optimized success-rate analysis, exhaustively screens every 1--5 amino
acid motif, curates natural and designed repeat units, and benchmarks
HURDLER-aware codon-diversified constructs.

The scientific rule set is frozen as `legacy-optimized-v1`. Reviewable code,
task files, summary tables, notebooks, and manifests stay in this directory.
Large indexes, exhaustive Parquet shards, downloaded structures, and raw run
outputs use the exact `/net/scratch` mirror recorded in `study.yaml`.

Raw module, structure, Foldseek, GA, and IDT audit shards are written to the
shared `/net/scratch` mirror. Compact finalized tables, FASTA files, figures,
executed notebooks, HTML reports, and run manifests are promoted to `/home`.

Completed active Stage-1/Stage-2 production shards are removed only after the
finalized-table integrity gate passes. The authoritative cleanup record is
[`DIGS_CLEANUP_MANIFEST.json`](DIGS_CLEANUP_MANIFEST.json): it preserves result
hashes, row counts, IDT audit hashes, Slurm job IDs, deleted sizes, and exact
regeneration commands. Stage-2 search traces are normalized in Parquet rather
than duplicated inside every maximum-result row; raw IDT responses are retained
as a validated `jsonl.gz` audit artifact.

Stages:

1. `step01_reference_lookup`: hash reference inputs and build the sparse index.
2. `step02_success_landscape`: exhaustive 1--5AA plus Monte Carlo 6--60AA.
3. `step03_module_corpus`: every retrievable RepeatsDB PDB/AlphaFoldDB
   annotation plus all recoverable public designed full sequences. Natural
   boundaries come directly from official RepeatsDB loci and can never be
   changed by inference. Designed boundaries require concordant Biotite/mkdssp
   eight-state DSSP and Foldseek 3Di/fragment structural evidence. Both routes
   select the earlier middle copy and retain all source mappings, exclusions,
   aligned copies, and fixed/variable positions.
4. `step04_module_optimization`: HURDLER queries, 1,800/3,000 bp constructs,
   repeated-RE-site GA refinement and IDT-gated adaptive maximum-copy search.
   Every locally acceptable candidate is scored through the IDT API; rejection
   reasons reweight the corresponding GA terms before retrying the same length.
   Repeated non-selected RE sites are a soft fitness term; selected HURDLER
   sites remain hard; GC is a high-weight fitness and live-IDT rule. This
   distinction is frozen as
   `nonselected-re-sites-soft-score-selected-sites-hard-v2`.
   The search uses a 10-generation binary stage followed by one-copy increments
   and generation escalation through 100; it emits the maximum explicitly
   orderable DNA for both natural and designed collections.
   Stage 1 screens every unique accepted module against all eight plasmids and
   stores all solutions in a normalized Parquet table. Stage 2 groups four
   independent modules per one-CPU recoverable task and adaptively increases
   Digs array concurrency only after observing zero live-IDT API failures.
   Serial-versus-parallel scientific equivalence and speed measurements are
   retained in `run102_stage2_parallel_benchmark`.
5. `step05_reproducibility`: notebooks, figures, archive matrix, and final QC.
6. `step06_repetitive_dna_assembly`: exact-DNA active/latent RE planning for
   regulatory arrays and other repetitive nucleic-acid targets. Public elements
   and the 900-case synthetic factorial remain separate, every purchase
   fragment is scored (not optimized) by IDT, and the exact final target hash is
   a hard route invariant.

## Final handoff

The active lineage is `expanded-middle-repeatsdb-foldseek-v1`:

1. runs 93/95/97 materialize every RepeatsDB annotation and select one direct
   earlier-middle unit per biological protein;
2. runs 91/92/96/99 acquire designed structures (including missing-only AF3),
   enforce strict DSSP/Foldseek evidence, and preserve exclusions;
3. run100 merges and sequence-deduplicates the two collections;
4. run101 screens all accepted modules against all eight HURDLER plasmids;
5. runs 102/103 benchmark concurrency and run the live-IDT adaptive copy search;
6. the three versioned authoritative notebooks regenerate the reference,
   Stage-1, and Stage-2 reports from finalized artifacts.

`FINAL_REPORT_PERIODIC_V4_MIDDLE.md`, `periodic_v4_middle_unit`, `periodic_v3`,
first-unit, and pre-run101 files remain immutable audit history and are never
inputs to the active tables.
