# HURDLER final handoff: periodic-v4 middle modules

This is the authoritative handoff for analysis version
`periodic_v4_middle_unit`, boundary method
`spectral-secondary-structure-v4-middle-unit`, and module-selection policy
`repeat-region-middle-unit-tie-earlier-v1`. Adaptive DNA optimization uses
`nonselected-re-sites-soft-score-selected-sites-hard-v2`: repeated
non-selected RE sites are a GA score term, while only selected HURDLER-site
counts and frozen GC bounds are local hard gates before IDT. It supersedes
first-repeat and
`periodic_v3` analyses without deleting their audit records.

## Main results

The final catalog contains 249 unique, exact repeat-module AA sequences: 100
natural modules and 149 designed modules. `designed_primary100` contains 100 of
the designed set for the equal-size comparison. Every downstream HURDLER,
GA, IDT, capacity, and scatter calculation uses a real copy from the middle of
the inferred repeat region, never the first copy. When there are two central
copies, the earlier one is selected deterministically.

| cohort | modules | HURDLER usable | usable fraction |
|---|---:|---:|---:|
| natural100 | 100 | 89 | 0.8900 |
| designed_all | 149 | 146 | 0.9799 |
| designed_primary100 | 100 | 98 | 0.9800 |

The selected middle sequence differs from the first repeat copy for 100/100
natural modules and 10/149 designed modules. The catalog retains the complete
source protein, inferred repeat region, every aligned copy, flanks,
fixed/variable positions, residue-level secondary structure, source hashes,
and original author/RepeatsDB boundaries.

Of the 235 HURDLER-compatible modules, 225 have an IDT-zero-violation,
boundary-proven construct at the 1,800 bp capacity (85 natural, 140 designed),
and 230 at 3,000 bp (87 natural, 143 designed). The actual repeat-construct
scatter contains only results with at least two copies: 223 at 1,800 bp and
227 at 3,000 bp. Two and three results, respectively, are single-module-only
and are classified separately because one module is not a repeat construct;
10 and 5 HURDLER-compatible modules have no IDT-orderable construct found.

Every 3,000 bp maximum is greater than or equal to its 1,800 bp maximum.
Thirty-two stochastic cross-cap decreases found after correcting the RE-site
gate were re-run from the exact proven 1,800 bp DNA lower bound; all 32
recovery searches have a fresh IDT score and an explicit maximum-boundary
proof.

The largest IDT-orderable designed construct is `designed_DHR23`: a 20AA
middle module repeated 13 times at 1,800 bp and 15 times at 3,000 bp. The
natural maximum is `natural_1hm0_B_324_341`: an 18AA middle module repeated 14
times at both capacities. These are
empirical IDT-gated maxima under the frozen GA/search rules, not merely
fragment-length arithmetic ceilings.

## Authoritative deliverables

- Complete per-module summary with AA sequence, length, HURDLER result, both
  maximum copy counts, final GA weights, HURDLER scheme, IDT status/rejection
  evidence, exact passed DNA, and source/boundary evidence:
  [`step04_module_optimization/tables/periodic_v4/module_final_summary.parquet`](step04_module_optimization/tables/periodic_v4/module_final_summary.parquet)
  and its
  [CSV mirror](step04_module_optimization/tables/periodic_v4/module_final_summary.csv).
- Exact passed DNA constructs:
  [`step04_module_optimization/tables/periodic_v4/maximum_passed_constructs.fasta`](step04_module_optimization/tables/periodic_v4/maximum_passed_constructs.fasta)
  with corresponding
  [Parquet](step04_module_optimization/tables/periodic_v4/maximum_passed_constructs.parquet).
- Final natural/designed scatter plot, with module length on the x-axis,
  maximum IDT-orderable copy count on the y-axis, separate 1,800/3,000 bp
  panels, and different colors for natural and designed modules:
  [PNG](step04_module_optimization/figures/periodic_v4/module_length_vs_max_orderable_copies.png)
  and [PDF](step04_module_optimization/figures/periodic_v4/module_length_vs_max_orderable_copies.pdf).
  The main scatter contains only IDT-orderable repeat constructs with at least
  two modules. Single-module-only and no-orderable-result outcomes remain
  explicitly classified in
  [`module_copy_orderability_status.csv`](step04_module_optimization/tables/periodic_v4/module_copy_orderability_status.csv)
  and are counted in each panel title; neither is plotted as a fictitious
  zero-module repeat protein.
- Clean-kernel module report:
  [`step05_reproducibility/html/05_repeat_module_benchmark_periodic_v4_middle.html`](step05_reproducibility/html/05_repeat_module_benchmark_periodic_v4_middle.html).
- Catalog:
  [`step03_module_corpus/tables/periodic_v4/module_catalog_periodic_v4_middle.parquet`](step03_module_corpus/tables/periodic_v4/module_catalog_periodic_v4_middle.parquet).
- Cross-cap proof:
  [`step04_module_optimization/tables/periodic_v4/soft_re_site_monotonic_capacity_validation.json`](step04_module_optimization/tables/periodic_v4/soft_re_site_monotonic_capacity_validation.json).
- Complete workflow status and exact rerun commands:
  [`step05_reproducibility/tables/execution_status.csv`](step05_reproducibility/tables/execution_status.csv).
- Figure manifest, visual report, and contact sheet:
  [`step05_reproducibility/figures/figure_manifest.csv`](step05_reproducibility/figures/figure_manifest.csv)
  and
  [`step05_reproducibility/figures/figure_report.html`](step05_reproducibility/figures/figure_report.html).

## IDT and adaptive-search semantics

IDT is used only to score the DNA produced by the repository optimizer. The
software calls the gBlocks complexity screening API and never adopts an
IDT-generated codon-optimized sequence. A candidate advances the repeat count
only when the matched IDT response contains zero violated rules. Named
rejection reasons and thresholds are stored and raise the corresponding GA
weights before retrying the same length. Repeated non-selected
restriction-enzyme sites are an explicit soft GA score component and never
block IDT scoring by themselves. Translation, selected HURDLER codons, selected
site counts, and GC limits remain hard constraints.

For each HURDLER-compatible module and capacity, the search first locates a
short-generation interval by binary search and then increases the construct by
exactly one module. Each count can escalate through 10, 20, 40, 60, 80, and
100 generations. A reported maximum is accepted only at the mathematical
fragment ceiling or when the next copy has an explicit failed local-plus-IDT
evaluation at 100 generations. The corrected soft-RE production and monotonic
recovery made 9,020 real IDT API calls; including the superseded hard-gate audit
runs, the project made 11,270 real calls. Credentials remain only in the
private external environment file and are absent from repository artifacts.

## Other completed analyses

The short-fragment analysis exhaustively evaluated exactly 20, 400, 8,000,
160,000, and 3,200,000 unique 1--5AA motifs: 3,368,420 total, including every
one of the 400 five-AA shards. Each motif was repeated to the shortest sequence
of at least 6AA before the frozen HURDLER matcher was applied. The combined
1--60AA result is
[`step02_success_landscape/tables/success_rate_1_60.csv`](step02_success_landscape/tables/success_rate_1_60.csv);
1--5AA are exhaustive and 6--60AA are the documented Monte Carlo route.

The final figure manifest contains 31 validated outputs. The final code
validation discovered and passed 59 tests, plus syntax, import, CLI, kernel,
and dependency checks.

## Digs lineage

| work | Digs job ID(s) |
|---|---|
| middle-module catalog | 17425257, 17425865, 17426952 |
| 249-shard HURDLER screen and merge | 17432987, 17436290 |
| superseded hard-gate adaptive audit | 17436340, 17436804, 17439446, 17439552, 17439584 |
| corrected soft-RE adaptive production and strict merge | 17442114, 17442605 |
| corrected cross-cap preparation, IDT recovery, strict merge | 17442627, 17442637, 17442743 |
| monotonic combination, promotion, deliverable synchronization | 17442759, 17442761, 17447046 |
| per-module summary and final visual fix | 17442801, 17442860 |
| clean-kernel module notebook | 17442866 |
| figure report | 17442876 |
| final code validation | 17447953 |
| reproducibility report and status notebook | 17451872, 17451900 |

The full resource requests, state files, task indexes, statuses, and exact
requeue commands are machine-readable in
[`step05_reproducibility/tables/job_registry.json`](step05_reproducibility/tables/job_registry.json)
and `execution_status.csv`.

## Reproduction status and accepted limitations

The final execution matrix contains 174 workflows: 159 `passed`, 7
`blocked_missing_input`, 8 `deferred_long`, and 0 `failed`. Six blocked
workflows require absent historical `.scn` gel files; the seventh requires the
missing `codon_opt_results_m1..m15` mutation-grid trajectories. The eight
deferred workflows are unshardable historical monolithic notebooks expected to
exceed the two-hour policy. Their exact recovery/rerun commands are recorded;
no input, sequence, image, or success result was fabricated.

Historical run directories that mention `periodic_v3`, `first-unit`, or
pre-run48 module optimization are retained only for traceability. They are not
valid inputs for the tables and figures listed above.
