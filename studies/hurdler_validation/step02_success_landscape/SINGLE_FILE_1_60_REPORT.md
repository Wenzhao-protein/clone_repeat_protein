# HURDLER 1–60AA historical success landscape

Status: **passed** (`historical-notebook-success-v1`)

The corrected production run used Digs job `17520910` with 16 CPU cores. It
completed in 14 seconds (1m34s aggregate CPU, 531 MB peak memory). On an
identical set of 160,000 motifs, 16 processes took 0.389538 s versus 1.896098 s
serial, a 4.868× speedup, with exact row equality.

## Corrected success criterion

The script reconstructs the plasmid-specific pattern population used by the
original `hurdler_success_rate_analysis.ipynb` from its tracked input tables.
Every module is scanned as `module + module`, including cross-boundary matches.
The two 3-mer start positions must satisfy `5 <= d < L`, where `L` is the
single effective module length. The historical helper's first-regex-match and
two-distinct-3-mer behavior are also retained.

All 432 historical 7–60AA length/plasmid success counts exactly match
`output/hurdler_success_rate_7_60aa_per_plasmid.csv`: **0 mismatches**. The
superseded refactored lookup produced 282,907 successes over these observations;
the restored historical population produces the committed 315,516 successes.

The remaining decrease from 5AA to 6AA is not caused by a different success
test. A 5AA motif is repeated twice before screening and is therefore evaluated
as a 10AA effective module, with valid distances 5–9. A 6AA random module is
evaluated directly as 6AA, with only distance 5 valid. Across the eight
plasmids, the 5AA-to-6AA decrease is 2.076–4.663 percentage points.

## Authoritative raw results

- `/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run03_single_file_16core/raw/short_motifs_1_5.parquet`
  - 3,368,420 unique ordered motifs.
  - Per-length counts: 20, 400, 8,000, 160,000, and 3,200,000.
  - SHA256: `70a51d0796f8774304b0afd3cbe31a835a31777f27a0c5a5c28733e290887808`.
- `/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run03_single_file_16core/raw/random_modules_6_60.parquet`
  - 440,000 observations: 55 lengths × 8 plasmids × 1,000 modules.
  - 6AA seed: 420006; 7–60AA seed: 42 in the original length→plasmid→test order.
  - SHA256: `54520893dbb4785a882b3b8532f558d5380a7ca2a4d78d4aaee069cd2bfd8766`.

## Figure and notebook

The PNG/PDF now reproduce the plotting cell in
`notebooks/utils/get_re_sites.ipynb`: default Matplotlib font/line/marker
settings, original plasmid/color order, `Sequence Length` and
`Probability (%)` axes, grid and legend, and the exact title
`3-mer Probability vs Sequence Length`, placed on the user-requested near-square
6 × 5 inch canvas. The raw table still contains 1–60AA,
but the figure displays only 1–50AA. The 1–5AA points are replaced by the
exhaustive results. There is one success-rate value per length/plasmid, with no
confidence band, vertical divider, or additional method annotation.

The authoritative source notebook is
`notebooks/tasks/02_success_rate_1_60.ipynb`. It reads only the two raw Parquet
files above and documents why the effective-length transition can produce the
5AA-to-6AA decrease. Digs job `17520969` executed it from a clean Papermill/SIF
kernel and exported HTML in 4.626 seconds with zero error outputs.
The PNG/PDF pair then passed the focused figure-manifest and contact-sheet
validation in Digs job `17521125`.

The superseded 404-file motif and hit shard directories were already removed
after the replacement files passed independent validation. They remain absent.

## Three-copy rescan and exact comparison

Status: **passed** (`historical-notebook-success-v1-three-copy-scan`)

The active rescan changes only the searched string from `module + module` to
`module + module + module`. The single-module distance condition remains
`5 <= d < L`, the original plasmid-specific pattern population and first-match
semantics are unchanged, and all exhaustive and seeded input sequences are
identical to the two-copy baseline. The comparison validates identity with
bidirectional `EXCEPT ALL` queries before calculating any difference.

The 1AA and 2AA results are confirmed exact zeros in both scans: all 20 1AA
motifs and all 400 ordered 2AA motifs fail on every one of the eight plasmids.
For 1AA, the expanded homopolymer cannot pass the historical requirement for
two distinct recognized 3-mers. Exhaustive evaluation also finds no compatible
directed 3-mer pair for any expanded 2AA motif. The third copy does not change
either result.

The fixed test set shows an increase at module lengths:

`6–39, 43, 46, 47, 48AA`.

There is no observed increase at `1–5, 40–42, 44–45, 49–60AA`, and no
length/plasmid result decreases. Lengths 6–16 improve on all eight plasmids.
The largest individual-plasmid increase is 1.6 percentage points at 6AA; the
eight-plasmid mean at 6AA rises from 7.8% to 9.0%. At longer lengths the added
hits become progressively rarer in the fixed 1,000-sequence samples.

The mechanism is a cyclic-window edge effect. A 3-mer that begins at the final
residue of the first module needs two residues from the next module. A second
candidate 3-mer up to `L-1` residues later can then require residues from the
third module. The doubled string truncates that valid cyclic window, whereas
the tripled string completes it without changing the HURDLER distance rule.

The production calculation used Digs job `17521527`, 16 CPU cores, 35 seconds
wall time, 4m49s aggregate CPU and 1.5 GB peak memory. The identical 160,000-row
benchmark ran in 0.910851 s on 16 processes versus 6.039270 s serial, a 6.63×
speedup with exact equality. Digs job `17521591` independently generated the
comparison, and job `17521659` executed the notebook from a clean SIF kernel
with zero error outputs. Job `17521851` passed focused PNG/PDF and contact-sheet
validation. Digs job `17527153` generated the maintained near-square 6 × 5 inch
figure, executed the updated notebook, passed all 8 focused success-landscape
tests and revalidated the PNG/PDF pair. The full suite passed all 89 tests in
Digs job `17526817`.

Authoritative artifacts:

- Three-copy exhaustive 1–5AA Parquet:
  `/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run06_three_copy_16core/raw/short_motifs_1_5.parquet`
  (`3,368,420` rows; SHA256
  `18ff7736c49e08d5905aa558c55a6363653557360cca1fc17e62778c4616f78d`).
- Three-copy random 6–60AA Parquet:
  `/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run06_three_copy_16core/raw/random_modules_6_60.parquet`
  (`440,000` rows; SHA256
  `bd741cb3e05262ded37ed0c1ee39c51be0d6573e304b23807dad0c5ab9639ac7`).
- Exact 480-row length/plasmid comparison:
  `tables/scan_copy_comparison_2x_vs_3x.csv`.
- 60-row length summary:
  `tables/scan_copy_improvement_by_length.csv`.
- Three-copy PNG/PDF:
  `figures/scan_3x/success_rate_1_60_scan_3x.{png,pdf}`.
  The maintained near-square PNG is exactly 3600 × 3000 px at 600 dpi.
- Executed notebook and HTML:
  `step05_reproducibility/notebooks/04_success_rate_1_60_three_copy_executed.ipynb`
  and `step05_reproducibility/html/04_success_rate_1_60_three_copy.html`.
