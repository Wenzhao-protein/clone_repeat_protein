# Handoff

Status: implementation in progress.

- Rule profile: `legacy-optimized-v1`
- Repository revision at start: `a0d6bb0febc3bf45d0f55ea95d877f1b17b17cf2`
- Missing accepted input class: historical agarose `.scn` files
- Long-task policy: shard when safe; otherwise defer any task expected to exceed 2 hours

The current run state, scheduler IDs, completed outputs, and blockers are kept
under each step's `runs/run01_*` directory and summarized in
`step05_reproducibility/tables/execution_status.csv`.

## 2026-08-13 exact-DNA component orderability

- Analysis: `complete-route-purchase-orderability-v1`.
- Input: all 512 completed shards from
  `run004_complete_route_v2_primer_fix`; 15,535 selected complete routes from
  3,129 unique public elements.
- Digs smoke `18226538`: completed; actual 209-bp production gBlock, live-IDT
  score 0.0.
- Digs production `18226641`: completed in 1m24s, 1 CPU, peak 1.5 GB.
- Result: 15,535/15,535 found routes have every component orderable under the
  explicit oligo-pair/gBlock policy. Eleven unique gBlocks passed live IDT
  (`max score=2.1`); every actual donor primer is 25--64 nt.
- Public compact artifacts:
  `data/results/exact_dna_complete_route_purchase_orderability.{csv,json}`.
- Full route, element and unique-purchase CSV/Parquet tables plus the 11-record
  credential-free IDT audit remain in the exact scratch mirror under
  `run005_purchase_orderability_v1/production`.
