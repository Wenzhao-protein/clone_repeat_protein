# Handoff

Status: implementation in progress.

- Rule profile: `legacy-optimized-v1`
- Repository revision at start: `a0d6bb0febc3bf45d0f55ea95d877f1b17b17cf2`
- Missing accepted input class: historical agarose `.scn` files
- Long-task policy: shard when safe; otherwise defer any task expected to exceed 2 hours

The current run state, scheduler IDs, completed outputs, and blockers are kept
under each step's `runs/run01_*` directory and summarized in
`step05_reproducibility/tables/execution_status.csv`.
