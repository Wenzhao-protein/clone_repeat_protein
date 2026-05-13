# src/hurdler

Maintained Python package for the HURDLER cloning strategy.

| Module | Purpose | Run as |
|--------|---------|--------|
| `pipeline.py` | Full data-generation pipeline (`df1` → `df2` → lookup → success-rate). | `python -m hurdler.pipeline` |
| `query.py` | CLI to look up valid three-site combinations for a given 3-mer AA pair + plasmid. | `python -m hurdler.query --site-i-aa … --site-ii-aa … --plasmid …` |
| `validate.py` | Sanity-checks the generated data and prints a statistics report. | `python -m hurdler.validate` |
| `success_rate.py` | Success-rate analysis for random sequences of varying lengths. | `python -m hurdler.success_rate` |
| `utils.py` | Shared helpers: NEB-quality, pair compatibility, plasmid loading, heatmaps. | (library only) |

## Environment variables

- `HURDLER_INPUT_DIR` — overrides the default `./data/hurdler_analysis_input`.
- `HURDLER_OUTPUT_DIR` — overrides the default `./output`.

## Documentation

- Workflow guide: [`../../docs/workflows/hurdler_site_combinations.md`](../../docs/workflows/hurdler_site_combinations.md)
- Success-rate workflow: [`../../docs/workflows/hurdler_success_rate.md`](../../docs/workflows/hurdler_success_rate.md)
- Glossary: [`../../docs/glossary.md`](../../docs/glossary.md)
