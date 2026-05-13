# scripts

CLI helpers for the maintained workflows. The Python modules under
`../src/hurdler/` are already runnable as modules:

```bash
PYTHONPATH=src python -m hurdler.pipeline
PYTHONPATH=src python -m hurdler.query    --site-i-aa NEQ --site-ii-aa IQA --plasmid "pET-28a(+)"
PYTHONPATH=src python -m hurdler.validate
PYTHONPATH=src python -m hurdler.success_rate
```

Add new thin shell or Python wrappers here only when they cannot live
inside the package (e.g. SLURM submission scripts). Historical helpers
have moved to [`../archive/scripts/`](../archive/scripts/).
