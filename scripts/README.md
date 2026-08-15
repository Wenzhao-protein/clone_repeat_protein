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

The V2 notebook suite is regenerated with
`python scripts/generate_notebook_suite_v2.py`. Notebook 07 calls the
`hurdler production` interface; its bundles use
`create_missing_v2_tasks.py` and `finalize_v2_workflow.py` for conservative,
non-deleting recovery and finalization. Workflow-specific maintained
finalizers remain authoritative where they provide stronger scientific
validation.
