# src

Maintained Python source code. Currently a single package, `hurdler`.

| Subpackage | Purpose |
|------------|---------|
| [`hurdler/`](hurdler/) | HURDLER cloning toolkit — pipeline, query CLI, validation, success-rate analysis, shared enzyme/plasmid helpers. |

Run any module from the repository root:

```bash
PYTHONPATH=src python -m hurdler.pipeline
PYTHONPATH=src python -m hurdler.query    --site-i-aa NEQ --site-ii-aa IQA --plasmid "pET-28a(+)"
PYTHONPATH=src python -m hurdler.validate
PYTHONPATH=src python -m hurdler.success_rate
```

New importable code should be added as a new subpackage under `src/`,
not at the repository root. Historical scripts are preserved under
[`../archive/scripts/`](../archive/scripts/).
