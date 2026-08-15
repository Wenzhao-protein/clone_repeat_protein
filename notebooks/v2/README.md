# HURDLER authoritative notebook suite V2

These 14 source notebooks are generated, output-free, Colab-compatible entry
points. Scientific computation lives in `src/hurdler/notebook_backends/`; the
notebooks only explain parameters, call one backend, display its result and
export a credential-free workspace.

Run them in this order for a complete rebuild:

1. `01_reference_database.ipynb`
2. `02_lookup_plasmid.ipynb`
3. `03_module_corpus.ipynb`
4. `04_query_batch.ipynb`, `05_repeat_designer.ipynb`, or `06_exact_dna_designer.ipynb`
5. `07_production_builder.ipynb` for work that must run on Digs
6. `08_success_landscape_analysis.ipynb`, `09_module_result_analysis.ipynb`, and `10_exact_dna_result_analysis.ipynb`
7. `11_reproducibility.ipynb`

The gel, SEC and plasmid-sequencing notebooks (`12`–`14`) are independent
experimental validation workflows.

## Execution modes

- `tutorial`: bundled fixtures only; every notebook is required to Run All.
- `colab_full`: frozen snapshots by default, with explicit open-source refresh.
- `production_bundle`: write Digs task/submission/recovery/finalization files;
  never submit from Colab.
- `analyze`: download or import compact finalized results.

Every run writes a `hurdler_workspace_<run_id>.zip`. Secrets are never accepted
as notebook parameters and are never included in the workspace or a production
bundle. Regenerate the suite with:

```bash
python scripts/generate_notebook_suite_v2.py
```
