# Interactive workflows

1. [`01_interactive_hurdler_designer.ipynb`](01_interactive_hurdler_designer.ipynb) — Jupyter interface for split or exact full-protein input, 776 protein-level RE pairs, eight annotated plasmid profiles/four cut schemes, restoration/silencing review, GA controls, live IDT scoring or Bulk Input export.
2. [`02_colab_hurdler_designer.ipynb`](02_colab_hurdler_designer.ipynb) — [open the current branch directly in Colab](https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/agent/vector-aware-designer-v2/notebooks/workflows/02_colab_hurdler_designer.ipynb); it clones and installs this exact review branch before presenting credential-safe form controls with code hidden by default.
3. [`../../apps/hurdler_designer.py`](../../apps/hurdler_designer.py) — run the browser app locally with `./scripts/start_hurdler_web.sh`; computation and credentials remain on the user's machine.

The entries are thin layers over `hurdler.vector_design`. They do not run DSSP,
Foldseek, or structure prediction. The plasmid-independent lookup is committed
at `data/artifacts/vector-aware-hurdler-v2`, so compatibility queries are
offline. Full-protein mode never homogenizes repeat variants. Compatibility is
determined by the frozen two-copy protein geometry followed by annotated
retained-backbone filtering, never by an arbitrary optimized-DNA scan.

For a credential-free end-to-end smoke test:

```bash
HURDLER_NOTEBOOK_SMOKE=1 \
HURDLER_NOTEBOOK_SMOKE_OUTPUT=output/interactive_designer_smoke \
jupyter execute notebooks/workflows/01_interactive_hurdler_designer.ipynb
```

That mode deliberately uses Batch export and makes no HTTP call; it must never
be interpreted as IDT acceptance. Interactive API mode uses the live IDT scorer
only. No entry point implements ordering.
