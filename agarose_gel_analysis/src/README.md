# agarose_gel_analysis/src

Python sources for the agarose-gel GUI and SCN-file readers.

| File | Purpose |
|------|---------|
| `interactive_gui.py` | Base Tkinter GUI for interactive band/lane picking. |
| `interactive_gui_enhanced.py` | Enhanced GUI with extra lane statistics. |
| `interactive_gui_new.py` | Newer GUI iteration. |
| `interactive_demo.py` | Headless CLI demo, useful for batch runs. |
| `scn_reader.py` | Bio-Rad `.scn` raw reader. |
| `scn_metadata_reader.py` | `.scn` metadata-only reader. |
| `scn_detailed_metadata.py` | Verbose metadata dumper, used for debugging. |

Run from this folder (after activating `codon_opt` env):

```bash
python interactive_gui_enhanced.py
```

Backup / experimental GUI variants are preserved under
[`../../archive/agarose_gel_analysis/`](../../archive/agarose_gel_analysis/).
