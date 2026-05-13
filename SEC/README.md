# SEC

Size-exclusion chromatography analysis subproject. Reads MATLAB `.mat`
chromatogram exports, processes them into tidy DataFrames, and plots
per-module / per-assembly chromatograms.

## Layout

| Folder | Contents |
|--------|----------|
| [`src/`](src/) | `sec_utils.py` — shared helpers for reading `.mat` files and building the SEC DataFrame. |
| [`input/`](input/) | Raw chromatogram exports (`.mat`, large; gitignored). |
| [`output/`](output/) | Per-figure PDFs/PNGs used in publications and lab reports. |

## Notebooks

[`../notebooks/sec/result_analysis_total.ipynb`](../notebooks/sec/result_analysis_total.ipynb)
is the canonical driver notebook.

## Status

Exploratory subproject. Historical chunk-plot scratch (~hundreds of
PNGs) moved to [`../archive/sec_temp_plots/`](../archive/sec_temp_plots/).
