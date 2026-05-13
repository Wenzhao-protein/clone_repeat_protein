# agarose_gel_analysis

Agarose-gel image quantification subproject. Reads `.scn` (Bio-Rad
ChemiDoc) and `.png` gel images, allows interactive band picking, and
writes band-area and lane-statistics tables.

## Layout

| Folder | Contents |
|--------|----------|
| [`src/`](src/) | Python sources: interactive GUIs (`interactive_gui*.py`), CLI demo (`interactive_demo.py`), and `.scn` readers (`scn_*reader*.py`). |
| [`data/`](data/) | Input gel images and SCN files, organised by experiment date. |
| [`output/`](output/) | Generated band/lane statistics (`band_analysis_results.csv`, `lane_statistics.csv`). |

## Notebooks

The driver notebooks live at
[`../notebooks/agarose_gel/`](../notebooks/agarose_gel/). The GUI source
under `src/` is launched from those notebooks.

## Status

Exploratory subproject — kept in active use but not part of the
HURDLER maintained path. Backups and historical GUI variants live in
[`../archive/agarose_gel_analysis/`](../archive/agarose_gel_analysis/).
