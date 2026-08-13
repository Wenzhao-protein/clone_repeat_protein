# Notebooks

Start here: **[`workflows/01_interactive_hurdler_designer.ipynb`](workflows/01_interactive_hurdler_designer.ipynb)** — the annotation-aware v2 Jupyter designer. The [current Colab preview](https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/agent/vector-aware-designer-v2/notebooks/workflows/02_colab_hurdler_designer.ipynb) is an English tutorial with two mutually exclusive protein-input buttons, live joint RE/plasmid route statistics, automatic or bounded secondary-copy search, opt-in Drive ZIP recovery, and 180-second checkpoints. Accepted runs export every insert/intermediate plasmid as annotated GenBank, full translations, and circular/linear DNA Features Viewer maps. The [main-branch Colab entry](https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/main/notebooks/workflows/02_colab_hurdler_designer.ipynb) and local `./scripts/start_hurdler_web.sh` page use the same scientific controller.

For immutable nucleotide targets, use the [exact-DNA Colab preview](https://colab.research.google.com/github/Wenzhao-protein/clone_repeat_protein/blob/agent/vector-aware-designer-v2/notebooks/workflows/03_colab_exact_dna_hurdler_designer.ipynb). It accepts exact DNA/FASTA or a repeat unit, optional spacer and copy count; preserves active+latent and latent+latent RE geometry (not active+active); requires a complete exact seed-to-target route; and defaults to the RF00059 TPP riboswitch four-copy regulatory-element array. **Run all** asks for a temporary `idt.env` immediately, then automatically searches and exports only a route whose actual 125--3,000 bp purchase gBlocks all pass live IDT complexity scoring. Vector cutting is MCS-first with an inclusive 100-bp default restoration cutoff; set the control to 0 for strictly MCS-internal cutting. Any permitted outside cut is restored in the same scored primary fragment around the exact target seed. Its main ZIP contains only `cloning_steps.csv` and `purchase_inserts.fasta`; detailed GenBank and audit evidence is an optional separate ZIP. The live result displays every attempted fragment's total and per-rule IDT score, a deduplicated final gBlock table, and an interactive circular/linear assembly viewer from `step00_plasmid` through the final plasmid. The optional validation ZIP preserves the score trajectory, all step GenBanks, and static PNG/SVG maps.

The canonical notebooks are deliberately thin.

| Directory | Responsibility |
|---|---|
| `reference/` | Reference manifests, sparse lookup QC, direct RepeatsDB acquisition, and strict designed DSSP/Foldseek boundary evidence |
| `tasks/` | Queries, 1--60AA rates, Stage-1 compatibility, Stage-2 adaptive IDT/GA capacity, exact long-DNA active/latent assembly, and run status |
| `workflows/` | Interactive end-user HURDLER construct design and export |

Each canonical notebook has a tagged parameter cell, frozen rule profile,
input hashes, row/filter summaries, limitations, and deterministic PDF/PNG
outputs where applicable. Scientific computation belongs in `src/hurdler`.

The older topic directories (`hurdler`, `enzyme_selection`, `utils`,
`codon_optimization`, `agarose_gel`, and `sec`) are retained for historical
reproduction. Their clean-kernel results and exact failures are recorded in
`studies/hurdler_validation/step05_reproducibility/tables/execution_status.csv`.
They do not override canonical package results.

The interactive workflow uses fast sequence periodicity only, requires a
person to confirm or modify full-input unit coordinates, and does not run DSSP,
Foldseek, or structure prediction. Its scientific computations and explicit v2
schemas live in `src/hurdler/vector_design.py`; notebooks contain parameters,
widgets, confirmation and display logic only. Protein matching is independent
of plasmids. Annotated retained-backbone filtering, four MCS schemes,
restoration, strict feature protection, IDT scoring and Bulk Input export occur
only after a protein RE pair exists.

`tasks/08_long_repetitive_dna_assembly.ipynb` is the separate exact-nucleotide
complete-route workflow (`arbitrary-dna-complete-route-v2`). It documents both
credential formats, rejects interactive manual input in headless runs, and
reads compact Digs-sharded summaries from `step06_repetitive_dna_assembly`.
Functional donor cores shorter than 90 bp are reported as complementary
sticky-end primer pairs and deliberately bypass IDT gBlocks scoring; the 90 bp
boundary is strict. The legacy final-replacement percentage is shown only as
QC and never as a reviewer-response result.

Execute the canonical production notebooks through the pinned universal SIF
and Papermill wrapper via
`studies/hurdler_validation/scripts/execute_sif_notebook.py`. Executed
notebooks, HTML and manifests are separate artifacts; source notebooks remain
output-free.
