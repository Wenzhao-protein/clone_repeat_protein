# `hurdler` package

All maintained scientific logic is import-safe and exposed through `hurdler`:

- `reference`: provenance and reference manifests;
- `rules`/`constants`: `legacy-optimized-v1` conventions;
- `index`: sparse pattern index and normalized solution catalog;
- `matching`: unified module queries and candidate materialization;
- `short_screen`/`rate`: exhaustive 1--5AA and sampled 6--60AA analysis;
- `periodicity`: Fourier/self-similarity inference of primitive units and
  complete repeat-region boundaries from full protein sequences;
- `secondary_structure`: DSSP/author residue annotations, H/E/C transition
  periodicity, sequence-to-structure chain validation, and joint boundary
  selection;
- `modules`: RepeatsDB/designed curation, exact RCSB author-to-label chain
  mapping, auditable source-versus-primitive boundaries, and fixed/variable
  residue ranges;
- `repeatsdb`: direct parsing of both official RepeatsDB annotation schemas,
  one longest region per biological protein, exact annotated earlier-middle
  unit selection, and source-map-preserving sequence deduplication;
- `structural_repeats`: strict designed-only boundary evidence from Biotite
  eight-state DSSP, Foldseek 3Di lag agreement, fragment TM/LDDT validation,
  and MAFFT fixed/variable residue calls;
- `module_experiments`: all-plasmid Stage-1 compatibility, Stage-2 selected-pair
  input freezing, complete adaptive trace validation, final tables/FASTA, and
  the requested compatibility and maximum-copy figures;
- `optimization`: HURDLER-aware synonymous construct optimization;
- `design`: strict interactive `DesignRequest`/`DesignResult`, sequence-only
  boundary confirmation, frozen-index route enumeration, fragment planning,
  and final-plasmid simulation;
- `notebook_ui`: credential-clearing widgets and explicit mock headless smoke;
- `ga_optimization`/`idt`: genetic refinement with an explicit repeated-RE-site
  score, mathematical-bound binary-to-linear maximum-copy search, and
  credential-safe live IDT SciTools score-sum-below-10 gating with
  positive-score rule feedback into the corresponding GA weights;
- `dna_assembly`: immutable-DNA active/one-base-latent RE scanning on both
  strands, explicit cut/overhang geometry, exact linear route simulation,
  IDT-feedback breakpoint retries, purchase-fragment deduplication, and final
  target SHA256 validation;
- `complete_route`: purchasable-seed copy-number state graphs, fixed-plasmid
  full-route validation, exact whole-unit gains, candidate-only live IDT
  scoring, five-target element matrices, and strict production finalization;
- `dna_assembly_visualization`: production-first public-element, scalability,
  failure/rescue, route-complexity and fixed RF00050 figures;
- `schemas`/`io`/`paths`: versioned contracts and run metadata;
- `cli`: public command-line interface.

`pipeline.py`, `query.py`, `success_rate.py`, and `validate.py` are compatibility
wrappers. New code should call the library API or the `hurdler` executable.

The package has no import-time computation or filesystem mutation.
