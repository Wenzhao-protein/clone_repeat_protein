# Glossary

Shared terminology used in source code, notebooks, docs, and filenames.

## HURDLER cloning method

A three-restriction-site strategy for cloning repeat proteins. Three
restriction sites are placed inside or around a coding sequence so that
the repeat unit can be inserted seamlessly while a Type IIS enzyme
provides direction and scarless ligation.

- **Site I — seamless insert.** Regular Type II enzyme whose overhang is
  used to insert the repeat unit without changing the encoded protein.
  Must be methylation-compatible with the cloning host (DH5α) and must
  not cut the chosen plasmid backbone.
- **Site II — silent mutation.** Regular Type II enzyme whose recognition
  site overlaps the protein-coding sequence and can be removed by a
  silent (synonymous) mutation after ligation. Same methylation /
  plasmid-compatibility constraints as Site I, and must share an
  overhang pattern with Site III.
- **Site III — Type IIS.** Enzyme whose cut site is outside its
  recognition sequence (e.g. BsaI, BbsI). Provides directionality and
  scarless joining when paired with Site II.

## Sequence / enzyme terminology

- **3-mer AA.** A 3-amino-acid window of the protein sequence used to
  match an enzyme’s recognition site against the encoded peptide.
- **9-mer bp.** The 9 base pairs corresponding to a 3-mer AA window.
- **Overhang (ovhg).** Length and orientation of the sticky end left by
  a restriction enzyme. Compatibility between two enzymes is determined
  by their overhang length and sticky-end sequence (including the
  reverse complement).
- **Orthogonality.** Two enzymes are orthogonal when their sticky ends
  cannot mis-ligate to each other.
- **Methylation compatibility.** An enzyme is DH5α-compatible if its
  activity is not blocked by 6mA / 5mC modifications introduced by the
  cloning host.
- **Star activity.** Relaxed sequence specificity under suboptimal
  conditions; HURDLER candidates must be free of significant star
  activity.

## Files and dataframes

- **df1** — All valid Site I × Site II × Site III combinations together
  with per-plasmid compatibility booleans.
- **df2** — `df1` expanded with concrete 9-mer DNA sequences and 3-mer
  AA windows for Sites I and II (including silent-mutation variants for
  Site II).
- **lookup** — A dictionary keyed by `(3mer_i, 3mer_ii)` returning the
  list of admissible `(site_i, site_ii, site_iii, plasmid_dict)` tuples,
  used for fast queries.
- **Plasmid compatibility.** A boolean per (enzyme, plasmid) indicating
  that the enzyme does not cut the plasmid backbone.

## Plasmids covered

`pGEX-4T-1`, `pMAL-c5X`, `pET-21a(+)`, `pET-28a(+)`,
`pET-28a(+)_start_codon`, `pCold_I`, `pUC18`, `pQE-3`.
