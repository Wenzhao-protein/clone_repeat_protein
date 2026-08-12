# Environments

`hurdler.yml` is the maintained portable environment definition.
`hurdler-linux-64.lock` is the explicit package lock used for the validation
campaign. `condarc.yml` records the channel policy without machine-specific
paths or credentials.

The maintained environment includes DSSP 4 (`mkdssp`) so residue-level
secondary structure is generated reproducibly rather than inferred from PDB
header records.

`hurdler-pip.lock` pins the pip-only packages. Kaleido is fixed at 0.2.1
because it bundles a headless renderer; Kaleido 1.x requires a system Chrome
binary that is absent from the Digs universal container.

```bash
/net/software/conda/bin/conda env create \
  --prefix /home/wendai/.conda/envs/hurdler \
  --file envs/hurdler.yml
/home/wendai/.conda/envs/hurdler/bin/pip install -e .
```

The older environment files remain historical snapshots for adjacent
subprojects and are not the HURDLER validation environment.
