# envs

Conda environment files for the workflows in this repository.

| File | Purpose | Used by |
|------|---------|---------|
| `codon_opt.yml` | Core analysis env: `pandas`, `numpy`, `biopython`, `matplotlib`, `seaborn`, `tqdm`, `scipy`, plus codon-optimisation deps. | HURDLER pipeline, enzyme-selection notebooks, codon-optimisation notebooks. |
| `visualization.yml` | Visualisation env (PyMOL / Schrödinger channel). | Structure / 3D-visualisation notebooks only. |

## Usage

```bash
conda env create -f envs/codon_opt.yml
conda activate codon_opt
```

Re-create after the YAML changes:

```bash
conda env update -f envs/codon_opt.yml --prune
```
