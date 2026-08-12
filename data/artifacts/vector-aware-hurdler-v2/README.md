# vector-aware-hurdler-v2 protein index

This committed sparse CSR index contains every protein-level pair allowed by
the frozen `legacy-optimized-v1` geometry and orthogonality rules. It contains
776 Site-I/Site-II enzyme pairs and deliberately has no plasmid mask. The old
`legacy-optimized-v1` artifact remains immutable and contains 512 pairs after
its historical whole-plasmid prefilter.

`hurdler design-query` first queries this artifact and then evaluates each hit
against `data/reference_output/plasmid_reference_v2.json`, including the
selected retained long backbone and MCS cut scheme. Rebuild deterministically
with:

```bash
hurdler lookup protein-build \
  --input-dir output \
  --orthogonality data/reference_output/orthogonality.csv \
  --output-dir data/artifacts/vector-aware-hurdler-v2
```
