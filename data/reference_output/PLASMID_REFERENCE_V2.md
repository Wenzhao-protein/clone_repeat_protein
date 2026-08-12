# Annotation-aware plasmid reference v2

`plasmid_reference_v2.json` contains seven physical circular sequences, eight
expression profiles, normalized feature annotations, source URLs/download
dates/SHA256 values, and four retained cut schemes per profile.

SnapGene `.dna` inputs are downloaded only during materialization and are not
stored in Git. NCBI GenBank annotations are conservatively unioned when a
matching accession is available. Every physical sequence is validated
base-for-base against `data/reference_input/plasmids/*.fa`.

Coordinates are zero-based, half-open. Cut geometry is calculated in the
profile's expression direction and mapped back to physical plasmid
coordinates. Rebuild and validate with:

```bash
hurdler plasmid-reference build
hurdler plasmid-reference validate
```

Sources include public SnapGene plasmid resources and NCBI accessions
AB186388 (pCold I), U13853 (pGEX-4T-1), and L09136 (pUC18). Source records and
hashes are embedded in the JSON; no proprietary/raw binary is committed.
