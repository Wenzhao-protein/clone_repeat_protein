# VALIDATION.md — `constructs_and_sequencing_result/`

Automated QC of the deposited HURDLER design and sequencing files.

## 1. Design files (fragments and plasmid constructs)

All 27 GenBank files parse cleanly with Biopython.

| File | Length (bp) | Topology | sha256 (16) |
|---|---|---|---|
| `extended_armrp/1_fragment/hurdler_extended_armrp_primary.gbk` | 1487 | linear | `559e617da8d53b54` |
| `extended_armrp/1_fragment/hurdler_extended_armrp_secondary_M12.gbk` | 1621 | linear | `ca7615df29b55271` |
| `extended_armrp/1_fragment/hurdler_extended_armrp_secondary_M6Q6.gbk` | 1649 | linear | `d7f206a1aeeea625` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M13A.gbk` | 8169 | circular | `3d5ff4b603ee6ea1` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M25A.gbk` | 9681 | circular | `927d73185d9ce5e7` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M37A.gbk` | 11193 | circular | `a371c9f800a0433c` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M49A.gbk` | 12705 | circular | `aeafa9122433b037` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M6Q6_1_MA.gbk` | 8214 | circular | `444d1d795d29d0bc` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M6Q6_2_MA.gbk` | 9771 | circular | `4285be63af2a75bb` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M6Q6_3_MA.gbk` | 11328 | circular | `c60f1c0dbd521c18` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M6Q6_4_MA.gbk` | 12885 | circular | `a77df7a86c353954` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M6Q6_5_MA.gbk` | 14442 | circular | `2c8921ec4a5fda84` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M6Q6_6_MA.gbk` | 15999 | circular | `a32cf555ce9b18af` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4M6Q6_7_MA.gbk` | 17556 | circular | `5e1543af0507c170` |
| `extended_armrp/2_plasmid_construct/hurdler_extended_armrp_plasmid_Na4MA.gbk` | 6657 | circular | `924dd7321f4cfb1e` |
| `single_chain_armrp_assembly/1_fragment/hurdler_single_chain_fragment_1.gbk` | 2210 | linear | `c352f5579d600bf5` |
| `single_chain_armrp_assembly/1_fragment/hurdler_single_chain_fragment_2.gbk` | 941 | linear | `210660a1b2b9b9c8` |
| `single_chain_armrp_assembly/1_fragment/hurdler_single_chain_fragment_3.gbk` | 863 | linear | `4abc529fa6526611` |
| `single_chain_armrp_assembly/1_fragment/hurdler_single_chain_fragment_4.gbk` | 1578 | circular | `4dd993b50ab94908` |
| `single_chain_armrp_assembly/1_fragment/hurdler_single_chain_fragment_5.gbk` | 821 | linear | `f71f32eb1fe763a7` |
| `single_chain_armrp_assembly/1_fragment/hurdler_single_chain_fragment_6.gbk` | 1578 | circular | `27d864c8ef3d7c11` |
| `single_chain_armrp_assembly/2_plasmid_construct/hurdler_extended_armrp_plasmid_triangle_step1.gbk` | 7380 | circular | `4ce0abb871342eed` |
| `single_chain_armrp_assembly/2_plasmid_construct/hurdler_extended_armrp_plasmid_triangle_step2.gbk` | 8277 | circular | `a1924ec1014f2aa0` |
| `single_chain_armrp_assembly/2_plasmid_construct/hurdler_extended_armrp_plasmid_triangle_step3.gbk` | 9048 | circular | `32f55081a66cb3eb` |
| `single_chain_armrp_assembly/2_plasmid_construct/hurdler_extended_armrp_plasmid_triangle_step4.gbk` | 10575 | circular | `ff3dca9f314d0e50` |
| `single_chain_armrp_assembly/2_plasmid_construct/hurdler_extended_armrp_plasmid_triangle_step5.gbk` | 11304 | circular | `1dcf0341dd922cc0` |
| `single_chain_armrp_assembly/2_plasmid_construct/hurdler_extended_armrp_plasmid_triangle_step6.gbk` | 12831 | circular | `e897978874ab16a1` |
| `single_chain_armrp_assembly/1_fragment/designed_model.pse` | — (PyMOL session) | — | `af936e2eec4e473d` |

## 2. Internal arithmetic consistency

The construct series show exactly constant increments, as expected for recursive
directional ligation of a fixed module:

- `Na4MA` 6,657 bp → `Na4M13A` 8,169 → `Na4M25A` 9,681 → `Na4M37A` 11,193 → `Na4M49A` 12,705 bp.
  Increment = **1,512 bp = 12 modules × 126 bp (42 aa) per ArmRP module**.
- `Na4M6Q6_1_MA` 8,214 bp → … → `Na4M6Q6_7_MA` 17,556 bp.
  Increment = **1,557 bp per round**, constant across all seven steps.
- `triangle_step1` 7,380 → `step2` 8,277 → `step3` 9,048 → `step4` 10,575 → `step5` 11,304 → `step6` **12,831 bp**.
  The final value matches the whole-plasmid sequencing file names
  (`dArmRP-triangle-step6-2_12501-12831.gbk`) exactly.

## 3. Design vs. whole-plasmid sequencing

Each designed construct was compared against the corresponding Nanopore consensus
(circular-permutation and reverse-complement aware, full-length comparison).
**Every construct matched its consensus in length exactly (Δ = 0 bp).**

| Construct | Length (bp) | Mismatches | Identity | Note |
|---|---|---|---|---|
| Na4M13A | 8,169 | 0 | 100% | exact |
| Na4M25A | 9,681 | 0 | 100% | exact |
| Na4M37A | 11,193 | 1 | 99.991% | see 3.1 |
| Na4M49A | 12,705 | 2 | 99.984% | see 3.1 |
| Na4_M6Q6_2MA | 9,771 | 0 | 100% | exact |
| Na4_M6Q6_4MA | 12,885 | 0 | 100% | exact |
| Na4_M6Q6_5MA | 14,442 | 0 | 100% | exact |
| Na4_M6Q6_7_MA | 17,556 | 8 | 99.954% | see 3.2 — **silent** |

### 3.1 Na4M37A and Na4M49A

Both carry the same single substitution, G→A, in one copy of the repeat module
(`insert_part` feature; positions 1695 and 3207 respectively). The design context
occurs 3–4× per plasmid and matches `hurdler_extended_armrp_secondary_M12.gbk`; the
sequenced clone differs in exactly one of those copies. In the coding frame this
changes codon GCC→ACC, i.e. **Ala→Thr** in a single module
(context `SGGNEQIQAVIDAGALPALVQLLS`). Na4M49A carries one further unique-context
substitution at position 10,633 (AAC→GAC, Asn→Asp).

These are single-base differences at ≥99.98% identity in a long repeat array
sequenced by Oxford Nanopore, where consensus errors within repeats are a known
failure mode. They should be reported as observed rather than silently corrected,
and the design file should not be edited to match the read.

### 3.2 Na4_M6Q6_7_MA

Eight differences, occurring as two identical clusters of four, spaced **1,557 bp
apart — exactly one ligation round**. In the coding frame both clusters are
**translationally silent** (design and read both give `…NEQILQEALWALSNIAGGP…`).
The read variant is the one present verbatim in
`hurdler_extended_armrp_secondary_M6Q6.gbk` and occurs 7× in the sequenced plasmid
versus 5× in the deposited design file. Interpretation: the deposited design file
records an alternative synonymous codon set for two of the seven modules; the
assembled plasmid is uniformly the secondary-fragment codon variant.
**Protein-level identity is 100%.**

## 4. Tree consistency with the Supplementary Information

The generated tree (`directory_structure.txt`) reproduces the listing at the end of
the SI **file-for-file** after the three private-use U+F026 filename characters were
replaced with ASCII underscores. No non-ASCII repository path remains in this
deposition.

Hidden CLC-native metadata (`.acl`, `.clcinfo`, `.orderlist2`), CLC project files
(`.clc`), navigation `README.md` files, and validation/control files are retained in
the repository but intentionally omitted from the SI tree. The first two comment
lines of `directory_structure.txt` state these exclusions explicitly, and the
repository test suite regenerates the tree to prevent drift.

## 5. Deposition provenance

The integrated source bundle contained 90 files (45,958,161 uncompressed bytes)
and had SHA256
`8c08b86b36a8b8eadde5b306da8f468367038a136b6aeffb4f03d92c8fb57605`.
Repository-only navigation, CLC project, and control files were retained during
integration. No scientific sequence file was renamed except for replacing three
private-use U+F026 filename characters with ASCII underscores; file contents were
not changed.
