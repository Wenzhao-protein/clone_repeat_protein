# Regulatory-array provenance and validation

Analysis version: `arbitrary-dna-complete-route-v2`.

## Public sequence sources

| Source | Unique elements | Exact source material |
|---|---:|---|
| CRISPRCasdb | 28,539 | Direct-repeat sequences from the official `dr_34.zip` download |
| Rfam | 15 | Earlier-middle exact member of each selected official SEED alignment |
| Ribocentre Aptamer | 488 | Curated aptamer sequence records with source metadata |

The source URLs, accessions, download timestamps and SHA256 hashes are stored
in `tables/public_elements/public_element_manifest.json`. The inventory has
29,052 accepted source rows, 29,042 within-source exact-sequence
representatives, 184 explicit exclusions, and complete original-to-
representative mappings.

Normalization removes whitespace and alignment gaps, uppercases, converts RNA
U to DNA T, rejects non-ACGT symbols, and hashes the normalized exact sequence.
Each representative is expanded without mutation to five independent exact
targets: 2, 4, 8, 16 and 32 copies. These are derived test arrays, not claims
that all copy numbers occur naturally.

## Complete-route validation

For every element, the planner:

1. finds the shortest purchasable exact seed;
2. searches active and one-base latent restriction sites on both strands;
3. requires at least one latent boundary and donor-derived restoration of its
   mismatch between the inward cuts;
4. accepts a graph edge only when removing the donor interval from the final
   state recovers the exact shorter integer-copy state;
5. checks signed overhangs, top/bottom cuts, orientation, all eight maintained
   plasmids, donor digestion, double-strand provenance and unintended cuts;
6. keeps the plasmid fixed across the full path while allowing pair changes;
7. tests 3,000/200/60/55-bp fragmentation ceilings and ranks complete paths by
   the frozen experimental objective;
8. calls IDT only for intact-target evidence and actual candidate purchases of
   at least 90 bp; every accepted long purchase requires a live response hash
   and score sum `<10`;
9. emits shorter donor cores as unscored complementary sticky-end primers;
10. requires the final sequence SHA256 to match `unit × target_copy_count`.

The route audit distinguishes `no_active_latent_pair`,
`no_exact_repeat_gain_pair`, vector/digest failure, purchase/IDT failure and a
verified exact route. API-unclassified rows remain in the tables but are
excluded from reviewer headline calculations.

## RF00050 fixed worked example

The notebook always reports the 125-bp RF00050 member
`Rfam:3cc020d8d3025df7`; it is never silently replaced by a positive control.
Under the earlier final-replacement calculation, a 3→4 interpretation had
been assigned to an EcoRI/SpeI route. Re-analysis with the full state graph
shows that this route does not establish a complete seed-to-target assembly:
RF00050 is currently `vector_or_digest_failure` at 2/4/8/16/32 copies.

Live IDT score-only validation on 2026-08-11 succeeded technically. The exact
target scores were 133.3, 308.3, 1625.4, 1265.1 and 1265.1 for 2/4/8/16/32
copies, respectively; the 125-bp seed scored 0.0. Every response has a stored
hash in the scratch audit. These scores do not change the molecular route
failure and are not a manufacturing guarantee.

If production contains a successful Rfam or Ribocentre case, the notebook
shows it separately as a labelled positive control while preserving the
RF00050 conclusion.
