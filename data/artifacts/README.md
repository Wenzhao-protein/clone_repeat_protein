# Bundled analysis artifacts

`legacy-optimized-v1/` is the complete frozen HURDLER pattern index used by
the interactive construct designer. It is committed intentionally so a fresh
clone can enumerate every supported Site-I/Site-II/Site-III candidate for all
eight maintained plasmids without a scratch mount or network connection.

The artifact contains 1,335,463 distinct pattern keys and 6,549,904 normalized
plasmid-specific candidate rows. `metadata.json` records the rule profile and
upstream source hashes. Runtime code validates the schema before any query.

This directory contains no credentials and no generated construct sequences.
