# Complete-route arbitrary-DNA HURDLER assembly

Current analysis: `arbitrary-dna-complete-route-v2`.

The immutable `arbitrary-dna-active-latent-v1` production run is retained as a
final-replacement baseline only. Its 53.67% success rate does not prove a path
from a purchasable seed and is not eligible for the reviewer response.

Version 2 represents exact repeat copy numbers as graph states. A target passes
only when a route starts from the shortest accepted exact seed, contains at
least one HURDLER growth cycle, uses one fixed maintained plasmid, preserves
integer-copy intermediates, has no unintended selected-enzyme cuts, and ends
with a DNA SHA256 identical to the target. Enzyme pairs may change between
cycles. Candidate routes are ranked by experimental steps, unique purchases,
purchased bp, pair changes, IDT score, and stable molecular identity.

Purchase policy:

- donor core `<90 bp`: two complementary 5′→3′ primers exposing the required
  sticky ends; no IDT complexity call;
- longer purchase DNA: a live IDT rule-score sum strictly below 10 and a stored
  response hash are required;
- IDT is used for scoring only and never supplies an optimized sequence.

The main corpus contains 29,042 public elements and five exact targets per
element (2/4/8/16/32 copies; 145,210 targets). Production is sharded by element
on Digs with 16-way array concurrency and no nested multiprocessing.

The authoritative, output-free notebook is
`notebooks/tasks/08_long_repetitive_dna_assembly.ipynb`. Executed notebooks,
HTML, figures, manifests, compact production summaries and exact rerun records
remain under this study directory; large route and IDT audit data use the
matching `/net/scratch` tree.
