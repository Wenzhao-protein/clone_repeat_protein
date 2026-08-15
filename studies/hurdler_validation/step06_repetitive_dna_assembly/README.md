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

## Component orderability production audit

`complete-route-purchase-orderability-v1` audits the selected routes from the
complete `run004_complete_route_v2_primer_fix` shard set. Digs smoke job
`18226538` tested the actual 209-bp production gBlock; final job `18226641`
checked 15,535 selected routes, 3,129 source elements, 129,564 component
occurrences and 9,916 unique logical purchases. Every selected route passed.

The 9,916 unique purchases comprise 6,787 sticky-end standard-primer pairs,
3,083 complementary standard-primer seed pairs, 35 complementary
90--124-nt Ultramer seed pairs and 11 gBlocks. All actual sticky-end primers
are 25--64 nt. Each gBlock is 125--3,000 bp and was evaluated through the live
IDT API under `idt-rule-score-sum-lt10-v1`; all 11 passed and the maximum score
was 2.1. The public compact result is
[`data/results/exact_dna_complete_route_purchase_orderability.csv`](../../../data/results/exact_dna_complete_route_purchase_orderability.csv).
The 35 Ultramer-dependent elements account for 175 routes. Excluding that
long-oligo product class gives a stricter conventional-primer-or-gBlock count
of 15,360/15,535 (98.87%).

This denominator is the set of already-found complete molecular routes. It
does not imply that all 145,210 catalog targets have a route, does not submit
an order, and does not replace a future full rerun under the newer
annotation-aware MCS/vector policy.

The authoritative, output-free notebook is
`notebooks/tasks/08_long_repetitive_dna_assembly.ipynb`. Executed notebooks,
HTML, figures, manifests, compact production summaries and exact rerun records
remain under this study directory; large route and IDT audit data use the
matching `/net/scratch` tree.
