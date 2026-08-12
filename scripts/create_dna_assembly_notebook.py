#!/usr/bin/env python3
"""Create the output-free complete-route regulatory-array notebook."""

from pathlib import Path

import nbformat as nbf


REPO = Path(__file__).resolve().parents[1]
OUTPUT = REPO / "notebooks/tasks/08_long_repetitive_dna_assembly.ipynb"


def md(text: str):
    return nbf.v4.new_markdown_cell(text)


def py(text: str, *, parameters: bool = False):
    metadata = {"tags": ["parameters"]} if parameters else {}
    return nbf.v4.new_code_cell(text, metadata=metadata)


def main() -> int:
    book = nbf.v4.new_notebook()
    book.metadata.kernelspec = {
        "display_name": "HURDLER",
        "language": "python",
        "name": "hurdler",
    }
    book.metadata.language_info = {"name": "python", "version": "3.11"}
    book.cells = [
        md("""# Complete-route HURDLER assembly of regulatory-element arrays

Analysis version: `arbitrary-dna-complete-route-v2`.

This notebook asks a stricter question than the immutable legacy
`arbitrary-dna-active-latent-v1` run. A target passes only when it can be built
from a purchasable exact seed through one or more verified HURDLER growth
cycles. One plasmid is fixed for the route, the restriction-enzyme pair may
change between cycles, every intermediate is an exact integer-copy array, and
the final SHA256 must equal the requested DNA.

The legacy 53.67% value tested mainly a final replacement against an assumed
recipient. It is retained only as a QC baseline and is never used as a
reviewer-response number here."""),
        py("""REPO_ROOT = "/home/wendai/projects/hurdler/clone_repeat_protein"
STUDY_ROOT = f"{REPO_ROOT}/studies/hurdler_validation"
SCRATCH_ROOT = "/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation"
ANALYSIS_VERSION = "arbitrary-dna-complete-route-v2"
PRODUCTION_RESULT_DIR = f"{STUDY_ROOT}/step06_repetitive_dna_assembly/tables/{ANALYSIS_VERSION}/production"
TARGET_CATALOG = f"{STUDY_ROOT}/step06_repetitive_dna_assembly/tables/{ANALYSIS_VERSION}/inputs/target_catalog.parquet"
IDT_CREDENTIAL_MODE = "path"       # "path" or "manual"
IDT_CREDENTIAL_PATH = "~/.config/hurdler/idt.env"
IDT_AUTH_METHOD = "password"       # "password" or "access_token"
RUN_LIVE_IDT = False
USE_PRODUCTION_RESULTS = True
RUN_SMOKE_IF_PRODUCTION_MISSING = True
HEADLESS_EXECUTION = False
SMOKE_ELEMENTS_PER_SOURCE = 1
RF00050_ELEMENT_ID = "Rfam:3cc020d8d3025df7"
""", parameters=True),
        md("""## IDT credentials

IDT supplies complexity scores only; its codon-optimization output is never
used. Path mode accepts an owner-only file outside the repository:

```dotenv
# Password-grant form
IDT_CLIENT_ID=your_client_id
IDT_CLIENT_SECRET=your_client_secret
IDT_USERNAME=your_idt_username
IDT_PASSWORD=your_idt_password
```

or the mutually exclusive token form:

```dotenv
IDT_ACCESS_TOKEN=your_current_access_token
```

```bash
mkdir -p ~/.config/hurdler
cp config/idt.env.example ~/.config/hurdler/idt.env
chmod 600 ~/.config/hurdler/idt.env
```

Set `IDT_CREDENTIAL_MODE="manual"` to enter values invisibly with `getpass`.
Manual credentials remain in kernel memory only and are rejected during
Papermill/Digs execution. Secrets are never notebook parameters, outputs,
hashes, manifests, or log fields."""),
        py("""from pathlib import Path
import json, sys
import pandas as pd
from IPython.display import Image, Markdown, display

repo = Path(REPO_ROOT)
sys.path.insert(0, str(repo / "src"))
from hurdler.complete_route import (
    TARGET_COPY_COUNTS, build_element_matrix, plan_complete_route_catalog,
)
from hurdler.dna_assembly import DNA_COMPLETE_ROUTE_VERSION
from hurdler.dna_assembly_visualization import (
    plot_complete_production_report, plot_production_qc,
    plot_synthetic_factorial_landscape, write_production_figure_manifest,
)
from hurdler.idt import IDTComplexityScorer, configure_idt_credentials
from hurdler.io import sha256_file

assert ANALYSIS_VERSION == DNA_COMPLETE_ROUTE_VERSION
step = Path(STUDY_ROOT) / "step06_repetitive_dna_assembly"
production = Path(PRODUCTION_RESULT_DIR)
figure_dir = step / "figures" / ANALYSIS_VERSION
analysis_dir = step / "tables" / ANALYSIS_VERSION / "notebook_analysis"
analysis_dir.mkdir(parents=True, exist_ok=True)
run_context = {
    "version": ANALYSIS_VERSION, "execution_mode": None,
    "input_hashes": {}, "row_counts": {}, "limitations": [],
}
"""),
        py("""if IDT_CREDENTIAL_MODE == "manual" and HEADLESS_EXECUTION:
    raise RuntimeError(
        "Manual IDT credentials are disabled in Papermill/Digs; use path mode"
    )
if RUN_LIVE_IDT:
    credential_status = configure_idt_credentials(
        mode=IDT_CREDENTIAL_MODE,
        path=IDT_CREDENTIAL_PATH,
        auth_method=IDT_AUTH_METHOD,
        headless=HEADLESS_EXECUTION,
        repository_root=repo,
        include_path_in_status=False,
    )
else:
    credential_status = {
        "credential_mode": IDT_CREDENTIAL_MODE,
        "auth_method": IDT_AUTH_METHOD,
        "required_fields_complete": "not_checked_live_idt_disabled",
        "authentication_verified": "not_requested",
    }
credential_status
"""),
        md("""## Public data and deterministic derivation

The main corpus contains every normalized public element currently acquired:

- CRISPRCasdb direct-repeat records;
- the earlier-middle exact member of each selected official Rfam SEED
  alignment, including RF00050;
- Ribocentre Aptamer exact sequences.

Normalization removes alignment gaps/whitespace, uppercases, converts RNA U
to DNA T, rejects symbols outside A/C/G/T, and deduplicates exact sequences
within each source while retaining source mappings. Each unique element is
expanded independently to exact 2-, 4-, 8-, 16-, and 32-copy targets. The
public inventory, exclusions, raw-download URLs, accessions, timestamps, and
SHA256 hashes are reported before any HURDLER calculation."""),
        py("""public_dir = step / "tables" / "public_elements"
inventory = pd.read_parquet(public_dir / "public_element_inventory.parquet")
exclusions = pd.read_parquet(public_dir / "public_element_exclusions.parquet")
source_mappings = pd.read_parquet(
    public_dir / "public_element_source_mappings.parquet"
)
public_manifest = json.loads(
    (public_dir / "public_element_manifest.json").read_text()
)
run_context["input_hashes"].update({
    "public_element_inventory.parquet": sha256_file(
        public_dir / "public_element_inventory.parquet"
    ),
    "target_catalog.parquet": sha256_file(TARGET_CATALOG),
})
source_summary = inventory.groupby("source_database").agg(
    unique_elements=("element_id", "size"),
    minimum_bp=("element_length_bp", "min"),
    median_bp=("element_length_bp", "median"),
    maximum_bp=("element_length_bp", "max"),
).reset_index()
display(source_summary)
display(pd.DataFrame(public_manifest["source_downloads"]))
{
    "unique_public_elements": len(inventory),
    "excluded_source_rows": len(exclusions),
    "retained_source_mappings": len(source_mappings),
    "derived_exact_targets": public_manifest["derived_target_rows"],
}
"""),
        md("""## Complete-route production tables

The production planner searches a copy-number state graph from 1 through 32.
Every edge deletes a whole number of repeat units from the final state to
recover the exact precursor, then verifies active/latent sites, top/bottom
cuts, sticky ends, donor provenance, fixed-plasmid compatibility, and absence
of unintended selected-enzyme cuts. Candidate routes are ranked by total
experimental steps, unique purchases, purchased bp, pair changes, IDT score,
and stable molecular identifiers.

For a donor core shorter than 90 bp, the output is two complementary 5′→3′
primers exposing the required sticky ends and no IDT request is made. Every
longer actual purchase must carry a live response hash and rule-score sum
strictly below 10."""),
        py("""required_files = {
    "targets": "production_target_analysis.parquet",
    "elements": "production_element_matrix.parquet",
    "routes": "production_selected_routes.parquet",
    "transitions": "production_transitions.parquet",
    "steps": "production_steps.parquet",
    "fragments": "production_fragments.parquet",
    "seeds": "production_seeds.parquet",
}
production_ready = all((production / name).is_file() for name in required_files.values())
if USE_PRODUCTION_RESULTS and production_ready:
    loaded = {
        key: pd.read_parquet(production / filename)
        for key, filename in required_files.items()
    }
    targets, elements, routes = loaded["targets"], loaded["elements"], loaded["routes"]
    transitions, steps = loaded["transitions"], loaded["steps"]
    fragments, seeds = loaded["fragments"], loaded["seeds"]
    headline = json.loads((production / "production_headline_summary.json").read_text())
    run_context["execution_mode"] = "complete_route_production"
else:
    if not RUN_SMOKE_IF_PRODUCTION_MISSING:
        raise FileNotFoundError("Complete-route production tables are incomplete")
    catalog = pd.read_parquet(TARGET_CATALOG)
    chosen = []
    for source in ("CRISPRCasdb", "Rfam", "Ribocentre_Aptamer"):
        source_rows = catalog.loc[catalog.source_database.eq(source)]
        chosen.extend(source_rows.element_id.drop_duplicates().head(
            SMOKE_ELEMENTS_PER_SOURCE
        ).tolist())
    chosen.append(RF00050_ELEMENT_ID)
    smoke_catalog = catalog.loc[catalog.element_id.isin(set(chosen))]
    smoke_input = Path(SCRATCH_ROOT) / (
        "step06_repetitive_dna_assembly/runs/run003_complete_route_v2/"
        "notebook_smoke/target_catalog.parquet"
    )
    smoke_input.parent.mkdir(parents=True, exist_ok=True)
    smoke_catalog.to_parquet(smoke_input, index=False)
    smoke_raw = smoke_input.parent / "raw"
    scorer = (
        IDTComplexityScorer(smoke_raw / "idt_audit.jsonl")
        if RUN_LIVE_IDT else None
    )
    loaded = plan_complete_route_catalog(
        smoke_input,
        repo / "data/reference_output",
        smoke_raw,
        artifact_dir=repo / "data/artifacts",
        idt_scorer=scorer,
        require_idt=RUN_LIVE_IDT,
    )
    targets, routes = loaded["targets"], loaded["selected_routes"]
    transitions, steps = loaded["transitions"], loaded["steps"]
    fragments, seeds = loaded["fragments"], loaded["seeds"]
    elements = build_element_matrix(targets)
    headline = {
        "status": "SMOKE_ONLY_NOT_REVIEWER_ELIGIBLE",
        "elements": len(elements), "targets": len(targets),
    }
    run_context["execution_mode"] = (
        "live_idt_smoke" if RUN_LIVE_IDT else "molecular_smoke_no_idt"
    )
run_context["row_counts"] = {
    "targets": len(targets), "elements": len(elements),
    "selected_routes": len(routes), "transitions": len(transitions),
    "steps": len(steps), "fragments": len(fragments),
}
headline
"""),
        py("""if run_context["execution_mode"] == "complete_route_production":
    assert len(elements) == 29_042
    assert len(targets) == 145_210
    assert set(targets.target_copy_count.astype(int)) == set(TARGET_COPY_COUNTS)
assert not targets.duplicated(
    ["source_database", "element_id", "target_copy_count"]
).any()
assert (
    targets.groupby(["source_database", "element_id"]).size() == 5
).all()
assert not (
    targets.fragment_rescued_by_hurdler & ~targets.complete_route_verified
).any()
if not routes.empty:
    assert routes.complete_route_id.is_unique
    assert routes.final_target_exact.all()
    assert routes.target_sequence_sha256.eq(routes.final_sequence_sha256).all()
if not steps.empty:
    assert steps.unintended_cut_count.fillna(0).eq(0).all()
    assert steps.double_strand_source_verified.all()
if RUN_LIVE_IDT and not fragments.empty:
    long_purchases = fragments.loc[
        fragments.product_type.ne("annealed_sticky_end_primer_pair")
    ]
    assert long_purchases.idt_response_sha256.ne("").all()
    assert long_purchases.idt_score.lt(10).all()
"""),
        md("""## Figure 1 — source coverage and five independent outcomes

Panel B keeps 0/1/2/3/4/5 passed target lengths separate for every element;
there is no “any copy number passed” collapse. Panel C reports one exact
percentage at each target copy number without intervals."""),
        md("""## Figure 2 — exact-target scalability landscape

Each heatmap cell shows the complete-route success percentage and successful
count/total for one shared unit-length bin and one exact target copy number.
Empty bins remain visible as `n=0`."""),
        md("""## Figure 3 — mutually exclusive failures and valid rescue

Rescue requires whole-target IDT failure, a fully verified route, and accepted
purchase fragments. API-unclassified records are retained but excluded from
reviewer headline calculations."""),
        md("""## Figure 4 — experimental complexity

This summarizes HURDLER cycles, maximum verified copies, fixed plasmids,
pair-change counts, Site-I/Site-II usage, and purchase-product count plus
length distribution."""),
        py("""regeneration = elements[[
    "source_database", "element_id", "unit_length_bp",
    "maximum_verified_copy_count", "successful_target_count",
    "all_five_complete",
]].copy()
regeneration.to_parquet(
    analysis_dir / "element_length_vs_maximum_copy_regeneration.parquet",
    index=False,
)
regeneration.to_csv(
    analysis_dir / "element_length_vs_maximum_copy_regeneration.csv",
    index=False,
)
figures = plot_complete_production_report(
    targets, elements, routes, transitions, fragments, seeds, figure_dir
)
if run_context["execution_mode"] == "complete_route_production":
    metrics_path = production / "production_run_metrics.parquet"
    benchmark_path = (
        step / "tables" / ANALYSIS_VERSION / "worker_benchmark.json"
    )
    if metrics_path.is_file() and benchmark_path.is_file():
        figures += plot_production_qc(
            pd.read_parquet(metrics_path),
            json.loads(benchmark_path.read_text()),
            figure_dir,
        )
    synthetic_path = production / "production_synthetic_target_analysis.parquet"
    if synthetic_path.is_file():
        figures += plot_synthetic_factorial_landscape(
            pd.read_parquet(synthetic_path), figure_dir
        )
manifest = write_production_figure_manifest(
    figures,
    figure_dir / "figure_manifest.csv",
    input_tables=[
        production / filename for filename in required_files.values()
        if (production / filename).is_file()
    ] or [TARGET_CATALOG],
)
for path in figures:
    if path.suffix == ".png":
        display(Image(filename=str(path)))
manifest
"""),
        md("""## Figure 5 — RF00050 worked example

RF00050 is never replaced silently. Its five exact target outcomes are shown
as calculated. If the four-copy target fails, the failure stage is reported
verbatim and a separately labelled Rfam/Ribocentre positive control is shown
only when production contains one."""),
        py("""rf = targets.loc[targets.element_id.eq(RF00050_ELEMENT_ID)].sort_values(
    "target_copy_count"
)
if len(rf) != 5:
    raise ValueError("RF00050 must have exactly five independent target rows")
display(rf[[
    "target_copy_count", "target_length_bp", "seed_copy_count",
    "complete_route_verified", "hurdler_step_count", "plasmid",
    "whole_target_idt_status", "whole_target_idt_score",
    "fragment_rescued_by_hurdler", "failure_reason",
]])
rf4 = rf.loc[rf.target_copy_count.eq(4)].iloc[0]
if bool(rf4.complete_route_verified):
    display(Markdown(
        "**RF00050 four-copy verdict:** complete route verified; final DNA is exact."
    ))
else:
    display(Markdown(
        f"**RF00050 four-copy verdict:** failed at `{rf4.failure_reason}`."
    ))
    positive = targets.loc[
        targets.source_database.isin(["Rfam", "Ribocentre_Aptamer"])
        & targets.target_copy_count.eq(4)
        & targets.complete_route_verified
    ].head(1)
    if not positive.empty:
        display(Markdown("**Separately labelled positive control:**"))
        display(positive[[
            "source_database", "element_id", "unit_length_bp", "plasmid",
            "hurdler_step_count", "complete_route_id",
        ]])
"""),
        md("""## Output interpretation

- `complete_route_verified`: a path exists from a purchasable exact seed.
- `final_target_exact`: the final sequence SHA equals the requested array.
- `hurdler_step_count`: digest–ligation growth cycles, excluding seed setup.
- `pair_change_count`: RE-pair changes within the fixed plasmid route.
- `whole_target_idt_status`: score-only evidence for purchasing the intact
  target, not a manufacturing guarantee.
- `fragment_rescued_by_hurdler`: intact target failed IDT, but the full exact
  route passed and every actual purchase was accepted.
- `primer_only_unscored`: all functional donor cores were under 90 bp and were
  emitted as complementary sticky-end primer pairs with no complexity call.
- `no_exact_repeat_gain_pair`: an active/latent replacement exists formally,
  but it does not add a complete repeat unit and therefore is not array growth.
"""),
        py("""if run_context["execution_mode"] == "complete_route_production":
    display(Markdown((production / "reviewer_response.md").read_text()))
else:
    display(Markdown(
        "**Smoke run only. No reviewer N/X is emitted until all 29,042 elements "
        "and 145,210 exact targets are finalized with live IDT evidence.**"
    ))
run_context["limitations"] = [
    "IDT complexity score is not a quote or wet-lab guarantee",
    "donor cores below 90 bp are unscored complementary sticky-end primers",
    "legacy final-step-only success is excluded from reviewer conclusions",
    "smoke-mode values are never reviewer-response numbers",
]
run_context["credentials"] = credential_status
run_context
"""),
    ]
    for cell in book.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(book, OUTPUT)
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
