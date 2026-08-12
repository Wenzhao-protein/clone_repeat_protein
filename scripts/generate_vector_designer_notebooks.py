#!/usr/bin/env python3
"""Generate output-free Jupyter and Colab entry notebooks for designer v2."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]


INTRO = """# Annotation-aware HURDLER designer v2

This interface first enumerates protein-level Site-I/Site-II RE pairs and only
then evaluates the annotated retained long backbone of seven physical plasmids
(eight selectable profiles) under four MCS cut schemes. The default input is
`N-cap + module × n + C-cap`; full-protein mode preserves every residue and
repeat variant. It designs files only and never submits an order.

Workflow: enter/confirm sequence → query all allowed RE pairs → inspect
`RE pair → plasmid profile → cut scheme → restoration/silencing` → select a
route → optionally optimize → use live IDT complexity scoring or export Bulk
Input files. IDT is used only for scoring, never for codon optimization.

**Colab users:** choose **Runtime → Run all** once.  The implementation cells
are shown as compact forms by default; double-click a form header or choose
**Show code** to inspect the underlying Python.  After initialization, use the
displayed fields, menus, checkboxes, and buttons from top to bottom.
"""


PARAMETERS = """# Papermill parameters
headless_smoke = False
smoke_output_dir = "output/vector_aware_designer_smoke"
default_input_mode = "split"
default_n_cap = "M"
default_repeat_module = "ACDEFGHIKLMNPQRSTVWY"
default_repeat_copies = 3
default_c_cap = "G"
"""


IMPORTS = """import json
import os
from pathlib import Path

import pandas as pd
import ipywidgets as widgets
from IPython.display import display, clear_output

from hurdler.idt import (
    IDTComplexityScorer,
    configure_idt_credentials,
    configure_idt_credentials_from_bytes,
    configure_idt_credentials_from_values,
)
from hurdler.vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    DesignSelection,
    design_construct_v2,
    design_query,
    write_design_outputs_v2,
)
"""


WIDGETS = r'''input_mode = widgets.ToggleButtons(options=[("N-cap + module + C-cap", "split"), ("Complete exact protein", "full")], value=default_input_mode)
n_cap = widgets.Text(value=default_n_cap, description="N-cap AA")
module = widgets.Textarea(value=default_repeat_module, description="Module AA")
copies = widgets.BoundedIntText(value=default_repeat_copies, min=2, max=10000, description="Copies")
c_cap = widgets.Text(value=default_c_cap, description="C-cap AA")
full_protein = widgets.Textarea(description="Full AA/FASTA", layout=widgets.Layout(width="95%", height="130px"))
repeat_start = widgets.IntText(value=1, description="Repeat start")
repeat_end = widgets.IntText(value=0, description="Repeat end")
repeat_period = widgets.IntText(value=0, description="Period")
site_i_allow = widgets.Text(description="Site I allow")
site_ii_allow = widgets.Text(description="Site II allow")
site_iii_allow = widgets.Text(description="Site III allow")
plasmid_allow = widgets.SelectMultiple(options=("pGEX-4T-1", "pMAL-c5X", "pET-21a(+)", "pET-28a(+)", "pET-28a(+)_start_codon", "pCold_I", "pUC18", "pQE-3"), description="Profiles")
allow_left = widgets.Checkbox(False, description="Allow left cutter as fallback")
allow_right = widgets.Checkbox(False, description="Allow right cutter as fallback")
query_button = widgets.Button(description="1. Query RE pairs and vectors", button_style="primary")
query_output = widgets.Output()
route_choice = widgets.Dropdown(description="Route", layout=widgets.Layout(width="95%"))

def _csv_tuple(value):
    return tuple(part.strip() for part in value.split(",") if part.strip())

state = {"query_result": None}

split_inputs = widgets.VBox([widgets.HBox([n_cap, copies, c_cap]), module])
full_inputs = widgets.VBox([
    full_protein,
    widgets.HTML("Confirm the inferred 1-based repeat coordinates before querying."),
    widgets.HBox([repeat_start, repeat_end, repeat_period]),
])

def _sync_input_mode(_=None):
    split_inputs.layout.display = "" if input_mode.value == "split" else "none"
    full_inputs.layout.display = "" if input_mode.value == "full" else "none"

input_mode.observe(_sync_input_mode, names="value")
_sync_input_mode()

def _current_query():
    return CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode=input_mode.value,
        sequence_id="interactive_design",
        n_cap=n_cap.value,
        repeat_module=module.value,
        c_cap=c_cap.value,
        repeat_copies=copies.value,
        full_protein_sequence=full_protein.value,
        repeat_region_start=repeat_start.value or None,
        repeat_region_end=repeat_end.value or None,
        repeat_period=repeat_period.value or None,
        site_i_allowlist=_csv_tuple(site_i_allow.value),
        site_ii_allowlist=_csv_tuple(site_ii_allow.value),
        site_iii_allowlist=_csv_tuple(site_iii_allow.value),
        plasmid_allowlist=tuple(plasmid_allow.value),
        allow_left_cutter_in_hurdler_pair=allow_left.value,
        allow_right_cutter_in_hurdler_pair=allow_right.value,
    )

def _run_query(_=None):
    with query_output:
        clear_output()
        result = design_query(_current_query())
        state["query_result"] = result
        print(result.status, "—", result.message)
        if result.boundary_analysis and result.status == "needs_boundary_confirmation":
            display(pd.DataFrame(result.boundary_analysis["candidates"]).head(20))
            print("Copy a proposed start/end/period above, verify it, then query again.")
        if result.protein_candidates:
            display(pd.DataFrame(result.protein_candidates).head(100))
        if result.vector_routes:
            frame = pd.DataFrame(result.vector_routes)
            display(frame[["rank", "site_i_enzyme", "site_ii_enzyme", "profile_id", "cut_scheme", "left_cutter", "right_cutter", "restoration_length_bp", "cutter_reuse"]].head(200))
            route_choice.options = [
                (f"#{row['rank']} {row['site_i_enzyme']}/{row['site_ii_enzyme']} → {row['profile_id']} {row['cut_scheme']}", index)
                for index, row in enumerate(result.vector_routes)
            ]

query_button.on_click(_run_query)
display(widgets.VBox([
    input_mode,
    split_inputs, full_inputs,
    widgets.HBox([site_i_allow, site_ii_allow, site_iii_allow]), plasmid_allow,
    widgets.HBox([allow_left, allow_right]), query_button, query_output,
    widgets.HTML("<h3>Confirm one HURDLER RE-pair/vector route</h3>"), route_choice,
]))
'''


OPTIMIZER = r'''validation_mode = widgets.ToggleButtons(options=[("No optimization", "none"), ("Live IDT API", "api"), ("IDT Bulk Input files", "batch")], value="none")
credential_mode = widgets.Dropdown(options=[("External mode-600 env file", "path"), ("Manual OAuth fields", "manual"), ("Temporary env upload", "upload")], value="path", description="Credentials")
credential_path = widgets.Text(placeholder="Path outside this repo", description="Env path")
client_id = widgets.Password(description="Client ID")
client_secret = widgets.Password(description="Client secret")
username = widgets.Password(description="Username")
password = widgets.Password(description="Password")
access_token = widgets.Password(description="Access token")
credential_upload = widgets.FileUpload(accept=".env,text/plain", multiple=False, description="Temporary env")
population = widgets.IntSlider(value=16, min=4, max=256, step=4, description="Population")
max_copies = widgets.BoundedIntText(value=20, min=2, max=10000, description="Max copies")
mutation = widgets.FloatSlider(value=0.08, min=0.001, max=0.5, step=0.001, description="Mutation")
crossover = widgets.FloatSlider(value=0.75, min=0, max=1, step=0.01, description="Crossover")
elite = widgets.FloatSlider(value=0.15, min=0.01, max=0.5, step=0.01, description="Elite")
seed = widgets.IntText(value=42, description="Seed")
auto_feedback = widgets.Checkbox(True, description="Auto-adjust weights from IDT rules")
weights = widgets.Textarea(value=json.dumps({
    "selected_re_site_excess": 1e9, "repeated_re_site_excess": 1e4,
    "gc_window_violation": 1e9, "gc_window_soft_violation": 100,
    "repeated_8mer": 5, "repeated_13mer": 100, "repeated_14mer": 250,
    "hairpin_10mer_proxy": 25, "homopolymer_excess": 250,
    "terminal_repeat_proxy": 100, "negative_log_cai": 50,
}, sort_keys=True), description="GA weights", layout=widgets.Layout(width="95%", height="130px"))
output_dir = widgets.Text(value="output/interactive_vector_aware_design", description="Output")
run_button = widgets.Button(description="2. Optimize / export", button_style="success")
download_button = widgets.Button(description="Download design bundle", icon="download")
download_button.layout.display = "none"
run_output = widgets.Output()
credential_controls = widgets.VBox([])
path_credentials = widgets.VBox([credential_path])
manual_credentials = widgets.VBox([
    widgets.HBox([client_id, client_secret]),
    widgets.HBox([username, password]), access_token,
])
upload_credentials = widgets.VBox([credential_upload])

def _sync_credential_mode(_=None):
    path_credentials.layout.display = "" if credential_mode.value == "path" else "none"
    manual_credentials.layout.display = "" if credential_mode.value == "manual" else "none"
    upload_credentials.layout.display = "" if credential_mode.value == "upload" else "none"

def _sync_validation_mode(_=None):
    credential_controls.layout.display = "" if validation_mode.value == "api" else "none"

credential_controls.children = (
    credential_mode, path_credentials, manual_credentials, upload_credentials,
)
credential_mode.observe(_sync_credential_mode, names="value")
validation_mode.observe(_sync_validation_mode, names="value")
_sync_credential_mode()
_sync_validation_mode()

def _configure_credentials():
    if validation_mode.value != "api":
        return None
    if credential_mode.value == "path":
        return configure_idt_credentials(mode="path", path=credential_path.value, include_path_in_status=False)
    if credential_mode.value == "manual":
        entered = {
            "IDT_CLIENT_ID": client_id.value, "IDT_CLIENT_SECRET": client_secret.value,
            "IDT_USERNAME": username.value, "IDT_PASSWORD": password.value,
            "IDT_ACCESS_TOKEN": access_token.value,
        }
        try:
            return configure_idt_credentials_from_values(entered)
        finally:
            entered.clear()
            for widget in (client_id, client_secret, username, password, access_token):
                widget.value = ""
    uploaded = credential_upload.value
    if not uploaded:
        raise ValueError("Choose one temporary env file")
    item = next(iter(uploaded.values())) if isinstance(uploaded, dict) else uploaded[0]
    payload = bytes(item["content"] if isinstance(item, dict) else item.content)
    try:
        return configure_idt_credentials_from_bytes(payload)
    finally:
        payload = b""
        credential_upload.value = ()

def _run_design(_=None):
    with run_output:
        clear_output()
        queried = state.get("query_result")
        if queried is None or not queried.vector_routes:
            raise RuntimeError("Run a successful query and choose a route first")
        route = queried.vector_routes[route_choice.value]
        _configure_credentials()
        scorer = IDTComplexityScorer(Path(output_dir.value) / "idt_audit.jsonl") if validation_mode.value == "api" else None
        request = DesignRequestV2(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            query=_current_query(),
            selection=DesignSelection(route["candidate_id"], route["profile_id"], route["scheme_id"], route["site_iii_options"][0]),
            validation_mode=validation_mode.value,
            max_repeat_copies=max_copies.value if input_mode.value == "split" else None,
            population_size=population.value,
            mutation_rate=mutation.value,
            crossover_rate=crossover.value,
            elite_fraction=elite.value,
            seed=seed.value,
            generation_schedule=(10, 20, 40, 60, 80, 100),
            score_weights=json.loads(weights.value),
            auto_adjust_weights_from_idt=auto_feedback.value,
        )
        result = design_construct_v2(request, idt_scorer=scorer)
        files = write_design_outputs_v2(result, output_dir.value)
        state["design_files"] = files
        print(result.status, "—", result.message)
        display(pd.DataFrame(result.primary_fragments))
        display(pd.DataFrame(result.cloning_steps))
        print(files)
        download_button.layout.display = ""

def _download_design(_=None):
    import shutil
    source = Path(output_dir.value)
    if not source.is_dir():
        raise FileNotFoundError("Run the design step before downloading files")
    archive = Path(shutil.make_archive(str(source), "zip", source))
    if "COLAB_RELEASE_TAG" in os.environ:
        from google.colab import files as colab_files
        colab_files.download(str(archive))
    else:
        with run_output:
            print(f"Design bundle: {archive.resolve()}")

run_button.on_click(_run_design)
download_button.on_click(_download_design)
display(widgets.VBox([
    validation_mode, credential_controls,
    widgets.HBox([population, max_copies, mutation]), widgets.HBox([crossover, elite]), widgets.HBox([seed, auto_feedback]),
    weights, output_dir, run_button, download_button, run_output,
]))
'''


SMOKE = r'''if headless_smoke or os.environ.get("HURDLER_NOTEBOOK_SMOKE") == "1":
    smoke_dir = Path(os.environ.get("HURDLER_NOTEBOOK_SMOKE_OUTPUT", smoke_output_dir))
    smoke_query = CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2, input_mode="split", sequence_id="vector_aware_smoke",
        n_cap="M", repeat_module="ACDEFGHIKLMNPQRSTVWY", repeat_copies=3, c_cap="G",
    )
    smoke_queried = design_query(smoke_query)
    smoke_route = smoke_queried.vector_routes[0]
    smoke_request = DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=smoke_query,
        selection=DesignSelection(smoke_route["candidate_id"], smoke_route["profile_id"], smoke_route["scheme_id"], smoke_route["site_iii_options"][0]),
        validation_mode="batch", max_repeat_copies=3, population_size=4, generation_schedule=(10, 100),
    )
    smoke_result = design_construct_v2(smoke_request)
    smoke_files = write_design_outputs_v2(smoke_result, smoke_dir)
    assert smoke_result.status == "optimized_unvalidated_batch"
    assert not smoke_result.idt_audit
    print({"status": smoke_result.status, "routes": len(smoke_result.vector_routes), "files": smoke_files})
'''


COLAB_BOOTSTRAP = r'''# Colab-only bootstrap: computation and credentials stay in this runtime.
import os, subprocess, sys
from pathlib import Path

repository_ref = os.environ.get("HURDLER_REPOSITORY_REF", "agent/vector-aware-designer-v2")
if "COLAB_RELEASE_TAG" in os.environ:
    if not Path("/content/clone_repeat_protein").exists():
        subprocess.run([
            "git", "clone", "--branch", repository_ref, "--single-branch",
            "https://github.com/Wenzhao-protein/clone_repeat_protein",
            "/content/clone_repeat_protein",
        ], check=True)
    os.chdir("/content/clone_repeat_protein")
    subprocess.run([sys.executable, "-m", "pip", "install", "-e", ".[notebooks,optimization]"], check=True)
'''


def _code_cell(
    source: str,
    *,
    colab: bool,
    cell_id: str,
    title: str,
    tags: list[str] | None = None,
):
    if colab:
        source = f"#@title {title}\n" + source
    cell = nbf.v4.new_code_cell(source)
    cell.metadata["id"] = cell_id
    if tags:
        cell.metadata["tags"] = tags
    if colab:
        cell.metadata["cellView"] = "form"
        cell.metadata["colab"] = {}
    return cell


def notebook(*, colab: bool = False):
    nb = nbf.v4.new_notebook()
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    }
    intro = nbf.v4.new_markdown_cell(INTRO)
    intro.metadata["id"] = "hurdler-introduction"
    cells = [intro]
    if colab:
        cells.append(
            _code_cell(
                COLAB_BOOTSTRAP,
                colab=True,
                cell_id="hurdler-initialize",
                title="0. Initialize HURDLER (run once)",
            )
        )
    cells.extend(
        [
            _code_cell(
                PARAMETERS,
                colab=colab,
                cell_id="hurdler-parameters",
                title="Notebook defaults",
                tags=["parameters"],
            ),
            _code_cell(
                IMPORTS,
                colab=colab,
                cell_id="hurdler-imports",
                title="Load the HURDLER design engine",
            ),
            nbf.v4.new_markdown_cell(
                "## 1. Enter a protein, query RE pairs, and confirm a route"
            ),
            _code_cell(
                WIDGETS,
                colab=colab,
                cell_id="hurdler-query-controls",
                title="Protein and RE-pair/vector controls",
            ),
            nbf.v4.new_markdown_cell(
                "## 2. Optional codon optimization and IDT scoring/export"
            ),
            _code_cell(
                OPTIMIZER,
                colab=colab,
                cell_id="hurdler-optimization-controls",
                title="GA, IDT, and file export controls",
            ),
            _code_cell(
                SMOKE,
                colab=colab,
                cell_id="hurdler-headless-smoke",
                title="Automated validation hook (normally inactive)",
            ),
        ]
    )
    for index, cell in enumerate(cells):
        cell.metadata.setdefault("id", f"hurdler-cell-{index:02d}")
    for cell in cells:
        if cell.cell_type == "code":
            cell.execution_count = None
            cell.outputs = []
    nb.cells = cells
    return nb


def main() -> None:
    nbf.write(notebook(), ROOT / "notebooks" / "workflows" / "01_interactive_hurdler_designer.ipynb")
    nbf.write(notebook(colab=True), ROOT / "notebooks" / "workflows" / "02_colab_hurdler_designer.ipynb")


if __name__ == "__main__":
    main()
