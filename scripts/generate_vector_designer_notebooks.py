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


COLAB_INTRO = """# HURDLER repeat-protein designer — interactive tutorial

This notebook designs a scarless, iterative cloning route for a repeat protein
and exports an auditable molecular record for every step. It **does not place an
order**. Run the notebook once with **Runtime → Run all**, then work through the
numbered panels from top to bottom.

Here **RDL** means the repeat-directional-ligation cycle: one optimized
secondary donor can be digested and inserted repeatedly at the growing array
boundary while Site II is re-silenced after ligation.

## What the three restriction-enzyme roles mean

- **Site I** is active in the coding sequence and defines one HURDLER boundary.
- **Site II** supplies the compatible second boundary but is silenced in the
  assembled sequence by a synonymous codon. The encoded protein is unchanged.
- **Site III** releases the reusable secondary fragment from disposable
  adapters. Those adapters are removed during digestion and are not present in
  the expressed construct.
- **Vector cutters** open an annotation-aware plasmid cut scheme. A route may
  restore bases removed around the MCS; `restore length` is the combined left
  and right restoration sequence, excluding temporary Type-IIS adapters.

## Tutorial workflow

1. Configure temporary or Google Drive storage and choose one credential-safe IDT mode.
   Computation always runs under `/content`; Drive receives only checkpoint/final ZIP files.
2. Choose exactly one protein input mode. The supplied split example is visible
   by default; the complete-protein/FASTA box stays hidden until selected.
3. Select individual RE enzymes and plasmids, then run the molecular query. The live cards count only routes
   jointly supported by the current RE, plasmid and restore-length filters.
4. Explicitly choose Site I/II, Site III,
   plasmid and cut scheme. Changing an upstream field invalidates confirmation.
5. Choose automatic secondary exploration (from one repeat to the physical
   limit) or a bounded repeat-copy range. The GA preserves translation and uses
   repeated RE sites, GC, repeated k-mers, hairpin proxies and codon usage in
   its score. In Live API mode every completed candidate is sent to IDT for
   **complexity scoring only**; HURDLER never uses an IDT optimization result.
   A live trajectory plots every IDT fragment score and all returned rule components.
   Either run GA inside Colab or export a reproducible Local + Slurm bundle.
   The external bundle freezes the complete request, exact Git commit, Conda
   YAML, 16-CPU/32-GB/24-hour defaults, checkpoint commands and an external
   mode-600 IDT credential path; it never copies credential values.
6. Inspect the independent plasmid/insert viewer after the GA panel. Route
   confirmation preloads the annotated circular step00 plasmid; completed GA
   runs add every insert and intermediate plasmid.
7. Review the final status and download the timestamped results ZIP.
7. Download the UTC-stamped result ZIP. It contains purchase FASTA/CSV, IDT and GA
   audits, `step00_plasmid.gb`, every `stepXX_insert.gb` and
   `stepXX_plasmid.gb`, translations, manifests, and static plasmid maps.

## Recovery and interpretation

The longest live-IDT-accepted secondary is checkpointed immediately and the
checkpoint is refreshed every 180 seconds. IDT score sum `<10` is this
notebook's complexity-screen criterion, not a quotation or wet-lab guarantee.
The interactive viewer can switch step/molecule, draw plasmids circularly or
linearly, focus on the cloning region, and show bases/codon translation when
zoomed. See the [DNA Features Viewer documentation](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/index.html).

Credentials used inside Colab are kept only in runtime memory/environment, cleared after use,
and never copied to Drive, GenBank, manifests, notebook output, or ZIP files.
The external runner independently reads its configured repo-external env file
on the target machine and performs a live 125-bp IDT API preflight before GA.
"""


PARAMETERS = """# Papermill parameters
headless_smoke = False
smoke_output_dir = "output/vector_aware_designer_smoke"
default_input_mode = "split"
default_n_cap = "MGSHHHHHHSSGIEGRSSGYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGGGGSGGGGSLEVLFQGPDLPKLVKLLKSSNEEILLKALRALAEIASGG"
default_repeat_module = "NEQIQAVIDAGALPALVQLLSSPNEQILQEALWALSNIASGG"
default_repeat_copies = 25
default_c_cap = "NEQIQAVIDAGALPALVQLLSSPNEQILQEALWALSNIASGGNEQKQAVKEAGALEKLEQLQSHENEKIQKEAQEALEKLQSHGGGLEVLFQGPSSGEFGGGGSMVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYKGGGSSGHHHHHH"
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
minimum_secondary = widgets.BoundedIntText(value=12, min=1, max=1000, description="Min secondary")
feedback_rounds = widgets.BoundedIntText(value=100, min=1, max=1000, description="IDT rounds")
generations_per_round = widgets.BoundedIntText(value=10, min=1, max=1000, description="GA/round")
elite_seed_count = widgets.BoundedIntText(value=10, min=1, max=256, description="Top seeds")
mutation = widgets.FloatSlider(value=0.08, min=0.001, max=0.5, step=0.001, description="Mutation")
crossover = widgets.FloatSlider(value=0.75, min=0, max=1, step=0.01, description="Crossover")
elite = widgets.FloatSlider(value=0.15, min=0.01, max=0.5, step=0.01, description="Elite")
seed = widgets.IntText(value=42, description="Seed")
auto_feedback = widgets.Checkbox(True, description="Auto-adjust weights from IDT rules")
auto_parameter_feedback = widgets.Checkbox(True, description="Auto-adjust GA parameters")
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
            assembly_strategy="exact_reused_secondary_rdl" if input_mode.value == "split" else "single_exact",
            population_size=population.value,
            mutation_rate=mutation.value,
            crossover_rate=crossover.value,
            elite_fraction=elite.value,
            seed=seed.value,
            generation_schedule=(10, 20, 40, 60, 80, 100),
            score_weights=json.loads(weights.value),
            auto_adjust_weights_from_idt=auto_feedback.value,
            minimum_secondary_copies=minimum_secondary.value,
            max_idt_feedback_rounds=feedback_rounds.value,
            generations_per_feedback_round=generations_per_round.value,
            elite_seed_count=elite_seed_count.value,
            auto_adjust_ga_parameters_from_idt=auto_parameter_feedback.value,
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
    widgets.HBox([minimum_secondary, feedback_rounds]),
    widgets.HBox([population, mutation]), widgets.HBox([crossover, elite]),
    widgets.HBox([generations_per_round, elite_seed_count]),
    widgets.HBox([seed, auto_feedback]), auto_parameter_feedback,
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
        validation_mode="batch", assembly_strategy="exact_reused_secondary_rdl",
        population_size=4, generation_schedule=(10, 100),
    )
    smoke_result = design_construct_v2(smoke_request)
    smoke_files = write_design_outputs_v2(smoke_result, smoke_dir)
    assert smoke_result.status == "optimized_unvalidated_batch"
    assert not smoke_result.idt_audit
    print({"status": smoke_result.status, "routes": len(smoke_result.vector_routes), "files": smoke_files})
'''


COLAB_BOOTSTRAP = r'''#@markdown Select the Git branch used by this notebook preview.
repository_ref = "agent/vector-aware-designer-v2" #@param ["agent/vector-aware-designer-v2", "main"] {allow-input: true}

# Colab-only bootstrap: computation and credentials stay in this runtime.
import importlib
import importlib.util
import os
import subprocess
import sys
from pathlib import Path


def _run_hidden(command):
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode:
        tail = (completed.stdout or "")[-4000:]
        raise RuntimeError(f"Initialization command failed ({completed.returncode}): {tail}")
    return completed


repository_dir = Path("/content/clone_repeat_protein")
try:
    import google.colab  # noqa: F401
except ImportError:
    running_in_colab = False
else:
    running_in_colab = True

if running_in_colab:
    if (repository_dir / ".git").is_dir():
        _run_hidden(
            ["git", "-C", str(repository_dir), "fetch", "--depth=1", "origin", repository_ref],
        )
        _run_hidden(
            ["git", "-C", str(repository_dir), "checkout", "--detach", "FETCH_HEAD"],
        )
    else:
        _run_hidden([
            "git", "clone", "--branch", repository_ref, "--single-branch",
            "https://github.com/Wenzhao-protein/clone_repeat_protein",
            str(repository_dir),
        ])
    os.chdir(repository_dir)
    _run_hidden([sys.executable, "-m", "pip", "install", "-e", ".[notebooks,optimization]"])
    source_dir = str(repository_dir / "src")
    if source_dir not in sys.path:
        sys.path.insert(0, source_dir)
    for module_name in tuple(sys.modules):
        if module_name == "hurdler" or module_name.startswith("hurdler."):
            del sys.modules[module_name]
    importlib.invalidate_caches()
elif importlib.util.find_spec("hurdler") is None:
    raise RuntimeError(
        "HURDLER is not installed. In Colab, run this initialization cell; "
        "locally, install the project with python -m pip install -e ."
    )

hurdler_package = importlib.import_module("hurdler")
hurdler_initialization_message = f"HURDLER ready: {Path(hurdler_package.__file__).resolve()}"
'''


COLAB_STORAGE_PANEL = r'''display(storage_panel)'''


COLAB_PROTEIN_FORM = r'''display(protein_input_panel)'''


COLAB_SELECTOR_POLICY_FORM = r'''display(cutter_policy_panel)'''


COLAB_IMPORTS = r'''import asyncio
import hashlib
import html
import json
import os
import queue
import shutil
import subprocess
import tempfile
import threading
import time
import traceback
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import ipywidgets as widgets
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
from IPython import get_ipython
from IPython.display import Javascript, Markdown, clear_output, display
import matplotlib.pyplot as plt

from hurdler.design import parse_protein_input
from hurdler.idt import (
    IDTComplexityScorer,
    clear_idt_secret_environment,
    configure_idt_credentials_from_bytes,
    configure_idt_credentials_from_values,
)
from hurdler.idt_trajectory import idt_score_history_rows, plot_idt_score_trajectory
from hurdler.design import role_enzyme_options
from hurdler.optimization import translate_dna
from hurdler.design_artifacts import (
    build_step00_plasmid_record,
    timestamped_results_archive,
    write_secondary_checkpoint,
)
from hurdler.external_ga import ExternalGAResources, create_external_ga_bundle
from hurdler.protein_index import ProteinPatternIndex
from hurdler.progress import DesignProgressEvent, DesignRunControl, DesignRunStopped
from hurdler.vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    DesignSelection,
    build_route_universe,
    design_construct_v2,
    design_query,
    filter_route_universe,
    bundled_protein_index_dir,
    _secondary_adapters,
    write_design_outputs_v2,
)

COLAB_WIDGET_MANAGER_ENABLED = False
try:
    # Colab does not render ipywidgets created by a background/application
    # bootstrap reliably until its custom widget manager is explicitly enabled.
    # Keep this in the initialization cell, before any tutorial module is shown.
    from google.colab import output as colab_output
    colab_output.enable_custom_widget_manager()
    COLAB_WIDGET_MANAGER_ENABLED = True
except ImportError:
    # A regular Jupyter kernel uses its native widget manager.
    pass

try:
    from google.colab import data_table
    data_table.enable_dataframe_formatter()
except ImportError:
    pass
'''


COLAB_CONTROLLER = r'''PLASMID_FORM_VALUES = (
    ("pGEX-4T-1", "use_pGEX_4T_1"),
    ("pMAL-c5X", "use_pMAL_c5X"),
    ("pET-21a(+)", "use_pET_21a"),
    ("pET-28a(+)", "use_pET_28a"),
    ("pET-28a(+)_start_codon", "use_pET_28a_start_codon"),
    ("pCold_I", "use_pCold_I"),
    ("pUC18", "use_pUC18"),
    ("pQE-3", "use_pQE_3"),
)

VALIDATION_FORM_VALUES = {
    "No optimization": "none",
    "Live IDT API": "api",
    "IDT Bulk Input files": "batch",
}


def _csv_tuple(value):
    return tuple(part.strip() for part in str(value).split(",") if part.strip())


def _selected_plasmids():
    selected = tuple(name for name, variable in PLASMID_FORM_VALUES if bool(globals()[variable]))
    if not selected:
        raise ValueError("Select at least one plasmid profile in Step 2")
    return selected


def _current_query():
    if input_mode == "N-cap + repeat module + C-cap":
        return CompatibilityQuery(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            input_mode="split",
            sequence_id=str(sequence_id).strip() or "interactive_design",
            n_cap=n_cap_aa,
            repeat_module=repeat_module_aa,
            c_cap=c_cap_aa,
            repeat_copies=int(initial_repeat_copies),
            site_i_allowlist=_csv_tuple(site_i_allowlist),
            site_ii_allowlist=_csv_tuple(site_ii_allowlist),
            site_iii_allowlist=_csv_tuple(site_iii_allowlist),
            plasmid_allowlist=_selected_plasmids(),
            allow_left_cutter_in_hurdler_pair=bool(allow_left_cutter_in_hurdler_pair),
            allow_right_cutter_in_hurdler_pair=bool(allow_right_cutter_in_hurdler_pair),
        )

    parsed_id, normalized = parse_protein_input(full_protein_or_fasta)
    coordinates = (
        int(repeat_region_start_1based),
        int(repeat_region_end_1based),
        int(repeat_period_aa),
    )
    if any(value > 0 for value in coordinates) and not all(value > 0 for value in coordinates):
        raise ValueError("Provide repeat start, end, and period together, or leave all three as 0")
    confirmed = all(value > 0 for value in coordinates)
    return CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode="full",
        sequence_id=str(sequence_id).strip() or parsed_id,
        full_protein_sequence=normalized,
        repeat_region_start=coordinates[0] if confirmed else None,
        repeat_region_end=coordinates[1] if confirmed else None,
        repeat_period=coordinates[2] if confirmed else None,
        site_i_allowlist=_csv_tuple(site_i_allowlist),
        site_ii_allowlist=_csv_tuple(site_ii_allowlist),
        site_iii_allowlist=_csv_tuple(site_iii_allowlist),
        plasmid_allowlist=_selected_plasmids(),
        allow_left_cutter_in_hurdler_pair=bool(allow_left_cutter_in_hurdler_pair),
        allow_right_cutter_in_hurdler_pair=bool(allow_right_cutter_in_hurdler_pair),
    )


def _query_fingerprint(query):
    payload = json.dumps(asdict(query), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def _generation_schedule():
    values = tuple(sorted({int(value.strip()) for value in str(generation_schedule).split(",") if value.strip()}))
    if not values or values[-1] != 100 or any(value <= 0 for value in values):
        raise ValueError("Generation schedule must contain positive integers and terminate at 100")
    return values


def _score_weights():
    return {
        "selected_re_site_excess": float(weight_selected_re_site_excess),
        "gc_window_violation": float(weight_gc_window_violation),
        "repeated_re_site_excess": float(weight_repeated_re_site_excess),
        "repeated_14mer": float(weight_repeated_14mer),
        "repeated_13mer": float(weight_repeated_13mer),
        "repeated_8mer": float(weight_repeated_8mer),
        "hairpin_10mer_proxy": float(weight_hairpin_10mer_proxy),
        "homopolymer_excess": float(weight_homopolymer_excess),
        "terminal_repeat_proxy": float(weight_terminal_repeat_proxy),
        "gc_window_soft_violation": float(weight_gc_window_soft_violation),
        "negative_log_cai": float(weight_negative_log_cai),
    }


state = {
    "query_result": None,
    "query_fingerprint": None,
    "confirmed_route": None,
    "confirmed_fingerprint": None,
    "design_files": None,
    "archive": None,
}

route_choice = widgets.Dropdown(
    options=[("Run the query before selecting a route", None)],
    value=None,
    description="Route",
    layout=widgets.Layout(width="98%"),
)
confirm_button = widgets.Button(description="Confirm selected route", button_style="warning", disabled=True)
design_button = widgets.Button(description="Optimize / export", button_style="success", disabled=True)
download_button = widgets.Button(description="Download design ZIP", icon="download", disabled=True)
query_button = widgets.Button(description="Re-run HURDLER query", button_style="primary")
credential_upload = widgets.FileUpload(accept=".env,text/plain", multiple=False, description="Temporary IDT env")
query_output = widgets.Output()
route_output = widgets.Output()
design_output = widgets.Output()


def _invalidate_confirmation():
    state["confirmed_route"] = None
    state["confirmed_fingerprint"] = None
    confirm_button.disabled = True
    design_button.disabled = True
    download_button.disabled = True


def _route_label(row):
    site_iii = ",".join(row.get("site_iii_options", ()))
    return (
        f"#{row['rank']} · {row['site_i_enzyme']} / {row['site_ii_enzyme']} / {site_iii}"
        f" → {row['profile_id']} · {row['cut_scheme']}"
    )


def _run_query(_button=None):
    with query_output:
        clear_output(wait=True)
        with route_output:
            clear_output(wait=True)
        _invalidate_confirmation()
        route_choice.options = [("Select a route after a successful query", None)]
        route_choice.value = None
        try:
            query = _current_query()
            result = design_query(query)
        except Exception as exc:
            state["query_result"] = None
            state["query_fingerprint"] = None
            display(Markdown(f"**Query input error:** `{type(exc).__name__}: {exc}`"))
            return
        state["query_result"] = result
        state["query_fingerprint"] = _query_fingerprint(query)
        display(Markdown(f"**Query status:** `{result.status}` — {result.message}"))
        if result.status == "needs_boundary_confirmation" and result.boundary_analysis:
            candidates = pd.DataFrame(result.boundary_analysis.get("candidates", []))
            display(candidates)
            display(Markdown(
                "Enter one verified start/end/period triplet in **Step 1**, then press "
                "**Re-run HURDLER query**. Optimization remains disabled."
            ))
            return
        if result.protein_candidates:
            display(Markdown(f"**Protein-level candidates:** {len(result.protein_candidates):,}"))
            display(pd.DataFrame(result.protein_candidates))
        if not result.vector_routes:
            display(Markdown("No annotation-safe vector route is available; optimization remains disabled."))
            return
        routes = pd.DataFrame(result.vector_routes)
        route_columns = [
            "rank", "site_i_enzyme", "site_ii_enzyme", "site_iii_options",
            "profile_id", "cut_scheme", "left_cutter", "right_cutter",
            "restoration_length_bp", "cutter_reuse", "silencing_decisions",
        ]
        display(Markdown(f"**Annotation-safe routes:** {len(routes):,}. Review the complete table, then select one below."))
        display(routes[route_columns])
        route_choice.options = [
            ("Select one route — no automatic rank-1 choice", None),
            *((_route_label(row), index) for index, row in enumerate(result.vector_routes)),
        ]
        route_choice.value = None
        confirm_button.disabled = False


def _confirm_route(_button=None):
    with route_output:
        clear_output(wait=True)
        result = state.get("query_result")
        if result is None or not result.vector_routes:
            display(Markdown("**Run a successful query first.**"))
            return
        try:
            current_fingerprint = _query_fingerprint(_current_query())
        except Exception as exc:
            display(Markdown(f"**Current form values are invalid:** `{type(exc).__name__}: {exc}`"))
            return
        if current_fingerprint != state.get("query_fingerprint"):
            _invalidate_confirmation()
            display(Markdown("**The protein/RE/vector form changed. Re-run the query before confirming a route.**"))
            return
        if route_choice.value is None:
            display(Markdown("**Select a route explicitly; rank 1 is never chosen automatically.**"))
            return
        route = result.vector_routes[int(route_choice.value)]
        state["confirmed_route"] = dict(route)
        state["confirmed_fingerprint"] = current_fingerprint
        design_button.disabled = False
        display(Markdown(
            f"**Confirmed route #{route['rank']}:** {route['site_i_enzyme']} / "
            f"{route['site_ii_enzyme']} → {route['profile_id']} · {route['cut_scheme']}"
        ))


def _secret_value(reader, name):
    try:
        return str(reader(name) or "").strip()
    except Exception:
        return ""


def _configure_colab_secrets(reader=None):
    if reader is None:
        try:
            from google.colab import userdata
        except ImportError as exc:
            raise RuntimeError("Colab Secrets are available only in a Colab runtime") from exc
        reader = userdata.get
    token = _secret_value(reader, "IDT_ACCESS_TOKEN")
    if token:
        return configure_idt_credentials_from_values(
            {"IDT_ACCESS_TOKEN": token}, auth_method="access_token"
        )
    values = {
        name: _secret_value(reader, name)
        for name in ("IDT_CLIENT_ID", "IDT_CLIENT_SECRET", "IDT_USERNAME", "IDT_PASSWORD")
    }
    try:
        return configure_idt_credentials_from_values(values, auth_method="password")
    finally:
        values.clear()


def _configure_api_credentials():
    if idt_credential_source == "Colab Secrets":
        return _configure_colab_secrets()
    if idt_credential_source == "Hidden runtime prompt":
        auth_method = "access_token" if idt_prompt_auth_method == "Access token" else "password"
        return configure_idt_credentials(mode="manual", auth_method=auth_method, prompt=getpass.getpass)
    uploaded = credential_upload.value
    if not uploaded:
        raise ValueError("Choose one temporary env file in the runtime control panel")
    item = next(iter(uploaded.values())) if isinstance(uploaded, dict) else uploaded[0]
    payload = bytes(item["content"] if isinstance(item, dict) else item.content)
    try:
        return configure_idt_credentials_from_bytes(payload)
    finally:
        payload = b""
        credential_upload.value = ()


def _run_design(_button=None):
    with design_output:
        clear_output(wait=True)
        route = state.get("confirmed_route")
        if route is None:
            display(Markdown("**Select and confirm one route before optimization.**"))
            return
        try:
            query = _current_query()
            current_fingerprint = _query_fingerprint(query)
            if current_fingerprint != state.get("confirmed_fingerprint"):
                _invalidate_confirmation()
                raise RuntimeError("Protein/RE/vector settings changed; re-run the query and confirm a route again")
            mode = VALIDATION_FORM_VALUES[validation_mode]
            scorer = None
            if mode == "api":
                _configure_api_credentials()
                scorer = IDTComplexityScorer(Path(output_directory) / "idt_audit.jsonl")
            request = DesignRequestV2(
                schema_version=DESIGN_SCHEMA_VERSION_V2,
                query=query,
                selection=DesignSelection(
                    route["candidate_id"], route["profile_id"], route["scheme_id"],
                    route["site_iii_options"][0],
                ),
                validation_mode=mode,
                max_repeat_copies=int(maximum_repeat_copies) if query.input_mode == "split" else None,
                population_size=int(population_size),
                mutation_rate=float(mutation_rate),
                crossover_rate=float(crossover_rate),
                elite_fraction=float(elite_fraction),
                seed=int(random_seed),
                generation_schedule=_generation_schedule(),
                score_weights=_score_weights(),
                auto_adjust_weights_from_idt=bool(auto_adjust_weights_from_idt),
            )
            result = design_construct_v2(request, idt_scorer=scorer)
            files = write_design_outputs_v2(result, output_directory)
            source = Path(output_directory).resolve()
            archive = Path(shutil.make_archive(str(source), "zip", root_dir=source))
            state["design_files"] = files
            state["archive"] = archive
            download_button.disabled = not bool(enable_zip_download)
            display(Markdown(f"**Design status:** `{result.status}` — {result.message}"))
            if result.selected_route:
                display(pd.DataFrame([result.selected_route]))
            if result.primary_fragments:
                display(Markdown("### Purchase fragments"))
                display(pd.DataFrame(result.primary_fragments))
            if result.cloning_steps:
                display(Markdown("### Cloning plan"))
                display(pd.DataFrame(result.cloning_steps))
            display(Markdown(f"Design bundle prepared: `{archive.name}`. No order was submitted."))
        except Exception as exc:
            display(Markdown(f"**Design failed safely:** `{type(exc).__name__}: {exc}`"))
        finally:
            clear_idt_secret_environment()


def _download_design(_button=None):
    archive = state.get("archive")
    if archive is None or not Path(archive).is_file():
        with design_output:
            display(Markdown("**Run the design first; no ZIP is available.**"))
        return
    try:
        from google.colab import files as colab_files
    except ImportError:
        with design_output:
            display(Markdown(f"ZIP path: `{Path(archive).resolve()}`"))
        return
    colab_files.download(str(archive))


query_button.on_click(_run_query)
confirm_button.on_click(_confirm_route)
design_button.on_click(_run_design)
download_button.on_click(_download_design)

credential_children = [
    widgets.HTML(
        "<b>IDT credentials:</b> Colab Secrets is the default. The upload control is used only "
        "when Step 3 selects Temporary env upload; values are never written to notebook output."
    ),
    credential_upload,
]
display(widgets.VBox([
    widgets.HTML("<h2>4. Query and manually confirm a route</h2>"),
    query_button,
    query_output,
    route_choice,
    confirm_button,
    route_output,
    widgets.HTML("<h2>5. Optional optimization and file export</h2>"),
    *credential_children,
    design_button,
    download_button,
    design_output,
]))

# Runtime → Run all reaches this call after every native form value has been assigned.
_run_query()
'''


COLAB_CONTROLLER_V2 = r'''PLASMID_OPTIONS = (
    "pGEX-4T-1", "pMAL-c5X", "pET-21a(+)", "pET-28a(+)",
    "pET-28a(+)_start_codon", "pCold_I", "pUC18", "pQE-3",
)


def _help_card(title, widget, *, unit, default, purpose, allowed, effect):
    return widgets.VBox([
        widgets.HTML(
            f"<b style='color:#111827'>{title}</b> <span style='color:#4b5563'>[{unit}]</span><br>"
            f"<small><b>Default:</b> {default} · <b>Allowed:</b> {allowed}<br>"
            f"<b>Purpose:</b> {purpose}<br><b>Effect:</b> {effect}</small>"
        ),
        widget,
    ], layout=widgets.Layout(border="1px solid #d1d5db", padding="7px", margin="3px 0"))


# Storage is opt-in. Merely choosing Drive does not authenticate or mount it.
storage_mode_widget = widgets.ToggleButtons(
    options=(("Colab temporary storage", "runtime"), ("Google Drive", "drive")),
    value="runtime", description="Storage",
)
drive_root_widget = widgets.Text(
    value="/content/drive/MyDrive/HURDLER",
    description="Drive folder", layout=widgets.Layout(width="98%"),
)
mount_drive_button = widgets.Button(description="Mount Google Drive", icon="folder-open")
storage_status = widgets.HTML(
    "<b>Runtime storage active.</b> A Colab restart removes local files; download the final ZIP or opt into Drive."
)
storage_state = {"drive_mounted": False}


def _mount_google_drive(_button=None):
    if storage_mode_widget.value != "drive":
        storage_status.value = "Select <b>Google Drive</b> first; no mount was attempted."
        return
    try:
        from google.colab import drive
    except ImportError:
        storage_status.value = "Google Drive mounting is available only in hosted Colab."
        return
    drive.mount("/content/drive")
    Path(drive_root_widget.value).mkdir(parents=True, exist_ok=True)
    storage_state["drive_mounted"] = True
    storage_status.value = (
        "<b>Drive mounted.</b> Only checkpoint and final ZIP files will be copied; computation remains under /content."
    )


mount_drive_button.on_click(_mount_google_drive)


def _storage_mode_changed(_change=None):
    if storage_mode_widget.value == "runtime":
        storage_status.value = (
            "<b>Runtime storage active.</b> A Colab restart removes local files; download the final ZIP or opt into Drive."
        )
    elif storage_state.get("drive_mounted"):
        storage_status.value = "<b>Drive already mounted.</b> Checkpoint and final ZIP files will be copied."
    else:
        storage_status.value = "<b>Drive selected but not mounted.</b> Click Mount Google Drive before optimization."


storage_mode_widget.observe(_storage_mode_changed, names="value")
storage_panel = widgets.VBox([
    widgets.HTML("<h3>Result storage</h3>"),
    _help_card(
        "Storage destination", storage_mode_widget, unit="mode", default="Colab temporary storage",
        purpose="Controls whether recovery/final ZIP files are also copied to Google Drive.",
        allowed="runtime or Drive", effect="Drive is mounted only after the button is clicked; credentials are never saved.",
    ),
    _help_card(
        "Drive output folder", drive_root_widget, unit="path", default="MyDrive/HURDLER",
        purpose="Destination for checkpoint and final ZIP archives.", allowed="a folder under /content/drive",
        effect="Does not change the /content computation directory.",
    ),
    mount_drive_button, storage_status,
])


# Exactly one protein-input panel is visible. The supplied tutorial example is
# intentionally visible at first; the complete-sequence box is not.
input_mode_widget = widgets.ToggleButtons(
    options=(("N-cap / module / C-cap", "split"), ("Complete protein / FASTA", "full")),
    value="split", description="Input mode",
)
sequence_id_widget = widgets.Text(
    value="", placeholder="interactive_design", description="Sequence ID",
    layout=widgets.Layout(width="98%"),
)
n_cap_widget = widgets.Textarea(
    value=default_n_cap, description="N-cap AA", layout=widgets.Layout(width="98%", height="80px"),
)
repeat_module_widget = widgets.Textarea(
    value=default_repeat_module, description="Module AA", layout=widgets.Layout(width="98%", height="75px"),
)
initial_copies_widget = widgets.BoundedIntText(
    value=default_repeat_copies, min=2, max=10000, description="Target copies",
)
c_cap_widget = widgets.Textarea(
    value=default_c_cap, description="C-cap AA", layout=widgets.Layout(width="98%", height="100px"),
)
full_protein_widget = widgets.Textarea(
    value="", placeholder=">protein_id\nFULL_AMINO_ACID_SEQUENCE",
    description="Full AA/FASTA", layout=widgets.Layout(width="98%", height="150px"),
)
repeat_start_widget = widgets.BoundedIntText(value=0, min=0, max=1_000_000, description="Start")
repeat_end_widget = widgets.BoundedIntText(value=0, min=0, max=1_000_000, description="End")
repeat_period_widget = widgets.BoundedIntText(value=0, min=0, max=100_000, description="Period")

split_input_panel = widgets.VBox([
    _help_card("N-terminal cap", n_cap_widget, unit="AA", default="tutorial construct N-cap", purpose="Fixed protein before the repeat array.", allowed="standard one-letter amino acids; may be empty", effect="Included once in the primary CDS and translation."),
    _help_card("Repeat module", repeat_module_widget, unit="AA", default="43-AA tutorial armadillo module", purpose="Exact protein unit repeated in the array.", allowed="standard one-letter amino acids", effect="Determines HURDLER sites, secondary core bp and all copy counts."),
    _help_card("Final target copies", initial_copies_widget, unit="repeat copies", default=str(default_repeat_copies), purpose="Exact repeat count required in the final protein.", allowed="integer 2–10,000", effect="Sets the final CDS and the primary + secondary reuse equation."),
    _help_card("C-terminal cap", c_cap_widget, unit="AA", default="tutorial construct C-cap", purpose="Fixed protein after the repeat array.", allowed="standard one-letter amino acids; may be empty", effect="Included once in the primary CDS and translation."),
])
full_input_panel = widgets.VBox([
    _help_card("Complete protein or FASTA", full_protein_widget, unit="AA", default="empty", purpose="Preserves every residue and repeat variant in one complete protein.", allowed="one raw sequence or one FASTA record", effect="HURDLER proposes boundaries when coordinates remain zero."),
    _help_card("Repeat-region coordinates", widgets.HBox([repeat_start_widget, repeat_end_widget, repeat_period_widget]), unit="1-based AA", default="0 / 0 / 0", purpose="Confirms repeat start, inclusive end and period.", allowed="all zero or all positive", effect="Required together before full-protein optimization."),
])


def _sync_input_panel(_change=None):
    split_input_panel.layout.display = "" if input_mode_widget.value == "split" else "none"
    full_input_panel.layout.display = "" if input_mode_widget.value == "full" else "none"


_sync_input_panel()
protein_input_panel = widgets.VBox([
    widgets.HTML("<p>Use the two buttons to switch input methods; only the active method is read.</p>"),
    _help_card("Protein input method", input_mode_widget, unit="mode", default="split tutorial example", purpose="Chooses explicit cap/module/copy input or one complete protein.", allowed="one of two buttons", effect="Switching invalidates cached routes and confirmed designs."),
    _help_card("Sequence identifier", sequence_id_widget, unit="text", default="interactive_design or FASTA header", purpose="Names result records and ZIP files.", allowed="short text", effect="Does not change molecular compatibility."),
    split_input_panel, full_input_panel,
])

allow_left_cutter_widget = widgets.Checkbox(value=False, description="Allow left vector cutter reuse")
allow_right_cutter_widget = widgets.Checkbox(value=False, description="Allow right vector cutter reuse")
cutter_policy_panel = widgets.VBox([
    widgets.HTML("<h3>1b. Optional cutter-reuse fallback</h3>"),
    _help_card("Left cutter reuse", allow_left_cutter_widget, unit="boolean", default="off", purpose="Allows the left vector cutter to also fill a HURDLER RE role only as a fallback.", allowed="on/off", effect="Can recover routes but is evaluated after restore filtering."),
    _help_card("Right cutter reuse", allow_right_cutter_widget, unit="boolean", default="off", purpose="Allows the right vector cutter to also fill a HURDLER RE role only as a fallback.", allowed="on/off", effect="Can recover routes but is evaluated after restore filtering."),
])

protein_index = ProteinPatternIndex.load(bundled_protein_index_dir())
enzyme_roles = role_enzyme_options(protein_index)
project_root = bundled_protein_index_dir().parents[2]
declared_site_iii = tuple(sorted(
    pd.read_csv(project_root / "output" / "selected_site_iii_enzymes.csv")["enzyme"].astype(str).unique()
))
all_enzyme_options = tuple(sorted(
    set(enzyme_roles["site_i"]) | set(enzyme_roles["site_ii"]) | set(declared_site_iii)
))

def _individual_checkbox_group(options, *, columns):
    boxes = {
        str(option): widgets.Checkbox(
            value=True,
            description=str(option),
            indent=False,
            layout=widgets.Layout(width="auto"),
        )
        for option in options
    }
    grid = widgets.GridBox(
        children=tuple(boxes.values()),
        layout=widgets.Layout(
            width="98%",
            grid_template_columns=f"repeat({columns}, minmax(0, 1fr))",
            grid_gap="4px 12px",
            border="1px solid #ddd",
            padding="8px",
        ),
    )
    return boxes, grid


enzyme_checkboxes, enzyme_checkbox_grid = _individual_checkbox_group(
    all_enzyme_options, columns=4
)
enzyme_bulk_control = widgets.ToggleButtons(
    options=(("Select all RE", "all"), ("Select none", "none"), ("Custom", "custom")),
    value="all",
    description="RE selection",
    button_style="",
)
enzyme_selection_status = widgets.HTML()

plasmid_checkboxes, plasmid_checkbox_grid = _individual_checkbox_group(
    PLASMID_OPTIONS, columns=2
)
plasmid_bulk_control = widgets.ToggleButtons(
    options=(("Select all plasmids", "all"), ("Select none", "none"), ("Custom", "custom")),
    value="all",
    description="Plasmid selection",
    button_style="",
)
plasmid_selection_status = widgets.HTML()
enzyme_route_support = widgets.HTML()
plasmid_route_support = widgets.HTML()

max_restoration_length_widget = widgets.BoundedIntText(
    value=100, min=0, max=10000,
    description="Max restore (bp)",
    layout=widgets.Layout(width="320px"),
)
advanced_route_filters = widgets.Accordion(children=[widgets.VBox([
    _help_card(
        "Maximum restoration length", max_restoration_length_widget,
        unit="bp", default="100", purpose="Caps left + right vector sequence restored by the donor.",
        allowed="integer 0–10,000; inclusive", effect="Longer routes disappear from live counts and final query; Type-IIS adapters are excluded.",
    ),
])])
advanced_route_filters.set_title(0, "Advanced route filters")
advanced_route_filters.selected_index = None

_selection_sync = {"enzymes": False, "plasmids": False}

state = {
    "query_result": None,
    "query_fingerprint": None,
    "confirmed_route": None,
    "confirmed_site_iii": None,
    "confirmed_fingerprint": None,
    "design_files": None,
    "archive": None,
    "progress_events": [],
    "route_universe": None,
    "route_universe_fingerprint": None,
    "design_result": None,
    "run_directory": None,
    "best_checkpoint": None,
    "last_checkpoint_write": 0.0,
    "checkpoint_archive": None,
    "external_bundle": None,
    "external_bundle_fingerprint": None,
    "external_bundle_error": None,
    "credential_payload": None,
    "credential_auth_method": None,
    "run_control": None,
    "run_thread": None,
    "run_active": False,
    "run_terminal_status": None,
    "run_id": 0,
    "progress_queue": queue.Queue(),
    "progress_lock": threading.Lock(),
    "ui_pump_task": None,
    "ui_asyncio_loop": None,
    "ui_io_loop": None,
    "ui_schedule_lock": threading.Lock(),
    "ui_drain_pending": False,
    "checkpoint_thread": None,
    "visible_log_lines": [],
    "idt_score_events": [],
    "run_started_monotonic": None,
    "last_progress_monotonic": None,
    "viewer_rows": [],
    "viewer_directory": None,
    "checkpoint_lock": threading.Lock(),
}


def _runtime_scratch_path(name, fallback):
    content = Path("/content")
    if content.is_dir() and os.access(content, os.W_OK):
        return content / name
    fallback = Path(fallback)
    try:
        fallback.relative_to(content)
    except ValueError:
        return fallback
    return Path(tempfile.gettempdir()) / name

query_button = widgets.Button(description="Run / re-run HURDLER query", button_style="primary")
query_output = widgets.Output()
pair_choice = widgets.Dropdown(
    options=[("Run the query first", None)], value=None,
    description="Site I / II", layout=widgets.Layout(width="98%"),
)
site_iii_choice = widgets.Dropdown(
    options=[("Choose Site I / II first", None)], value=None,
    description="Site III", layout=widgets.Layout(width="98%"),
)
profile_choice = widgets.Dropdown(
    options=[("Choose all three RE first", None)], value=None,
    description="Plasmid", layout=widgets.Layout(width="98%"),
)
scheme_choice = widgets.Dropdown(
    options=[("Choose a plasmid first", None)], value=None,
    description="Cut scheme", layout=widgets.Layout(width="98%"),
)
confirm_button = widgets.Button(description="Confirm RE / plasmid route", button_style="warning", disabled=True)
route_output = widgets.Output()


def _chosen_role_enzymes(role):
    selected = {
        name for name, checkbox in enzyme_checkboxes.items() if checkbox.value
    }
    return tuple(name for name in enzyme_roles[role] if name in selected)


def _selected_role_enzymes(role):
    values = _chosen_role_enzymes(role)
    if not values:
        raise ValueError(f"Select at least one enzyme eligible for {role.replace('_', ' ').title()}")
    return values


def _selected_plasmids():
    selected = tuple(
        name for name in PLASMID_OPTIONS if plasmid_checkboxes[name].value
    )
    if not selected:
        raise ValueError("Select at least one plasmid profile")
    return selected


def _form_query(
    *,
    site_i_allowlist,
    site_ii_allowlist,
    site_iii_allowlist,
    plasmid_allowlist,
    max_restoration_length_bp,
):
    common = dict(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        sequence_id=str(sequence_id_widget.value).strip() or "interactive_design",
        site_i_allowlist=site_i_allowlist,
        site_ii_allowlist=site_ii_allowlist,
        site_iii_allowlist=site_iii_allowlist,
        plasmid_allowlist=plasmid_allowlist,
        allow_left_cutter_in_hurdler_pair=bool(allow_left_cutter_widget.value),
        allow_right_cutter_in_hurdler_pair=bool(allow_right_cutter_widget.value),
        max_restoration_length_bp=max_restoration_length_bp,
    )
    if input_mode_widget.value == "split":
        return CompatibilityQuery(
            input_mode="split", n_cap=n_cap_widget.value, repeat_module=repeat_module_widget.value,
            c_cap=c_cap_widget.value, repeat_copies=int(initial_copies_widget.value), **common,
        )
    parsed_id, normalized = parse_protein_input(full_protein_widget.value)
    coordinates = (
        int(repeat_start_widget.value), int(repeat_end_widget.value), int(repeat_period_widget.value),
    )
    if any(value > 0 for value in coordinates) and not all(value > 0 for value in coordinates):
        raise ValueError("Provide repeat start, end, and period together, or leave all three as 0")
    confirmed = all(value > 0 for value in coordinates)
    common["sequence_id"] = str(sequence_id_widget.value).strip() or parsed_id
    return CompatibilityQuery(
        input_mode="full", full_protein_sequence=normalized,
        repeat_region_start=coordinates[0] if confirmed else None,
        repeat_region_end=coordinates[1] if confirmed else None,
        repeat_period=coordinates[2] if confirmed else None,
        **common,
    )


def _current_query():
    return _form_query(
        site_i_allowlist=_selected_role_enzymes("site_i"),
        site_ii_allowlist=_selected_role_enzymes("site_ii"),
        site_iii_allowlist=_selected_role_enzymes("site_iii"),
        plasmid_allowlist=_selected_plasmids(),
        max_restoration_length_bp=int(max_restoration_length_widget.value),
    )


def _universe_query():
    return _form_query(
        site_i_allowlist=(),
        site_ii_allowlist=(),
        site_iii_allowlist=(),
        plasmid_allowlist=(),
        max_restoration_length_bp=None,
    )


def _query_fingerprint(query):
    payload = json.dumps(asdict(query), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def _ensure_route_universe():
    query = _universe_query()
    fingerprint = _query_fingerprint(query)
    if (
        state.get("route_universe") is None
        or state.get("route_universe_fingerprint") != fingerprint
    ):
        state["route_universe"] = build_route_universe(query, protein_index=protein_index)
        state["route_universe_fingerprint"] = fingerprint
    return state["route_universe"]


def _invalidate_confirmation(message=""):
    state["confirmed_route"] = None
    state["confirmed_site_iii"] = None
    state["confirmed_fingerprint"] = None
    confirm_button.disabled = True
    design_button.disabled = True
    download_button.disabled = True
    state["design_result"] = None
    state["external_bundle"] = None
    state["external_bundle_fingerprint"] = None
    if "export_bundle_button" in globals():
        export_bundle_button.disabled = True
    if "viewer_panel" in globals():
        _reset_viewer_placeholder()
    if "results_status" in globals():
        results_status.value = (
            "<div style='border:2px dashed #4b2e83;background:#ffffff;color:#111827;"
            "border-radius:8px;padding:12px'>No current result. Confirm a route and start a new run.</div>"
        )
    if "design_output" in globals():
        with design_output:
            clear_output(wait=True)
    if "idt_plot_status" in globals() and not state.get("run_active"):
        _set_idt_plot_placeholder(
            "No IDT scores yet. The trajectory updates after each scored purchase fragment."
        )
    if message:
        with route_output:
            clear_output(wait=True)
            display(Markdown(message))


def _input_changed(_change=None):
    _sync_input_panel()
    state["route_universe"] = None
    state["route_universe_fingerprint"] = None
    state["query_result"] = None
    state["query_fingerprint"] = None
    _populate_pair_choices(None)
    _set_support_summary([])
    _invalidate_confirmation("**Protein input changed. Run the HURDLER query again.**")


def _selection_changed(_change=None):
    _refresh_live_support()


def _selection_mode(boxes):
    selected_count = sum(bool(checkbox.value) for checkbox in boxes.values())
    if selected_count == len(boxes):
        return "all", selected_count
    if selected_count == 0:
        return "none", selected_count
    return "custom", selected_count


def _refresh_selection_group(group):
    if group == "enzymes":
        boxes = enzyme_checkboxes
        control = enzyme_bulk_control
        status = enzyme_selection_status
        label = "RE enzymes"
    else:
        boxes = plasmid_checkboxes
        control = plasmid_bulk_control
        status = plasmid_selection_status
        label = "plasmids"
    mode, selected_count = _selection_mode(boxes)
    _selection_sync[group] = True
    try:
        control.value = mode
    finally:
        _selection_sync[group] = False
    status.value = (
        f"<b>{selected_count}/{len(boxes)}</b> {label} selected"
        + ("" if mode != "custom" else " · individual selection")
    )


def _set_selection_group(group, selected):
    boxes = enzyme_checkboxes if group == "enzymes" else plasmid_checkboxes
    _selection_sync[group] = True
    try:
        for checkbox in boxes.values():
            checkbox.value = bool(selected)
    finally:
        _selection_sync[group] = False
    _refresh_selection_group(group)
    _selection_changed()


def _individual_selection_changed(group, _change=None):
    if _selection_sync[group]:
        return
    _refresh_selection_group(group)
    _selection_changed()


def _bulk_selection_changed(group, change):
    if _selection_sync[group]:
        return
    if change.get("name") != "value":
        return
    if change.get("new") == "custom":
        _refresh_selection_group(group)
        return
    _set_selection_group(group, change["new"] == "all")


def _set_all_enzymes(_button=None):
    _set_selection_group("enzymes", True)


def _set_no_enzymes(_button=None):
    _set_selection_group("enzymes", False)


def _set_all_plasmids(_button=None):
    _set_selection_group("plasmids", True)


def _set_no_plasmids(_button=None):
    _set_selection_group("plasmids", False)


def _reset_route_selectors():
    pair_choice.options = [("No supported pair under the current filters", None)]
    pair_choice.value = None
    site_iii_choice.options = [("Choose Site I / II first", None)]
    site_iii_choice.value = None
    profile_choice.options = [("Choose all three RE first", None)]
    profile_choice.value = None
    scheme_choice.options = [("Choose a plasmid first", None)]
    scheme_choice.value = None


def _support_summary(routes):
    pairs = {
        (row["site_i_enzyme"], row["site_ii_enzyme"])
        for row in routes
    }
    site_iii = {
        enzyme
        for row in routes
        for enzyme in row.get("site_iii_options", ())
    }
    plasmids = {row["profile_id"] for row in routes}
    minimum = min(
        (int(row["restoration_length_bp"]) for row in routes),
        default=None,
    )
    minimum_text = "—" if minimum is None else f"{minimum} bp"
    zero = not routes
    border = "#b31b1b" if zero else "#4b2e83"
    background = "#fff2f2" if zero else "#f4f0fa"
    color = "#b31b1b" if zero else "#4b2e83"
    cards = (
        ("Supported RE pairs", f"{len(pairs):,}"),
        ("Available Site III", f"{len(site_iii):,}"),
        ("Supported plasmids", f"{len(plasmids):,}"),
        ("Minimum restore", minimum_text),
    )
    body = "".join(
        "<div style='min-width:145px;padding:10px 14px'>"
        f"<div style='font-size:12px;color:#555'><b>{label}:</b></div>"
        f"<div style='font-size:27px;line-height:1.15;font-weight:800;color:{color}'>{value}</div>"
        "</div>"
        for label, value in cards
    )
    return (
        f"<div style='display:flex;flex-wrap:wrap;border:2px solid {border};"
        f"background:{background};border-radius:8px;margin:8px 0'>{body}</div>"
    )


def _set_support_summary(routes):
    summary = _support_summary(routes)
    enzyme_route_support.value = summary
    plasmid_route_support.value = summary


def _set_support_error(message):
    summary = (
        "<div style='border:2px solid #b31b1b;background:#fff2f2;border-radius:8px;"
        "padding:12px;margin:8px 0;font-size:16px'>"
        f"<b>Route filter error:</b> {message}"
        "</div>"
    )
    enzyme_route_support.value = summary
    plasmid_route_support.value = summary


def _populate_pair_choices(result):
    if result is None or not result.vector_routes:
        _reset_route_selectors()
        return []
    routes = pd.DataFrame(result.vector_routes)
    pair_rows = routes.sort_values("rank").drop_duplicates(
        ["site_i_enzyme", "site_ii_enzyme"]
    )
    pairs = [
        ((row.site_i_enzyme, row.site_ii_enzyme), int(row.rank))
        for row in pair_rows.itertuples(index=False)
    ]
    pair_choice.options = [
        ("Select Site I / Site II explicitly", None),
        *((f"#{rank} {pair[0]} / {pair[1]}", pair) for pair, rank in pairs),
    ]
    pair_choice.value = None
    site_iii_choice.options = [("Choose Site I / II first", None)]
    profile_choice.options = [("Choose all three RE first", None)]
    scheme_choice.options = [("Choose a plasmid first", None)]
    return pairs


def _refresh_live_support(_change=None):
    _invalidate_confirmation()
    selected_roles = {
        role: _chosen_role_enzymes(role)
        for role in ("site_i", "site_ii", "site_iii")
    }
    selected_plasmids = tuple(
        name for name in PLASMID_OPTIONS if plasmid_checkboxes[name].value
    )
    if not all(selected_roles.values()) or not selected_plasmids:
        state["query_result"] = None
        state["query_fingerprint"] = None
        _populate_pair_choices(None)
        _set_support_summary([])
        return None
    try:
        query = _form_query(
            site_i_allowlist=selected_roles["site_i"],
            site_ii_allowlist=selected_roles["site_ii"],
            site_iii_allowlist=selected_roles["site_iii"],
            plasmid_allowlist=selected_plasmids,
            max_restoration_length_bp=int(max_restoration_length_widget.value),
        )
        result = filter_route_universe(_ensure_route_universe(), query)
    except Exception as exc:
        state["query_result"] = None
        state["query_fingerprint"] = None
        _populate_pair_choices(None)
        _set_support_error(f"{type(exc).__name__}: {exc}")
        return None
    state["query_result"] = result
    state["query_fingerprint"] = _query_fingerprint(query)
    _populate_pair_choices(result)
    _set_support_summary(result.vector_routes)
    return result


def _route_filter_changed(change):
    if change.get("name") == "value":
        _refresh_live_support()


def _routes_for_pair():
    result = state.get("query_result")
    if result is None or pair_choice.value is None:
        return []
    site_i, site_ii = pair_choice.value
    return [
        row for row in result.vector_routes
        if row["site_i_enzyme"] == site_i and row["site_ii_enzyme"] == site_ii
    ]


def _update_site_iii(_change=None):
    _invalidate_confirmation()
    options = sorted({value for row in _routes_for_pair() for value in row["site_iii_options"]})
    site_iii_choice.options = [("Select Site III explicitly", None), *((value, value) for value in options)]
    site_iii_choice.value = None
    profile_choice.options = [("Choose all three RE first", None)]
    scheme_choice.options = [("Choose a plasmid first", None)]


def _routes_for_three_enzymes():
    if site_iii_choice.value is None:
        return []
    return [row for row in _routes_for_pair() if site_iii_choice.value in row["site_iii_options"]]


def _update_profiles(_change=None):
    _invalidate_confirmation()
    profiles = sorted({row["profile_id"] for row in _routes_for_three_enzymes()})
    profile_choice.options = [("Select a supporting plasmid", None), *((value, value) for value in profiles)]
    profile_choice.value = None
    scheme_choice.options = [("Choose a plasmid first", None)]


def _update_schemes(_change=None):
    _invalidate_confirmation()
    rows = [row for row in _routes_for_three_enzymes() if row["profile_id"] == profile_choice.value]
    rows = sorted(rows, key=lambda row: int(row["rank"]))
    scheme_choice.options = [
        ("Select a cut scheme explicitly", None),
        *((
            f"#{row['rank']} {row['cut_scheme']} · {row['left_cutter']}/{row['right_cutter']}"
            f" · restore {row['restoration_length_bp']} bp",
            int(row["rank"]),
        ) for row in rows),
    ]
    scheme_choice.value = None
    confirm_button.disabled = not bool(rows)


def _run_query(_button=None):
    with query_output:
        clear_output(wait=True)
        try:
            query = _current_query()
            result = filter_route_universe(_ensure_route_universe(), query)
        except Exception as exc:
            state["query_result"] = None
            state["query_fingerprint"] = None
            _populate_pair_choices(None)
            _set_support_summary([])
            display(Markdown(f"**Query input error:** `{type(exc).__name__}: {exc}`"))
            return
        _invalidate_confirmation()
        state["query_result"] = result
        state["query_fingerprint"] = _query_fingerprint(query)
        pairs = _populate_pair_choices(result)
        _set_support_summary(result.vector_routes)
        display(Markdown(f"**Query status:** `{result.status}` — {result.message}"))
        if result.status == "needs_boundary_confirmation" and result.boundary_analysis:
            display(pd.DataFrame(result.boundary_analysis.get("candidates", [])))
            display(Markdown("Confirm start/end/period in the protein form, then re-run this query."))
            return
        if not result.vector_routes:
            display(Markdown("No annotation-safe vector route is available."))
            return
        routes = pd.DataFrame(result.vector_routes)
        display(Markdown(
            f"**{len(result.protein_candidates):,} protein candidates; "
            f"{len(pairs):,} RE pairs; {len(result.vector_routes):,} vector routes.**"
        ))
        display(routes[[
            "rank", "site_i_enzyme", "site_ii_enzyme", "site_iii_options",
            "profile_id", "cut_scheme", "left_cutter", "right_cutter",
            "restoration_length_bp", "cutter_reuse",
        ]])


def _confirm_route(_button=None):
    with route_output:
        clear_output(wait=True)
        result = state.get("query_result")
        if result is None or scheme_choice.value is None or site_iii_choice.value is None:
            display(Markdown("**Select all three enzymes, plasmid, and cut scheme first.**"))
            return
        current = _query_fingerprint(_current_query())
        if current != state.get("query_fingerprint"):
            _invalidate_confirmation()
            display(Markdown("**Inputs changed. Re-run the query before confirming.**"))
            return
        matches = [
            row for row in _routes_for_three_enzymes()
            if row["profile_id"] == profile_choice.value and int(row["rank"]) == int(scheme_choice.value)
        ]
        if len(matches) != 1:
            display(Markdown("**The selected route is stale or ambiguous; re-run the query.**"))
            return
        route = dict(matches[0])
        state["confirmed_route"] = route
        state["confirmed_site_iii"] = str(site_iii_choice.value)
        state["confirmed_fingerprint"] = current
        _sync_execution_target()
        _update_secondary_lengths()
        try:
            _prepare_route_preview()
        except Exception as exc:
            viewer_status.value = (
                f"<div style='color:#b31b1b'><b>Preview unavailable:</b> "
                f"{type(exc).__name__}: {str(exc)[:300]}. GA controls remain usable.</div>"
            )
        display(Markdown(
            f"**Confirmed:** {route['site_i_enzyme']} / {route['site_ii_enzyme']} / "
            f"{site_iii_choice.value} → {route['profile_id']} · {route['cut_scheme']}. "
            "GA/IDT controls are now enabled below."
        ))


for checkbox in enzyme_checkboxes.values():
    checkbox.observe(
        lambda change: _individual_selection_changed("enzymes", change),
        names="value",
    )
for checkbox in plasmid_checkboxes.values():
    checkbox.observe(
        lambda change: _individual_selection_changed("plasmids", change),
        names="value",
    )
enzyme_bulk_control.observe(
    lambda change: _bulk_selection_changed("enzymes", change), names="value"
)
plasmid_bulk_control.observe(
    lambda change: _bulk_selection_changed("plasmids", change), names="value"
)
max_restoration_length_widget.observe(_route_filter_changed, names="value")
for input_widget in (
    input_mode_widget, sequence_id_widget, n_cap_widget, repeat_module_widget,
    initial_copies_widget, c_cap_widget, full_protein_widget,
    repeat_start_widget, repeat_end_widget, repeat_period_widget,
    allow_left_cutter_widget, allow_right_cutter_widget,
):
    input_widget.observe(_input_changed, names="value")
_refresh_selection_group("enzymes")
_refresh_selection_group("plasmids")
pair_choice.observe(_update_site_iii, names="value")
site_iii_choice.observe(_update_profiles, names="value")
profile_choice.observe(_update_schemes, names="value")
query_button.on_click(_run_query)
confirm_button.on_click(_confirm_route)


def _numeric_control(description, value, minimum, maximum, step, *, integer=False, help_text=""):
    slider_type = widgets.IntSlider if integer else widgets.FloatSlider
    text_type = widgets.BoundedIntText if integer else widgets.BoundedFloatText
    slider = slider_type(value=value, min=minimum, max=maximum, step=step, readout=False,
                         layout=widgets.Layout(width="62%"))
    number = text_type(value=value, min=minimum, max=maximum, step=step,
                       layout=widgets.Layout(width="34%"))
    widgets.link((slider, "value"), (number, "value"))
    card = widgets.VBox([
        widgets.HTML(
            f"<b>{description}</b><br><small>Default {value}; allowed {minimum}–{maximum}. "
            f"{help_text}</small>"
        ), widgets.HBox([slider, number]),
    ], layout=widgets.Layout(width="49%", border="1px solid #ddd", padding="6px"))
    return card, number


settings_mode = widgets.ToggleButtons(
    options=[("Keep recommended defaults", "basic"), ("Advanced settings", "advanced")],
    value="basic", description="Settings",
)
idt_setup_mode_widget = widgets.ToggleButtons(
    options=(
        ("Create Credentials", "create"),
        ("Upload idt.env", "upload"),
        ("Do not use IDT API — export batch input", "batch"),
    ),
    value="create",
    description="IDT setup",
    layout=widgets.Layout(width="100%"),
)
idt_auth_method_widget = widgets.ToggleButtons(
    options=(("Client credentials", "password"), ("Access token", "access_token")),
    value="password", description="Authentication",
)
idt_client_id_widget = widgets.Password(description="Client ID", layout=widgets.Layout(width="49%"))
idt_client_secret_widget = widgets.Password(description="Client secret", layout=widgets.Layout(width="49%"))
idt_username_widget = widgets.Password(description="IDT username", layout=widgets.Layout(width="49%"))
idt_password_widget = widgets.Password(description="IDT password", layout=widgets.Layout(width="49%"))
idt_access_token_widget = widgets.Password(description="Access token", layout=widgets.Layout(width="98%"))
credential_upload = widgets.FileUpload(
    accept=".env,text/plain", multiple=False, description="Upload idt.env"
)
credential_test_button = widgets.Button(description="Test credentials", icon="check", button_style="info")
credential_upload_test_button = widgets.Button(description="Test uploaded credentials", icon="check", button_style="info")
credential_download_button = widgets.Button(description="Download idt.env", icon="download")
credential_status = widgets.HTML("<b>IDT status:</b> credentials not configured")
output_directory_widget = widgets.Text(value="/content/hurdler_runs/current", description="Runtime work folder", layout=widgets.Layout(width="98%"))
auto_download_widget = widgets.Checkbox(value=True, description="Auto-download ZIP after success")
verbose_generations = widgets.Checkbox(value=False, description="Show every GA generation in Advanced log")

external_worker_cpus = widgets.BoundedIntText(value=16, min=1, max=1024, description="GA worker CPUs")
external_memory_gb = widgets.BoundedIntText(value=32, min=1, max=1_048_576, description="Total memory (GB)")
external_walltime = widgets.Text(value="24:00:00", description="Walltime")
external_partition = widgets.Text(value="cpu", description="Partition")
external_account = widgets.Text(value="", description="Account")
external_qos = widgets.Text(value="", description="QoS")
external_constraint = widgets.Text(value="", description="Constraint")
external_conda_environment = widgets.Text(value="hurdler", description="Conda env")
external_result_directory = widgets.Text(value="results", description="Results folder")
EXTERNAL_IDT_CREDENTIAL_PATH = "~/.config/hurdler/idt.env"

population_card, population_number = _numeric_control("Population", 16, 4, 256, 4, integer=True, help_text="Candidates per GA generation; larger values improve exploration but cost time.")
mutation_card, mutation_number = _numeric_control("Mutation rate", 0.08, 0.001, 0.5, 0.001, help_text="Probability of synonymous codon mutation; higher values explore more aggressively.")
crossover_card, crossover_number = _numeric_control("Crossover rate", 0.75, 0.0, 1.0, 0.01, help_text="Probability of recombining parent DNA candidates while preserving translation.")
elite_card, elite_number = _numeric_control("Elite fraction", 0.15, 0.01, 0.5, 0.01, help_text="Best fraction retained each generation; high values reduce diversity.")
secondary_copy_range_widget = widgets.IntRangeSlider(
    value=(12, 20), min=1, max=50, step=1, description="Copy range",
    continuous_update=False, readout=True,
    layout=widgets.Layout(width="98%"),
)
minimum_secondary_number = widgets.BoundedIntText(
    value=12, min=1, max=50, step=1, description="Minimum",
    layout=widgets.Layout(width="48%"),
)
maximum_secondary_number = widgets.BoundedIntText(
    value=20, min=1, max=50, step=1, description="Maximum",
    layout=widgets.Layout(width="48%"),
)
secondary_range_card = widgets.VBox([
    widgets.HTML(
        "<b>Bounded secondary-copy range</b> <span style='color:#666'>[repeat copies]</span><br>"
        "<small>Default 12–20; allowed 1–50. Drag either handle on the shared axis or type both bounds manually. "
        "Automatic mode ignores this range and explores to the physical/route limit.</small>"
    ),
    secondary_copy_range_widget,
    widgets.HBox([minimum_secondary_number, maximum_secondary_number]),
], layout=widgets.Layout(border="1px solid #ddd", padding="7px", margin="3px 0"))
feedback_round_card, feedback_round_number = _numeric_control(
    "Maximum GA→IDT feedback rounds", 100, 1, 1000, 1, integer=True,
    help_text="Retries per copy count. Positive IDT rules adjust weights before the next retry."
)
generations_per_round_card, generations_per_round_number = _numeric_control(
    "GA generations per feedback round", 10, 1, 1000, 1, integer=True,
    help_text="Generations added before each exact-DNA IDT score request."
)
elite_seed_card, elite_seed_number = _numeric_control(
    "Warm-start top candidates", 10, 1, 256, 1, integer=True,
    help_text="Number of prior elite sequences carried into the next feedback round."
)
max_population_card, max_population_number = _numeric_control(
    "Adaptive population cap", 256, 4, 2048, 4, integer=True,
    help_text="Upper bound when IDT rejection increases the GA population."
)
max_mutation_card, max_mutation_number = _numeric_control(
    "Adaptive mutation cap", 0.35, 0.001, 1.0, 0.001,
    help_text="Maximum mutation rate reached during IDT feedback."
)
max_crossover_card, max_crossover_number = _numeric_control(
    "Adaptive crossover cap", 0.95, 0.0, 1.0, 0.01,
    help_text="Maximum crossover rate reached during IDT feedback."
)
secondary_search_mode_widget = widgets.ToggleButtons(
    options=(("Automatic to limit", "automatic"), ("Bounded copy range", "bounded")),
    value="bounded", description="Secondary search",
)
secondary_length_status = widgets.HTML()
seed_number = widgets.IntText(value=42, description="Random seed")
generation_schedule_widget = widgets.Text(value="10,20,40,60,80,100", description="Generations")
auto_weight_feedback = widgets.Checkbox(value=True, description="Adjust weights from IDT positive rules")
auto_parameter_feedback = widgets.Checkbox(value=True, description="Adapt population / mutation / crossover from IDT score")

_secondary_range_sync = {"active": False}


def _range_slider_changed(change):
    if _secondary_range_sync["active"]:
        return
    _secondary_range_sync["active"] = True
    try:
        minimum_secondary_number.value, maximum_secondary_number.value = map(int, change["new"])
    finally:
        _secondary_range_sync["active"] = False
    _update_secondary_lengths()


def _range_number_changed(change):
    if _secondary_range_sync["active"]:
        return
    minimum = int(minimum_secondary_number.value)
    maximum = int(maximum_secondary_number.value)
    _secondary_range_sync["active"] = True
    try:
        if minimum > maximum:
            if change.get("owner") is minimum_secondary_number:
                maximum = minimum
                maximum_secondary_number.value = maximum
            else:
                minimum = maximum
                minimum_secondary_number.value = minimum
        secondary_copy_range_widget.value = (minimum, maximum)
    finally:
        _secondary_range_sync["active"] = False
    _update_secondary_lengths()


def _secondary_bounds():
    if secondary_search_mode_widget.value == "automatic":
        return 1, None
    minimum = int(minimum_secondary_number.value)
    maximum = int(maximum_secondary_number.value)
    if maximum < minimum:
        raise ValueError("Maximum secondary copies cannot be smaller than minimum secondary copies")
    return minimum, maximum


def _update_secondary_lengths(_change=None):
    bounded = secondary_search_mode_widget.value == "bounded"
    secondary_copy_range_widget.disabled = not bounded
    minimum_secondary_number.disabled = not bounded
    maximum_secondary_number.disabled = not bounded
    module_bp = len("".join(str(repeat_module_widget.value).split())) * 3
    try:
        minimum, maximum = _secondary_bounds()
    except ValueError as exc:
        secondary_length_status.value = f"<span style='color:#b31b1b'><b>{exc}</b></span>"
        return
    adapter_bp = 0
    route = state.get("confirmed_route")
    if route is not None and state.get("confirmed_site_iii"):
        try:
            selected = {**route, "site_iii_enzyme": state["confirmed_site_iii"]}
            left_adapter, right_adapter, _evidence = _secondary_adapters(selected, project_root=project_root)
            adapter_bp = len(left_adapter) + len(right_adapter)
        except Exception:
            adapter_bp = 0
    if maximum is None:
        physical = max(0, (3000 - adapter_bp) // max(1, module_bp))
        text = f"Automatic: starts at 1 copy ({module_bp:,} core bp) and explores up to about {physical:,} copies under the 3,000-bp purchase cap"
    else:
        text = (
            f"Bounded: {minimum}–{maximum} copies; core {minimum * module_bp:,}–{maximum * module_bp:,} bp; "
            f"purchase {minimum * module_bp + adapter_bp:,}–{maximum * module_bp + adapter_bp:,} bp including {adapter_bp} adapter bp"
        )
    secondary_length_status.value = f"<b>{text}.</b>"

weight_defaults = {
    "selected_re_site_excess": 1_000_000_000.0,
    "gc_window_violation": 1_000_000_000.0,
    "repeated_re_site_excess": 10_000.0,
    "repeated_14mer": 250.0,
    "repeated_13mer": 100.0,
    "repeated_8mer": 5.0,
    "hairpin_10mer_proxy": 25.0,
    "homopolymer_excess": 250.0,
    "terminal_repeat_proxy": 100.0,
    "gc_window_soft_violation": 100.0,
    "negative_log_cai": 50.0,
}
weight_widgets = {
    name: widgets.FloatText(value=value, description=name, layout=widgets.Layout(width="49%"))
    for name, value in weight_defaults.items()
}


def _two_per_row(items):
    rows = []
    for index in range(0, len(items), 2):
        rows.append(widgets.HBox(items[index:index + 2], layout=widgets.Layout(width="100%")))
    return widgets.VBox(rows)


advanced_panel = widgets.VBox([
    widgets.HTML(
        "<p><b>Advanced GA settings.</b> All sequence changes are synonymous. "
        "The selected Site-I/Site-II excess count is a hard constraint; repeated non-selected RE sites and IDT rule feedback are weighted score terms.</p>"
    ),
    _two_per_row([population_card, mutation_card, crossover_card, elite_card]),
    _two_per_row([
        generations_per_round_card, elite_seed_card,
        max_population_card, max_mutation_card, max_crossover_card,
    ]),
    _help_card("Random seed", seed_number, unit="integer", default="42", purpose="Makes stochastic GA choices reproducible.", allowed="any integer", effect="Changing it explores a different synonymous-DNA trajectory."),
    _help_card("Generation schedule", generation_schedule_widget, unit="generations", default="10,20,40,60,80,100", purpose="Legacy schedule used by applicable optimization routes.", allowed="positive comma-separated integers ending at 100", effect="Longer schedules cost more CPU."),
    widgets.HBox([auto_weight_feedback, auto_parameter_feedback]),
    verbose_generations,
    widgets.HTML("<b>GA score weights</b><br><small>Each value multiplies its named violation. Larger values make that defect less tolerable; IDT positive-rule feedback may raise mapped weights but never changes the protein.</small>"),
    _two_per_row(list(weight_widgets.values())),
])

external_resource_panel = widgets.Accordion(children=[widgets.VBox([
    widgets.HTML(
        "<p><b>Portable Local + Slurm run.</b> These resources apply to the exported bundle. "
        "Fitness calculations use the requested worker processes; random operations and IDT calls stay serial. "
        "Account, QoS and constraint are optional and accept scheduler-safe names only.</p>"
    ),
    _two_per_row([
        _help_card("GA worker CPUs", external_worker_cpus, unit="CPUs", default="16", purpose="Parallel fitness workers and Slurm cpus-per-task.", allowed="integer 1–1,024", effect="The external preflight fails if fewer CPUs are available."),
        _help_card("Total memory", external_memory_gb, unit="GB", default="32", purpose="Total Slurm memory request.", allowed="positive integer", effect="Local mode reports the request; Slurm enforces it."),
        _help_card("Walltime", external_walltime, unit="HH:MM:SS", default="24:00:00", purpose="Slurm job limit.", allowed="HH:MM:SS or D-HH:MM:SS", effect="The scheduler stops jobs that exceed it."),
        _help_card("Partition", external_partition, unit="Slurm name", default="cpu", purpose="CPU partition used by sbatch.", allowed="letters, digits, dot, underscore, hyphen", effect="Must exist on the target cluster."),
    ]),
    _two_per_row([
        _help_card("Account", external_account, unit="optional", default="empty", purpose="Optional Slurm billing account.", allowed="safe scheduler name", effect="Adds a structured --account directive."),
        _help_card("QoS", external_qos, unit="optional", default="empty", purpose="Optional Slurm quality of service.", allowed="safe scheduler name", effect="Adds a structured --qos directive."),
        _help_card("Constraint", external_constraint, unit="optional", default="empty", purpose="Optional node feature constraint.", allowed="safe scheduler name", effect="Adds a structured --constraint directive."),
        _help_card("Conda environment", external_conda_environment, unit="name", default="hurdler", purpose="Environment created/activated by run_ga.sh.", allowed="safe environment name", effect="setup uses the bundled YAML."),
    ]),
    _help_card("External results folder", external_result_directory, unit="path", default="results", purpose="Stores progress, checkpoint and final archives.", allowed="relative bundle path or absolute shared path", effect="Compute nodes must be able to write it."),
    widgets.HTML(
        "<div style='border-left:5px solid #b7a57a;background:#fffaf0;padding:10px'>"
        "<b>External IDT credentials:</b> Live-API bundles always read "
        "<code>~/.config/hurdler/idt.env</code> on the target machine and auto-detect its format. "
        "Colab credentials are never copied into the bundle. Batch bundles omit all IDT arguments.</div>"
    ),
])])
external_resource_panel.set_title(0, "Resource request")
external_resource_panel.selected_index = None


def _sync_settings(_change=None):
    advanced_panel.layout.display = "" if settings_mode.value == "advanced" else "none"


settings_mode.observe(_sync_settings, names="value")
_sync_settings()

credential_registration_help = widgets.HTML(
    "<div style='border:2px solid #4b2e83;border-left-width:6px;background:#ffffff;"
    "color:#111827;padding:14px;border-radius:8px;line-height:1.45'>"
    "<b style='color:#3b1f69;font-size:16px'>Create an IDT SciTools API client</b><ol style='color:#111827'>"
    "<li><a href='https://www.idtdna.com/page/tools/scitools-plus-api-overview' target='_blank'>"
    "Sign in or create an IDT account</a>.</li>"
    "<li>Open <b>My Account → API access → Request new API key</b>.</li>"
    "<li>Choose a unique Client ID, accept the API terms, and securely copy the generated Client secret.</li>"
    "<li>Enter the four password-grant fields below, or choose Access token.</li></ol>"
    "<b style='color:#111827'>Password-grant file</b><pre style='background:#f9fafb;color:#111827;"
    "border:1px solid #d1d5db;padding:8px'>IDT_CLIENT_ID=your_client_id\nIDT_CLIENT_SECRET=your_client_secret\n"
    "IDT_USERNAME=your_idt_username\nIDT_PASSWORD=your_idt_password</pre>"
    "<b style='color:#111827'>Access-token file</b><pre style='background:#f9fafb;color:#111827;"
    "border:1px solid #d1d5db;padding:8px'>IDT_ACCESS_TOKEN=your_current_access_token</pre></div>"
)
credential_security_notice = widgets.HTML(
    "<div style='border:2px solid #2d6a4f;background:#effaf4;color:#111827;border-radius:8px;padding:10px'>"
    "<b>Credential handling:</b> this notebook does not write secrets to notebook output, the repository, "
    "logs, Drive, checkpoints, or result bundles. Values remain only in this Colab kernel and are sent only "
    "to IDT OAuth/API endpoints. Colab is still a third-party runtime; use a temporary access token if that is preferred."
    "</div>"
)
idt_password_fields_panel = widgets.VBox([
    widgets.HBox([idt_client_id_widget, idt_client_secret_widget]),
    widgets.HBox([idt_username_widget, idt_password_widget]),
])
idt_token_field_panel = widgets.VBox([idt_access_token_widget])
credential_create_panel = widgets.VBox([
    credential_registration_help,
    idt_auth_method_widget,
    idt_password_fields_panel,
    idt_token_field_panel,
    widgets.HBox([credential_test_button, credential_download_button]),
])
credential_upload_panel = widgets.VBox([
    widgets.HTML(
        "<b>Upload one UTF-8 <code>idt.env</code>.</b> It must contain exactly one of the two formats shown above. "
        "The uploaded bytes are parsed in memory and are never included in any output archive."
    ),
    credential_registration_help,
    widgets.HBox([credential_upload, credential_upload_test_button]),
])
back_to_ga_button = widgets.Button(description="Back to GA settings", icon="arrow-up")
credential_batch_panel = widgets.VBox([
    widgets.HTML(
        "<div style='border-left:5px solid #b7a57a;background:#fffaf0;padding:12px'>"
        "<b>No live API calls will be made.</b> GA exports IDT Bulk Input CSV, TSV and FASTA plus elite candidates. "
        "This mode is unvalidated and never claims IDT acceptance. If IDT finds no orderable candidate, return here, "
        "adjust GA settings, and run again.</div>"
    ),
    back_to_ga_button,
])
idt_credential_panel = widgets.VBox([
    _help_card(
        "IDT scoring setup", idt_setup_mode_widget, unit="mode", default="Create Credentials",
        purpose="Chooses live complexity scoring from an in-memory credential or an offline Bulk Input export.",
        allowed="create, upload, or no API", effect="Only the two live modes may report IDT score-sum <10 acceptance.",
    ),
    credential_security_notice,
    credential_create_panel,
    credential_upload_panel,
    credential_batch_panel,
    credential_status,
])

stage_html = widgets.HTML(
    "<div style='border:2px solid #4b2e83;background:#f4f0fa;color:#111827;border-radius:8px;padding:10px'>"
    "<b>Status:</b> waiting for route confirmation</div>"
)
generation_progress = widgets.IntProgress(value=0, min=0, max=1, description="GA")
candidate_progress = widgets.IntProgress(
    value=0, min=0, max=1, description="Fitness", bar_style="info"
)
current_html = widgets.HTML("")
attempt_log_html = widgets.HTML(
    "<pre style='color:#111827;background:#ffffff'>Waiting for a confirmed route and GA start.</pre>",
    layout=widgets.Layout(height="230px", overflow="auto", border="1px solid #ccc", padding="8px"),
)
idt_plot_status = widgets.HTML(
    "<div style='border:2px dashed #4b2e83;background:#ffffff;color:#111827;"
    "border-radius:8px;padding:12px'>No IDT scores yet. The trajectory updates after each scored purchase fragment.</div>"
)
idt_plot_output = widgets.Output()
idt_score_table_output = widgets.Output()
design_output = widgets.Output()
results_status = widgets.HTML(
    "<div style='border:2px dashed #4b2e83;background:#ffffff;color:#111827;"
    "border-radius:8px;padding:12px'>No optimization result yet. Confirm a route, then run in Colab or export the external bundle.</div>"
)
execution_target_widget = widgets.ToggleButtons(
    options=(("Run in Colab", "colab"), ("Local / Slurm bundle", "external")),
    value="colab",
    description="Execution target",
    layout=widgets.Layout(width="100%"),
)
design_button = widgets.Button(
    description="Run GA in Colab", icon="play", disabled=True,
    layout=widgets.Layout(width="230px", height="44px"),
)
design_button.style.button_color = "#4b2e83"
design_button.style.font_weight = "bold"
pause_button = widgets.Button(description="Pause GA", icon="pause", disabled=True)
pause_button.style.button_color = "#b7a57a"
stop_button = widgets.Button(description="Stop GA", icon="stop", disabled=True, button_style="danger")
export_bundle_button = widgets.Button(
    description="Download Local / Slurm bundle", button_style="info", icon="archive", disabled=True
)
download_button = widgets.Button(description="Download results ZIP", icon="download", disabled=True)
external_bundle_output = widgets.Output()


def _generation_schedule():
    values = tuple(sorted({
        int(value.strip()) for value in generation_schedule_widget.value.split(",") if value.strip()
    }))
    if not values or values[-1] != 100 or any(value <= 0 for value in values):
        raise ValueError("Generation schedule must contain positive integers and terminate at 100")
    return values


def _set_idt_plot_placeholder(message):
    idt_plot_status.value = (
        "<div style='border:2px dashed #4b2e83;background:#ffffff;color:#111827;"
        f"border-radius:8px;padding:12px'>{html.escape(str(message))}</div>"
    )
    with idt_plot_output:
        clear_output(wait=True)
    with idt_score_table_output:
        clear_output(wait=True)


def _render_idt_trajectory():
    rows = idt_score_history_rows(state["idt_score_events"])
    if not rows:
        _set_idt_plot_placeholder(
            "No IDT scores yet. The trajectory updates after each scored purchase fragment."
        )
        return
    idt_plot_status.value = (
        "<div style='border:2px solid #2d6a4f;background:#effaf4;color:#111827;"
        f"border-radius:8px;padding:10px'><b>{len(rows)} IDT evaluations received.</b> "
        "Green points pass, red points fail, and grey points are unclassified.</div>"
    )
    with idt_plot_output:
        clear_output(wait=True)
        figure = plot_idt_score_trajectory(rows, title="Live IDT complexity score trajectory")
        display(figure)
        plt.close(figure)
    with idt_score_table_output:
        clear_output(wait=True)
        columns = [
            "evaluation_index", "fragment_id", "fragment_kind", "repeat_copies",
            "feedback_round", "idt_total_score", "idt_classification",
            "idt_cache_hit", "positive_rule_names_json",
        ]
        display(pd.DataFrame(rows)[columns].tail(10))


def _should_log_progress(event):
    if event.stage == "ga" and event.status == "fitness_running":
        index = int(event.candidate_index or 0)
        total = int(event.candidate_total or 0)
        return bool(index == 1 or index == total or index % 4 == 0)
    if event.stage == "ga" and event.status == "running":
        generation = int(event.generation or 0)
        final_generation = int(event.generations or 0)
        return bool(
            verbose_generations.value
            or generation == 1
            or generation == final_generation
            or generation % 5 == 0
        )
    return event.status in {
        "started", "worker_entered", "preparing", "population_initializing", "fitness_started",
        "baseline_fitness_started", "baseline_fitness_completed",
        "attempt_started", "attempt_completed", "request_started",
        "fragment_scored", "request_completed", "completed", "failed",
        "parameters_adjusted", "no_novel_candidate",
    }


def _progress_line(event):
    return (
        f"{event.stage:<12} {event.status:<18} {event.fragment_kind or '-':<12} "
        f"copies={event.copies if event.copies is not None else '-'} "
        f"feedback={event.feedback_round if event.feedback_round is not None else '-'}/"
        f"{event.max_feedback_rounds if event.max_feedback_rounds is not None else '-'} "
        f"gen={event.generation if event.generation is not None else '-'}/"
        f"{event.generations if event.generations is not None else '-'} "
        f"candidate={event.candidate_index if event.candidate_index is not None else '-'}/"
        f"{event.candidate_total if event.candidate_total is not None else '-'} "
        f"ga={event.ga_score if event.ga_score is not None else '-'} "
        f"idt={event.idt_score if event.idt_score is not None else '-'}"
    )


def _render_progress(event: DesignProgressEvent):
    stage_html.value = (
        "<div style='border:2px solid #4b2e83;background:#f4f0fa;color:#111827;border-radius:8px;padding:10px'>"
        f"<b>Status:</b> {html.escape(event.stage)} · {html.escape(event.status)}</div>"
    )
    if event.generations:
        generation_progress.max = max(1, int(event.generations))
        generation_progress.value = min(generation_progress.max, int(event.generation or 0))
    if event.candidate_total:
        candidate_progress.max = max(1, int(event.candidate_total))
        candidate_progress.value = min(
            candidate_progress.max, int(event.candidate_index or 0)
        )
    current_html.value = (
        f"<b>{event.fragment_kind or 'design'}</b> · copies={event.copies if event.copies is not None else '—'} "
        f"· feedback={event.feedback_round if event.feedback_round is not None else '—'}/"
        f"{event.max_feedback_rounds if event.max_feedback_rounds is not None else '—'} "
        f"· generation={event.generation if event.generation is not None else '—'}/"
        f"{event.generations if event.generations is not None else '—'} "
        f"· fitness candidate={event.candidate_index if event.candidate_index is not None else '—'}/"
        f"{event.candidate_total if event.candidate_total is not None else '—'} "
        f"· best score={event.ga_score if event.ga_score is not None else '—'} "
        f"· IDT={event.idt_score if event.idt_score is not None else '—'} "
        f"· pop/mut/xover={event.population_size or '—'}/"
        f"{event.mutation_rate if event.mutation_rate is not None else '—'}/"
        f"{event.crossover_rate if event.crossover_rate is not None else '—'} "
        f"· elapsed={event.elapsed_seconds or 0:.1f}s"
    )
    if _should_log_progress(event):
        state["visible_log_lines"].append(_progress_line(event))
        attempt_log_html.value = (
            "<pre style='color:#111827;background:#ffffff'>"
            + html.escape("\n".join(state["visible_log_lines"][-18:]))
            + "</pre>"
        )
    if event.stage == "idt" and event.status == "fragment_scored":
        state["idt_score_events"].append(event.to_dict())
        _render_idt_trajectory()


def _progress_update(event: DesignProgressEvent, run_id=None):
    active_run_id = int(state.get("run_id", 0))
    run_id = active_run_id if run_id is None else int(run_id)
    if run_id != active_run_id:
        return
    with state["progress_lock"]:
        state["progress_events"].append(event.to_dict())
        state["last_progress_monotonic"] = time.monotonic()
    _enqueue_ui_event("progress", run_id, event)


def _drain_ui_events(max_items=100):
    handled = 0
    while handled < int(max_items):
        try:
            kind, run_id, payload = state["progress_queue"].get_nowait()
        except queue.Empty:
            break
        handled += 1
        if int(run_id) != int(state.get("run_id", 0)):
            continue
        if kind == "progress":
            _render_progress(payload)
        elif kind == "success":
            _finish_design_success(*payload)
        elif kind == "stopped":
            _finish_design_stopped()
        elif kind == "error":
            _finish_design_error(*payload)
        elif kind == "stage_html":
            stage_html.value = str(payload)
    return handled


def _capture_ui_dispatcher():
    """Capture Colab/ipykernel's UI loop before the GA worker starts."""
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None
    shell = get_ipython()
    io_loop = getattr(getattr(shell, "kernel", None), "io_loop", None)
    state["ui_asyncio_loop"] = loop
    state["ui_io_loop"] = io_loop
    return loop, io_loop


def _schedule_ui_callback(callback, *args):
    """Post one callback to the kernel UI loop from any worker thread."""
    loop = state.get("ui_asyncio_loop")
    if loop is not None and loop.is_running():
        loop.call_soon_threadsafe(callback, *args)
        return True
    io_loop = state.get("ui_io_loop")
    if io_loop is not None:
        io_loop.add_callback(callback, *args)
        return True
    return False


def _request_ui_drain(run_id):
    """Coalesce worker events while guaranteeing a UI-thread drain request."""
    with state["ui_schedule_lock"]:
        if state.get("ui_drain_pending"):
            return True
        state["ui_drain_pending"] = True

    def _scheduled_drain():
        try:
            if int(run_id) == int(state.get("run_id", 0)):
                _drain_ui_events()
        finally:
            with state["ui_schedule_lock"]:
                state["ui_drain_pending"] = False
        if not state["progress_queue"].empty():
            _request_ui_drain(int(state.get("run_id", 0)))

    scheduled = _schedule_ui_callback(_scheduled_drain)
    if not scheduled:
        with state["ui_schedule_lock"]:
            state["ui_drain_pending"] = False
    return scheduled


def _enqueue_ui_event(kind, run_id, payload):
    state["progress_queue"].put((str(kind), int(run_id), payload))
    return _request_ui_drain(int(run_id))


async def _ui_event_pump(run_id):
    try:
        while int(run_id) == int(state.get("run_id", 0)):
            _drain_ui_events()
            if not state.get("run_active") and state["progress_queue"].empty():
                break
            last = state.get("last_progress_monotonic") or state.get("run_started_monotonic")
            if state.get("run_active") and last is not None:
                idle = time.monotonic() - float(last)
                elapsed = time.monotonic() - float(state["run_started_monotonic"])
                if idle >= 2.0:
                    stage_html.value = (
                        "<div style='border:2px solid #4b2e83;background:#f4f0fa;color:#111827;"
                        "border-radius:8px;padding:10px'><b>Status:</b> GA worker active inside the current "
                        f"fitness/API unit · elapsed {elapsed:.1f}s · last event {idle:.1f}s ago</div>"
                    )
                if (
                    time.monotonic() - float(state.get("last_checkpoint_write") or 0.0)
                    >= 180.0
                ):
                    _start_periodic_checkpoint()
            worker = state.get("run_thread")
            if state.get("run_active") and worker is not None and not worker.is_alive() and state["progress_queue"].empty():
                _finish_design_error("WorkerExited", "GA worker ended without a terminal event")
                break
            await asyncio.sleep(0.2)
    except asyncio.CancelledError:
        return
    except Exception as exc:
        if int(run_id) == int(state.get("run_id", 0)) and state.get("run_active"):
            _finish_design_error(type(exc).__name__, f"UI event pump failed: {exc}")


def _start_ui_event_pump(run_id):
    previous = state.get("ui_pump_task")
    if previous is not None and not previous.done():
        previous.cancel()
    loop, io_loop = _capture_ui_dispatcher()
    if loop is None:
        # Some ipykernel/Colab releases execute widget callbacks outside the
        # asyncio task while their Tornado IOLoop is still the authoritative
        # UI thread. Schedule creation there instead of letting the worker
        # thread touch widget state.
        if io_loop is None:
            state["ui_pump_task"] = None
            return None

        def _schedule_on_kernel_loop():
            if int(run_id) != int(state.get("run_id", 0)):
                return
            state["ui_pump_task"] = asyncio.ensure_future(_ui_event_pump(int(run_id)))

        io_loop.add_callback(_schedule_on_kernel_loop)
        return "scheduled_on_kernel_loop"
    task = loop.create_task(_ui_event_pump(int(run_id)))
    state["ui_pump_task"] = task
    return task


def _checkpoint_local_path():
    run_directory = Path(state.get("run_directory") or output_directory_widget.value)
    root = _runtime_scratch_path("hurdler_checkpoints", run_directory.parent / "checkpoints")
    root.mkdir(parents=True, exist_ok=True)
    safe_id = "".join(character if character.isalnum() or character in "._-" else "_" for character in str(sequence_id_widget.value or "interactive_design"))
    return root / f"hurdler_{safe_id}_checkpoint_latest.zip"


def _start_periodic_checkpoint():
    """Write the 180-second heartbeat checkpoint without blocking the UI/GA."""
    existing = state.get("checkpoint_thread")
    if existing is not None and existing.is_alive():
        return existing
    # Reserve the interval before starting so a slow Drive mount cannot spawn
    # duplicate writers while the previous copy is still in progress.
    state["last_checkpoint_write"] = time.monotonic()

    def _write():
        try:
            _persist_checkpoint(force=True)
        except Exception as exc:
            _enqueue_ui_event(
                "stage_html", int(state.get("run_id", 0)),
                "<div style='border:2px solid #b7a57a;background:#fffaf0;color:#111827;"
                f"border-radius:8px;padding:10px'><b>Checkpoint warning:</b> {html.escape(str(exc)[:300])}</div>",
            )

    worker = threading.Thread(
        target=_write, name="hurdler-checkpoint-writer", daemon=True
    )
    state["checkpoint_thread"] = worker
    worker.start()
    return worker


def _copy_archive_to_drive(path):
    if storage_mode_widget.value != "drive":
        return None
    if not storage_state.get("drive_mounted"):
        raise RuntimeError("Google Drive was selected but is not mounted; click Mount Google Drive in Step 0")
    destination = Path(drive_root_widget.value)
    destination.mkdir(parents=True, exist_ok=True)
    copied = destination / Path(path).name
    shutil.copy2(path, copied)
    return copied


def _persist_checkpoint_unlocked(payload=None, *, force=False):
    now = time.monotonic()
    if payload is not None:
        score = payload.get("idt_complexity_score")
        accepted = (
            payload.get("event") == "accepted_secondary"
            and payload.get("validation_mode") == "api"
            and isinstance(score, (int, float)) and not isinstance(score, bool)
            and float(score) < 10
        )
        previous = state.get("best_checkpoint")
        if accepted and (previous is None or int(payload["repeat_copies"]) > int(previous["repeat_copies"])):
            state["best_checkpoint"] = dict(payload)
            force = True
    if not force and now - float(state.get("last_checkpoint_write") or 0.0) < 180:
        return None
    best = state.get("best_checkpoint")
    public_payload = best or {
        "event": "heartbeat",
        "sequence_id": str(sequence_id_widget.value or "interactive_design"),
        "accepted_secondary_available": False,
        "status": stage_html.value,
        "query_fingerprint": state.get("confirmed_fingerprint"),
        "route_fingerprint": hashlib.sha256(
            json.dumps(state.get("confirmed_route") or {}, sort_keys=True).encode()
        ).hexdigest(),
        "tested_lengths": sorted({
            row.get("copies") for row in state.get("progress_events", []) if row.get("copies") is not None
        }),
        "failure_reason": "No live-IDT-accepted secondary has been obtained yet",
    }
    public_payload = dict(public_payload)
    public_payload["run_status"] = state.get("run_terminal_status") or "running"
    checkpoint = write_secondary_checkpoint(public_payload, _checkpoint_local_path())
    state["checkpoint_archive"] = checkpoint
    state["last_checkpoint_write"] = now
    if storage_mode_widget.value == "drive" and storage_state.get("drive_mounted"):
        _copy_archive_to_drive(checkpoint)
    return checkpoint


def _persist_checkpoint(payload=None, *, force=False):
    with state["checkpoint_lock"]:
        return _persist_checkpoint_unlocked(payload, force=force)


def _checkpoint_update(payload):
    _persist_checkpoint(payload, force=False)
    checkpoint = state.get("best_checkpoint")
    if checkpoint:
        message = (
            "<div style='border:2px solid #4b2e83;background:#f4f0fa;border-radius:8px;padding:10px'>"
            f"<b>Checkpoint saved:</b> {checkpoint['repeat_copies']} secondary copies · "
            f"IDT {checkpoint['idt_complexity_score']} · {Path(state['checkpoint_archive']).name}</div>"
        )
        _enqueue_ui_event(
            "stage_html", int(state.get("run_id", 0)), message
        )


def _validation_mode():
    return "batch" if idt_setup_mode_widget.value == "batch" else "api"


def _wipe_bytearray(payload):
    if isinstance(payload, bytearray):
        for index in range(len(payload)):
            payload[index] = 0


def _wipe_cached_credentials():
    _wipe_bytearray(state.get("credential_payload"))
    state["credential_payload"] = None
    state["credential_auth_method"] = None
    clear_idt_secret_environment()


def _clear_create_secret_widgets():
    idt_client_secret_widget.value = ""
    idt_password_widget.value = ""
    idt_access_token_widget.value = ""


def _create_credential_values():
    if idt_auth_method_widget.value == "access_token":
        return {"IDT_ACCESS_TOKEN": str(idt_access_token_widget.value).strip()}
    return {
        "IDT_CLIENT_ID": str(idt_client_id_widget.value).strip(),
        "IDT_CLIENT_SECRET": str(idt_client_secret_widget.value).strip(),
        "IDT_USERNAME": str(idt_username_widget.value).strip(),
        "IDT_PASSWORD": str(idt_password_widget.value).strip(),
    }


def _dotenv_payload(values):
    ordered = (
        ("IDT_ACCESS_TOKEN",)
        if "IDT_ACCESS_TOKEN" in values
        else ("IDT_CLIENT_ID", "IDT_CLIENT_SECRET", "IDT_USERNAME", "IDT_PASSWORD")
    )
    return ("\n".join(f"{name}={values.get(name, '')}" for name in ordered) + "\n").encode()


def _uploaded_payload():
    uploaded = credential_upload.value
    if not uploaded:
        return None
    item = next(iter(uploaded.values())) if isinstance(uploaded, dict) else uploaded[0]
    return bytes(item["content"] if isinstance(item, dict) else item.content)


def _clear_credential_upload():
    """Drop uploaded bytes without assigning Colab's read-only value trait."""
    global credential_upload
    previous = credential_upload
    credential_upload = widgets.FileUpload(
        accept=".env,text/plain", multiple=False,
        description="Upload idt.env",
    )
    credential_upload_panel.children = (
        credential_upload_panel.children[0],
        credential_upload_panel.children[1],
        widgets.HBox([credential_upload, credential_upload_test_button]),
    )
    try:
        previous.close()
    except Exception:
        # Replacement already removed the only live UI reference.  Some
        # hosted widget backends do not implement close() completely.
        pass


def _configure_api_credentials():
    if _validation_mode() != "api":
        raise RuntimeError("IDT credentials are not used in Bulk Input mode")
    cached = state.get("credential_payload")
    if isinstance(cached, bytearray) and cached:
        return configure_idt_credentials_from_bytes(bytes(cached))
    if idt_setup_mode_widget.value == "create":
        values = _create_credential_values()
        try:
            status = configure_idt_credentials_from_values(
                values, auth_method=str(idt_auth_method_widget.value)
            )
            payload = _dotenv_payload(values)
            state["credential_payload"] = bytearray(payload)
            state["credential_auth_method"] = status["auth_method"]
            return status
        finally:
            values.clear()
            _clear_create_secret_widgets()
    payload = _uploaded_payload()
    if payload is None:
        raise FileNotFoundError("Upload one idt.env file before testing or running Live IDT scoring")
    try:
        status = configure_idt_credentials_from_bytes(payload)
        state["credential_payload"] = bytearray(payload)
        state["credential_auth_method"] = status["auth_method"]
        return status
    finally:
        payload = b""
        _clear_credential_upload()


def _sync_idt_auth_fields(_change=None):
    password = idt_auth_method_widget.value == "password"
    idt_password_fields_panel.layout.display = "" if password else "none"
    idt_token_field_panel.layout.display = "none" if password else ""


def _sync_idt_setup(change=None):
    if change is not None:
        _wipe_cached_credentials()
        _clear_create_secret_widgets()
    mode = idt_setup_mode_widget.value
    credential_create_panel.layout.display = "" if mode == "create" else "none"
    credential_upload_panel.layout.display = "" if mode == "upload" else "none"
    credential_batch_panel.layout.display = "" if mode == "batch" else "none"
    credential_status.value = (
        "<b>IDT status:</b> offline Bulk Input mode; no credentials or API calls"
        if mode == "batch"
        else "<b>IDT status:</b> credentials remain only in this kernel and have not been tested"
    )
    if "idt_plot_status" in globals() and not state.get("run_active"):
        if mode == "batch":
            _set_idt_plot_placeholder(
                "IDT API disabled. GA will export Bulk Input files; no score trajectory is claimed."
            )
        else:
            _set_idt_plot_placeholder(
                "No IDT scores yet. The trajectory updates after each scored purchase fragment."
            )
    if "_invalidate_external_bundle" in globals():
        _invalidate_external_bundle()


def _test_idt_credentials(_button=None):
    credential_test_button.disabled = True
    credential_upload_test_button.disabled = True
    credential_status.value = "<b>IDT status:</b> testing OAuth and one 125-bp complexity request…"
    try:
        status = _configure_api_credentials()
        with tempfile.TemporaryDirectory(prefix="hurdler_idt_test_") as temporary:
            scorer = IDTComplexityScorer(Path(temporary) / "audit.jsonl")
            summary = scorer.score("hurdler_credential_test", "ACGT" * 31 + "A")
        score = summary.get("idt_complexity_score")
        if not isinstance(score, (int, float)) or isinstance(score, bool):
            raise RuntimeError("IDT returned no finite numeric complexity score")
        credential_status.value = (
            f"<div style='color:#2d6a4f'><b>IDT status: verified.</b> "
            f"Authentication method: {status['auth_method']}; numeric response parsed successfully.</div>"
        )
    except Exception as exc:
        credential_status.value = (
            f"<div style='color:#b31b1b'><b>IDT test failed safely:</b> {type(exc).__name__}: "
            f"{str(exc)[:300]}</div>"
        )
    finally:
        clear_idt_secret_environment()
        credential_test_button.disabled = False
        credential_upload_test_button.disabled = False


def _download_credential_env(_button=None):
    try:
        if state.get("credential_payload") is None:
            if idt_setup_mode_widget.value != "create":
                raise RuntimeError("Create credentials in the hidden form before downloading idt.env")
            values = _create_credential_values()
            try:
                configure_idt_credentials_from_values(
                    values, auth_method=str(idt_auth_method_widget.value)
                )
                state["credential_payload"] = bytearray(_dotenv_payload(values))
                state["credential_auth_method"] = str(idt_auth_method_widget.value)
            finally:
                values.clear()
                clear_idt_secret_environment()
                _clear_create_secret_widgets()
        from google.colab import files as colab_files
        handle = tempfile.NamedTemporaryFile(
            prefix="hurdler_idt_", suffix=".env", dir="/tmp", delete=False
        )
        temporary_path = Path(handle.name)
        try:
            handle.write(bytes(state["credential_payload"]))
            handle.close()
            temporary_path.chmod(0o600)
            colab_files.download(str(temporary_path))
        finally:
            temporary_path.unlink(missing_ok=True)
        credential_status.value = "<b>IDT status:</b> idt.env sent to your browser; no copy was retained on the runtime filesystem"
    except ImportError:
        credential_status.value = "<b>IDT status:</b> browser download is available in hosted Colab only"
    except Exception as exc:
        credential_status.value = f"<b>IDT download failed safely:</b> {type(exc).__name__}: {str(exc)[:300]}"


idt_auth_method_widget.observe(_sync_idt_auth_fields, names="value")
idt_setup_mode_widget.observe(_sync_idt_setup, names="value")
credential_test_button.on_click(_test_idt_credentials)
credential_upload_test_button.on_click(_test_idt_credentials)
credential_download_button.on_click(_download_credential_env)
_sync_idt_auth_fields()
_sync_idt_setup()


def _build_design_request(*, ga_workers):
    route = state.get("confirmed_route")
    if route is None:
        raise RuntimeError("Confirm the RE/plasmid route before creating a GA request")
    query = _current_query()
    current = _query_fingerprint(query)
    if current != state.get("confirmed_fingerprint"):
        _invalidate_confirmation()
        raise RuntimeError("Protein/RE/plasmid settings changed; re-run the query and confirm again")
    minimum_secondary, maximum_secondary = _secondary_bounds()
    return DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection(
            route["candidate_id"], route["profile_id"], route["scheme_id"],
            str(state["confirmed_site_iii"]),
        ),
        validation_mode=_validation_mode(),
        assembly_strategy="exact_reused_secondary_rdl",
        population_size=int(population_number.value),
        mutation_rate=float(mutation_number.value),
        crossover_rate=float(crossover_number.value),
        elite_fraction=float(elite_number.value),
        seed=int(seed_number.value),
        generation_schedule=_generation_schedule(),
        score_weights={name: float(widget.value) for name, widget in weight_widgets.items()},
        auto_adjust_weights_from_idt=bool(auto_weight_feedback.value),
        minimum_secondary_copies=minimum_secondary,
        maximum_secondary_copies=maximum_secondary,
        max_idt_feedback_rounds=int(feedback_round_number.value),
        generations_per_feedback_round=int(generations_per_round_number.value),
        elite_seed_count=int(elite_seed_number.value),
        auto_adjust_ga_parameters_from_idt=bool(auto_parameter_feedback.value),
        max_population_size=int(max_population_number.value),
        max_mutation_rate=float(max_mutation_number.value),
        max_crossover_rate=float(max_crossover_number.value),
        ga_workers=int(ga_workers),
    )


def _bundle_fingerprint(request, resources):
    payload = {
        "request": asdict(request),
        "resources": asdict(resources),
        "credential_path": EXTERNAL_IDT_CREDENTIAL_PATH if request.validation_mode == "api" else None,
        "auth_method": "auto" if request.validation_mode == "api" else None,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _invalidate_external_bundle(_change=None):
    state["external_bundle"] = None
    state["external_bundle_fingerprint"] = None
    state["external_bundle_error"] = None


def _export_external_bundle(_button=None):
    with external_bundle_output:
        clear_output(wait=True)
        state["external_bundle_error"] = None
        try:
            resources = ExternalGAResources(
                worker_cpus=int(external_worker_cpus.value),
                memory_gb=int(external_memory_gb.value),
                walltime=str(external_walltime.value).strip(),
                partition=str(external_partition.value).strip(),
                account=str(external_account.value).strip(),
                qos=str(external_qos.value).strip(),
                constraint=str(external_constraint.value).strip(),
                conda_environment=str(external_conda_environment.value).strip(),
                result_directory=str(external_result_directory.value).strip(),
            )
            request = _build_design_request(ga_workers=resources.worker_cpus)
            commit = subprocess.check_output(
                ["git", "-C", str(project_root), "rev-parse", "HEAD"], text=True
            ).strip()
            bundle_root = _runtime_scratch_path(
                "hurdler_bundles", Path(output_directory_widget.value).parent / "bundles"
            )
            bundle = create_external_ga_bundle(
                request,
                bundle_root,
                repository_commit=commit,
                environment_file=project_root / "envs" / "hurdler.yml",
                resources=resources,
                idt_credential_path=EXTERNAL_IDT_CREDENTIAL_PATH,
                auth_method="auto",
            )
            state["external_bundle"] = bundle
            state["external_bundle_fingerprint"] = _bundle_fingerprint(request, resources)
            display(Markdown(
                f"**External GA bundle prepared:** `{bundle.name}`  \n"
                f"Frozen commit: `{commit}` · workers: `{resources.worker_cpus}` · memory: `{resources.memory_gb} GB`.  \n"
                "The ZIP contains only the external IDT path/auth method—no credential values."
            ))
            try:
                from google.colab import files as colab_files
            except ImportError:
                display(Markdown(f"Bundle path: `{bundle.resolve()}`"))
            else:
                colab_files.download(str(bundle))
        except Exception as exc:
            state["external_bundle_error"] = f"{type(exc).__name__}: {exc}"
            display(Markdown(f"**Bundle export failed safely:** `{type(exc).__name__}: {exc}`"))


def _download_design(_button=None):
    archive = state.get("archive")
    if archive is None or not Path(archive).is_file():
        with design_output:
            display(Markdown("**No ZIP is available yet.**"))
        return
    try:
        from google.colab import files as colab_files
    except ImportError:
        with design_output:
            display(Markdown(f"ZIP path: `{Path(archive).resolve()}`"))
        return
    colab_files.download(str(archive))


def _set_run_active(active):
    state["run_active"] = bool(active)
    for widget in ga_request_widgets:
        widget.disabled = bool(active)
    credential_upload.disabled = bool(active)
    confirmed = state.get("confirmed_route") is not None
    design_button.disabled = bool(active) or not confirmed or execution_target_widget.value != "colab"
    export_bundle_button.disabled = bool(active) or not confirmed or execution_target_widget.value != "external"
    pause_button.disabled = not active or execution_target_widget.value != "colab"
    stop_button.disabled = not active or execution_target_widget.value != "colab"
    if not active:
        pause_button.description = "Pause GA"
        pause_button.icon = "pause"
        design_button.description = "Run GA in Colab"
        design_button.icon = "play"


def _finish_design_success(result, files, archive, output_directory, drive_archive):
    state["design_files"] = files
    state["archive"] = archive
    state["design_result"] = result
    state["run_terminal_status"] = result.status
    download_button.disabled = False
    _prepare_viewer(result, output_directory)
    results_status.value = (
        "<div style='border:2px solid #2d6a4f;background:#effaf4;color:#111827;"
        f"border-radius:8px;padding:10px'><b>Result ready:</b> {html.escape(result.status)} · "
        f"{html.escape(result.message)}</div>"
    )
    stage_html.value = (
        "<div style='border:2px solid #2d6a4f;background:#effaf4;border-radius:8px;padding:10px'>"
        f"<b>Status:</b> {result.status}</div>"
    )
    with design_output:
        display(Markdown(f"**Design status:** `{result.status}` — {result.message}"))
        if result.rdl_plan:
            display(Markdown("### Exact-copy RDL equation"))
            display(pd.DataFrame([result.rdl_plan]))
        fragments = [*result.primary_fragments, *result.secondary_fragments]
        if fragments:
            display(Markdown("### Unique purchase fragments"))
            display(pd.DataFrame(fragments))
        if result.cloning_steps:
            display(Markdown("### Cloning plan"))
            display(pd.DataFrame(result.cloning_steps))
        if result.status == "optimized_unvalidated_batch":
            display(Markdown(
                "### IDT Bulk Input ready\nThis result has **not** been scored by the IDT API. "
                "Use `idt_bulk_input.csv`, `.tsv`, or `.fasta` from the ZIP, then return to the GA settings if no candidate is orderable."
            ))
            display(back_to_ga_button)
        drive_text = f" Drive copy: `{drive_archive}`." if drive_archive else ""
        display(Markdown(f"ZIP prepared: `{archive.name}`.{drive_text} No order was submitted."))
    _set_run_active(False)
    state["run_control"] = None
    if auto_download_widget.value and result.status in {"idt_accepted", "optimized_unvalidated_batch"}:
        _download_design()


def _finish_design_stopped():
    state["run_terminal_status"] = "stopped_by_user"
    stage_html.value = (
        "<div style='border:2px solid #b31b1b;background:#fff2f2;border-radius:8px;padding:10px'>"
        "<b>Status:</b> stopped_by_user · checkpoint preserved</div>"
    )
    results_status.value = (
        "<div style='border:2px solid #b7a57a;background:#fffaf0;color:#111827;"
        "border-radius:8px;padding:10px'><b>Stopped by user.</b> The latest checkpoint and audit were preserved.</div>"
    )
    with design_output:
        display(Markdown(
            "**GA stopped by the user at a safe point.** The current audit and checkpoint were preserved; "
            "no unfinished candidate is reported as accepted."
        ))
    _set_run_active(False)
    state["run_control"] = None


def _finish_design_error(error_type, message):
    state["run_terminal_status"] = "failed"
    stage_html.value = (
        "<div style='border:2px solid #b31b1b;background:#fff2f2;border-radius:8px;padding:10px'>"
        f"<b>Status:</b> failed · {error_type}</div>"
    )
    results_status.value = (
        "<div style='border:2px solid #b31b1b;background:#fff2f2;color:#111827;"
        f"border-radius:8px;padding:10px'><b>Run failed:</b> {html.escape(str(error_type))} · "
        f"{html.escape(str(message)[:500])}</div>"
    )
    with design_output:
        display(Markdown(f"**Design failed safely:** `{error_type}: {message[:500]}`"))
    _set_run_active(False)
    state["run_control"] = None


def _run_design_worker(request, query, output_directory, control, run_id):
    try:
        _enqueue_ui_event(
            "progress", int(run_id),
            DesignProgressEvent(
                stage="design",
                status="worker_entered",
                message="GA worker entered the design engine",
                elapsed_seconds=0.0,
            ),
        )
        scorer = (
            IDTComplexityScorer(output_directory / "idt_audit.jsonl")
            if request.validation_mode == "api"
            else None
        )
        result = design_construct_v2(
            request,
            idt_scorer=scorer,
            progress_callback=lambda event: _progress_update(event, run_id),
            checkpoint_callback=_checkpoint_update,
            run_control=control,
        )
        files = write_design_outputs_v2(result, output_directory)
        archive_root = _runtime_scratch_path(
            "hurdler_archives", output_directory.parent / "archives"
        )
        archive = timestamped_results_archive(
            output_directory, archive_root, sequence_id=query.sequence_id,
        )
        drive_archive = (
            _copy_archive_to_drive(archive)
            if storage_mode_widget.value == "drive"
            else None
        )
        _enqueue_ui_event(
            "success", int(run_id),
            (result, files, archive, output_directory, drive_archive),
        )
    except DesignRunStopped:
        state["run_terminal_status"] = "stopped_by_user"
        _persist_checkpoint(force=True)
        _enqueue_ui_event("stopped", int(run_id), None)
    except Exception as exc:
        _enqueue_ui_event(
            "error", int(run_id), (type(exc).__name__, str(exc))
        )
    finally:
        clear_idt_secret_environment()


def _run_design(_button=None):
    if state.get("run_active"):
        return
    route = state.get("confirmed_route")
    if route is None:
        with design_output:
            clear_output(wait=True)
            display(Markdown("**Confirm the RE/plasmid route before optimization.**"))
        return
    state["run_id"] = int(state.get("run_id", 0)) + 1
    run_id = int(state["run_id"])
    while True:
        try:
            state["progress_queue"].get_nowait()
        except queue.Empty:
            break
    state["progress_events"] = []
    state["visible_log_lines"] = []
    state["idt_score_events"] = []
    state["archive"] = None
    state["design_result"] = None
    state["best_checkpoint"] = None
    state["last_checkpoint_write"] = 0.0
    state["run_terminal_status"] = "running"
    state["run_started_monotonic"] = time.monotonic()
    state["last_progress_monotonic"] = state["run_started_monotonic"]
    generation_progress.value = 0
    candidate_progress.value = 0
    stage_html.value = (
        "<div style='border:2px solid #4b2e83;background:#f4f0fa;color:#111827;border-radius:8px;padding:10px'>"
        "<b>Status:</b> starting GA worker</div>"
    )
    request_line = f"run_requested  run_id={run_id}  preparing credentials and GA request"
    state["visible_log_lines"].append(request_line)
    attempt_log_html.value = (
        "<pre style='color:#111827;background:#ffffff'>"
        + html.escape(request_line)
        + "</pre>"
    )
    if _validation_mode() == "batch":
        _set_idt_plot_placeholder(
            "IDT API disabled. GA will export Bulk Input files; no score trajectory is claimed."
        )
    else:
        _set_idt_plot_placeholder(
            "GA started. Waiting for the first locally valid fragment to be scored by IDT."
        )
    results_status.value = (
        "<div style='border:2px solid #4b2e83;background:#f4f0fa;color:#111827;"
        f"border-radius:8px;padding:10px'><b>Run {run_id} started.</b> Live progress appears above.</div>"
    )
    with design_output:
        clear_output(wait=True)
    try:
        query = _current_query()
        output_directory = Path(output_directory_widget.value)
        output_directory.mkdir(parents=True, exist_ok=True)
        state["run_directory"] = output_directory
        state["last_checkpoint_write"] = time.monotonic()
        if storage_mode_widget.value == "drive" and not storage_state.get("drive_mounted"):
            raise RuntimeError("Google Drive was selected but is not mounted; click Mount Google Drive in Step 0")
        if _validation_mode() == "api":
            status = _configure_api_credentials()
            credential_status.value = (
                f"<b>IDT status:</b> configured in memory via {status['auth_method']}; values are hidden"
            )
        request = _build_design_request(ga_workers=1)
        control = DesignRunControl()
        state["run_control"] = control
        _set_run_active(True)
        worker = threading.Thread(
            target=_run_design_worker,
            args=(request, query, output_directory, control, run_id),
            name="hurdler-ga-worker",
            daemon=True,
        )
        state["run_thread"] = worker
        _start_ui_event_pump(run_id)
        worker.start()
        state["visible_log_lines"].append(f"worker_started run_id={run_id} thread={worker.name}")
        attempt_log_html.value = (
            "<pre style='color:#111827;background:#ffffff'>"
            + html.escape("\n".join(state["visible_log_lines"]))
            + "</pre>"
        )
    except Exception as exc:
        clear_idt_secret_environment()
        _finish_design_error(type(exc).__name__, str(exc))


def _pause_or_resume(_button=None):
    control = state.get("run_control")
    if control is None or not state.get("run_active"):
        return
    if control.paused:
        control.resume()
        pause_button.description = "Pause GA"
        pause_button.icon = "pause"
        stage_html.value = (
            "<div style='background:#f4f0fa;color:#111827;padding:8px'>"
            "<b>Status:</b> resume requested; continuing from the same GA population</div>"
        )
    else:
        control.pause()
        _persist_checkpoint(force=True)
        pause_button.description = "Resume GA"
        pause_button.icon = "play"
        stage_html.value = (
            "<div style='background:#fffaf0;color:#111827;padding:8px'>"
            "<b>Status:</b> pause requested; waiting for the next candidate safe point</div>"
        )


def _stop_design(_button=None):
    control = state.get("run_control")
    if control is None or not state.get("run_active"):
        return
    state["run_terminal_status"] = "stopped_by_user"
    control.stop()
    _persist_checkpoint(force=True)
    pause_button.disabled = True
    stop_button.disabled = True
    stage_html.value = (
        "<div style='background:#fff2f2;color:#111827;padding:8px'>"
        "<b>Status:</b> stop requested; finishing the current candidate/API unit</div>"
    )


def _back_to_ga_settings(_button=None):
    settings_mode.value = "advanced"
    display(Javascript(
        "document.getElementById('hurdler-ga-settings')?.scrollIntoView({behavior:'smooth', block:'start'});"
    ))


pause_button.on_click(_pause_or_resume)
stop_button.on_click(_stop_design)
back_to_ga_button.on_click(_back_to_ga_settings)


viewer_step_widget = widgets.Dropdown(description="Assembly step", layout=widgets.Layout(width="48%"))
viewer_molecule_widget = widgets.ToggleButtons(
    options=(("Plasmid", "plasmid"), ("Insert", "insert")), value="plasmid", description="Molecule",
)
viewer_view_widget = widgets.ToggleButtons(
    options=(("Circular", "circular"), ("Linear", "linear")), value="circular", description="View",
)
viewer_range_widget = widgets.IntRangeSlider(
    value=(0, 1), min=0, max=1, step=1, description="Range (bp)",
    continuous_update=False, layout=widgets.Layout(width="98%"),
)
viewer_focus_button = widgets.Button(description="Focus cloning region", icon="search-plus")
viewer_reset_button = widgets.Button(description="Reset full view", icon="expand")
viewer_render_button = widgets.Button(description="Render selected molecule", button_style="info")
viewer_status = widgets.HTML(
    "<div style='border:2px dashed #4b2e83;border-radius:8px;padding:12px'>"
    "Confirm an RE/plasmid route to preview its annotated step00 plasmid before GA.</div>"
)
viewer_output = widgets.Output()
viewer_translation_output = widgets.Output()
viewer_panel = widgets.VBox([
    widgets.HTML(
        "<h3>Stepwise plasmid / insert viewer</h3>"
        "<p>Select a step and molecule. Plasmids support circular and linear maps; inserts are linear. "
        "Use Focus cloning region or the range slider to zoom. At ≤300 bp the exact bases and per-codon translation are shown.</p>"
    ),
    viewer_status,
    widgets.HBox([viewer_step_widget, viewer_molecule_widget, viewer_view_widget]),
    viewer_range_widget,
    widgets.HBox([viewer_focus_button, viewer_reset_button, viewer_render_button]),
    viewer_output, viewer_translation_output,
])
viewer_panel.layout.display = ""


def _viewer_rows():
    return list(state.get("viewer_rows") or [])


def _viewer_record_path(row):
    directory = state.get("viewer_directory")
    if directory is None:
        raise ValueError("No plasmid preview or assembly output is available")
    return Path(directory) / row["file"]


def _reset_viewer_placeholder():
    state["viewer_rows"] = []
    state["viewer_directory"] = None
    viewer_step_widget.options = (("Confirm a route first", None),)
    viewer_step_widget.value = None
    viewer_status.value = (
        "<div style='border:2px dashed #4b2e83;border-radius:8px;padding:12px'>"
        "Confirm an RE/plasmid route to preview its annotated step00 plasmid before GA.</div>"
    )
    with viewer_output:
        clear_output(wait=True)
    with viewer_translation_output:
        clear_output(wait=True)


def _viewer_selected_row():
    matches = [
        row for row in _viewer_rows()
        if int(row["step"]) == int(viewer_step_widget.value)
        and row["molecule"] == viewer_molecule_widget.value
    ]
    if len(matches) != 1:
        raise ValueError("The selected step/molecule is unavailable")
    return matches[0]


def _viewer_update_molecules(_change=None):
    if viewer_step_widget.value is None:
        return
    available = [
        row["molecule"] for row in _viewer_rows()
        if int(row["step"]) == int(viewer_step_widget.value)
    ]
    options = tuple((value.title(), value) for value in ("plasmid", "insert") if value in available)
    viewer_molecule_widget.options = options
    if options and viewer_molecule_widget.value not in {value for _label, value in options}:
        viewer_molecule_widget.value = options[0][1]
    _viewer_reset()


def _viewer_reset(_button=None):
    try:
        row = _viewer_selected_row()
    except (ValueError, TypeError):
        return
    record = SeqIO.read(_viewer_record_path(row), "genbank")
    viewer_range_widget.max = max(1, len(record))
    viewer_range_widget.value = (0, len(record))
    viewer_view_widget.disabled = row["molecule"] == "insert"
    if row["molecule"] == "insert":
        viewer_view_widget.value = "linear"


def _viewer_focus(_button=None):
    row = _viewer_selected_row()
    start = int(row.get("cloning_region_start_0based", 0))
    end = int(row.get("cloning_region_end_0based_exclusive", row["length_bp"]))
    padding = max(20, min(200, (end - start) // 10))
    viewer_range_widget.value = (
        max(0, start - padding), min(int(row["length_bp"]), end + padding)
    )
    viewer_view_widget.value = "linear"
    _render_viewer()


def _render_viewer(_button=None):
    with viewer_output:
        clear_output(wait=True)
        try:
            row = _viewer_selected_row()
            record = SeqIO.read(_viewer_record_path(row), "genbank")
            start, end = map(int, viewer_range_widget.value)
            circular = row["molecule"] == "plasmid" and viewer_view_widget.value == "circular" and (start, end) == (0, len(record))
            def feature_filter(feature):
                if feature.type == "source":
                    return False
                if not circular:
                    return True
                qualifiers = feature.qualifiers
                if qualifiers.get("feature_kind", [""])[0] == "restriction_site":
                    return True
                if feature.type == "CDS":
                    return True
                feature_class = qualifiers.get("feature_class", [""])[0]
                return feature_class in {
                    "antibiotic_resistance", "origin", "replication_origin",
                    "promoter", "terminator", "operator",
                }
            translator = BiopythonTranslator(features_filters=(feature_filter,))
            if circular:
                graphic = translator.translate_record(record, record_class=CircularGraphicRecord)
                figure, axis = plt.subplots(figsize=(8, 8))
                graphic.plot(ax=axis)
            else:
                graphic = translator.translate_record(record)
                if (start, end) != (0, len(record)):
                    graphic = graphic.crop((start, end))
                figure, axis = plt.subplots(figsize=(13, 3.6))
                graphic.plot(ax=axis, figure_width=13)
            axis.set_title(f"{row['file']} · {start:,}–{end:,} bp")
            figure.tight_layout()
            display(figure)
            plt.close(figure)
        except Exception as exc:
            display(Markdown(f"**Viewer error:** `{type(exc).__name__}: {exc}`"))
            return
    with viewer_translation_output:
        clear_output(wait=True)
        window = str(record.seq[start:end])
        display(Markdown(
            f"**Molecule:** `{row['file']}` · **length:** {len(record):,} bp · "
            f"**SHA256:** `{row['sequence_sha256']}` · **copy count:** {row.get('copy_count', '—')}"
        ))
        target_cds = next((
            feature for feature in record.features
            if feature.type == "CDS" and "repeat-protein" in feature.qualifiers.get("label", [""])[0]
        ), None)
        if target_cds is not None:
            display(Markdown(f"**Expressed translation ({len(target_cds.qualifiers['translation'][0]):,} AA):** `{target_cds.qualifiers['translation'][0]}`"))
        if end - start <= 300:
            display(Markdown(f"**Bases {start + 1}–{end}:** `{window}`"))
            if target_cds is not None:
                cds_start, cds_end = int(target_cds.location.start), int(target_cds.location.end)
                aligned_start = max(start, cds_start)
                aligned_start += (3 - (aligned_start - cds_start) % 3) % 3
                aligned_end = min(end, cds_end)
                aligned_end -= (aligned_end - cds_start) % 3
                codon_rows = []
                for position in range(aligned_start, aligned_end, 3):
                    codon = str(record.seq[position:position + 3])
                    codon_rows.append({
                        "bp": f"{position + 1}-{position + 3}",
                        "codon": codon,
                        "AA": translate_dna(codon),
                        "protein_position": (position - cds_start) // 3 + 1,
                    })
                if codon_rows:
                    display(pd.DataFrame(codon_rows))


def _prepare_viewer(result, output_directory):
    rows = list(result.assembly_steps)
    if not rows:
        _reset_viewer_placeholder()
        return
    state["viewer_rows"] = rows
    state["viewer_directory"] = Path(output_directory)
    steps = sorted({int(row["step"]) for row in rows})
    viewer_step_widget.options = tuple((f"Step {step:02d}", step) for step in steps)
    viewer_step_widget.value = steps[-1]
    viewer_status.value = "<b>Complete assembly timeline loaded.</b> The final plasmid is selected by default."
    _viewer_update_molecules()
    _render_viewer()


def _prepare_route_preview():
    result = state.get("query_result")
    route = state.get("confirmed_route")
    site_iii = state.get("confirmed_site_iii")
    if result is None or route is None or not site_iii:
        _reset_viewer_placeholder()
        return
    candidates = [
        row for row in result.protein_candidates
        if row["candidate_id"] == route["candidate_id"]
    ]
    if len(candidates) != 1:
        raise ValueError("Confirmed route does not resolve to one protein candidate")
    preview_root = _runtime_scratch_path(
        "hurdler_route_preview", Path(output_directory_widget.value).parent / "route_preview"
    )
    preview_root.mkdir(parents=True, exist_ok=True)
    record, row = build_step00_plasmid_record(route, candidates[0], str(site_iii))
    SeqIO.write(record, preview_root / row["file"], "genbank")
    state["viewer_rows"] = [row]
    state["viewer_directory"] = preview_root
    viewer_step_widget.options = (("Step 00", 0),)
    viewer_step_widget.value = 0
    viewer_molecule_widget.options = (("Plasmid", "plasmid"),)
    viewer_molecule_widget.value = "plasmid"
    viewer_view_widget.disabled = False
    viewer_view_widget.value = "circular"
    viewer_status.value = (
        f"<div style='border:2px solid #2d6a4f;background:#effaf4;border-radius:8px;padding:10px'>"
        f"<b>Route preview ready:</b> {route['profile_id']} · {route['site_i_enzyme']} / "
        f"{route['site_ii_enzyme']} / {site_iii}. Full annotations remain in GenBank; the circular map "
        "shows only key vector functions and selected RE sites to prevent label overlap.</div>"
    )
    _viewer_reset()
    _render_viewer()


viewer_step_widget.observe(_viewer_update_molecules, names="value")
viewer_molecule_widget.observe(lambda _change: _viewer_reset(), names="value")
viewer_focus_button.on_click(_viewer_focus)
viewer_reset_button.on_click(_viewer_reset)
viewer_render_button.on_click(_render_viewer)


design_button.on_click(_run_design)
export_bundle_button.on_click(_export_external_bundle)
download_button.on_click(_download_design)

for ga_bundle_widget in (
    idt_setup_mode_widget, idt_auth_method_widget, secondary_search_mode_widget,
    secondary_copy_range_widget, minimum_secondary_number, maximum_secondary_number,
    population_number, mutation_number, crossover_number, elite_number,
    feedback_round_number, generations_per_round_number, elite_seed_number,
    max_population_number, max_mutation_number, max_crossover_number,
    seed_number, generation_schedule_widget, auto_weight_feedback,
    auto_parameter_feedback, external_worker_cpus, external_memory_gb,
    external_walltime, external_partition, external_account, external_qos,
    external_constraint, external_conda_environment, external_result_directory,
    *weight_widgets.values(),
):
    ga_bundle_widget.observe(_invalidate_external_bundle, names="value")

basic_panel = widgets.VBox([
    widgets.HTML("<div id='hurdler-ga-settings'></div>"),
    _help_card("Settings level", settings_mode, unit="mode", default="recommended", purpose="Shows or hides low-level GA controls.", allowed="recommended or advanced", effect="Recommended mode still uses the displayed frozen defaults."),
    auto_download_widget,
    _help_card("Secondary-copy search", secondary_search_mode_widget, unit="mode", default="bounded copy range", purpose="Explores secondary donor repeat count.", allowed="bounded 1–50 or automatic to physical limit", effect="Bounded is the tutorial default; automatic proves the physical/route limit."),
    secondary_range_card,
    secondary_length_status,
    feedback_round_card,
    _help_card("Runtime work folder", output_directory_widget, unit="path", default="/content/hurdler_runs/current", purpose="Stores active computation and complete uncompressed outputs.", allowed="runtime path", effect="Drive mode still computes here and copies only ZIP archives."),
])
ga_request_widgets = (
    input_mode_widget, sequence_id_widget, n_cap_widget, repeat_module_widget,
    initial_copies_widget, c_cap_widget, full_protein_widget,
    repeat_start_widget, repeat_end_widget, repeat_period_widget,
    allow_left_cutter_widget, allow_right_cutter_widget,
    enzyme_bulk_control, plasmid_bulk_control, max_restoration_length_widget,
    query_button, pair_choice, site_iii_choice, profile_choice, scheme_choice, confirm_button,
    storage_mode_widget, drive_root_widget, mount_drive_button,
    execution_target_widget,
    settings_mode, idt_setup_mode_widget, idt_auth_method_widget,
    idt_client_id_widget, idt_client_secret_widget, idt_username_widget,
    idt_password_widget, idt_access_token_widget, credential_upload,
    credential_test_button, credential_upload_test_button, credential_download_button,
    secondary_search_mode_widget, secondary_copy_range_widget,
    minimum_secondary_number, maximum_secondary_number,
    population_number, mutation_number, crossover_number, elite_number,
    feedback_round_number, generations_per_round_number, elite_seed_number,
    max_population_number, max_mutation_number, max_crossover_number,
    seed_number, generation_schedule_widget, auto_weight_feedback,
    auto_parameter_feedback, output_directory_widget, auto_download_widget,
    external_worker_cpus, external_memory_gb, external_walltime,
    external_partition, external_account, external_qos, external_constraint,
    external_conda_environment, external_result_directory,
    *enzyme_checkboxes.values(), *plasmid_checkboxes.values(),
    *weight_widgets.values(),
)
colab_execution_panel = widgets.VBox([
    widgets.HTML(
        "<div style='border-left:6px solid #4b2e83;background:#f4f0fa;color:#111827;"
        "padding:10px'><b>Run inside this Colab runtime.</b> Pause and Stop act at safe points.</div>"
    ),
    widgets.HBox([design_button, pause_button, stop_button]),
])
external_execution_panel = widgets.VBox([
    external_resource_panel,
    export_bundle_button,
    external_bundle_output,
])


def _sync_execution_target(_change=None):
    colab = execution_target_widget.value == "colab"
    colab_execution_panel.layout.display = "" if colab else "none"
    external_execution_panel.layout.display = "none" if colab else ""
    confirmed = state.get("confirmed_route") is not None
    active = bool(state.get("run_active"))
    design_button.disabled = active or not confirmed or not colab
    export_bundle_button.disabled = active or not confirmed or colab
    pause_button.disabled = not active or not colab
    stop_button.disabled = not active or not colab


execution_target_widget.observe(_sync_execution_target, names="value")
ga_panel = widgets.VBox([
    basic_panel,
    advanced_panel,
    _help_card(
        "Execution target", execution_target_widget, unit="mode", default="Run in Colab",
        purpose="Chooses either the interactive Colab worker or one portable Local/Slurm bundle.",
        allowed="exactly one target", effect="Only controls for the selected target are shown.",
    ),
    colab_execution_panel,
    external_execution_panel,
    widgets.HTML("<h3 style='color:#3b1f69'>Live GA / IDT log</h3>"),
    stage_html, generation_progress, candidate_progress, current_html, attempt_log_html,
    widgets.HTML("<h3 style='color:#3b1f69'>IDT score trajectory</h3>"),
    idt_plot_status, idt_plot_output, idt_score_table_output,
])
results_panel = widgets.VBox([
    results_status,
    download_button,
    design_output,
])
secondary_search_mode_widget.observe(_update_secondary_lengths, names="value")
secondary_copy_range_widget.observe(_range_slider_changed, names="value")
minimum_secondary_number.observe(_range_number_changed, names="value")
maximum_secondary_number.observe(_range_number_changed, names="value")
repeat_module_widget.observe(_update_secondary_lengths, names="value")
_update_secondary_lengths()
_reset_viewer_placeholder()
_ = _refresh_live_support()
_sync_execution_target()

setup_module = widgets.VBox([
    widgets.HTML("<h2>1. Storage and IDT setup</h2>"),
    storage_panel,
    idt_credential_panel,
])
protein_module = widgets.VBox([
    widgets.HTML("<h2>2. Protein sequence and repeat definition</h2>"),
    protein_input_panel,
])
route_advanced_panel = widgets.Accordion(children=[widgets.VBox([
    cutter_policy_panel,
    advanced_route_filters.children[0],
])])
route_advanced_panel.set_title(0, "Advanced route filters and cutter fallback")
route_advanced_panel.selected_index = None
route_filter_module = widgets.VBox([
    widgets.HTML(
        "<h2>3. Enzyme/plasmid filters and HURDLER query</h2>"
        "<p>RE and plasmid selections share the same filtered route universe. The cards report "
        "jointly usable routes, not checkbox totals.</p>"
    ),
    widgets.HBox([enzyme_bulk_control, enzyme_selection_status]),
    enzyme_route_support,
    enzyme_checkbox_grid,
    widgets.HBox([plasmid_bulk_control, plasmid_selection_status]),
    plasmid_route_support,
    plasmid_checkbox_grid,
    route_advanced_panel,
    query_button,
    query_output,
])
route_selection_module = widgets.VBox([
    widgets.HTML(
        "<h2>4. Route selection and confirmation</h2>"
        "<p>Select Site I/II, Site III, plasmid and cut scheme explicitly. Placeholder menus remain "
        "visible until the preceding choice is available.</p>"
    ),
    pair_choice, site_iii_choice, profile_choice, scheme_choice,
    confirm_button, route_output,
])
ga_module = widgets.VBox([
    widgets.HTML("<h2>5. GA optimization and execution</h2>"),
    ga_panel,
])
viewer_module = widgets.VBox([
    widgets.HTML("<h2>6. Interactive plasmid and insert viewer</h2>"),
    viewer_panel,
])
result_module = widgets.VBox([
    widgets.HTML("<h2>7. Results and downloads</h2>"),
    results_panel,
])
None
'''


COLAB_TUTORIAL_STEPS = (
    (
        "hurdler-step-1-setup",
        "1. Storage and IDT setup",
        "display(setup_module)\nNone",
    ),
    (
        "hurdler-step-2-protein",
        "2. Protein sequence and repeat definition",
        "display(protein_module)\nNone",
    ),
    (
        "hurdler-step-3-query",
        "3. Enzyme/plasmid filters and HURDLER query",
        "display(route_filter_module)\nNone",
    ),
    (
        "hurdler-step-4-route",
        "4. Route selection and confirmation",
        "display(route_selection_module)\nNone",
    ),
    (
        "hurdler-step-5-ga",
        "5. GA optimization and execution",
        "display(ga_module)\nNone",
    ),
    (
        "hurdler-step-6-viewer",
        "6. Interactive plasmid and insert viewer",
        "display(viewer_module)\nNone",
    ),
    (
        "hurdler-step-7-results",
        "7. Results and downloads",
        "display(result_module)\nNone",
    ),
)


COLAB_ENZYME_SELECTOR = r'''display(widgets.VBox([
    widgets.HTML(
        f"<h2>2. Select individual HURDLER enzymes</h2>"
        f"<p>One shared pool: {len(enzyme_roles['site_i'])} Site-I/II enzymes and "
        f"{len(declared_site_iii)} maintained Site-III enzymes "
        f"({len(enzyme_roles['site_iii'])} occur in the current protein-pair index). "
        f"Each selected enzyme is used only in legal roles. Site I is retained active, Site II is synonymously silenced, "
        f"and Site III releases disposable secondary adapters. The card below reports jointly usable routes—not checkbox counts.</p>"
    ),
    widgets.HBox([enzyme_bulk_control, enzyme_selection_status]),
    enzyme_route_support,
    enzyme_checkbox_grid,
]))'''


COLAB_PLASMID_SELECTOR = r'''display(widgets.VBox([
    widgets.HTML("<h2>3. Select plasmid profiles</h2><p>Profiles are annotation-aware expression orientations. Selecting a plasmid does not guarantee a route: its cut scheme must also support the chosen RE pair and restore threshold.</p>"),
    widgets.HBox([plasmid_bulk_control, plasmid_selection_status]),
    plasmid_route_support,
    plasmid_checkbox_grid,
]))'''


COLAB_QUERY_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>4. Query protein patterns and annotation-safe vector routes</h2><p>The query enumerates protein-compatible Site-I/Site-II codon windows, then intersects them with selected Site III enzymes, plasmid profiles, protected annotations and the inclusive restore-length cap.</p>"),
    advanced_route_filters,
    query_button, query_output,
]))
_run_query()'''


COLAB_RE_ROUTE_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>5. Select the three RE enzymes</h2><p>Choose a ranked Site-I/Site-II pair, then one compatible Site III. No first-ranked route is auto-confirmed.</p>"),
    pair_choice, site_iii_choice,
]))'''


COLAB_VECTOR_ROUTE_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>6. Select plasmid and cut scheme, then confirm</h2><p>The cut scheme specifies the two vector cutters, retained long backbone, removed MCS arc and exact restoration bases. Confirmation freezes its fingerprint for optimization.</p>"),
    profile_choice, scheme_choice, confirm_button, route_output,
]))'''


COLAB_GA_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>7. Run in Colab or export a Local / Slurm GA bundle</h2>"),
    widgets.HTML("This panel remains hidden until an RE/plasmid route is confirmed."),
    ga_panel, viewer_panel,
]))'''


def _code_cell(
    source: str,
    *,
    colab: bool,
    cell_id: str,
    title: str,
    tags: list[str] | None = None,
    form_view: bool = False,
):
    if colab:
        display_mode = ' { display-mode: "form" }' if form_view else ""
        source = f"#@title {title}{display_mode}\n" + source
    cell = nbf.v4.new_code_cell(source)
    cell["id"] = cell_id
    cell.metadata["id"] = cell_id
    if tags:
        cell.metadata["tags"] = tags
    if colab:
        # Keep the large bootstrap hidden, but let each tiny display cell use the
        # normal Colab execution surface.  Empty parameter forms can otherwise be
        # mistaken for a code-snippet launcher and show no widget output.
        cell.metadata["cellView"] = "form" if form_view else "both"
        cell.metadata["colab"] = {}
    return cell


def notebook(*, colab: bool = False):
    nb = nbf.v4.new_notebook()
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    }
    intro = nbf.v4.new_markdown_cell(COLAB_INTRO if colab else INTRO)
    intro["id"] = "hurdler-introduction"
    intro.metadata["id"] = "hurdler-introduction"
    cells = [intro]
    if colab:
        application_source = "\n\n".join(
            (COLAB_BOOTSTRAP, PARAMETERS, COLAB_IMPORTS, COLAB_CONTROLLER_V2, SMOKE)
        )
        cells.append(
            _code_cell(
                application_source,
                colab=True,
                cell_id="hurdler-initialize",
                title="Initialize HURDLER tutorial",
                form_view=True,
            )
        )
        for cell_id, title, source in COLAB_TUTORIAL_STEPS:
            cells.append(
                _code_cell(
                    source,
                    colab=True,
                    cell_id=cell_id,
                    title=title,
                )
            )
    else:
        cells.extend([
            _code_cell(
                PARAMETERS,
                colab=False,
                cell_id="hurdler-parameters",
                title="Notebook defaults",
                tags=["parameters"],
            ),
            _code_cell(
                IMPORTS,
                colab=False,
                cell_id="hurdler-imports",
                title="Load the HURDLER design engine",
            ),
            nbf.v4.new_markdown_cell(
                "## 1. Enter a protein, query RE pairs, and confirm a route"
            ),
            _code_cell(
                WIDGETS,
                colab=False,
                cell_id="hurdler-query-controls",
                title="Protein and RE-pair/vector controls",
            ),
            nbf.v4.new_markdown_cell(
                "## 2. Optional codon optimization and IDT scoring/export"
            ),
            _code_cell(
                OPTIMIZER,
                colab=False,
                cell_id="hurdler-optimization-controls",
                title="GA, IDT, and file export controls",
            ),
            _code_cell(
                SMOKE,
                colab=False,
                cell_id="hurdler-headless-smoke",
                title="Automated validation hook (normally inactive)",
            ),
        ])
    for index, cell in enumerate(cells):
        stable_id = cell.metadata.setdefault("id", f"hurdler-cell-{index:02d}")
        cell["id"] = stable_id
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
