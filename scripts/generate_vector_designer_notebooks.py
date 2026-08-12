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


COLAB_INTRO = """# Annotation-aware HURDLER designer v2 · exact-repeat RDL workflow

The protein defaults are visible as native Colab Forms before code runs. Choose
**Runtime → Run all** once to install HURDLER and reveal separate single-enzyme,
plasmid, RE-solution, and cut-scheme selectors. Source code stays folded.

HURDLER never silently chooses route 1 or a Site-III enzyme. After the query,
select the three enzymes, plasmid, and cut scheme, then confirm the route. Only
then does the GA/IDT panel appear. In
full-protein mode, a missing boundary intentionally stops after showing
candidate start/end/period values; enter the confirmed values in Step 1 and
press **Re-run HURDLER query**. No source code needs to be opened.

IDT is used only for complexity scoring, never for codon optimization or
ordering. Locally the default credential file is
`~/.config/hurdler/idt.env`; hosted Colab requests a temporary upload. Values
are read only when optimization starts and cleared from the process afterward.
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
            assembly_strategy="legacy_adaptive_max" if input_mode.value == "split" else "single_exact",
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

repository_dir = Path("/content/clone_repeat_protein")
try:
    import google.colab  # noqa: F401
except ImportError:
    running_in_colab = repository_dir.parent.is_dir()
else:
    running_in_colab = True

if running_in_colab:
    if (repository_dir / ".git").is_dir():
        subprocess.run(
            ["git", "-C", str(repository_dir), "fetch", "--depth=1", "origin", repository_ref],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(repository_dir), "checkout", "--detach", "FETCH_HEAD"],
            check=True,
        )
    else:
        subprocess.run([
            "git", "clone", "--branch", repository_ref, "--single-branch",
            "https://github.com/Wenzhao-protein/clone_repeat_protein",
            str(repository_dir),
        ], check=True)
    os.chdir(repository_dir)
    subprocess.run([sys.executable, "-m", "pip", "install", "-e", ".[notebooks,optimization]"], check=True)
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
print(f"HURDLER ready: {Path(hurdler_package.__file__).resolve()}")
'''


COLAB_PROTEIN_FORM = r'''#@markdown ### Choose one protein-input mode
input_mode = "N-cap + repeat module + C-cap" #@param ["N-cap + repeat module + C-cap", "Complete exact protein / FASTA"]
sequence_id = "" #@param {type:"string", placeholder:"optional; FASTA header is used when blank"}

#@markdown ### Split input (used only in N-cap/module/C-cap mode)
n_cap_aa = "MGSHHHHHHSSGIEGRSSGYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGGGGSGGGGSLEVLFQGPDLPKLVKLLKSSNEEILLKALRALAEIASGG" #@param {type:"string"}
repeat_module_aa = "NEQIQAVIDAGALPALVQLLSSPNEQILQEALWALSNIASGG" #@param {type:"string"}
initial_repeat_copies = 25 #@param {type:"integer"}
c_cap_aa = "NEQIQAVIDAGALPALVQLLSSPNEQILQEALWALSNIASGGNEQKQAVKEAGALEKLEQLQSHENEKIQKEAQEALEKLQSHGGGLEVLFQGPSSGEFGGGGSMVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYKGGGSSGHHHHHH" #@param {type:"string"}

#@markdown ### Complete-protein input (used only in complete-protein mode)
#@markdown Paste one raw AA sequence or one FASTA record. Every residue and repeat variant is preserved.
full_protein_or_fasta = "" #@param {type:"string", placeholder:"raw AA or one FASTA record"}
#@markdown Enter all three confirmed 1-based values, or leave all as 0 to request boundary candidates.
repeat_region_start_1based = 0 #@param {type:"integer"}
repeat_region_end_1based = 0 #@param {type:"integer"}
repeat_period_aa = 0 #@param {type:"integer"}
'''


COLAB_SELECTOR_POLICY_FORM = r'''#@markdown ### Cutter reuse policy
#@markdown Individual RE and plasmid selectors appear in the next two cells after initialization.
allow_left_cutter_in_hurdler_pair = False #@param {type:"boolean"}
allow_right_cutter_in_hurdler_pair = False #@param {type:"boolean"}
'''


COLAB_IMPORTS = r'''import getpass
import hashlib
import json
import os
import shutil
import time
import traceback
from dataclasses import asdict
from pathlib import Path

import pandas as pd
import ipywidgets as widgets
from IPython.display import Markdown, clear_output, display

from hurdler.design import parse_protein_input
from hurdler.idt import (
    IDT_CREDENTIAL_PATH,
    IDTComplexityScorer,
    clear_idt_secret_environment,
    configure_idt_credentials,
    configure_idt_credentials_from_bytes,
    configure_idt_credentials_from_values,
)
from hurdler.design import role_enzyme_options
from hurdler.protein_index import ProteinPatternIndex
from hurdler.progress import DesignProgressEvent
from hurdler.vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    DesignSelection,
    design_construct_v2,
    design_query,
    bundled_protein_index_dir,
    write_design_outputs_v2,
)

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

protein_index = ProteinPatternIndex.load(bundled_protein_index_dir())
enzyme_roles = role_enzyme_options(protein_index)
project_root = bundled_protein_index_dir().parents[2]
declared_site_iii = tuple(sorted(
    pd.read_csv(project_root / "output" / "selected_site_iii_enzymes.csv")["enzyme"].astype(str).unique()
))
all_enzyme_options = tuple(sorted(
    set(enzyme_roles["site_i"]) | set(enzyme_roles["site_ii"]) | set(declared_site_iii)
))

enzyme_selector = widgets.SelectMultiple(
    options=all_enzyme_options, value=all_enzyme_options, rows=18,
    description="Allowed RE", layout=widgets.Layout(width="98%"),
)
enzyme_all_button = widgets.Button(description="Select all RE", icon="check-square")
enzyme_none_button = widgets.Button(description="Select no RE", icon="square-o")
plasmid_selector = widgets.SelectMultiple(
    options=PLASMID_OPTIONS, value=PLASMID_OPTIONS, rows=8,
    description="Plasmids", layout=widgets.Layout(width="98%"),
)
plasmid_all_button = widgets.Button(description="Select all plasmids", icon="check-square")
plasmid_none_button = widgets.Button(description="Select no plasmids", icon="square-o")

state = {
    "query_result": None,
    "query_fingerprint": None,
    "confirmed_route": None,
    "confirmed_site_iii": None,
    "confirmed_fingerprint": None,
    "design_files": None,
    "archive": None,
    "progress_events": [],
}

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


def _selected_role_enzymes(role):
    selected = set(enzyme_selector.value)
    values = tuple(name for name in enzyme_roles[role] if name in selected)
    if not values:
        raise ValueError(f"Select at least one enzyme eligible for {role.replace('_', ' ').title()}")
    return values


def _selected_plasmids():
    selected = tuple(str(value) for value in plasmid_selector.value)
    if not selected:
        raise ValueError("Select at least one plasmid profile")
    return selected


def _current_query():
    common = dict(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        sequence_id=str(sequence_id).strip() or "interactive_design",
        site_i_allowlist=_selected_role_enzymes("site_i"),
        site_ii_allowlist=_selected_role_enzymes("site_ii"),
        site_iii_allowlist=_selected_role_enzymes("site_iii"),
        plasmid_allowlist=_selected_plasmids(),
        allow_left_cutter_in_hurdler_pair=bool(allow_left_cutter_in_hurdler_pair),
        allow_right_cutter_in_hurdler_pair=bool(allow_right_cutter_in_hurdler_pair),
    )
    if input_mode == "N-cap + repeat module + C-cap":
        return CompatibilityQuery(
            input_mode="split", n_cap=n_cap_aa, repeat_module=repeat_module_aa,
            c_cap=c_cap_aa, repeat_copies=int(initial_repeat_copies), **common,
        )
    parsed_id, normalized = parse_protein_input(full_protein_or_fasta)
    coordinates = (
        int(repeat_region_start_1based), int(repeat_region_end_1based), int(repeat_period_aa),
    )
    if any(value > 0 for value in coordinates) and not all(value > 0 for value in coordinates):
        raise ValueError("Provide repeat start, end, and period together, or leave all three as 0")
    confirmed = all(value > 0 for value in coordinates)
    common["sequence_id"] = str(sequence_id).strip() or parsed_id
    return CompatibilityQuery(
        input_mode="full", full_protein_sequence=normalized,
        repeat_region_start=coordinates[0] if confirmed else None,
        repeat_region_end=coordinates[1] if confirmed else None,
        repeat_period=coordinates[2] if confirmed else None,
        **common,
    )


def _query_fingerprint(query):
    payload = json.dumps(asdict(query), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def _invalidate_confirmation(message=""):
    state["confirmed_route"] = None
    state["confirmed_site_iii"] = None
    state["confirmed_fingerprint"] = None
    confirm_button.disabled = True
    ga_panel.layout.display = "none"
    design_button.disabled = True
    download_button.disabled = True
    if message:
        with route_output:
            clear_output(wait=True)
            display(Markdown(message))


def _selection_changed(_change=None):
    if state.get("confirmed_route") is not None:
        _invalidate_confirmation("**Selection changed. Re-run the query and confirm a new route.**")


def _set_all_enzymes(_button=None):
    enzyme_selector.value = all_enzyme_options


def _set_no_enzymes(_button=None):
    enzyme_selector.value = ()


def _set_all_plasmids(_button=None):
    plasmid_selector.value = PLASMID_OPTIONS


def _set_no_plasmids(_button=None):
    plasmid_selector.value = ()


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
        _invalidate_confirmation()
        pair_choice.options = [("Run the query first", None)]
        site_iii_choice.options = [("Choose Site I / II first", None)]
        profile_choice.options = [("Choose all three RE first", None)]
        scheme_choice.options = [("Choose a plasmid first", None)]
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
            display(pd.DataFrame(result.boundary_analysis.get("candidates", [])))
            display(Markdown("Confirm start/end/period in the protein form, then re-run this query."))
            return
        if not result.vector_routes:
            display(Markdown("No annotation-safe vector route is available."))
            return
        routes = pd.DataFrame(result.vector_routes)
        pair_rows = (
            routes.sort_values("rank").drop_duplicates(["site_i_enzyme", "site_ii_enzyme"])
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
        ga_panel.layout.display = ""
        design_button.disabled = False
        display(Markdown(
            f"**Confirmed:** {route['site_i_enzyme']} / {route['site_ii_enzyme']} / "
            f"{site_iii_choice.value} → {route['profile_id']} · {route['cut_scheme']}. "
            "GA/IDT controls are now enabled below."
        ))


enzyme_all_button.on_click(_set_all_enzymes)
enzyme_none_button.on_click(_set_no_enzymes)
plasmid_all_button.on_click(_set_all_plasmids)
plasmid_none_button.on_click(_set_no_plasmids)
enzyme_selector.observe(_selection_changed, names="value")
plasmid_selector.observe(_selection_changed, names="value")
pair_choice.observe(_update_site_iii, names="value")
site_iii_choice.observe(_update_profiles, names="value")
profile_choice.observe(_update_schemes, names="value")
query_button.on_click(_run_query)
confirm_button.on_click(_confirm_route)


def _numeric_control(description, value, minimum, maximum, step, *, integer=False):
    slider_type = widgets.IntSlider if integer else widgets.FloatSlider
    text_type = widgets.BoundedIntText if integer else widgets.BoundedFloatText
    slider = slider_type(value=value, min=minimum, max=maximum, step=step, readout=False,
                         layout=widgets.Layout(width="62%"))
    number = text_type(value=value, min=minimum, max=maximum, step=step,
                       layout=widgets.Layout(width="34%"))
    widgets.link((slider, "value"), (number, "value"))
    card = widgets.VBox([
        widgets.HTML(f"<b>{description}</b>"), widgets.HBox([slider, number]),
    ], layout=widgets.Layout(width="49%", border="1px solid #ddd", padding="6px"))
    return card, number


settings_mode = widgets.ToggleButtons(
    options=[("Keep recommended defaults", "basic"), ("Advanced settings", "advanced")],
    value="basic", description="Settings",
)
validation_mode_widget = widgets.ToggleButtons(
    options=[("Live IDT API", "api"), ("IDT Bulk files (unvalidated)", "batch"), ("Compatibility only", "none")],
    value="api", description="Validation",
)
credential_source = widgets.Dropdown(
    options=[("Automatic local env / Colab upload", "auto"), ("Colab Secrets", "secrets"), ("Hidden runtime prompt", "prompt")],
    value="auto", description="Credentials", layout=widgets.Layout(width="98%"),
)
credential_path = widgets.Text(
    value=str(IDT_CREDENTIAL_PATH), description="Local env", layout=widgets.Layout(width="98%"),
)
credential_upload = widgets.FileUpload(accept=".env,text/plain", multiple=False, description="Upload temporary idt.env")
output_directory_widget = widgets.Text(value="/content/hurdler_design", description="Output", layout=widgets.Layout(width="98%"))
auto_download_widget = widgets.Checkbox(value=True, description="Auto-download ZIP after success")
verbose_generations = widgets.Checkbox(value=False, description="Show every GA generation in Advanced log")

population_card, population_number = _numeric_control("Population", 16, 4, 256, 4, integer=True)
mutation_card, mutation_number = _numeric_control("Mutation rate", 0.08, 0.001, 0.5, 0.001)
crossover_card, crossover_number = _numeric_control("Crossover rate", 0.75, 0.0, 1.0, 0.01)
elite_card, elite_number = _numeric_control("Elite fraction", 0.15, 0.01, 0.5, 0.01)
seed_number = widgets.IntText(value=42, description="Random seed")
generation_schedule_widget = widgets.Text(value="10,20,40,60,80,100", description="Generations")
auto_weight_feedback = widgets.Checkbox(value=True, description="Adjust weights from IDT positive rules")

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
    _two_per_row([population_card, mutation_card, crossover_card, elite_card]),
    widgets.HBox([seed_number, generation_schedule_widget]),
    widgets.HBox([auto_weight_feedback, verbose_generations]),
    widgets.HTML("<b>GA score weights</b>"),
    _two_per_row(list(weight_widgets.values())),
])


def _sync_settings(_change=None):
    advanced_panel.layout.display = "" if settings_mode.value == "advanced" else "none"


settings_mode.observe(_sync_settings, names="value")
_sync_settings()

credential_help = widgets.HTML(
    "<b>IDT env format</b> (choose one method; never mix them):<br>"
    "<code>IDT_ACCESS_TOKEN=...</code><br>or<br>"
    "<code>IDT_CLIENT_ID=...</code><br><code>IDT_CLIENT_SECRET=...</code><br>"
    "<code>IDT_USERNAME=...</code><br><code>IDT_PASSWORD=...</code><br>"
    "Local setup: <code>mkdir -p ~/.config/hurdler &amp;&amp; chmod 700 ~/.config/hurdler</code>, "
    "save as <code>~/.config/hurdler/idt.env</code>, then <code>chmod 600</code>. "
    "Hosted Colab cannot read your local home directory; upload the env file temporarily."
)

stage_html = widgets.HTML("<b>Status:</b> waiting for route confirmation")
generation_progress = widgets.IntProgress(value=0, min=0, max=1, description="GA")
current_html = widgets.HTML("")
attempt_log_html = widgets.HTML("<pre>No attempts yet.</pre>")
design_output = widgets.Output()
design_button = widgets.Button(description="Optimize exact target / export", button_style="success", disabled=True)
download_button = widgets.Button(description="Download design ZIP", icon="download", disabled=True)


def _generation_schedule():
    values = tuple(sorted({
        int(value.strip()) for value in generation_schedule_widget.value.split(",") if value.strip()
    }))
    if not values or values[-1] != 100 or any(value <= 0 for value in values):
        raise ValueError("Generation schedule must contain positive integers and terminate at 100")
    return values


def _progress_update(event: DesignProgressEvent):
    state["progress_events"].append(event.to_dict())
    stage_html.value = f"<b>Status:</b> {event.stage} · {event.status}"
    if event.generations:
        generation_progress.max = max(1, int(event.generations))
        generation_progress.value = min(generation_progress.max, int(event.generation or 0))
    current_html.value = (
        f"<b>{event.fragment_kind or 'design'}</b> · copies={event.copies if event.copies is not None else '—'} "
        f"· generation={event.generation if event.generation is not None else '—'}/"
        f"{event.generations if event.generations is not None else '—'} "
        f"· best score={event.ga_score if event.ga_score is not None else '—'} "
        f"· elapsed={event.elapsed_seconds or 0:.1f}s"
    )
    keep = event.status in {"attempt_completed", "request_completed", "completed", "failed"}
    keep = keep or (verbose_generations.value and event.stage == "ga")
    if keep:
        lines = [
            f"{row['stage']:<12} {row['status']:<18} {row.get('fragment_kind') or '-':<10} "
            f"copies={row.get('copies')} gen={row.get('generation')}/{row.get('generations')} "
            f"score={row.get('ga_score')}"
            for row in state["progress_events"][-16:]
            if row["status"] in {"attempt_completed", "request_completed", "completed", "failed"}
            or (verbose_generations.value and row["stage"] == "ga")
        ]
        attempt_log_html.value = "<pre>" + "\n".join(lines[-12:]) + "</pre>"


def _secret_value(reader, name):
    try:
        return str(reader(name) or "").strip()
    except Exception:
        return ""


def _configure_colab_secrets(reader=None):
    if reader is None:
        from google.colab import userdata
        reader = userdata.get
    token = _secret_value(reader, "IDT_ACCESS_TOKEN")
    if token:
        return configure_idt_credentials_from_values({"IDT_ACCESS_TOKEN": token}, auth_method="access_token")
    values = {name: _secret_value(reader, name) for name in ("IDT_CLIENT_ID", "IDT_CLIENT_SECRET", "IDT_USERNAME", "IDT_PASSWORD")}
    try:
        return configure_idt_credentials_from_values(values, auth_method="password")
    finally:
        values.clear()


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
        description="Upload temporary idt.env",
    )
    credential_upload_row.children = (credential_path, credential_upload)
    try:
        previous.close()
    except Exception:
        # Replacement already removed the only live UI reference.  Some
        # hosted widget backends do not implement close() completely.
        pass


def _configure_api_credentials():
    if credential_source.value == "secrets":
        return _configure_colab_secrets()
    if credential_source.value == "prompt":
        return configure_idt_credentials(mode="manual", auth_method="access_token", prompt=getpass.getpass)
    local_path = Path(credential_path.value).expanduser()
    if local_path.is_file():
        return configure_idt_credentials(mode="path", path=local_path, include_path_in_status=False)
    payload = _uploaded_payload()
    if payload is None:
        raise FileNotFoundError(
            "No external IDT env file is available. Hosted Colab cannot access ~/.config on your computer; "
            "upload idt.env above or choose Colab Secrets."
        )
    try:
        return configure_idt_credentials_from_bytes(payload)
    finally:
        payload = b""
        _clear_credential_upload()


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


def _run_design(_button=None):
    route = state.get("confirmed_route")
    if route is None:
        with design_output:
            clear_output(wait=True)
            display(Markdown("**Confirm the RE/plasmid route before optimization.**"))
        return
    design_button.disabled = True
    design_button.description = "Running GA / IDT…"
    download_button.disabled = True
    state["progress_events"] = []
    state["archive"] = None
    generation_progress.value = 0
    stage_html.value = "<b>Status:</b> starting"
    with design_output:
        clear_output(wait=True)
    try:
        query = _current_query()
        current = _query_fingerprint(query)
        if current != state.get("confirmed_fingerprint"):
            _invalidate_confirmation()
            raise RuntimeError("Protein/RE/plasmid settings changed; re-run the query and confirm again")
        mode = validation_mode_widget.value
        scorer = None
        output_directory = Path(output_directory_widget.value)
        if mode == "api":
            _configure_api_credentials()
            scorer = IDTComplexityScorer(output_directory / "idt_audit.jsonl")
        request = DesignRequestV2(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            query=query,
            selection=DesignSelection(
                route["candidate_id"], route["profile_id"], route["scheme_id"],
                str(state["confirmed_site_iii"]),
            ),
            validation_mode=mode,
            assembly_strategy="exact_reused_secondary_rdl",
            population_size=int(population_number.value),
            mutation_rate=float(mutation_number.value),
            crossover_rate=float(crossover_number.value),
            elite_fraction=float(elite_number.value),
            seed=int(seed_number.value),
            generation_schedule=_generation_schedule(),
            score_weights={name: float(widget.value) for name, widget in weight_widgets.items()},
            auto_adjust_weights_from_idt=bool(auto_weight_feedback.value),
        )
        result = design_construct_v2(request, idt_scorer=scorer, progress_callback=_progress_update)
        files = write_design_outputs_v2(result, output_directory)
        archive = Path(shutil.make_archive(str(output_directory.resolve()), "zip", root_dir=output_directory.resolve()))
        state["design_files"] = files
        state["archive"] = archive
        download_button.disabled = False
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
            display(Markdown(f"ZIP prepared: `{archive.name}`. No order was submitted."))
        if auto_download_widget.value and result.status in {"idt_accepted", "optimized_unvalidated_batch"}:
            _download_design()
    except Exception as exc:
        stage_html.value = f"<b>Status:</b> failed · {type(exc).__name__}"
        with design_output:
            display(Markdown(f"**Design failed safely:** `{type(exc).__name__}: {exc}`"))
            if settings_mode.value == "advanced":
                display(widgets.HTML("<details><summary>Sanitized traceback</summary><pre>" + traceback.format_exc() + "</pre></details>"))
    finally:
        clear_idt_secret_environment()
        design_button.description = "Optimize exact target / export"
        design_button.disabled = state.get("confirmed_route") is None


design_button.on_click(_run_design)
download_button.on_click(_download_design)

credential_upload_row = widgets.HBox([credential_path, credential_upload])
basic_panel = widgets.VBox([
    settings_mode,
    widgets.HBox([validation_mode_widget, auto_download_widget]),
    credential_source,
    credential_upload_row,
    credential_help,
    output_directory_widget,
])
ga_panel = widgets.VBox([
    basic_panel, advanced_panel,
    widgets.HTML("<h3>Live progress</h3>"), stage_html, generation_progress,
    current_html, attempt_log_html,
    widgets.HBox([design_button, download_button]), design_output,
])
ga_panel.layout.display = "none"
'''


COLAB_ENZYME_SELECTOR = r'''display(widgets.VBox([
    widgets.HTML(
        f"<h2>2. Select individual HURDLER enzymes</h2>"
        f"<p>One shared pool: {len(enzyme_roles['site_i'])} Site-I/II enzymes and "
        f"{len(declared_site_iii)} maintained Site-III enzymes "
        f"({len(enzyme_roles['site_iii'])} occur in the current protein-pair index). "
        f"Each selected enzyme is used only in legal roles.</p>"
    ),
    widgets.HBox([enzyme_all_button, enzyme_none_button]), enzyme_selector,
]))'''


COLAB_PLASMID_SELECTOR = r'''display(widgets.VBox([
    widgets.HTML("<h2>3. Select plasmid profiles</h2>"),
    widgets.HBox([plasmid_all_button, plasmid_none_button]), plasmid_selector,
]))'''


COLAB_QUERY_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>4. Query protein patterns and annotation-safe vector routes</h2>"),
    query_button, query_output,
]))
_run_query()'''


COLAB_RE_ROUTE_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>5. Select the three RE enzymes</h2>"),
    pair_choice, site_iii_choice,
]))'''


COLAB_VECTOR_ROUTE_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>6. Select plasmid and cut scheme, then confirm</h2>"),
    profile_choice, scheme_choice, confirm_button, route_output,
]))'''


COLAB_GA_PANEL = r'''display(widgets.VBox([
    widgets.HTML("<h2>7. GA, IDT scoring, progress, and export</h2>"),
    widgets.HTML("This panel remains hidden until an RE/plasmid route is confirmed."),
    ga_panel,
]))'''


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
    cell["id"] = cell_id
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
    intro = nbf.v4.new_markdown_cell(COLAB_INTRO if colab else INTRO)
    intro["id"] = "hurdler-introduction"
    intro.metadata["id"] = "hurdler-introduction"
    cells = [intro]
    if colab:
        cells.extend(
            [
                _code_cell(
                    COLAB_BOOTSTRAP,
                    colab=True,
                    cell_id="hurdler-initialize",
                    title='0. Initialize HURDLER { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_PROTEIN_FORM,
                    colab=True,
                    cell_id="hurdler-protein-form",
                    title='1. Protein and repeat boundary { display-mode: "form", run: "auto" }',
                    tags=["parameters", "colab-native-form"],
                ),
                _code_cell(
                    COLAB_SELECTOR_POLICY_FORM,
                    colab=True,
                    cell_id="hurdler-selector-policy-form",
                    title='1b. Cutter fallback policy { display-mode: "form", run: "auto" }',
                    tags=["parameters", "colab-native-form"],
                ),
                _code_cell(
                    PARAMETERS,
                    colab=True,
                    cell_id="hurdler-test-defaults",
                    title='Internal smoke-test defaults { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_IMPORTS,
                    colab=True,
                    cell_id="hurdler-imports",
                    title='Load the HURDLER design engine { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_CONTROLLER_V2,
                    colab=True,
                    cell_id="hurdler-controller-v2",
                    title='Prepare interactive controllers { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_ENZYME_SELECTOR,
                    colab=True,
                    cell_id="hurdler-enzyme-selector",
                    title='2. Select individual RE enzymes { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_PLASMID_SELECTOR,
                    colab=True,
                    cell_id="hurdler-plasmid-selector",
                    title='3. Select plasmids { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_QUERY_PANEL,
                    colab=True,
                    cell_id="hurdler-query-panel",
                    title='4. Run HURDLER query { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_RE_ROUTE_PANEL,
                    colab=True,
                    cell_id="hurdler-re-route-panel",
                    title='5. Select Site I, II, and III { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_VECTOR_ROUTE_PANEL,
                    colab=True,
                    cell_id="hurdler-vector-route-panel",
                    title='6. Select plasmid and cut scheme { display-mode: "form" }',
                ),
                _code_cell(
                    COLAB_GA_PANEL,
                    colab=True,
                    cell_id="hurdler-ga-panel",
                    title='7. Optimize exact target and export { display-mode: "form" }',
                ),
                _code_cell(
                    SMOKE,
                    colab=True,
                    cell_id="hurdler-headless-smoke",
                    title='Automated validation hook (normally inactive) { display-mode: "form" }',
                ),
            ]
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
