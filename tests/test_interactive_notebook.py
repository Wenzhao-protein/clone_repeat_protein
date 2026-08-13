from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import threading
import time
import zipfile
from types import ModuleType, SimpleNamespace
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

from hurdler.cli import main
from hurdler.design import bundled_index_dir
from hurdler.idt import clear_idt_secret_environment


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "workflows" / "01_interactive_hurdler_designer.ipynb"
COLAB_NOTEBOOK = ROOT / "notebooks" / "workflows" / "02_colab_hurdler_designer.ipynb"


def _colab_payload():
    return json.loads(COLAB_NOTEBOOK.read_text())


def _colab_cell_source(cell_id: str) -> str:
    cell = next(
        item for item in _colab_payload()["cells"]
        if item.get("metadata", {}).get("id") == cell_id
    )
    return "".join(cell["source"])


def _colab_runtime_sources() -> tuple[str, str, str]:
    """Split the single generated application cell for fast controller tests."""
    source = _colab_cell_source("hurdler-initialize")
    parameters_start = source.index("# Papermill parameters")
    imports_start = source.index("import asyncio", parameters_start)
    controller_start = source.index("PLASMID_OPTIONS = (", imports_start)
    smoke_start = source.index("if headless_smoke or", controller_start)
    return (
        source[parameters_start:imports_start],
        source[imports_start:controller_start],
        source[controller_start:smoke_start],
    )


@pytest.fixture
def colab_runtime_namespace():
    namespace: dict[str, object] = {}
    parameters, imports, controller = _colab_runtime_sources()
    exec(compile(parameters, "hurdler-parameters", "exec"), namespace)
    exec(compile(imports, "hurdler-imports", "exec"), namespace)
    namespace["display"] = lambda *_args, **_kwargs: None
    namespace["clear_output"] = lambda *_args, **_kwargs: None
    exec(compile(controller, "hurdler-controller-v2", "exec"), namespace)
    yield namespace
    clear_idt_secret_environment()


def test_notebook_source_is_output_free_and_has_no_widget_state():
    payload = json.loads(NOTEBOOK.read_text())
    assert payload["nbformat"] == 4
    assert "widgets" not in payload.get("metadata", {})
    assert any("parameters" in cell.get("metadata", {}).get("tags", []) for cell in payload["cells"])
    for cell in payload["cells"]:
        if cell["cell_type"] == "code":
            assert cell.get("outputs", []) == []
            assert cell.get("execution_count") is None
    text = NOTEBOOK.read_text()
    assert "IDT_CLIENT_SECRET=" not in text
    assert "IDT_PASSWORD=" not in text
    assert "/home/wendai/.config/hurdler/idt.env" not in text


def test_colab_cells_are_named_hidden_forms_without_duplicate_visible_titles():
    payload = _colab_payload()
    ids = []
    for cell in payload["cells"]:
        cell_id = cell.get("metadata", {}).get("id")
        assert cell_id
        assert cell.get("id") == cell_id
        ids.append(cell_id)
        if cell["cell_type"] != "code":
            continue
        assert cell["metadata"]["cellView"] == "form"
        assert cell["metadata"]["colab"] == {}
        assert "".join(cell["source"]).startswith("#@title ")
        assert '{ display-mode: "form" }' in "".join(cell["source"]).splitlines()[0]
        assert cell.get("outputs", []) == []
        assert cell.get("execution_count") is None
    assert len(ids) == len(set(ids))


def test_colab_tutorial_uses_two_mutually_exclusive_input_panels():
    payload = _colab_payload()
    intro = "".join(next(cell for cell in payload["cells"] if cell["id"] == "hurdler-introduction")["source"])
    controller = _colab_runtime_sources()[2]
    assert "Tutorial workflow" in intro
    assert "Site I" in intro and "Site II" in intro and "Site III" in intro
    assert "Google Drive" in intro and "RDL" in intro and "IDT" in intro
    assert 'options=(("N-cap / module / C-cap", "split"), ("Complete protein / FASTA", "full"))' in controller
    assert 'value="split"' in controller
    assert 'full_input_panel.layout.display = "" if input_mode_widget.value == "full" else "none"' in controller
    assert "_help_card" in controller
    assert [cell["id"] for cell in payload["cells"]] == [
        "hurdler-introduction",
        "hurdler-initialize",
        "hurdler-step-1-setup",
        "hurdler-step-2-protein",
        "hurdler-step-3-query",
        "hurdler-step-4-route",
        "hurdler-step-5-ga",
        "hurdler-step-6-viewer",
        "hurdler-step-7-results",
        "hurdler-reconnect-interface",
    ]
    assert "widgets.Password(" in COLAB_NOTEBOOK.read_text()
    assert "Create Credentials" in controller
    assert "Upload idt.env" in controller
    assert "Manual credential input" in controller
    assert "Do not use IDT API — export batch input" in controller
    assert 'description="Test uploaded credentials"' not in controller
    assert "credential_upload_test_button" not in controller
    assert "credential_test_button" not in controller
    assert "widgets.FileUpload(" not in controller
    for secret in (
        "IDT_ACCESS_TOKEN", "IDT_CLIENT_ID", "IDT_CLIENT_SECRET",
        "IDT_USERNAME", "IDT_PASSWORD",
    ):
        assert not re.search(rf"^{secret}\s*=.*#@param", COLAB_NOTEBOOK.read_text(), re.MULTILINE)


def test_colab_has_separate_tutorial_steps_with_viewer_after_ga():
    sources = {
        cell["metadata"]["id"]: "".join(cell["source"])
        for cell in _colab_payload()["cells"]
    }
    assert "display(setup_module)" in sources["hurdler-step-1-setup"]
    assert "display(protein_module)" in sources["hurdler-step-2-protein"]
    assert "display(route_filter_module)" in sources["hurdler-step-3-query"]
    assert "display(route_selection_module)" in sources["hurdler-step-4-route"]
    assert "display(ga_module)" in sources["hurdler-step-5-ga"]
    assert "display(viewer_module)" in sources["hurdler-step-6-viewer"]
    assert "display(result_module)" in sources["hurdler-step-7-results"]
    assert "display(reconnect_module)" in sources["hurdler-reconnect-interface"]
    assert "_install_colab_ui_bridge()" in sources["hurdler-step-5-ga"]
    assert "_install_colab_ui_bridge()" in sources["hurdler-reconnect-interface"]
    controller = _colab_runtime_sources()[2]
    assert "Select all RE" in controller
    assert "Select none" in controller
    assert "Select all plasmids" in controller
    assert "enzyme_checkboxes" in controller
    assert "plasmid_checkboxes" in controller
    assert "Advanced route filters" in controller
    assert 'value=100' in controller
    assert "max_restoration_length_bp=int(max_restoration_length_widget.value)" in controller
    assert "widgets.GridBox" in controller
    for module in (
        "setup_module", "protein_module", "route_filter_module",
        "route_selection_module", "viewer_module", "ga_module", "result_module",
    ):
        assert module in controller
    assert 'value="create"' in controller
    assert 'assembly_strategy="exact_reused_secondary_rdl"' in controller
    assert "widgets.BoundedIntText" in controller
    assert "widgets.BoundedFloatText" in controller
    assert "widgets.IntRangeSlider" in controller
    assert 'value=(12, 20), min=1, max=50' in controller
    assert 'options=(("Automatic to limit", "automatic"), ("Bounded copy range", "bounded"))' in controller
    assert '"Maximum GA→IDT feedback rounds", 100' in controller
    assert '"GA generations per feedback round", 10' in controller
    assert '"Warm-start top candidates", 10' in controller
    assert "minimum_secondary_copies=minimum_secondary" in controller
    assert "maximum_secondary_copies=maximum_secondary" in controller
    assert "max_idt_feedback_rounds=int(feedback_round_number.value)" in controller
    assert "auto_adjust_ga_parameters_from_idt=bool(auto_parameter_feedback.value)" in controller
    assert "feedback={event.feedback_round" in controller
    assert 'ga_panel.layout.display = "none"' not in controller
    assert "credential_upload.value = ()" not in controller
    assert "_choose_colab_credential_payload()" in controller
    imports = _colab_runtime_sources()[1]
    assert "colab_output.enable_custom_widget_manager()" in imports
    assert "write_secondary_checkpoint" in imports
    assert "timestamped_results_archive" in imports
    assert "create_external_ga_bundle" in imports
    assert 'description="Run GA in Colab"' in controller
    assert 'description="Download Local / Slurm bundle"' in controller
    assert controller.count('description="Download Local / Slurm bundle"') == 1
    assert 'options=(("Run in Colab", "colab"), ("Local / Slurm bundle", "external"))' in controller
    assert 'value=16, min=1, max=1024, description="GA worker CPUs"' in controller
    assert 'value=32, min=1, max=1_048_576, description="Total memory (GB)"' in controller
    assert 'value="24:00:00", description="Walltime"' in controller
    assert 'value="cpu", description="Partition"' in controller
    assert 'EXTERNAL_IDT_CREDENTIAL_PATH = "~/.config/hurdler/idt.env"' in controller
    assert "external_idt_credential_path" not in controller
    assert 'description="Pause GA"' in controller
    assert 'description="Stop GA"' in controller
    assert "CircularGraphicRecord" in imports
    assert "idt_plot_output" in controller
    assert 'status="fragment_scored"' not in controller
    assert "_ui_event_pump" in controller
    assert "_enqueue_ui_event" in controller
    assert 'register_callback("hurdler.ui_drain"' in controller
    assert "google.colab.kernel.invokeFunction('hurdler.ui_drain'" in controller
    assert "ui_bridge_output" in controller
    assert "reconnect_tabs" in controller
    assert "call_soon_threadsafe" in controller
    assert "tutorial_app" not in controller


def test_colab_bootstrap_does_not_depend_on_optional_release_tag():
    source = _colab_cell_source("hurdler-initialize")
    assert "COLAB_RELEASE_TAG" not in source
    assert "import google.colab" in source
    assert '"pip", "install", "-e"' in source
    assert 'sys.path.insert(0, source_dir)' in source
    namespace: dict[str, object] = {}
    exec(compile(source, "hurdler-initialize", "exec"), namespace)
    assert namespace["hurdler_package"].__name__ == "hurdler"


def test_colab_run_all_queries_but_never_auto_confirms_rank_one(colab_runtime_namespace):
    colab_runtime_namespace["_run_query"]()
    state = colab_runtime_namespace["state"]
    assert state["query_result"].status == "compatible_unoptimized"
    assert len(state["query_result"].protein_candidates) == 770
    assert len(state["query_result"].vector_routes) == 1001
    assert all(
        row["restoration_length_bp"] <= 100
        for row in state["query_result"].vector_routes
    )
    assert state["query_result"].request["max_restoration_length_bp"] == 100
    assert state["confirmed_route"] is None
    assert colab_runtime_namespace["pair_choice"].value is None
    assert colab_runtime_namespace["site_iii_choice"].value is None
    assert colab_runtime_namespace["design_button"].disabled is True
    assert colab_runtime_namespace["storage_mode_widget"].value == "runtime"
    assert colab_runtime_namespace["storage_state"]["drive_mounted"] is False
    assert colab_runtime_namespace["input_mode_widget"].value == "split"
    assert colab_runtime_namespace["split_input_panel"].layout.display == ""
    assert colab_runtime_namespace["full_input_panel"].layout.display == "none"


def test_colab_input_buttons_switch_visibility_and_invalidate_routes(colab_runtime_namespace):
    namespace = colab_runtime_namespace
    assert namespace["state"]["route_universe"] is not None
    namespace["input_mode_widget"].value = "full"
    assert namespace["split_input_panel"].layout.display == "none"
    assert namespace["full_input_panel"].layout.display == ""
    assert namespace["state"]["route_universe"] is None
    assert namespace["state"]["confirmed_route"] is None


def test_colab_secondary_search_supports_automatic_and_bounded_copy_ranges(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    assert namespace["secondary_search_mode_widget"].value == "bounded"
    assert namespace["_secondary_bounds"]() == (12, 20)
    assert namespace["secondary_copy_range_widget"].value == (12, 20)
    assert namespace["secondary_copy_range_widget"].min == 1
    assert namespace["secondary_copy_range_widget"].max == 50
    namespace["minimum_secondary_number"].value = 4
    namespace["maximum_secondary_number"].value = 9
    assert namespace["_secondary_bounds"]() == (4, 9)
    assert namespace["secondary_copy_range_widget"].value == (4, 9)
    namespace["secondary_copy_range_widget"].value = (6, 11)
    assert namespace["minimum_secondary_number"].value == 6
    assert namespace["maximum_secondary_number"].value == 11
    assert "core" in namespace["secondary_length_status"].value
    assert "purchase" in namespace["secondary_length_status"].value
    namespace["minimum_secondary_number"].value = 15
    assert namespace["maximum_secondary_number"].value == 15
    namespace["secondary_search_mode_widget"].value = "automatic"
    assert namespace["_secondary_bounds"]() == (1, None)
    assert namespace["secondary_copy_range_widget"].disabled is True


def test_colab_complete_fasta_preserves_header_and_requires_boundary_confirmation(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    unit = "ACDEFGHIKLMNPQRSTVWY"
    namespace["input_mode_widget"].value = "full"
    namespace["sequence_id_widget"].value = ""
    namespace["full_protein_widget"].value = f">header_from_fasta details\nM{unit * 3}G"
    namespace["repeat_start_widget"].value = 0
    namespace["repeat_end_widget"].value = 0
    namespace["repeat_period_widget"].value = 0
    query = namespace["_current_query"]()
    assert query.sequence_id == "header_from_fasta"
    namespace["_run_query"]()
    assert namespace["state"]["query_result"].status == "needs_boundary_confirmation"
    assert namespace["design_button"].disabled is True


def test_colab_form_edit_invalidates_confirmed_route(colab_runtime_namespace):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    assert namespace["state"]["confirmed_route"] is not None
    assert namespace["design_button"].disabled is False
    assert namespace["ga_module"].layout.display in (None, "")
    namespace["n_cap_widget"].value = "MM"
    assert namespace["state"]["confirmed_route"] is None
    assert namespace["design_button"].disabled is True
    assert namespace["state"]["viewer_rows"] == []


def test_colab_route_confirmation_previews_step00_before_ga(colab_runtime_namespace):
    namespace = colab_runtime_namespace
    assert namespace["viewer_panel"].layout.display == ""
    assert namespace["state"]["viewer_rows"] == []
    _confirm_first_colab_route(namespace)
    rows = namespace["state"]["viewer_rows"]
    assert len(rows) == 1
    assert rows[0]["file"] == "step00_plasmid.gb"
    assert Path(namespace["state"]["viewer_directory"], rows[0]["file"]).is_file()
    assert namespace["viewer_step_widget"].value == 0
    assert namespace["viewer_molecule_widget"].value == "plasmid"
    assert namespace["viewer_view_widget"].value == "circular"
    assert "Route preview ready" in namespace["viewer_status"].value


def test_colab_external_bundle_freezes_request_and_is_invalidated_by_ga_edits(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    assert namespace["export_bundle_button"].disabled is True
    _confirm_first_colab_route(namespace)
    namespace["execution_target_widget"].value = "external"
    assert namespace["export_bundle_button"].disabled is False
    namespace["idt_setup_mode_widget"].value = "batch"
    namespace["output_directory_widget"].value = str(tmp_path / "runtime")
    namespace["external_worker_cpus"].value = 4
    namespace["external_memory_gb"].value = 12
    namespace["external_walltime"].value = "02:30:00"
    namespace["external_partition"].value = "cpu"
    namespace["secondary_search_mode_widget"].value = "bounded"
    namespace["minimum_secondary_number"].value = 3
    namespace["maximum_secondary_number"].value = 7
    namespace["_export_external_bundle"]()
    assert namespace["state"]["external_bundle_error"] is None, namespace["state"]["external_bundle_error"]
    bundle = Path(namespace["state"]["external_bundle"])
    assert bundle.is_file()
    with zipfile.ZipFile(bundle) as archive:
        request = json.loads(archive.read("request.json"))
        script = archive.read("run_ga.sh").decode()
    assert request["ga_workers"] == 4
    assert request["minimum_secondary_copies"] == 3
    assert request["maximum_secondary_copies"] == 7
    assert request["query"]["max_restoration_length_bp"] == 100
    assert "#SBATCH --cpus-per-task=4" in script
    assert "#SBATCH --mem=12G" in script
    assert "#SBATCH --time=02:30:00" in script
    assert "--idt-credential-file" not in script
    namespace["population_number"].value = int(namespace["population_number"].value) + 4
    assert namespace["state"]["external_bundle"] is None


def test_colab_manual_credentials_are_hidden_auto_tested_and_cleared(colab_runtime_namespace):
    namespace = colab_runtime_namespace

    class PassingScorer:
        def __init__(self, _audit_path):
            pass

        def score(self, _sequence_id, _sequence):
            return {"idt_complexity_score": 0.0}

    namespace["IDTComplexityScorer"] = PassingScorer
    namespace["idt_manual_mode_button"].click()
    namespace["idt_auth_method_widget"].value = "access_token"
    namespace["idt_access_token_widget"].value = "temporary-token"
    assert namespace["idt_setup_mode_widget"].value == "manual"
    assert namespace["credential_create_panel"].layout.display == "none"
    assert namespace["credential_manual_panel"].layout.display == ""
    assert "IDT status: verified" in namespace["credential_status"].value
    assert namespace["idt_access_token_widget"].value == ""
    assert isinstance(namespace["state"]["credential_payload"], bytearray)
    assert b"temporary-token" in bytes(namespace["state"]["credential_payload"])
    assert "IDT_ACCESS_TOKEN" not in os.environ
    clear_idt_secret_environment()


def test_colab_uploaded_env_is_selected_and_tested_automatically(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace

    class PassingScorer:
        def __init__(self, _audit_path):
            pass

        def score(self, _sequence_id, _sequence):
            return {"idt_complexity_score": 0.0}

    namespace["IDTComplexityScorer"] = PassingScorer
    payload = b"IDT_ACCESS_TOKEN=temporary-upload-token\n"
    namespace["_choose_colab_credential_payload"] = lambda: payload
    namespace["idt_upload_mode_button"].click()
    assert namespace["idt_setup_mode_widget"].value == "upload"
    assert namespace["credential_create_panel"].layout.display == "none"
    assert namespace["credential_upload_panel"].layout.display == ""
    assert "IDT status: verified" in namespace["credential_status"].value
    assert isinstance(namespace["state"]["credential_payload"], bytearray)
    assert bytes(namespace["state"]["credential_payload"]) == payload
    assert namespace["idt_mode_action_row"].children[1] is namespace["idt_upload_mode_button"]
    assert namespace["idt_upload_mode_button"].icon == "upload"
    assert namespace["state"]["pending_credential_upload"] is None
    assert "IDT_ACCESS_TOKEN" not in os.environ
    clear_idt_secret_environment()


def test_colab_native_chooser_returns_bytes_without_public_disk_upload(
    monkeypatch, colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    payload = b"IDT_ACCESS_TOKEN=native-picker-token\n"
    calls = []
    fake_files = SimpleNamespace(
        _upload_files=lambda *, multiple: calls.append(multiple) or {"idt.env": payload}
    )
    google_module = ModuleType("google")
    colab_module = ModuleType("google.colab")
    colab_module.files = fake_files
    google_module.colab = colab_module
    monkeypatch.setitem(sys.modules, "google", google_module)
    monkeypatch.setitem(sys.modules, "google.colab", colab_module)
    assert namespace["_choose_colab_credential_payload"]() == payload
    assert calls == [False]


def test_colab_manual_password_grant_auto_tests_after_fourth_field(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace

    class PassingScorer:
        def __init__(self, _audit_path):
            pass

        def score(self, _sequence_id, _sequence):
            return {"idt_complexity_score": 0.0}

    namespace["IDTComplexityScorer"] = PassingScorer
    namespace["idt_manual_mode_button"].click()
    namespace["idt_client_id_widget"].value = "client-id"
    namespace["idt_client_secret_widget"].value = "client-secret"
    namespace["idt_username_widget"].value = "username"
    assert "waiting for all required" in namespace["credential_status"].value
    namespace["idt_password_widget"].value = "password"
    assert "IDT status: verified" in namespace["credential_status"].value
    assert namespace["idt_client_secret_widget"].value == ""
    assert namespace["idt_password_widget"].value == ""
    assert namespace["idt_client_id_widget"].value == "client-id"
    assert namespace["idt_username_widget"].value == "username"
    assert "IDT_CLIENT_SECRET" in bytes(namespace["state"]["credential_payload"]).decode()
    assert "IDT_PASSWORD" in bytes(namespace["state"]["credential_payload"]).decode()
    assert "IDT_PASSWORD" not in os.environ
    clear_idt_secret_environment()


@pytest.mark.parametrize(
    "uploaded",
    [
        ({"name": "idt.env", "content": memoryview(b"IDT_ACCESS_TOKEN=test-token\n")},),
        {"idt.env": {"content": b"IDT_ACCESS_TOKEN=test-token\n"}},
        {"idt.env": b"IDT_ACCESS_TOKEN=test-token\n"},
        {"content": "IDT_ACCESS_TOKEN=test-token\n"},
    ],
)
def test_colab_upload_parser_accepts_hosted_widget_value_shapes(
    colab_runtime_namespace, uploaded
):
    namespace = colab_runtime_namespace
    assert namespace["_payload_from_upload_value"](uploaded) == b"IDT_ACCESS_TOKEN=test-token\n"


def test_colab_native_upload_cancel_stays_in_upload_mode_without_testing(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    namespace["_choose_colab_credential_payload"] = lambda: (_ for _ in ()).throw(
        FileNotFoundError("selection cancelled")
    )
    namespace["idt_upload_mode_button"].click()
    assert namespace["idt_setup_mode_widget"].value == "upload"
    assert "selection cancelled" in namespace["credential_status"].value
    assert namespace["state"]["pending_credential_upload"] is None


def test_colab_idt_mode_row_has_four_actions_and_no_duplicate_upload_or_test_button(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    row = namespace["idt_mode_action_row"]
    assert row.children == (
        namespace["idt_create_mode_button"],
        namespace["idt_upload_mode_button"],
        namespace["idt_manual_mode_button"],
        namespace["idt_batch_mode_button"],
    )
    assert namespace["idt_upload_mode_button"].description == "Upload idt.env"
    assert namespace["idt_upload_mode_button"].icon == "upload"
    assert len(namespace["credential_upload_panel"].children) == 2
    assert namespace["credential_registration_help"] not in namespace["credential_upload_panel"].children
    assert "credential_upload_test_button" not in namespace
    assert "credential_test_button" not in namespace


def test_colab_shared_enzyme_and_plasmid_select_all_none_controls(colab_runtime_namespace):
    namespace = colab_runtime_namespace
    assert len(namespace["all_enzyme_options"]) == 47
    assert len(namespace["enzyme_checkboxes"]) == 47
    assert all(box.value for box in namespace["enzyme_checkboxes"].values())
    assert namespace["enzyme_bulk_control"].value == "all"

    first_enzyme = namespace["all_enzyme_options"][0]
    namespace["enzyme_checkboxes"][first_enzyme].value = False
    assert namespace["enzyme_bulk_control"].value == "custom"
    assert first_enzyme not in namespace["_selected_role_enzymes"]("site_i") + namespace["_selected_role_enzymes"]("site_ii") + namespace["_selected_role_enzymes"]("site_iii")

    namespace["_set_no_enzymes"]()
    assert not any(box.value for box in namespace["enzyme_checkboxes"].values())
    assert namespace["enzyme_bulk_control"].value == "none"
    with pytest.raises(ValueError, match="Site I"):
        namespace["_current_query"]()

    namespace["enzyme_bulk_control"].value = "all"
    assert all(box.value for box in namespace["enzyme_checkboxes"].values())
    assert namespace["enzyme_bulk_control"].value == "all"
    namespace["enzyme_bulk_control"].value = "custom"
    assert namespace["enzyme_bulk_control"].value == "all"

    first_plasmid = namespace["PLASMID_OPTIONS"][0]
    namespace["plasmid_checkboxes"][first_plasmid].value = False
    assert namespace["plasmid_bulk_control"].value == "custom"
    assert first_plasmid not in namespace["_selected_plasmids"]()

    namespace["plasmid_bulk_control"].value = "none"
    assert not any(box.value for box in namespace["plasmid_checkboxes"].values())
    assert namespace["plasmid_bulk_control"].value == "none"
    with pytest.raises(ValueError, match="plasmid"):
        namespace["_current_query"]()

    namespace["_set_all_plasmids"]()
    assert all(box.value for box in namespace["plasmid_checkboxes"].values())
    assert namespace["_selected_plasmids"]() == namespace["PLASMID_OPTIONS"]
    assert namespace["plasmid_bulk_control"].value == "all"


def test_colab_live_support_cards_match_joint_filtered_routes(colab_runtime_namespace):
    namespace = colab_runtime_namespace
    result = namespace["state"]["query_result"]
    routes = result.vector_routes
    summary = namespace["enzyme_route_support"].value
    assert summary == namespace["plasmid_route_support"].value
    assert "Supported RE pairs:" in summary
    assert f">{len({(row['site_i_enzyme'], row['site_ii_enzyme']) for row in routes}):,}</div>" in summary
    assert "Available Site III:" in summary
    assert f">{len({enzyme for row in routes for enzyme in row['site_iii_options']}):,}</div>" in summary
    assert "Supported plasmids:" in summary
    assert f">{len({row['profile_id'] for row in routes}):,}</div>" in summary
    assert "Minimum restore:" in summary
    assert f">{min(row['restoration_length_bp'] for row in routes)} bp</div>" in summary


def test_colab_restore_limit_refilters_cache_and_invalidates_confirmation(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    universe = namespace["state"]["route_universe"]
    assert namespace["state"]["confirmed_route"] is not None

    namespace["max_restoration_length_widget"].value = 0
    result = namespace["state"]["query_result"]
    assert namespace["state"]["route_universe"] is universe
    assert namespace["state"]["confirmed_route"] is None
    assert namespace["ga_module"].layout.display in (None, "")
    assert namespace["design_button"].disabled is True
    assert result.request["max_restoration_length_bp"] == 0
    assert result.vector_routes
    assert all(row["restoration_length_bp"] == 0 for row in result.vector_routes)


def test_colab_joint_re_pair_and_unsupported_plasmid_reports_zero(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    namespace["_set_no_enzymes"]()
    for enzyme in ("AflII", "BglII", "BbsI"):
        namespace["enzyme_checkboxes"][enzyme].value = True
    namespace["_set_no_plasmids"]()
    namespace["plasmid_checkboxes"]["pET-21a(+)"].value = True

    assert namespace["state"]["query_result"].status == "no_vector_route"
    summary = namespace["enzyme_route_support"].value
    assert summary == namespace["plasmid_route_support"].value
    assert "Supported RE pairs:" in summary and ">0</div>" in summary
    assert "Available Site III:" in summary
    assert "Supported plasmids:" in summary
    assert "Minimum restore:" in summary and ">—</div>" in summary


@pytest.mark.parametrize(
    ("module", "expected_status"),
    [("WWWWWW", "no_hurdler_pair_match"), ("YSPTSPS", "no_vector_route")],
)
def test_colab_incompatible_queries_never_enable_optimization(
    colab_runtime_namespace, module, expected_status
):
    namespace = colab_runtime_namespace
    namespace["repeat_module_widget"].value = module
    namespace["_run_query"]()
    assert namespace["state"]["query_result"].status == expected_status
    assert namespace["state"]["confirmed_route"] is None
    assert namespace["confirm_button"].disabled is True
    assert namespace["design_button"].disabled is True


def _confirm_first_colab_route(namespace):
    namespace["_run_query"]()
    namespace["pair_choice"].value = namespace["pair_choice"].options[1][1]
    namespace["site_iii_choice"].value = namespace["site_iii_choice"].options[1][1]
    namespace["profile_choice"].value = namespace["profile_choice"].options[1][1]
    namespace["scheme_choice"].value = namespace["scheme_choice"].options[1][1]
    namespace["_confirm_route"]()
    assert namespace["state"]["confirmed_route"] is not None


def _wait_for_colab_worker(namespace, timeout: float = 30.0):
    worker = namespace["state"].get("run_thread")
    assert worker is not None
    deadline = time.monotonic() + timeout
    while worker.is_alive() and time.monotonic() < deadline:
        namespace["_drain_ui_events"]()
        worker.join(timeout=0.05)
    namespace["_drain_ui_events"]()
    assert not worker.is_alive(), "Colab GA worker did not finish"


def test_colab_manual_route_batch_export_never_builds_an_idt_client(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    namespace["idt_setup_mode_widget"].value = "batch"
    namespace["population_number"].value = 4
    namespace["generation_schedule_widget"].value = "10,100"
    namespace["output_directory_widget"].value = str(tmp_path / "batch")
    namespace["auto_download_widget"].value = False
    namespace["storage_mode_widget"].value = "drive"
    namespace["drive_root_widget"].value = str(tmp_path / "mock_drive")
    namespace["storage_state"]["drive_mounted"] = True

    class MustNotBeConstructed:
        def __init__(self, *_args, **_kwargs):
            raise AssertionError("Batch mode constructed an IDT API client")

    namespace["IDTComplexityScorer"] = MustNotBeConstructed
    namespace["_run_design"]()
    assert "worker_started" in namespace["attempt_log_html"].value
    assert "No attempts yet" not in namespace["attempt_log_html"].value
    _wait_for_colab_worker(namespace)
    assert namespace["generation_progress"].value > 0
    summary = json.loads((tmp_path / "batch" / "design_summary.json").read_text())
    assert summary["status"] == "optimized_unvalidated_batch"
    assert summary["idt_audit"] == []
    assert summary["rdl_plan"]["secondary_repeat_copies"] >= 12
    assert summary["rdl_plan"]["minimum_secondary_satisfied"] is True
    assert summary["request"]["query"]["max_restoration_length_bp"] == 100
    manifest = json.loads((tmp_path / "batch" / "run_manifest.json").read_text())
    assert manifest["max_restoration_length_bp"] == 100
    assert (tmp_path / "batch" / "rdl_plan.json").is_file()
    assert (tmp_path / "batch" / "secondary_fragments.csv").is_file()
    assert (tmp_path / "batch" / "idt_bulk_input.csv").is_file()
    assert (tmp_path / "batch" / "idt_bulk_input.tsv").is_file()
    assert (tmp_path / "batch" / "idt_bulk_input.fasta").is_file()
    assert (tmp_path / "batch" / "ga_elite_candidates.csv").is_file()
    assert (tmp_path / "batch" / "ga_elite_candidates.fasta").is_file()
    assert (tmp_path / "batch" / "ga_parameter_history.csv").is_file()
    assert (tmp_path / "batch" / "idt_feedback_history.csv").is_file()
    assert (tmp_path / "batch" / "step00_plasmid.gb").is_file()
    assert (tmp_path / "batch" / "step01_insert.gb").is_file()
    assert (tmp_path / "batch" / "step01_plasmid.gb").is_file()
    assert (tmp_path / "batch" / "assembly_step_manifest.json").is_file()
    assert (tmp_path / "batch" / "step_translations.csv").is_file()
    assert Path(namespace["state"]["archive"]).is_file()
    assert re.fullmatch(
        r"hurdler_interactive_design_\d{8}T\d{6}Z_results\.zip",
        Path(namespace["state"]["archive"]).name,
    )
    assert (tmp_path / "mock_drive" / Path(namespace["state"]["archive"]).name).is_file()
    assert namespace["viewer_panel"].layout.display == ""
    final_step = max(value for _label, value in namespace["viewer_step_widget"].options)
    namespace["viewer_step_widget"].value = final_step
    namespace["viewer_molecule_widget"].value = "insert"
    assert namespace["viewer_view_widget"].value == "linear"
    assert namespace["viewer_view_widget"].disabled is True
    namespace["viewer_molecule_widget"].value = "plasmid"
    assert namespace["viewer_view_widget"].disabled is False
    namespace["_viewer_focus"]()
    assert namespace["viewer_view_widget"].value == "linear"
    assert namespace["viewer_range_widget"].value != (
        0, namespace["viewer_range_widget"].max
    )


def test_colab_mock_secrets_api_flow_clears_credentials(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    namespace["idt_setup_mode_widget"].value = "manual"
    namespace["idt_auth_method_widget"].value = "access_token"
    namespace["population_number"].value = 4
    namespace["generation_schedule_widget"].value = "10,100"
    namespace["output_directory_widget"].value = str(tmp_path / "api")
    namespace["auto_download_widget"].value = False

    class PassingScorer:
        def score(self, name: str, sequence: str):
            sequence_sha = hashlib.sha256(sequence.encode()).hexdigest()
            return {
                "idt_status": "passed",
                "idt_explicit_pass": True,
                "idt_complexity_score": 0.0,
                "idt_score_complete": True,
                "idt_score_policy": "idt-rule-score-sum-lt10-v1",
                "idt_rule_details_json": "[]",
                "idt_positive_score_names_json": "[]",
                "idt_violation_names_json": "[]",
                "idt_scored_sequence_sha256": sequence_sha,
                "idt_response_sha256": hashlib.sha256(("response:" + sequence_sha).encode()).hexdigest(),
            }

    def configure_mock_secret():
        os.environ["IDT_ACCESS_TOKEN"] = "temporary-token"
        return {
            "credential_mode": "manual",
            "auth_method": "access_token",
            "required_fields_complete": True,
        }

    namespace["_configure_api_credentials"] = configure_mock_secret
    namespace["IDTComplexityScorer"] = lambda _audit_path: PassingScorer()
    namespace["_run_design"]()
    _wait_for_colab_worker(namespace)
    summary = json.loads((tmp_path / "api" / "design_summary.json").read_text())
    assert summary["status"] == "idt_accepted"
    assert summary["idt_audit"]
    assert summary["rdl_plan"]["final_copy_count_exact"] is True
    assert (tmp_path / "api" / "rdl_plan.json").is_file()
    assert (tmp_path / "api" / "idt_bulk_input.csv").is_file()
    assert (tmp_path / "api" / "idt_bulk_input.tsv").is_file()
    assert (tmp_path / "api" / "idt_bulk_input.fasta").is_file()
    checkpoint = Path(namespace["state"]["checkpoint_archive"])
    assert checkpoint.is_file()
    with zipfile.ZipFile(checkpoint) as archive:
        assert "checkpoint.json" in archive.namelist()
        assert "best_secondary_core.fasta" in archive.namelist()
        assert "best_secondary_purchase.fasta" in archive.namelist()
        assert "temporary-token" not in b"".join(archive.read(name) for name in archive.namelist()).decode()
    assert "IDT_ACCESS_TOKEN" not in os.environ


def test_colab_pause_resume_and_stop_use_cooperative_background_control(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    namespace["idt_setup_mode_widget"].value = "batch"
    namespace["output_directory_widget"].value = str(tmp_path / "stopped")
    namespace["auto_download_widget"].value = False
    started = threading.Event()
    iterations: list[int] = []

    def slow_design(_request, *, run_control, **_kwargs):
        started.set()
        for index in range(10_000):
            run_control.safe_point()
            iterations.append(index)
            time.sleep(0.002)
        raise AssertionError("test run should have been stopped")

    namespace["design_construct_v2"] = slow_design
    namespace["_run_design"]()
    assert started.wait(2)
    namespace["_pause_or_resume"]()
    time.sleep(0.08)
    paused_count = len(iterations)
    time.sleep(0.08)
    assert len(iterations) == paused_count
    assert namespace["pause_button"].description == "Resume GA"
    namespace["_pause_or_resume"]()
    deadline = time.time() + 2
    while len(iterations) == paused_count and time.time() < deadline:
        time.sleep(0.01)
    assert len(iterations) > paused_count
    namespace["_stop_design"]()
    _wait_for_colab_worker(namespace)
    assert namespace["state"]["run_terminal_status"] == "stopped_by_user"
    assert namespace["state"]["run_active"] is False
    checkpoint = Path(namespace["state"]["checkpoint_archive"])
    assert checkpoint.is_file()
    with zipfile.ZipFile(checkpoint) as archive:
        payload = json.loads(archive.read("checkpoint.json"))
    assert payload["run_status"] == "stopped_by_user"
    assert "best_secondary_core.fasta" not in zipfile.ZipFile(checkpoint).namelist()


def test_colab_kernel_loop_automatically_drains_ga_feedback(tmp_path):
    """Exercise the real ipykernel loop without manually draining the UI queue."""
    notebook = nbformat.read(COLAB_NOTEBOOK, as_version=4)
    output_directory = json.dumps(str(tmp_path / "kernel_feedback"))
    notebook.cells.append(
        nbformat.v4.new_code_cell(
            f'''_run_query()
pair_choice.value = pair_choice.options[1][1]
site_iii_choice.value = site_iii_choice.options[1][1]
profile_choice.value = profile_choice.options[1][1]
scheme_choice.value = scheme_choice.options[1][1]
_confirm_route()
idt_setup_mode_widget.value = "batch"
output_directory_widget.value = {output_directory}
auto_download_widget.value = False

def _unexpected_early_checkpoint(*_args, **_kwargs):
    raise AssertionError("ordinary progress must not write an immediate checkpoint")

write_secondary_checkpoint = _unexpected_early_checkpoint

def _slow_feedback_design(_request, *, progress_callback, run_control, **_kwargs):
    for generation in range(1, 4):
        run_control.safe_point()
        progress_callback(DesignProgressEvent(
            stage="ga", status="running", fragment_kind="kernel_feedback",
            copies=12, generation=generation, generations=3,
            ga_score=float(4 - generation), elapsed_seconds=generation * 0.1,
        ))
        time.sleep(0.15)
    raise RuntimeError("intentional kernel feedback terminator")

design_construct_v2 = _slow_feedback_design
_run_design()
assert "run_requested" in attempt_log_html.value
observed_generations = []
for _index in range(200):
    observed_generations.append(generation_progress.value)
    if state.get("run_terminal_status") == "failed":
        break
    await asyncio.sleep(0.05)
assert max(observed_generations) == 3
assert "worker_entered" in attempt_log_html.value
assert "gen=3/3" in attempt_log_html.value
assert state["run_terminal_status"] == "failed"
assert state["run_active"] is False
'''
        )
    )
    executed = NotebookClient(
        notebook,
        timeout=180,
        kernel_name="python3",
        resources={"metadata": {"path": str(ROOT)}},
    ).execute(cwd=str(ROOT))
    assert all(
        output.get("output_type") != "error"
        for cell in executed.cells
        for output in cell.get("outputs", [])
    )


def test_colab_frontend_callback_bridge_drains_worker_events(colab_runtime_namespace):
    namespace = colab_runtime_namespace

    class FakeColabOutput:
        def __init__(self):
            self.callbacks = {}

        def register_callback(self, name, callback):
            self.callbacks[name] = callback

    bridge = FakeColabOutput()
    namespace["COLAB_WIDGET_MANAGER_ENABLED"] = True
    namespace["colab_output"] = bridge
    assert namespace["_install_colab_ui_bridge"]() is True
    assert "hurdler.ui_drain" in bridge.callbacks
    assert namespace["state"]["ui_bridge_installed"] is True

    namespace["state"]["run_id"] = 17
    namespace["state"]["run_active"] = True
    event = namespace["DesignProgressEvent"](
        stage="ga", status="fitness_running", fragment_kind="bridge_test",
        candidate_index=1, candidate_total=4, elapsed_seconds=0.1,
    )
    namespace["state"]["progress_queue"].put(("progress", 17, event))
    result = bridge.callbacks["hurdler.ui_drain"]()
    assert result.data["handled"] == 1
    assert result.data["run_id"] == 17
    assert namespace["candidate_progress"].value == 1
    assert "fitness_running" in namespace["stage_html"].value
    assert namespace["state"]["ui_bridge_poll_count"] == 1


def test_readme_notebook_links_resolve():
    assert NOTEBOOK.is_file()
    assert COLAB_NOTEBOOK.is_file()
    assert (ROOT / "apps" / "hurdler_designer.py").is_file()
    assert (ROOT / "scripts" / "start_hurdler_web.sh").is_file()
    assert (ROOT / "notebooks" / "README.md").is_file()
    readme = (ROOT / "README.md").read_text()
    assert "notebooks/workflows/01_interactive_hurdler_designer.ipynb" in readme
    assert "./scripts/start_hurdler_web.sh" in readme
    assert "codespaces" not in readme.lower()


def test_bundled_index_checksum_manifest():
    root = bundled_index_dir()
    assert (root / "pattern_index.npz").stat().st_size > 3_000_000
    entries = [line.split(None, 1) for line in (root / "SHA256SUMS").read_text().splitlines()]
    assert len(entries) == 13
    for expected, relative in entries:
        path = root / relative
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


def test_design_construct_cli_writes_topology_bundle(tmp_path, capsys):
    output = tmp_path / "cli"
    assert main(
        [
            "design-construct",
            "--request",
            str(ROOT / "data" / "example_design_request.json"),
            "--output-dir",
            str(output),
            "--legacy-v1",
        ]
    ) == 0
    response = json.loads(capsys.readouterr().out)
    assert response["status"] == "hurdler_incompatible"
    assert (output / "design_summary.json").is_file()
    assert not (output / "optimized_construct.fasta").exists()


def test_unversioned_v1_request_requires_explicit_compatibility_flag(tmp_path):
    with pytest.raises(SystemExit, match="--legacy-v1"):
        main(
            [
                "design-construct",
                "--request",
                str(ROOT / "data" / "example_design_request.json"),
                "--output-dir",
                str(tmp_path / "rejected"),
            ]
        )


@pytest.mark.parametrize("notebook_path", [NOTEBOOK, COLAB_NOTEBOOK])
def test_notebook_headless_mock_executes_clean_kernel(
    tmp_path, monkeypatch, notebook_path
):
    monkeypatch.setenv("HURDLER_NOTEBOOK_SMOKE", "1")
    monkeypatch.setenv("HURDLER_NOTEBOOK_SMOKE_OUTPUT", str(tmp_path / "design"))
    notebook = nbformat.read(notebook_path, as_version=4)
    executed = NotebookClient(
        notebook,
        timeout=180,
        kernel_name="python3",
        resources={"metadata": {"path": str(ROOT)}},
    ).execute(cwd=str(ROOT))
    assert all(
        output.get("output_type") != "error"
        for cell in executed.cells
        for output in cell.get("outputs", [])
    )
    if notebook_path == COLAB_NOTEBOOK:
        tutorial_cells = [
            cell for cell in executed.cells
            if str(cell.get("id", "")).startswith("hurdler-step-")
        ]
        assert len(tutorial_cells) == 7
        for cell in tutorial_cells:
            assert any(
                "application/vnd.jupyter.widget-view+json" in output.get("data", {})
                for output in cell.get("outputs", [])
            ), f"{cell['id']} did not render an interactive widget"
        reconnect_cell = next(
            cell for cell in executed.cells
            if cell.get("id") == "hurdler-reconnect-interface"
        )
        assert any(
            "application/vnd.jupyter.widget-view+json" in output.get("data", {})
            for output in reconnect_cell.get("outputs", [])
        )
    assert (tmp_path / "design" / "optimized_construct.fasta").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.csv").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.tsv").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.fasta").stat().st_size > 0
    summary = json.loads((tmp_path / "design" / "design_summary.json").read_text())
    assert summary["status"] == "optimized_unvalidated_batch"
    assert summary["idt_audit"] == []
