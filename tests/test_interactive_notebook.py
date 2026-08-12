from __future__ import annotations

import hashlib
import json
import os
import re
from datetime import datetime, timezone
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


@pytest.fixture
def colab_runtime_namespace():
    namespace: dict[str, object] = {}
    for cell_id in (
        "hurdler-protein-form",
        "hurdler-selector-policy-form",
        "hurdler-test-defaults",
        "hurdler-imports",
    ):
        exec(compile(_colab_cell_source(cell_id), cell_id, "exec"), namespace)
    namespace["display"] = lambda *_args, **_kwargs: None
    namespace["clear_output"] = lambda *_args, **_kwargs: None
    exec(
        compile(
            _colab_cell_source("hurdler-controller-v2"),
            "hurdler-controller-v2",
            "exec",
        ),
        namespace,
    )
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


def test_colab_cells_are_named_hidden_forms():
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
        assert cell.get("outputs", []) == []
        assert cell.get("execution_count") is None
    assert len(ids) == len(set(ids))


def test_colab_user_inputs_are_native_params_visible_before_execution():
    payload = _colab_payload()
    form_ids = {
        "hurdler-protein-form",
        "hurdler-selector-policy-form",
    }
    sources = {
        cell["metadata"]["id"]: "".join(cell["source"])
        for cell in payload["cells"]
        if cell.get("metadata", {}).get("id") in form_ids
    }
    assert set(sources) == form_ids
    expected_fields = {
        "input_mode", "sequence_id", "n_cap_aa", "repeat_module_aa",
        "initial_repeat_copies", "c_cap_aa", "full_protein_or_fasta",
        "repeat_region_start_1based", "repeat_region_end_1based", "repeat_period_aa",
        "allow_left_cutter_in_hurdler_pair", "allow_right_cutter_in_hurdler_pair",
    }
    observed = set()
    for source in sources.values():
        assert '{ display-mode: "form", run: "auto" }' in source.splitlines()[0]
        assert "#@markdown" in source
        observed.update(
            match.group(1)
            for match in re.finditer(r"^([A-Za-z_]\w*)\s*=.*#@param", source, re.MULTILINE)
        )
    assert observed == expected_fields
    assert "widgets.Text(" not in "\n".join(sources.values())
    assert "widgets.Password(" not in COLAB_NOTEBOOK.read_text()
    for secret in (
        "IDT_ACCESS_TOKEN", "IDT_CLIENT_ID", "IDT_CLIENT_SECRET",
        "IDT_USERNAME", "IDT_PASSWORD",
    ):
        assert not re.search(rf"^{secret}\s*=.*#@param", COLAB_NOTEBOOK.read_text(), re.MULTILINE)


def test_colab_has_separate_individual_re_plasmid_route_and_ga_cells():
    sources = {
        cell["metadata"]["id"]: "".join(cell["source"])
        for cell in _colab_payload()["cells"]
    }
    assert "Select all RE" in sources["hurdler-controller-v2"]
    assert "Select no RE" in sources["hurdler-controller-v2"]
    assert "Select all plasmids" in sources["hurdler-controller-v2"]
    assert "Select no plasmids" in sources["hurdler-controller-v2"]
    assert "hurdler-enzyme-selector" in sources
    assert "hurdler-plasmid-selector" in sources
    assert "hurdler-re-route-panel" in sources
    assert "hurdler-vector-route-panel" in sources
    assert "hurdler-ga-panel" in sources
    controller = sources["hurdler-controller-v2"]
    assert 'value="api"' in controller
    assert 'assembly_strategy="exact_reused_secondary_rdl"' in controller
    assert "widgets.BoundedIntText" in controller
    assert "widgets.BoundedFloatText" in controller
    assert '"Minimum secondary modules (N)", 12' in controller
    assert '"Maximum GA→IDT feedback rounds", 100' in controller
    assert '"GA generations per feedback round", 10' in controller
    assert '"Warm-start top candidates", 10' in controller
    assert "minimum_secondary_copies=int(minimum_secondary_number.value)" in controller
    assert "max_idt_feedback_rounds=int(feedback_round_number.value)" in controller
    assert "auto_adjust_ga_parameters_from_idt=bool(auto_parameter_feedback.value)" in controller
    assert "feedback={event.feedback_round" in controller
    assert 'ga_panel.layout.display = "none"' in controller
    assert "credential_upload.value = ()" not in controller
    assert "_clear_credential_upload()" in controller


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
    assert len(state["query_result"].vector_routes) == 3072
    assert state["confirmed_route"] is None
    assert colab_runtime_namespace["pair_choice"].value is None
    assert colab_runtime_namespace["site_iii_choice"].value is None
    assert colab_runtime_namespace["design_button"].disabled is True


def test_colab_complete_fasta_preserves_header_and_requires_boundary_confirmation(
    colab_runtime_namespace,
):
    namespace = colab_runtime_namespace
    unit = "ACDEFGHIKLMNPQRSTVWY"
    namespace["input_mode"] = "Complete exact protein / FASTA"
    namespace["sequence_id"] = ""
    namespace["full_protein_or_fasta"] = f">header_from_fasta details\nM{unit * 3}G"
    namespace["repeat_region_start_1based"] = 0
    namespace["repeat_region_end_1based"] = 0
    namespace["repeat_period_aa"] = 0
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
    assert namespace["ga_panel"].layout.display == ""
    namespace["n_cap_aa"] = "MM"
    namespace["_confirm_route"]()
    assert namespace["state"]["confirmed_route"] is None
    assert namespace["design_button"].disabled is True


def test_colab_secrets_prefers_access_token_and_clears_environment(colab_runtime_namespace):
    namespace = colab_runtime_namespace
    requested: list[str] = []

    def reader(name: str) -> str:
        requested.append(name)
        return "temporary-token" if name == "IDT_ACCESS_TOKEN" else "must-not-be-read"

    status = namespace["_configure_colab_secrets"](reader)
    assert status["auth_method"] == "access_token"
    assert requested == ["IDT_ACCESS_TOKEN"]
    assert os.environ["IDT_ACCESS_TOKEN"] == "temporary-token"
    clear_idt_secret_environment()
    assert "IDT_ACCESS_TOKEN" not in os.environ


def test_colab_uploaded_env_clears_read_only_value_by_replacing_widget(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    upload = namespace["credential_upload"]
    payload = b"IDT_ACCESS_TOKEN=temporary-upload-token\n"
    upload.value = ({
        "name": "idt.env",
        "type": "text/plain",
        "size": len(payload),
        "content": memoryview(payload),
        "last_modified": datetime.now(timezone.utc),
    },)
    value_trait = upload.traits()["value"]
    original_read_only = value_trait.read_only
    value_trait.read_only = True
    namespace["credential_source"].value = "auto"
    namespace["credential_path"].value = str(tmp_path / "not-present.env")
    try:
        status = namespace["_configure_api_credentials"]()
    finally:
        value_trait.read_only = original_read_only
    assert status["credential_mode"] == "upload"
    assert status["upload_retained"] is False
    assert namespace["credential_upload"] is not upload
    assert namespace["credential_upload"].value == ()
    assert namespace["credential_upload_row"].children[1] is namespace["credential_upload"]
    clear_idt_secret_environment()


def test_colab_shared_enzyme_and_plasmid_select_all_none_controls(colab_runtime_namespace):
    namespace = colab_runtime_namespace
    assert len(namespace["all_enzyme_options"]) == 47
    namespace["_set_no_enzymes"]()
    assert namespace["enzyme_selector"].value == ()
    with pytest.raises(ValueError, match="Site I"):
        namespace["_current_query"]()
    namespace["_set_all_enzymes"]()
    assert len(namespace["enzyme_selector"].value) == 47
    namespace["_set_no_plasmids"]()
    with pytest.raises(ValueError, match="plasmid"):
        namespace["_current_query"]()
    namespace["_set_all_plasmids"]()
    assert namespace["plasmid_selector"].value == namespace["PLASMID_OPTIONS"]


@pytest.mark.parametrize(
    ("module", "expected_status"),
    [("WWWWWW", "no_hurdler_pair_match"), ("YSPTSPS", "no_vector_route")],
)
def test_colab_incompatible_queries_never_enable_optimization(
    colab_runtime_namespace, module, expected_status
):
    namespace = colab_runtime_namespace
    namespace["repeat_module_aa"] = module
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


def test_colab_manual_route_batch_export_never_builds_an_idt_client(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    namespace["validation_mode_widget"].value = "batch"
    namespace["population_number"].value = 4
    namespace["generation_schedule_widget"].value = "10,100"
    namespace["output_directory_widget"].value = str(tmp_path / "batch")
    namespace["auto_download_widget"].value = False

    class MustNotBeConstructed:
        def __init__(self, *_args, **_kwargs):
            raise AssertionError("Batch mode constructed an IDT API client")

    namespace["IDTComplexityScorer"] = MustNotBeConstructed
    namespace["_run_design"]()
    summary = json.loads((tmp_path / "batch" / "design_summary.json").read_text())
    assert summary["status"] == "optimized_unvalidated_batch"
    assert summary["idt_audit"] == []
    assert summary["rdl_plan"]["secondary_repeat_copies"] >= 12
    assert summary["rdl_plan"]["minimum_secondary_satisfied"] is True
    assert (tmp_path / "batch" / "rdl_plan.json").is_file()
    assert (tmp_path / "batch" / "secondary_fragments.csv").is_file()
    assert (tmp_path / "batch" / "idt_bulk_input.csv").is_file()
    assert (tmp_path / "batch" / "idt_bulk_input.tsv").is_file()
    assert (tmp_path / "batch" / "idt_bulk_input.fasta").is_file()
    assert (tmp_path / "batch" / "ga_elite_candidates.csv").is_file()
    assert (tmp_path / "batch" / "ga_elite_candidates.fasta").is_file()
    assert (tmp_path / "batch" / "ga_parameter_history.csv").is_file()
    assert (tmp_path / "batch" / "idt_feedback_history.csv").is_file()
    assert Path(namespace["state"]["archive"]).is_file()


def test_colab_mock_secrets_api_flow_clears_credentials(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    namespace["validation_mode_widget"].value = "api"
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
        return namespace["_configure_colab_secrets"](
            lambda name: "temporary-token" if name == "IDT_ACCESS_TOKEN" else ""
        )

    namespace["_configure_api_credentials"] = configure_mock_secret
    namespace["IDTComplexityScorer"] = lambda _audit_path: PassingScorer()
    namespace["_run_design"]()
    summary = json.loads((tmp_path / "api" / "design_summary.json").read_text())
    assert summary["status"] == "idt_accepted"
    assert summary["idt_audit"]
    assert summary["rdl_plan"]["final_copy_count_exact"] is True
    assert (tmp_path / "api" / "rdl_plan.json").is_file()
    assert (tmp_path / "api" / "idt_bulk_input.csv").is_file()
    assert (tmp_path / "api" / "idt_bulk_input.tsv").is_file()
    assert (tmp_path / "api" / "idt_bulk_input.fasta").is_file()
    assert "IDT_ACCESS_TOKEN" not in os.environ


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
    assert (tmp_path / "design" / "optimized_construct.fasta").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.csv").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.tsv").stat().st_size > 0
    assert (tmp_path / "design" / "idt_bulk_input.fasta").stat().st_size > 0
    summary = json.loads((tmp_path / "design" / "design_summary.json").read_text())
    assert summary["status"] == "optimized_unvalidated_batch"
    assert summary["idt_audit"] == []
