from __future__ import annotations

import hashlib
import json
import os
import re
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
        "hurdler-re-vector-form",
        "hurdler-optimization-form",
        "hurdler-test-defaults",
        "hurdler-imports",
    ):
        exec(compile(_colab_cell_source(cell_id), cell_id, "exec"), namespace)
    namespace["display"] = lambda *_args, **_kwargs: None
    namespace["clear_output"] = lambda *_args, **_kwargs: None
    exec(
        compile(
            _colab_cell_source("hurdler-query-and-design-controls"),
            "hurdler-query-and-design-controls",
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
    assert ".config/hurdler/idt.env" not in text


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
        "hurdler-re-vector-form",
        "hurdler-optimization-form",
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
        "site_i_allowlist", "site_ii_allowlist", "site_iii_allowlist",
        "use_pGEX_4T_1", "use_pMAL_c5X", "use_pET_21a", "use_pET_28a",
        "use_pET_28a_start_codon", "use_pCold_I", "use_pUC18", "use_pQE_3",
        "allow_left_cutter_in_hurdler_pair", "allow_right_cutter_in_hurdler_pair",
        "validation_mode", "idt_credential_source", "idt_prompt_auth_method",
        "output_directory", "enable_zip_download", "maximum_repeat_copies",
        "population_size", "mutation_rate", "crossover_rate", "elite_fraction",
        "random_seed", "generation_schedule", "auto_adjust_weights_from_idt",
        "weight_selected_re_site_excess", "weight_gc_window_violation",
        "weight_repeated_re_site_excess", "weight_repeated_14mer",
        "weight_repeated_13mer", "weight_repeated_8mer",
        "weight_hairpin_10mer_proxy", "weight_homopolymer_excess",
        "weight_terminal_repeat_proxy", "weight_gc_window_soft_violation",
        "weight_negative_log_cai",
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


def test_colab_run_all_queries_but_never_auto_confirms_rank_one(colab_runtime_namespace):
    state = colab_runtime_namespace["state"]
    assert state["query_result"].status == "compatible_unoptimized"
    assert state["query_result"].vector_routes
    assert state["confirmed_route"] is None
    assert colab_runtime_namespace["route_choice"].value is None
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
    namespace["route_choice"].value = 0
    namespace["_confirm_route"]()
    assert namespace["state"]["confirmed_route"] is not None
    assert namespace["design_button"].disabled is False
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
    namespace["route_choice"].value = 0
    namespace["_confirm_route"]()
    assert namespace["state"]["confirmed_route"] is not None


def test_colab_manual_route_batch_export_never_builds_an_idt_client(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    namespace["validation_mode"] = "IDT Bulk Input files"
    namespace["maximum_repeat_copies"] = 3
    namespace["population_size"] = 4
    namespace["generation_schedule"] = "10,100"
    namespace["output_directory"] = str(tmp_path / "batch")

    class MustNotBeConstructed:
        def __init__(self, *_args, **_kwargs):
            raise AssertionError("Batch mode constructed an IDT API client")

    namespace["IDTComplexityScorer"] = MustNotBeConstructed
    namespace["_run_design"]()
    summary = json.loads((tmp_path / "batch" / "design_summary.json").read_text())
    assert summary["status"] == "optimized_unvalidated_batch"
    assert summary["idt_audit"] == []
    assert Path(namespace["state"]["archive"]).is_file()


def test_colab_mock_secrets_api_flow_clears_credentials(
    tmp_path, colab_runtime_namespace
):
    namespace = colab_runtime_namespace
    _confirm_first_colab_route(namespace)
    namespace["validation_mode"] = "Live IDT API"
    namespace["maximum_repeat_copies"] = 3
    namespace["population_size"] = 4
    namespace["generation_schedule"] = "10,100"
    namespace["output_directory"] = str(tmp_path / "api")

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
