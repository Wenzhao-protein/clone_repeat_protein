from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks/workflows/03_colab_exact_dna_hurdler_designer.ipynb"


def _payload():
    return json.loads(NOTEBOOK.read_text())


def test_exact_dna_colab_is_output_free_hidden_and_stably_named():
    payload = _payload()
    assert payload["nbformat"] == 4
    assert "widgets" not in payload.get("metadata", {})
    identifiers = []
    for cell in payload["cells"]:
        if cell["cell_type"] != "code":
            continue
        identifier = cell["metadata"]["id"]
        identifiers.append(identifier)
        assert cell["id"] == identifier
        assert cell["metadata"]["cellView"] == "form"
        assert cell["metadata"]["jupyter"]["source_hidden"] is True
        assert "".join(cell["source"]).startswith("#@title ")
        assert cell["outputs"] == []
        assert cell["execution_count"] is None
    assert len(identifiers) == len(set(identifiers))


def test_exact_dna_colab_native_forms_and_runtime_selectors_are_present():
    text = NOTEBOOK.read_text()
    sources = {
        cell["metadata"]["id"]: "".join(cell["source"])
        for cell in _payload()["cells"]
        if cell["cell_type"] == "code"
    }
    parameters = set(
        re.findall(
            r"^([A-Za-z_]\w*)\s*=.*#@param",
            sources["exact-dna-input-form"] + sources["exact-dna-search-form"],
            re.MULTILINE,
        )
    )
    assert parameters >= {
        "input_mode",
        "sequence_id",
        "repeat_unit",
        "optional_spacer",
        "repeat_copies",
        "complete_exact_dna_or_fasta",
        "max_purchase_bp",
        "max_search_states",
        "search_timeout_seconds",
    }
    selector = sources["exact-dna-re-plasmid-selection"]
    assert "Select all" in selector and "Select none" in selector
    assert "site_i_boxes" in selector and "site_ii_boxes" in selector
    assert "plasmid_boxes" in selector
    assert '"AflII"' in selector and '"ApaI"' in selector
    query = sources["exact-dna-query-and-confirm"]
    assert "run_query()" in query
    assert "Confirm selected route" in query
    assert "query_exact_dna" in query
    assert "pair_dropdown" in query and "plasmid_dropdown" in query
    assert "scheme_dropdown" in query and "route_dropdown" in query
    assert "Verified active/latent transitions" in query
    export = sources["exact-dna-idt-export"]
    assert 'value="none"' in export
    assert "Bulk Input" in export and "Live IDT API" in export
    assert "tempfile.TemporaryDirectory" in export
    assert 'destination / "idt_raw_audit.jsonl"' not in export
    assert "clear_idt_secret_environment()" in export
    assert "files.download" in export
    assert "IDT_CLIENT_SECRET=" not in text
    assert "IDT_PASSWORD=" not in text


def test_exact_dna_colab_generator_is_the_single_source_of_truth():
    generator = (ROOT / "scripts/generate_exact_dna_colab.py").read_text()
    assert "03_colab_exact_dna_hurdler_designer.ipynb" in generator
    assert "RF00059" in generator
    assert "cellView" in generator
    assert "source_hidden" in generator


def test_readme_links_exact_dna_colab_and_notebook_index():
    readme = (ROOT / "README.md").read_text()
    index = (ROOT / "notebooks/README.md").read_text()
    relative = "notebooks/workflows/03_colab_exact_dna_hurdler_designer.ipynb"
    assert relative in readme
    assert relative in index
    assert "colab.research.google.com" in readme
