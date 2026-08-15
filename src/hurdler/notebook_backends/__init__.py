"""Authoritative Python backends for the HURDLER V2 notebook suite."""

from __future__ import annotations

import importlib


BACKEND_MODULES = {
    "01_reference_database": "reference_database",
    "02_lookup_plasmid": "lookup_plasmid",
    "03_module_corpus": "module_corpus",
    "04_query_batch": "query_batch",
    "05_repeat_designer": "repeat_designer",
    "06_exact_dna_designer": "exact_dna_designer",
    "07_production_builder": "production_builder",
    "08_success_landscape_analysis": "success_landscape_analysis",
    "09_module_result_analysis": "module_result_analysis",
    "10_exact_dna_result_analysis": "exact_dna_result_analysis",
    "11_reproducibility": "reproducibility",
    "12_agarose": "agarose",
    "13_sec": "sec",
    "14_plasmid_sequencing": "plasmid_sequencing",
}


def load_backend(notebook_id: str):
    try:
        name = BACKEND_MODULES[notebook_id]
    except KeyError as exc:
        raise KeyError(f"Unknown V2 notebook backend: {notebook_id}") from exc
    return importlib.import_module(f"hurdler.notebook_backends.{name}")


__all__ = ["BACKEND_MODULES", "load_backend"]
