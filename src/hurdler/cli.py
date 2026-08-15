"""Unified command-line interface for maintained HURDLER workflows."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

from .artifacts import ArtifactRegistry, registry_rows
from .design import (
    DesignRequest,
    bundled_index_dir,
    design_construct,
    write_design_outputs,
)
from .designed_inventory import build_designed_structure_inventory, validate_af3_outputs
from .ga_optimization import refine_construct_table
from .dna_assembly import (
    build_target_corpus,
    finalize_target_plans,
    plan_target_catalog,
    plot_dna_assembly_report,
    plot_graphical_abstract_panel,
)
from .complete_route import (
    finalize_complete_route_shards,
    plan_complete_route_catalog,
)
from .purchase_orderability import audit_complete_route_purchase_orderability
from .exact_dna_design import (
    ExactDNAQuery,
    ExactDNASelection,
    confirm_best_exact_dna_route,
    confirm_exact_dna_route,
    query_exact_dna,
    write_exact_dna_outputs,
)
from .dna_assembly_visualization import (
    plot_complete_production_report,
    write_production_figure_manifest,
)
from .idt import (
    IDT_CREDENTIAL_PATH,
    IDT_SCORE_POLICY,
    IDTComplexityScorer,
    clear_idt_secret_environment,
    configure_idt_credentials,
    load_idt_credentials,
)
from .design_artifacts import timestamped_results_archive, write_secondary_checkpoint
from .progress import DesignProgressEvent
from .index import PatternIndex, build_pattern_index
from .io import write_json_atomic
from .protein_index import ProteinPatternIndex, build_protein_pattern_index
from .production_bundle import (
    WORKFLOWS,
    ProductionBundleRequest,
    build_production_bundle,
    validate_production_bundle,
)
from .plasmid_reference import (
    build_plasmid_reference,
    bundled_plasmid_reference_path,
    load_plasmid_reference,
    validate_plasmid_reference,
)
from .vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    bundled_protein_index_dir,
    design_construct_v2,
    design_query,
    write_design_outputs_v2,
)
from .matching import materialize_best_solution, query_all_plasmids
from .modules import (
    fetch_natural_modules,
    merge_module_catalogs,
    parse_dhr_supplement,
    parse_fasta_modules,
    parse_module_manifest,
    parse_pdb_exact_modules,
    refine_module_boundaries,
)
from .optimization import optimize_module_catalog
from .module_experiments import (
    finalize_adaptive_copy_results,
    finalize_module_compatibility,
    plot_maximum_copy_scatter,
    prepare_adaptive_copy_inputs,
    run_module_compatibility,
)
from .module_results import PUBLIC_RESULT_FILENAME, export_module_results
from .paths import ProjectPaths
from .qc import legacy_qc
from .rate import run_success_rate
from .reference import build_reference_manifest
from .repeatsdb import (
    build_natural_corpus,
    finalize_natural_corpus,
    write_annotation_inventory,
)
from .short_screen import finalize_short_results, screen_short_shard, summarize_short_results
from .structural_repeats import (
    DEFAULT_FOLDSEEK,
    finalize_designed_catalog,
    infer_designed_catalog,
)


def _print(payload: object) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True, default=str))


def _default_paths() -> ProjectPaths:
    return ProjectPaths.discover()


def build_parser() -> argparse.ArgumentParser:
    paths = _default_paths()
    parser = argparse.ArgumentParser(prog="hurdler", description="HURDLER repeat-protein cloning toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    reference = subparsers.add_parser("reference", help="Reference-data operations")
    reference_sub = reference.add_subparsers(dest="reference_command", required=True)
    reference_build = reference_sub.add_parser("build", help="Validate and hash canonical reference data")
    reference_build.add_argument("--reference-dir", type=Path, default=paths.reference_output)
    reference_build.add_argument("--output", type=Path, default=paths.output / "artifacts" / "reference_manifest.json")

    lookup = subparsers.add_parser("lookup", help="Pattern-index operations")
    lookup_sub = lookup.add_subparsers(dest="lookup_command", required=True)
    lookup_build = lookup_sub.add_parser("build", help="Build the sparse legacy-compatible pattern index")
    lookup_build.add_argument("--rules", default="legacy-optimized-v1", choices=["legacy-optimized-v1"])
    lookup_build.add_argument("--input-dir", type=Path, default=paths.output)
    lookup_build.add_argument("--output-dir", type=Path, default=paths.output / "artifacts" / "legacy-optimized-v1")
    lookup_build.add_argument("--orthogonality", type=Path, default=paths.reference_output / "orthogonality.csv")
    lookup_protein = lookup_sub.add_parser(
        "protein-build", help="Build the plasmid-independent vector-aware protein index"
    )
    lookup_protein.add_argument("--input-dir", type=Path, default=paths.output)
    lookup_protein.add_argument("--output-dir", type=Path, default=bundled_protein_index_dir())
    lookup_protein.add_argument("--orthogonality", type=Path, default=paths.reference_output / "orthogonality.csv")

    plasmid_reference = subparsers.add_parser(
        "plasmid-reference", help="Build or validate the annotated seven-vector/eight-profile database"
    )
    plasmid_reference_sub = plasmid_reference.add_subparsers(dest="plasmid_reference_command", required=True)
    plasmid_reference_build = plasmid_reference_sub.add_parser("build")
    plasmid_reference_build.add_argument("--output", type=Path, default=bundled_plasmid_reference_path())
    plasmid_reference_build.add_argument("--without-ncbi", action="store_true")
    plasmid_reference_validate = plasmid_reference_sub.add_parser("validate")
    plasmid_reference_validate.add_argument("--input", type=Path, default=bundled_plasmid_reference_path())

    query = subparsers.add_parser("query", help="Query a repeat module against all plasmids")
    query.add_argument("--module", required=True)
    query.add_argument("--index-dir", type=Path, default=paths.output / "artifacts" / "legacy-optimized-v1")
    query.add_argument("--successful-only", action="store_true")

    short = subparsers.add_parser("screen-short", help="Exhaustively screen one short-motif shard")
    short.add_argument("--index-dir", type=Path, default=paths.output / "artifacts" / "legacy-optimized-v1")
    short.add_argument("--output-dir", type=Path, required=True)
    short.add_argument("--length", type=int, required=True)
    short.add_argument("--prefix", default="")
    short.add_argument("--summarize", nargs="*", type=Path)
    short.add_argument("--summary-output", type=Path)
    short.add_argument("--finalize", action="store_true", help="Validate and combine all 404 shards")

    rate = subparsers.add_parser("success-rate", help="Run deterministic legacy Monte Carlo analysis")
    rate.add_argument("--index-dir", type=Path, default=paths.output / "artifacts" / "legacy-optimized-v1")
    rate.add_argument("--output", type=Path, required=True)
    rate.add_argument("--min-length", type=int, default=7)
    rate.add_argument("--max-length", type=int, default=60)
    rate.add_argument("--tests", type=int, default=1000)
    rate.add_argument("--seed", type=int, default=42)

    curate = subparsers.add_parser("curate-modules", help="Curate natural or designed module sequences")
    curate.add_argument("--natural-output", type=Path)
    curate.add_argument("--natural-per-class", type=int, default=20)
    curate.add_argument("--all-repeatsdb", action="store_true")
    curate.add_argument("--one-per-protein", action="store_true")
    curate.add_argument("--natural-mappings-output", type=Path)
    curate.add_argument("--natural-exclusions-output", type=Path)
    curate.add_argument("--natural-cache-dir", type=Path)
    curate.add_argument("--natural-workers", type=int, default=12)
    curate.add_argument("--max-annotations", type=int, help="Smoke-only API row cap")
    curate.add_argument("--natural-shard-index", type=int, default=0)
    curate.add_argument("--natural-shard-count", type=int, default=1)
    curate.add_argument("--finalize-natural-mappings", nargs="*", type=Path)
    curate.add_argument("--finalize-natural-inventories", nargs="*", type=Path, default=[])
    curate.add_argument("--finalize-natural-exclusions", nargs="*", type=Path, default=[])
    curate.add_argument("--annotation-inventory", type=Path)
    curate.add_argument("--write-annotation-inventory", type=Path)
    curate.add_argument("--boundary-workers", type=int, default=1)
    curate.add_argument("--designed-fasta", nargs="*", type=Path)
    curate.add_argument("--designed-manifest", nargs="*", type=Path)
    curate.add_argument("--dhr-supplement", type=Path)
    curate.add_argument("--designed-pdb", nargs="*", type=Path)
    curate.add_argument("--designed-family", default="designed_repeat")
    curate.add_argument("--designed-evidence", default="C", choices=["A", "B", "C"])
    curate.add_argument("--designed-source-url", default="")
    curate.add_argument("--catalog-output", type=Path)

    merge_catalogs = subparsers.add_parser(
        "merge-module-catalogs",
        help="Merge finalized natural and strict-designed module catalogs",
    )
    merge_catalogs.add_argument("--input", nargs="+", type=Path, required=True)
    merge_catalogs.add_argument("--output", type=Path, required=True)

    boundaries = subparsers.add_parser(
        "infer-boundaries", help="Infer primitive repeat units from complete protein sequences"
    )
    boundaries.add_argument("--input", type=Path, required=True)
    boundaries.add_argument("--output", type=Path, required=True)
    boundaries.add_argument("--candidates-output", type=Path)
    boundaries.add_argument("--units-output", type=Path)
    boundaries.add_argument("--positions-output", type=Path)
    boundaries.add_argument("--workers", type=int, default=1)
    boundaries.add_argument("--fixed-threshold", type=float, default=0.8)

    designed_boundaries = subparsers.add_parser(
        "infer-designed-boundaries",
        help="Infer strict DSSP/Foldseek-supported designed repeat modules",
    )
    designed_boundaries.add_argument("--input", type=Path)
    designed_boundaries.add_argument("--output", type=Path, required=True)
    designed_boundaries.add_argument("--candidates-output", type=Path)
    designed_boundaries.add_argument("--units-output", type=Path)
    designed_boundaries.add_argument("--positions-output", type=Path)
    designed_boundaries.add_argument("--exclusions-output", type=Path)
    designed_boundaries.add_argument("--dssp-engine", default="biotite", choices=["biotite"])
    designed_boundaries.add_argument("--mkdssp", type=Path, default=Path("/home/wendai/.conda/envs/hurdler/bin/mkdssp"))
    designed_boundaries.add_argument("--foldseek", type=Path, default=DEFAULT_FOLDSEEK)
    designed_boundaries.add_argument("--mafft", default="mafft")
    designed_boundaries.add_argument("--shard-index", type=int, default=0)
    designed_boundaries.add_argument("--shard-count", type=int, default=1)
    designed_boundaries.add_argument("--finalize-mappings", nargs="*", type=Path)
    designed_boundaries.add_argument("--finalize-exclusions", nargs="*", type=Path, default=[])
    designed_boundaries.add_argument("--finalize-candidate-tables", nargs="*", type=Path, default=[])
    designed_boundaries.add_argument("--finalize-unit-tables", nargs="*", type=Path, default=[])
    designed_boundaries.add_argument("--finalize-position-tables", nargs="*", type=Path, default=[])

    designed_inventory = subparsers.add_parser(
        "designed-inventory",
        help="Match designed structures and prepare missing-only AF3 GPU tasks",
    )
    designed_inventory.add_argument("--catalog", type=Path, required=True)
    designed_inventory.add_argument("--structure-roots", nargs="+", type=Path, required=True)
    designed_inventory.add_argument("--output", type=Path, required=True)
    designed_inventory.add_argument("--af3-output-root", type=Path, required=True)
    designed_inventory.add_argument("--af3-task-file", type=Path, required=True)

    validate_structures = subparsers.add_parser(
        "validate-designed-structures",
        help="Validate exact AF3 sequences and confidence provenance",
    )
    validate_structures.add_argument("--inventory", type=Path, required=True)
    validate_structures.add_argument("--output", type=Path, required=True)

    compatibility = subparsers.add_parser(
        "module-compatibility",
        help="Run or finalize Stage-1 compatibility against all eight plasmids",
    )
    compatibility.add_argument("--catalog", type=Path)
    compatibility.add_argument("--index-dir", type=Path, default=paths.output / "artifacts" / "legacy-optimized-v1")
    compatibility.add_argument("--output-dir", type=Path, required=True)
    compatibility.add_argument("--shard-index", type=int, default=0)
    compatibility.add_argument("--shard-count", type=int, default=1)
    compatibility.add_argument("--finalize-summaries", nargs="*", type=Path)
    compatibility.add_argument("--finalize-candidates", nargs="*", type=Path, default=[])

    optimize = subparsers.add_parser("optimize-modules", help="Query and codon-diversify a module catalog")
    optimize.add_argument("--catalog", type=Path, required=True)
    optimize.add_argument("--index-dir", type=Path, default=paths.output / "artifacts" / "legacy-optimized-v1")
    optimize.add_argument("--output-dir", type=Path, required=True)
    optimize.add_argument("--fragment-limits", type=int, nargs="+", default=[1800, 3000])
    optimize.add_argument("--external-deduction-bp", type=int, default=0)
    optimize.add_argument("--codon-usage", type=Path, default=paths.reference_output / "codon_usage.csv")
    optimize.add_argument("--shard-index", type=int, default=0)
    optimize.add_argument("--shard-count", type=int, default=1)
    optimize.add_argument("--limit", type=int)
    optimize.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Candidate-optimization worker processes; request the same CPU count from Digs",
    )

    refine = subparsers.add_parser(
        "refine-ga",
        help="Genetically refine synonymous constructs with optional IDT orderability feedback",
    )
    refine.add_argument("--constructs", type=Path, required=True)
    refine.add_argument("--output-dir", type=Path, required=True)
    refine.add_argument("--codon-usage", type=Path, default=paths.reference_output / "codon_usage.csv")
    refine.add_argument(
        "--restriction-sites",
        type=Path,
        default=paths.reference_output / "restriction_enzyme.csv",
    )
    refine.add_argument("--shard-index", type=int, default=0)
    refine.add_argument("--shard-count", type=int, default=1)
    refine.add_argument("--population-size", type=int, default=16)
    refine.add_argument("--generations", type=int, default=20)
    refine.add_argument("--seed", type=int, default=42)
    refine.add_argument(
        "--adaptive-copy-search",
        action="store_true",
        help="Use a local+IDT binary probe, then add one module at a time and escalate feedback-guided GA through 100 generations",
    )
    refine.add_argument("--short-generations", type=int, default=10)
    refine.add_argument(
        "--generation-schedule",
        type=int,
        nargs="+",
        default=[10, 20, 40, 60, 80, 100],
        help="Escalating generation budgets for each one-module extension; 100 is always enforced",
    )
    refine.add_argument(
        "--use-idt",
        action="store_true",
        help="Score exact GA DNA with IDT; adaptive search uses violations to reweight GA but never adopts IDT-generated DNA",
    )

    adaptive = subparsers.add_parser(
        "adaptive-copy-search",
        help="Run Stage-2 selected-pair-clean GA search with live IDT scoring",
    )
    adaptive.add_argument("--constructs", type=Path)
    adaptive.add_argument(
        "--compatibility",
        type=Path,
        help="Stage-1 per-module table; prepares both capacities with its frozen selected pair",
    )
    adaptive.add_argument("--output-dir", type=Path, required=True)
    adaptive.add_argument("--codon-usage", type=Path, default=paths.reference_output / "codon_usage.csv")
    adaptive.add_argument("--restriction-sites", type=Path, default=paths.reference_output / "restriction_enzyme.csv")
    adaptive.add_argument("--shard-index", type=int, default=0)
    adaptive.add_argument("--shard-count", type=int, default=1)
    adaptive.add_argument("--population-size", type=int, default=16)
    adaptive.add_argument("--seed", type=int, default=42)
    adaptive.add_argument("--short-generations", type=int, default=10)
    adaptive.add_argument("--generation-schedule", type=int, nargs="+", default=[10, 20, 40, 60, 80, 100])
    adaptive.add_argument("--idt-policy", default=IDT_SCORE_POLICY, choices=[IDT_SCORE_POLICY])
    adaptive.add_argument(
        "--credential-path",
        type=Path,
        default=IDT_CREDENTIAL_PATH,
        help="External mode-600 IDT env file used by live adaptive scoring",
    )
    adaptive.add_argument(
        "--auth-method", choices=["password", "access_token"], default=None
    )
    adaptive.add_argument(
        "--idt-batch",
        action="store_true",
        help="Do not call IDT; emit optimized candidates that require later batch validation",
    )
    adaptive.add_argument("--scatter-output", type=Path)
    adaptive.add_argument("--fragment-limits", type=int, nargs="+", default=[1800, 3000])
    adaptive.add_argument("--external-deduction-bp", type=int, default=0)
    adaptive.add_argument(
        "--finalize-results",
        nargs="*",
        type=Path,
        help="Merge completed adaptive Parquet shards without making new API calls",
    )
    adaptive.add_argument(
        "--finalize-result-list",
        type=Path,
        help="Newline-delimited adaptive Parquet paths (avoids OS argument limits)",
    )
    adaptive.add_argument(
        "--idt-audits",
        nargs="*",
        type=Path,
        default=[],
        help="Per-shard IDT JSONL audits required to prove accepted constructs",
    )
    adaptive.add_argument(
        "--idt-audit-list",
        type=Path,
        help="Newline-delimited IDT JSONL paths (avoids OS argument limits)",
    )

    export_results = subparsers.add_parser(
        "export-module-results",
        help="Build the public one-row-per-module HURDLER/IDT result CSV",
    )
    export_results.add_argument("--catalog", type=Path, required=True)
    export_results.add_argument(
        "--source-mappings", nargs="+", type=Path, required=True
    )
    export_results.add_argument("--compatibility", type=Path, required=True)
    export_results.add_argument("--maximum-results", type=Path, required=True)
    export_results.add_argument(
        "--output",
        type=Path,
        default=paths.root / "data" / "results" / PUBLIC_RESULT_FILENAME,
    )
    export_results.add_argument(
        "--generated-at-utc",
        help="Optional fixed ISO-8601 timestamp for byte-reproducible exports",
    )

    dna_assembly = subparsers.add_parser(
        "dna-assembly",
        help="Plan exact arbitrary-DNA assembly through active/latent RE sites",
    )
    dna_sub = dna_assembly.add_subparsers(dest="dna_assembly_command", required=True)
    dna_corpus = dna_sub.add_parser(
        "build-corpus",
        help="Build the versioned real/synthetic exact-DNA target catalog",
    )
    dna_corpus.add_argument("--source-table", nargs="*", type=Path, default=[])
    dna_corpus.add_argument("--output", type=Path, required=True)
    dna_corpus.add_argument("--no-synthetic", action="store_true")
    dna_corpus.add_argument("--seed", type=int, default=42)

    dna_plan = dna_sub.add_parser(
        "plan",
        help="Plan one resumable target-catalog shard",
    )
    dna_plan.add_argument("--catalog", type=Path, required=True)
    dna_plan.add_argument("--reference-dir", type=Path, default=paths.reference_output)
    dna_plan.add_argument("--artifact-dir", type=Path, default=paths.output)
    dna_plan.add_argument("--output-dir", type=Path, required=True)
    dna_plan.add_argument("--shard-index", type=int, default=0)
    dna_plan.add_argument("--shard-count", type=int, default=1)
    dna_plan.add_argument("--use-idt", action="store_true")
    dna_plan.add_argument("--credential-mode", choices=["path", "manual"], default="path")
    dna_plan.add_argument("--credential-path", type=Path, default=IDT_CREDENTIAL_PATH)
    dna_plan.add_argument("--auth-method", choices=["password", "access_token"], default=None)

    dna_finalize = dna_sub.add_parser(
        "finalize",
        help="Merge shards, validate cardinality, and generate reviewer artifacts",
    )
    dna_finalize.add_argument("--shard-dir", nargs="*", type=Path, default=[])
    dna_finalize.add_argument(
        "--shard-dir-list",
        type=Path,
        help="Newline-delimited shard directories, for large production runs",
    )
    dna_finalize.add_argument("--output-dir", type=Path, required=True)
    dna_finalize.add_argument("--figure-dir", type=Path)

    dna_complete = dna_sub.add_parser(
        "plan-complete",
        help="Plan one element shard from a purchasable seed to exact repeat targets",
    )
    dna_complete.add_argument("--catalog", type=Path, required=True)
    dna_complete.add_argument("--reference-dir", type=Path, default=paths.reference_output)
    dna_complete.add_argument("--artifact-dir", type=Path, default=paths.output)
    dna_complete.add_argument("--output-dir", type=Path, required=True)
    dna_complete.add_argument("--shard-index", type=int, default=0)
    dna_complete.add_argument("--shard-count", type=int, default=1)
    dna_complete.add_argument("--limit-elements", type=int)
    dna_complete.add_argument("--use-idt", action="store_true")
    dna_complete.add_argument("--credential-mode", choices=["path", "manual"], default="path")
    dna_complete.add_argument("--credential-path", type=Path, default=IDT_CREDENTIAL_PATH)
    dna_complete.add_argument("--auth-method", choices=["password", "access_token"], default=None)

    dna_complete_finalize = dna_sub.add_parser(
        "finalize-complete",
        help="Finalize complete-route element shards and reviewer analysis tables",
    )
    dna_complete_finalize.add_argument("--shard-dir", nargs="*", type=Path, default=[])
    dna_complete_finalize.add_argument("--shard-dir-list", type=Path)
    dna_complete_finalize.add_argument("--output-dir", type=Path, required=True)
    dna_complete_finalize.add_argument("--expected-elements", type=int)
    dna_complete_finalize.add_argument("--expected-real-targets", type=int)
    dna_complete_finalize.add_argument("--figure-dir", type=Path)

    dna_purchase_audit = dna_sub.add_parser(
        "audit-purchases",
        help="Audit whether every component of each selected route is an orderable oligo pair or gBlock",
    )
    dna_purchase_audit.add_argument("--raw-root", type=Path, required=True)
    dna_purchase_audit.add_argument("--output-dir", type=Path, required=True)
    dna_purchase_audit.add_argument("--expected-shards", type=int)
    dna_purchase_audit.add_argument("--expected-routes", type=int)
    dna_purchase_audit.add_argument("--expected-elements", type=int)
    dna_purchase_audit.add_argument("--use-idt", action="store_true")
    dna_purchase_audit.add_argument(
        "--credential-mode", choices=["path", "manual"], default="path"
    )
    dna_purchase_audit.add_argument(
        "--credential-path", type=Path, default=IDT_CREDENTIAL_PATH
    )
    dna_purchase_audit.add_argument(
        "--auth-method", choices=["password", "access_token"], default=None
    )

    dna_interactive = dna_sub.add_parser(
        "interactive-design",
        help="Query or confirm one exact arbitrary-DNA/array HURDLER route",
    )
    dna_interactive.add_argument("--request", type=Path, required=True)
    dna_interactive.add_argument("--output-dir", type=Path, required=True)
    dna_interactive.add_argument(
        "--plasmid-reference", type=Path, default=bundled_plasmid_reference_path()
    )
    dna_interactive.add_argument(
        "--idt-credential-file", type=Path, default=IDT_CREDENTIAL_PATH
    )
    dna_interactive.add_argument(
        "--auth-method", choices=["password", "access_token"], default=None
    )

    design = subparsers.add_parser(
        "design-construct",
        help="Design one confirmed repeat-protein HURDLER construct",
    )
    design.add_argument("--request", type=Path, required=True)
    design.add_argument("--output-dir", type=Path, required=True)
    design.add_argument("--index-dir", type=Path, default=bundled_index_dir())
    design.add_argument(
        "--idt-credential-file",
        type=Path,
        default=IDT_CREDENTIAL_PATH,
        help="Repo-external mode-600 env file; used only for validation_mode=api",
    )
    design.add_argument("--auth-method", choices=["password", "access_token"], default=None)
    design.add_argument("--progress-jsonl", type=Path)
    design.add_argument("--checkpoint-zip", type=Path)
    design.add_argument("--checkpoint-interval-seconds", type=float, default=180.0)
    design.add_argument("--final-archive-dir", type=Path)
    design.add_argument("--fail-on-nonaccepted", action="store_true")
    design.add_argument(
        "--legacy-v1",
        action="store_true",
        help="Explicitly read the historical unversioned DesignRequest schema",
    )
    design.add_argument("--protein-index-dir", type=Path, default=bundled_protein_index_dir())
    design.add_argument("--plasmid-reference", type=Path, default=bundled_plasmid_reference_path())

    idt_preflight = subparsers.add_parser(
        "idt-preflight",
        help="Validate external IDT credentials and one fixed 125-bp complexity call",
    )
    idt_preflight.add_argument(
        "--idt-credential-file", type=Path, default=IDT_CREDENTIAL_PATH
    )
    idt_preflight.add_argument(
        "--auth-method", choices=["password", "access_token"], default=None
    )

    design_query_parser = subparsers.add_parser(
        "design-query", help="Enumerate protein RE pairs and annotation-safe vector cut routes"
    )
    design_query_parser.add_argument("--request", type=Path, required=True)
    design_query_parser.add_argument("--output", type=Path)
    design_query_parser.add_argument("--protein-index-dir", type=Path, default=bundled_protein_index_dir())
    design_query_parser.add_argument("--plasmid-reference", type=Path, default=bundled_plasmid_reference_path())

    web = subparsers.add_parser("web", help="Run the local Marimo HURDLER designer")
    web.add_argument("--host", default="127.0.0.1")
    web.add_argument("--port", type=int, default=2718)
    web.add_argument("--no-browser", action="store_true")

    validate = subparsers.add_parser("validate-run", help="Validate index and audit legacy artifacts")
    validate.add_argument("--source-dir", type=Path, default=paths.output)
    validate.add_argument("--index-dir", type=Path, default=paths.output / "artifacts" / "legacy-optimized-v1")
    validate.add_argument("--output", type=Path, default=paths.output / "artifacts" / "legacy_qc.json")

    artifacts = subparsers.add_parser("artifacts", help="List, fetch and verify V2 notebook artifacts")
    artifacts_sub = artifacts.add_subparsers(dest="artifacts_command", required=True)
    artifacts_list = artifacts_sub.add_parser("list")
    artifacts_list.add_argument("--registry", type=Path, default=paths.root / "data/artifact_registry_v2.json")
    artifacts_list.add_argument("--level", choices=["fixture", "snapshot", "compact_result", "production_raw"])
    artifacts_fetch = artifacts_sub.add_parser("fetch")
    artifacts_fetch.add_argument("artifact_id")
    artifacts_fetch.add_argument("--registry", type=Path, default=paths.root / "data/artifact_registry_v2.json")
    artifacts_fetch.add_argument("--output", type=Path)
    artifacts_fetch.add_argument("--allow-production-raw", action="store_true")
    artifacts_verify = artifacts_sub.add_parser("verify")
    artifacts_verify.add_argument("artifact_id")
    artifacts_verify.add_argument("--registry", type=Path, default=paths.root / "data/artifact_registry_v2.json")
    artifacts_verify.add_argument("--path", type=Path)

    production = subparsers.add_parser("production", help="Generate portable Digs production bundles")
    production_sub = production.add_subparsers(dest="production_command", required=True)
    production_sub.add_parser("list")
    production_build = production_sub.add_parser("bundle")
    production_build.add_argument("--request", type=Path, required=True)
    production_build.add_argument("--output-dir", type=Path, required=True)
    production_validate = production_sub.add_parser("validate-bundle")
    production_validate.add_argument("bundle", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    paths = _default_paths()
    if args.command == "reference" and args.reference_command == "build":
        _print(build_reference_manifest(args.reference_dir, args.output))
    elif args.command == "lookup" and args.lookup_command == "build":
        _print(build_pattern_index(args.input_dir, args.output_dir, orthogonality_path=args.orthogonality))
    elif args.command == "lookup" and args.lookup_command == "protein-build":
        _print(build_protein_pattern_index(args.input_dir, args.output_dir, orthogonality_path=args.orthogonality))
    elif args.command == "plasmid-reference":
        if args.plasmid_reference_command == "build":
            database = build_plasmid_reference(args.output, include_ncbi=not args.without_ncbi)
            _print(validate_plasmid_reference(database))
        else:
            _print(validate_plasmid_reference(load_plasmid_reference(args.input)))
    elif args.command == "query":
        index = PatternIndex.load(args.index_dir)
        rows = [materialize_best_solution(result, index) for result in query_all_plasmids(args.module, index)]
        if args.successful_only:
            rows = [row for row in rows if row["success"]]
        _print(rows)
    elif args.command == "screen-short":
        if args.finalize:
            _print(finalize_short_results(args.output_dir))
        elif args.summarize:
            if args.summary_output is None:
                raise SystemExit("--summary-output is required with --summarize")
            frame = summarize_short_results(args.summarize, args.summary_output)
            _print({"rows": len(frame), "output": str(args.summary_output)})
        else:
            _print(screen_short_shard(args.index_dir, args.output_dir, length=args.length, prefix=args.prefix))
    elif args.command == "success-rate":
        frame = run_success_rate(
            args.index_dir,
            args.output,
            min_length=args.min_length,
            max_length=args.max_length,
            tests_per_plasmid=args.tests,
            seed=args.seed,
        )
        _print({"rows": len(frame), "output": str(args.output)})
    elif args.command == "curate-modules":
        frames = []
        if args.write_annotation_inventory:
            frame = write_annotation_inventory(
                args.write_annotation_inventory,
                max_annotations=args.max_annotations,
            )
            _print(
                {
                    "annotation_rows": len(frame),
                    "output": str(args.write_annotation_inventory),
                }
            )
            return 0
        if args.finalize_natural_mappings:
            if args.natural_output is None:
                raise SystemExit("--natural-output is required when finalizing natural shards")
            frames.append(
                finalize_natural_corpus(
                    args.finalize_natural_mappings,
                    args.natural_output,
                    region_inventory_paths=args.finalize_natural_inventories,
                    exclusion_paths=args.finalize_natural_exclusions,
                    annotation_inventory_path=args.annotation_inventory,
                )
            )
        elif args.all_repeatsdb:
            if not args.one_per_protein:
                raise SystemExit("--all-repeatsdb requires --one-per-protein for this corpus")
            if args.natural_output is None:
                raise SystemExit("--natural-output is required with --all-repeatsdb")
            frames.append(
                build_natural_corpus(
                    args.natural_output,
                    mappings_path=args.natural_mappings_output,
                    exclusions_path=args.natural_exclusions_output,
                    cache_dir=args.natural_cache_dir,
                    workers=args.natural_workers,
                    max_annotations=args.max_annotations,
                    shard_index=args.natural_shard_index,
                    shard_count=args.natural_shard_count,
                    annotation_inventory_path=args.annotation_inventory,
                )
            )
        elif args.natural_output:
            frames.append(
                fetch_natural_modules(
                    args.natural_output,
                    per_class=args.natural_per_class,
                    boundary_workers=args.boundary_workers,
                )
            )
        if args.designed_fasta:
            frames.append(
                parse_fasta_modules(
                    args.designed_fasta,
                    family=args.designed_family,
                    evidence_tier=args.designed_evidence,
                    source_url=args.designed_source_url,
                )
            )
        if args.designed_manifest:
            frames.extend(parse_module_manifest(path) for path in args.designed_manifest)
        if args.dhr_supplement:
            frames.append(parse_dhr_supplement(args.dhr_supplement))
        if args.designed_pdb:
            frames.append(
                parse_pdb_exact_modules(
                    args.designed_pdb,
                    family=args.designed_family,
                    source_url=args.designed_source_url,
                    experimental_accessions=[f"THR{number}" for number in range(1, 13)],
                    structure_accessions=["THR1", "THR2", "THR3", "THR5", "THR6"],
                )
            )
        if args.catalog_output:
            if not frames:
                raise SystemExit("No module source was selected")
            refined = refine_module_boundaries(
                pd.concat(frames, ignore_index=True), workers=args.boundary_workers
            )
            frame = merge_module_catalogs([refined], args.catalog_output)
            _print({"rows": len(frame), "output": str(args.catalog_output)})
        else:
            _print({"collections": [len(frame) for frame in frames]})
    elif args.command == "merge-module-catalogs":
        frames = [
            pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
            for path in args.input
        ]
        frame = merge_module_catalogs(frames, args.output)
        _print(
            {
                "rows": len(frame),
                "collection_counts": frame.collection.value_counts().to_dict(),
                "output": str(args.output),
            }
        )
    elif args.command == "infer-boundaries":
        source = (
            pd.read_parquet(args.input)
            if args.input.suffix == ".parquet"
            else pd.read_csv(args.input)
        )
        frame = refine_module_boundaries(
            source,
            audit_path=args.output,
            candidates_path=args.candidates_output,
            unit_alignment_path=args.units_output,
            position_variability_path=args.positions_output,
            workers=args.workers,
            fixed_threshold=args.fixed_threshold,
        )
        _print(
            {
                "rows": len(frame),
                "refined": int(frame.boundary_refinement_status.eq("refined").sum()),
                "output": str(args.output),
            }
        )
    elif args.command == "infer-designed-boundaries":
        if args.finalize_mappings:
            frame = finalize_designed_catalog(
                args.finalize_mappings,
                args.finalize_exclusions,
                args.output,
                candidate_paths=args.finalize_candidate_tables,
                unit_paths=args.finalize_unit_tables,
                position_paths=args.finalize_position_tables,
            )
        else:
            missing = [
                name
                for name, value in (
                    ("--input", args.input),
                    ("--candidates-output", args.candidates_output),
                    ("--units-output", args.units_output),
                    ("--positions-output", args.positions_output),
                    ("--exclusions-output", args.exclusions_output),
                )
                if value is None
            ]
            if missing:
                raise SystemExit(
                    "Designed inference requires: " + ", ".join(missing)
                )
            frame = infer_designed_catalog(
                args.input,
                args.output,
                candidates_path=args.candidates_output,
                units_path=args.units_output,
                positions_path=args.positions_output,
                exclusions_path=args.exclusions_output,
                dssp_executable=args.mkdssp,
                foldseek_binary=args.foldseek,
                mafft_binary=args.mafft,
                shard_index=args.shard_index,
                shard_count=args.shard_count,
            )
        _print({"accepted_rows": len(frame), "output": str(args.output)})
    elif args.command == "designed-inventory":
        frame = build_designed_structure_inventory(
            args.catalog,
            args.structure_roots,
            args.output,
            af3_output_root=args.af3_output_root,
            af3_task_path=args.af3_task_file,
        )
        _print(
            {
                "rows": len(frame),
                "matched_structures": int(
                    frame.structure_inventory_status.eq(
                        "author_or_pdb_structure_available"
                    ).sum()
                ),
                "af3_missing": int(
                    frame.structure_inventory_status.eq(
                        "missing_structure_af3_requested"
                    ).sum()
                ),
                "output": str(args.output),
            }
        )
    elif args.command == "validate-designed-structures":
        frame = validate_af3_outputs(args.inventory, args.output)
        _print(
            {
                "rows": len(frame),
                "status_counts": frame.af3_validation_status.value_counts().to_dict(),
                "output": str(args.output),
            }
        )
    elif args.command == "module-compatibility":
        if args.finalize_summaries:
            per_module, candidates, binned = finalize_module_compatibility(
                args.finalize_summaries,
                args.finalize_candidates,
                args.output_dir,
            )
            _print({"module_rows": len(per_module), "candidate_rows": int(candidates.attrs.get("candidate_row_count", len(candidates))), "bin_rows": len(binned), "output_dir": str(args.output_dir)})
        else:
            if args.catalog is None:
                raise SystemExit("--catalog is required unless --finalize-summaries is used")
            summary, candidates = run_module_compatibility(
                args.catalog,
                args.index_dir,
                args.output_dir,
                shard_index=args.shard_index,
                shard_count=args.shard_count,
            )
            _print({"module_rows": len(summary), "candidate_rows": len(candidates), "output_dir": str(args.output_dir)})
    elif args.command == "optimize-modules":
        results, constructs = optimize_module_catalog(
            args.catalog,
            args.index_dir,
            args.output_dir,
            fragment_limits=tuple(args.fragment_limits),
            external_deduction_bp=args.external_deduction_bp,
            codon_usage_path=args.codon_usage,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            limit=args.limit,
            workers=args.workers,
        )
        _print({"module_result_rows": len(results), "construct_rows": len(constructs), "output_dir": str(args.output_dir)})
    elif args.command == "refine-ga":
        if args.use_idt:
            load_idt_credentials()
        refined = refine_construct_table(
            args.constructs,
            args.codon_usage,
            args.restriction_sites,
            args.output_dir,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            population_size=args.population_size,
            generations=args.generations,
            seed=args.seed,
            use_idt=args.use_idt,
            adaptive_copy_search_enabled=args.adaptive_copy_search,
            short_generations=args.short_generations,
            generation_schedule=tuple(args.generation_schedule),
        )
        _print(
            {
                "rows": len(refined),
                "ga_passed": int(refined.ga_status.eq("passed").sum()),
                "final_passed": int(refined.final_passed.sum()),
                "output_dir": str(args.output_dir),
            }
        )
    elif args.command == "adaptive-copy-search":
        if args.finalize_results or args.finalize_result_list:
            if args.compatibility is None:
                raise SystemExit(
                    "--compatibility is required with --finalize-results"
                )
            result_paths = list(args.finalize_results or [])
            if args.finalize_result_list:
                result_paths.extend(
                    Path(line.strip())
                    for line in args.finalize_result_list.read_text().splitlines()
                    if line.strip()
                )
            idt_audit_paths = list(args.idt_audits or [])
            if args.idt_audit_list:
                idt_audit_paths.extend(
                    Path(line.strip())
                    for line in args.idt_audit_list.read_text().splitlines()
                    if line.strip()
                )
            results, traces, summary = finalize_adaptive_copy_results(
                result_paths,
                args.compatibility,
                args.output_dir,
                idt_audit_paths=idt_audit_paths,
            )
            _print(
                {
                    "result_rows": len(results),
                    "trace_rows": len(traces),
                    "summary_rows": len(summary),
                    "output_dir": str(args.output_dir),
                }
            )
            return 0
        if not args.idt_batch:
            configure_idt_credentials(
                mode="path",
                path=args.credential_path,
                auth_method=args.auth_method,
                headless=True,
                include_path_in_status=False,
            )
        constructs_path = args.constructs
        if args.compatibility is not None:
            if constructs_path is not None:
                raise SystemExit("Use either --constructs or --compatibility, not both")
            constructs_path = args.output_dir / "adaptive_copy_inputs.parquet"
            prepare_adaptive_copy_inputs(
                args.compatibility,
                constructs_path,
                codon_usage_path=args.codon_usage,
                fragment_limits=tuple(args.fragment_limits),
                external_deduction_bp=args.external_deduction_bp,
            )
        if constructs_path is None:
            raise SystemExit("One of --constructs or --compatibility is required")
        refined = refine_construct_table(
            constructs_path,
            args.codon_usage,
            args.restriction_sites,
            args.output_dir,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            population_size=args.population_size,
            seed=args.seed,
            use_idt=not args.idt_batch,
            adaptive_copy_search_enabled=True,
            short_generations=args.short_generations,
            generation_schedule=tuple(args.generation_schedule),
        )
        if args.scatter_output:
            plot_maximum_copy_scatter(refined, args.scatter_output)
        _print({"rows": len(refined), "passed": int(refined.final_passed.sum()), "output_dir": str(args.output_dir)})
    elif args.command == "export-module-results":
        exported = export_module_results(
            args.catalog,
            args.source_mappings,
            args.compatibility,
            args.maximum_results,
            args.output,
            repository_root=paths.root,
            generated_at_utc=args.generated_at_utc,
        )
        _print(
            {
                "rows": len(exported),
                "natural_rows": int(exported.collection.eq("Natural").sum()),
                "designed_rows": int(exported.collection.eq("Designed").sum()),
                "hurdler_compatible": int(
                    exported.hurdler_compatible.astype(bool).sum()
                ),
                "output": str(args.output),
                "size_bytes": args.output.stat().st_size,
            }
        )
    elif args.command == "design-query":
        query_request = CompatibilityQuery.from_dict(json.loads(args.request.read_text()))
        query_result = design_query(
            query_request,
            protein_index_dir=args.protein_index_dir,
            plasmid_reference_path=args.plasmid_reference,
        )
        if args.output:
            write_json_atomic(query_result.to_dict(), args.output)
        _print(query_result.to_dict())
    elif args.command == "idt-preflight":
        credential_status = configure_idt_credentials(
            mode="path",
            path=args.idt_credential_file,
            auth_method=args.auth_method,
            headless=True,
            include_path_in_status=False,
        )
        with tempfile.TemporaryDirectory(prefix="hurdler_idt_preflight_") as temporary:
            scorer = IDTComplexityScorer(Path(temporary) / "idt_preflight_audit.jsonl")
            summary = scorer.score(
                "hurdler_preflight_125bp",
                "ACGT" * 31 + "A",
            )
        score = summary.get("idt_complexity_score")
        if (
            not isinstance(score, (int, float))
            or isinstance(score, bool)
            or not math.isfinite(float(score))
        ):
            raise RuntimeError(
                "IDT preflight response did not contain a finite numeric rule-score sum"
            )
        _print(
            {
                "status": "passed",
                "credential_mode": credential_status["credential_mode"],
                "auth_method": credential_status["auth_method"],
                "control_length_bp": 125,
                "numeric_score_received": True,
                "score_policy": IDT_SCORE_POLICY,
                "orderability_not_required_for_preflight": True,
            }
        )
        return 0
    elif args.command == "design-construct":
        payload = json.loads(args.request.read_text())
        if payload.get("schema_version") == DESIGN_SCHEMA_VERSION_V2:
            request_v2 = DesignRequestV2.from_dict(payload)
            scorer = None
            credential_status = None
            progress_path = args.progress_jsonl
            checkpoint_path = args.checkpoint_zip
            if args.checkpoint_interval_seconds <= 0:
                raise ValueError("--checkpoint-interval-seconds must be positive")
            if progress_path is not None:
                progress_path.parent.mkdir(parents=True, exist_ok=True)
                progress_path.touch(exist_ok=True)
            checkpoint_state: dict[str, object] = {
                "best": None,
                "last_write": 0.0,
                "tested_copies": set(),
                "status": "starting",
            }

            def persist_checkpoint(*, force: bool = False) -> None:
                if checkpoint_path is None:
                    return
                now = time.monotonic()
                if (
                    not force
                    and now - float(checkpoint_state["last_write"])
                    < float(args.checkpoint_interval_seconds)
                ):
                    return
                best = checkpoint_state["best"]
                payload = dict(best) if isinstance(best, dict) else {
                    "event": "heartbeat",
                    "sequence_id": request_v2.query.sequence_id,
                    "validation_mode": request_v2.validation_mode,
                    "accepted_secondary_available": False,
                    "status": checkpoint_state["status"],
                    "tested_lengths": sorted(checkpoint_state["tested_copies"]),
                    "failure_reason": "No live-IDT-accepted secondary has been obtained yet",
                }
                write_secondary_checkpoint(payload, checkpoint_path)
                checkpoint_state["last_write"] = now

            def progress_callback(event: DesignProgressEvent) -> None:
                checkpoint_state["status"] = f"{event.stage}:{event.status}"
                if event.copies is not None:
                    checkpoint_state["tested_copies"].add(int(event.copies))
                if progress_path is not None:
                    with progress_path.open("a") as handle:
                        handle.write(json.dumps(event.to_dict(), sort_keys=True) + "\n")
                persist_checkpoint(force=False)

            def checkpoint_callback(payload: dict) -> None:
                previous = checkpoint_state["best"]
                improved = previous is None or int(payload.get("repeat_copies", 0)) > int(
                    previous.get("repeat_copies", 0)
                )
                if improved:
                    checkpoint_state["best"] = dict(payload)
                persist_checkpoint(force=improved)

            if request_v2.validation_mode == "api":
                credential_status = configure_idt_credentials(
                    mode="path",
                    path=args.idt_credential_file,
                    auth_method=args.auth_method,
                    headless=True,
                    include_path_in_status=False,
                )
                scorer = IDTComplexityScorer(args.output_dir / "idt_audit.jsonl")
            result_v2 = design_construct_v2(
                request_v2,
                protein_index_dir=args.protein_index_dir,
                plasmid_reference_path=args.plasmid_reference,
                idt_scorer=scorer,
                progress_callback=progress_callback,
                checkpoint_callback=checkpoint_callback,
            )
            files_v2 = write_design_outputs_v2(result_v2, args.output_dir)
            persist_checkpoint(force=True)
            final_archive = None
            if args.final_archive_dir is not None:
                final_archive = timestamped_results_archive(
                    args.output_dir,
                    args.final_archive_dir,
                    sequence_id=request_v2.query.sequence_id,
                )
            _print(
                {
                    "schema_version": result_v2.schema_version,
                    "status": result_v2.status,
                    "message": result_v2.message,
                    "protein_candidate_count": len(result_v2.protein_candidates),
                    "vector_route_count": len(result_v2.vector_routes),
                    "output_files": files_v2,
                    "credential_status": credential_status,
                    "progress_jsonl": str(progress_path) if progress_path else None,
                    "checkpoint_zip": str(checkpoint_path) if checkpoint_path else None,
                    "final_archive": str(final_archive) if final_archive else None,
                }
            )
            accepted_statuses = {
                "idt_accepted",
                "optimized_unvalidated_batch",
                "compatible_unoptimized",
            }
            return int(args.fail_on_nonaccepted and result_v2.status not in accepted_statuses)
        if not args.legacy_v1:
            raise SystemExit(
                f"A v2 request must declare schema_version={DESIGN_SCHEMA_VERSION_V2!r}; "
                "use --legacy-v1 only for an intentional historical request"
            )
        request = DesignRequest.from_dict(payload)
        scorer = None
        credential_status = None
        if request.optimize:
            credential_status = configure_idt_credentials(
                mode="path",
                path=args.idt_credential_file,
                auth_method=args.auth_method,
                headless=True,
                include_path_in_status=False,
            )
            scorer = IDTComplexityScorer(args.output_dir / "idt_audit.jsonl")
        result = design_construct(
            request,
            index_dir=args.index_dir,
            idt_scorer=scorer,
        )
        files = write_design_outputs(result, args.output_dir)
        _print(
            {
                "status": result.status,
                "message": result.message,
                "candidate_count": len(result.candidates),
                "output_files": files,
                "credential_status": credential_status,
            }
        )
    elif args.command == "web":
        app = paths.root / "apps" / "hurdler_designer.py"
        command = [sys.executable, "-m", "marimo", "run", str(app), "--host", args.host, "--port", str(args.port)]
        if args.no_browser:
            command.append("--headless")
        return subprocess.call(command)
    elif args.command == "dna-assembly":
        if args.dna_assembly_command == "build-corpus":
            frame = build_target_corpus(
                args.output,
                source_tables=args.source_table,
                include_synthetic=not args.no_synthetic,
                seed=args.seed,
            )
            companion = args.output.with_suffix(
                ".csv" if args.output.suffix.lower() in {".parquet", ".pq"} else ".parquet"
            )
            if companion.suffix == ".csv":
                frame.to_csv(companion, index=False)
            else:
                frame.to_parquet(companion, index=False)
            _print({"rows": len(frame), "output": str(args.output), "companion": str(companion)})
        elif args.dna_assembly_command == "plan":
            scorer = None
            credential_status = None
            if args.use_idt:
                credential_status = configure_idt_credentials(
                    mode=args.credential_mode,
                    path=args.credential_path,
                    auth_method=args.auth_method,
                    headless=not sys.stdin.isatty(),
                )
                scorer = IDTComplexityScorer(args.output_dir / "idt_audit.jsonl")
            tables = plan_target_catalog(
                args.catalog,
                args.reference_dir,
                args.output_dir,
                artifact_dir=args.artifact_dir,
                idt_scorer=scorer,
                require_idt=args.use_idt,
                shard_index=args.shard_index,
                shard_count=args.shard_count,
            )
            _print(
                {
                    "rows": {name: len(frame) for name, frame in tables.items()},
                    "output_dir": str(args.output_dir),
                    "credential_status": credential_status,
                }
            )
        elif args.dna_assembly_command == "finalize":
            shard_dirs = list(args.shard_dir)
            if args.shard_dir_list:
                shard_dirs.extend(
                    Path(line.strip())
                    for line in args.shard_dir_list.read_text().splitlines()
                    if line.strip()
                )
            if not shard_dirs:
                raise SystemExit("At least one --shard-dir or --shard-dir-list entry is required")
            tables = finalize_target_plans(shard_dirs, args.output_dir)
            figure_dir = args.figure_dir or args.output_dir / "figures"
            figures = plot_dna_assembly_report(tables["summary"], figure_dir)
            figures.extend(plot_graphical_abstract_panel(figure_dir))
            _print(
                {
                    "rows": {name: len(frame) for name, frame in tables.items()},
                    "figures": [str(path) for path in figures],
                    "output_dir": str(args.output_dir),
                }
            )
        elif args.dna_assembly_command == "plan-complete":
            scorer = None
            credential_status = None
            if args.use_idt:
                credential_status = configure_idt_credentials(
                    mode=args.credential_mode,
                    path=args.credential_path,
                    auth_method=args.auth_method,
                    headless=not sys.stdin.isatty(),
                )
                scorer = IDTComplexityScorer(args.output_dir / "idt_audit.jsonl")
            tables = plan_complete_route_catalog(
                args.catalog,
                args.reference_dir,
                args.output_dir,
                artifact_dir=args.artifact_dir,
                idt_scorer=scorer,
                require_idt=args.use_idt,
                shard_index=args.shard_index,
                shard_count=args.shard_count,
                limit_elements=args.limit_elements,
            )
            _print(
                {
                    "rows": {name: len(frame) for name, frame in tables.items()},
                    "output_dir": str(args.output_dir),
                    "credential_status": credential_status,
                }
            )
        elif args.dna_assembly_command == "finalize-complete":
            shard_dirs = list(args.shard_dir)
            if args.shard_dir_list:
                shard_dirs.extend(
                    Path(line.strip())
                    for line in args.shard_dir_list.read_text().splitlines()
                    if line.strip()
                )
            tables = finalize_complete_route_shards(
                shard_dirs,
                args.output_dir,
                expected_public_elements=args.expected_elements,
                expected_real_targets=args.expected_real_targets,
            )
            figures = []
            if args.figure_dir:
                figures = plot_complete_production_report(
                    tables["targets"],
                    tables["element_matrix"],
                    tables["selected_routes"],
                    tables["transitions"],
                    tables["fragments"],
                    tables["seeds"],
                    args.figure_dir,
                )
                write_production_figure_manifest(
                    figures,
                    args.figure_dir / "figure_manifest.csv",
                    input_tables=[
                        args.output_dir / "production_target_analysis.parquet",
                        args.output_dir / "production_element_matrix.parquet",
                        args.output_dir / "production_selected_routes.parquet",
                    ],
                )
            _print(
                {
                    "rows": {name: len(frame) for name, frame in tables.items()},
                    "figures": [str(path) for path in figures],
                    "output_dir": str(args.output_dir),
                }
            )
        elif args.dna_assembly_command == "audit-purchases":
            scorer = None
            credential_status = None
            try:
                if args.use_idt:
                    credential_status = configure_idt_credentials(
                        mode=args.credential_mode,
                        path=args.credential_path,
                        auth_method=args.auth_method,
                        headless=not sys.stdin.isatty(),
                        include_path_in_status=False,
                    )
                    scorer = IDTComplexityScorer(args.output_dir / "idt_audit.jsonl")
                tables = audit_complete_route_purchase_orderability(
                    args.raw_root,
                    args.output_dir,
                    idt_scorer=scorer,
                    expected_shards=args.expected_shards,
                    expected_routes=args.expected_routes,
                    expected_elements=args.expected_elements,
                )
            finally:
                clear_idt_secret_environment()
            _print(
                {
                    "rows": {name: len(frame) for name, frame in tables.items()},
                    "output_dir": str(args.output_dir),
                    "credential_status": credential_status,
                }
            )
        elif args.dna_assembly_command == "interactive-design":
            payload = json.loads(args.request.read_text())
            query_payload = payload.get("query", payload)
            query = ExactDNAQuery.from_dict(query_payload)
            result = query_exact_dna(
                query,
                plasmid_reference_path=args.plasmid_reference,
            )
            credential_status = None
            try:
                if payload.get("selection"):
                    selection = ExactDNASelection(**payload["selection"])
                    scorer = None
                    with tempfile.TemporaryDirectory(prefix="hurdler-idt-audit-") as temporary:
                        if selection.validation_mode == "api":
                            credential_status = configure_idt_credentials(
                                mode="path",
                                path=args.idt_credential_file,
                                auth_method=args.auth_method,
                                headless=True,
                                include_path_in_status=False,
                            )
                            scorer = IDTComplexityScorer(Path(temporary) / "raw.jsonl")
                        confirmer = (
                            confirm_best_exact_dna_route
                            if selection.validation_mode == "api"
                            else confirm_exact_dna_route
                        )
                        result = confirmer(
                            result,
                            selection,
                            idt_scorer=scorer,
                            plasmid_reference_path=args.plasmid_reference,
                        )
            finally:
                clear_idt_secret_environment()
            files = write_exact_dna_outputs(result, args.output_dir)
            _print(
                {
                    "schema_version": result.schema_version,
                    "status": result.status,
                    "message": result.message,
                    "route_count": len(result.route_candidates),
                    "output_files": files,
                    "credential_status": credential_status,
                }
            )
        else:
            raise AssertionError(f"Unhandled dna-assembly command: {args.dna_assembly_command}")
    elif args.command == "artifacts":
        registry = ArtifactRegistry(args.registry)
        if args.artifacts_command == "list":
            _print(registry_rows(registry.list(level=args.level)))
        elif args.artifacts_command == "fetch":
            path = registry.fetch(
                args.artifact_id,
                args.output,
                allow_production_raw=args.allow_production_raw,
            )
            _print({"artifact_id": args.artifact_id, "path": str(path), "status": "verified"})
        elif args.artifacts_command == "verify":
            path = registry.verify(args.artifact_id, args.path)
            _print({"artifact_id": args.artifact_id, "path": str(path), "status": "verified"})
        else:  # pragma: no cover - argparse guards this
            raise AssertionError(f"Unhandled artifacts command: {args.artifacts_command}")
    elif args.command == "production":
        if args.production_command == "list":
            _print({key: value.__dict__ for key, value in WORKFLOWS.items()})
        elif args.production_command == "bundle":
            request = ProductionBundleRequest.from_dict(json.loads(args.request.read_text()))
            output = build_production_bundle(request, args.output_dir)
            _print({"bundle": str(output), **validate_production_bundle(output)})
        elif args.production_command == "validate-bundle":
            _print(validate_production_bundle(args.bundle))
        else:  # pragma: no cover - argparse guards this
            raise AssertionError(f"Unhandled production command: {args.production_command}")
    elif args.command == "validate-run":
        _print(legacy_qc(args.source_dir, args.index_dir, args.output))
    else:
        raise AssertionError(f"Unhandled command: {args.command}")
    return 0
