#!/usr/bin/env python3
"""Build the one-row-per-module final summary and requested scatter plots."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from hurdler.ga_optimization import GA_RE_SITE_POLICY
from hurdler.periodicity import MODULE_SELECTION_POLICY


CAPS = (1800, 3000)


def dna_sha256(sequence: object) -> str:
    return (
        hashlib.sha256(sequence.encode()).hexdigest()
        if isinstance(sequence, str) and sequence
        else ""
    )


def selected_weights(row: pd.Series) -> str:
    if bool(row.get("final_passed")) and isinstance(
        row.get("ga_score_profile_json"), str
    ):
        return str(row["ga_score_profile_json"])
    raw = row.get("adaptive_search_trace_json")
    if not isinstance(raw, str):
        return ""
    trace = json.loads(raw)
    if not trace:
        return ""
    return str(trace[-1].get("ga_score_profile_after_idt_json", ""))


def terminal_rejection_reasons(row: pd.Series) -> str:
    raw = row.get("adaptive_search_trace_json")
    if not isinstance(raw, str):
        return "[]"
    trace = json.loads(raw)
    if not trace:
        return "[]"
    return str(trace[-1].get("idt_violation_names_json", "[]"))


def cap_record(row: pd.Series, cap: int) -> dict[str, object]:
    passed = bool(row.get("final_passed"))
    dna = row.get("dna_sequence") if passed else ""
    dna = dna if isinstance(dna, str) else ""
    return {
        f"cap{cap}_mathematical_max_copies": int(row["mathematical_max_copies"]),
        f"cap{cap}_max_orderable_module_copies": int(row["verified_max_copies"]) if passed else 0,
        f"cap{cap}_orderable_construct_found": passed,
        f"cap{cap}_boundary_proven": bool(row.get("adaptive_boundary_proven")),
        f"cap{cap}_boundary_evidence": str(row.get("adaptive_boundary_evidence", "")),
        f"cap{cap}_stop_reason": str(row.get("adaptive_stop_reason", "")),
        f"cap{cap}_selected_ga_weights_json": selected_weights(row),
        f"cap{cap}_ga_re_site_policy": str(row.get("ga_re_site_policy", "")),
        f"cap{cap}_idt_status": str(row.get("idt_status", "")),
        f"cap{cap}_idt_passed_dna_sequence": dna,
        f"cap{cap}_idt_passed_dna_length_bp": len(dna),
        f"cap{cap}_idt_passed_dna_sha256": dna_sha256(dna),
        f"cap{cap}_terminal_rejection_reasons_json": terminal_rejection_reasons(row),
        f"cap{cap}_plasmid": str(row.get("plasmid", "")),
        f"cap{cap}_direction": str(row.get("direction", "")),
        f"cap{cap}_site_i_position": row.get("site_i_position"),
        f"cap{cap}_site_ii_position": row.get("site_ii_position"),
        f"cap{cap}_site_i_enzyme": str(row.get("site_i_enzyme", "")),
        f"cap{cap}_site_ii_enzyme": str(row.get("site_ii_enzyme", "")),
        f"cap{cap}_site_iii_enzymes": str(row.get("site_iii_enzymes", "")),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo", type=Path, default=Path("/home/wendai/projects/hurdler/clone_repeat_protein")
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--hurdler-results", type=Path)
    parser.add_argument("--optimized-constructs", type=Path)
    args = parser.parse_args()
    repo = args.repo.resolve()
    catalog_path = args.catalog or (
        repo / "studies/hurdler_validation/step03_module_corpus/tables/module_catalog.parquet"
    )
    hurdler_path = args.hurdler_results or (
        repo / "studies/hurdler_validation/step04_module_optimization/tables/module_hurdler_results.parquet"
    )
    optimized_path = args.optimized_constructs or (
        repo / "studies/hurdler_validation/step04_module_optimization/tables/optimized_constructs.parquet"
    )
    catalog = pd.read_parquet(catalog_path)
    hurdler = pd.read_parquet(hurdler_path)
    optimized = pd.read_parquet(optimized_path)

    eligibility = (
        hurdler.groupby("module_id")
        .agg(
            hurdler_usable=("success", "any"),
            hurdler_successful_plasmid_count=("success", "sum"),
            hurdler_tested_plasmid_count=("plasmid", "nunique"),
            hurdler_successful_plasmids=(
                "plasmid",
                lambda values: ",".join(
                    sorted(
                        hurdler.loc[
                            values.index[hurdler.loc[values.index, "success"]], "plasmid"
                        ].astype(str)
                    )
                ),
            ),
        )
        .reset_index()
    )
    catalog = catalog.merge(eligibility, on="module_id", how="left", validate="one_to_one")
    catalog["hurdler_usable"] = catalog.hurdler_usable.fillna(False).astype(bool)
    catalog["hurdler_successful_plasmid_count"] = (
        catalog.hurdler_successful_plasmid_count.fillna(0).astype(int)
    )
    catalog["hurdler_tested_plasmid_count"] = (
        catalog.hurdler_tested_plasmid_count.fillna(0).astype(int)
    )
    catalog["hurdler_successful_plasmids"] = catalog.hurdler_successful_plasmids.fillna("")

    by_key = {
        (str(row.module_id), int(row.fragment_limit_bp)): row._asdict()
        for row in optimized.itertuples(index=False)
    }
    rows: list[dict[str, object]] = []
    for module in catalog.to_dict(orient="records"):
        record = dict(module)
        record["unit_length_aa"] = int(record.get("unit_length", len(record["unit_sequence"])))
        for cap in CAPS:
            value = by_key[(str(record["module_id"]), cap)]
            record.update(cap_record(pd.Series(value), cap))
        rows.append(record)
    summary = pd.DataFrame(rows).sort_values(["collection", "module_id"]).reset_index(drop=True)
    # Keep the canonical policy name unambiguous in the final deliverable.
    # Some historical HURDLER rows predate this field, so a merge can otherwise
    # leave a partially populated legacy column beside selected_module_policy.
    summary["module_selection_policy"] = summary["selected_module_policy"]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    summary.to_parquet(args.output_dir / "module_final_summary.parquet", index=False)
    summary.to_csv(args.output_dir / "module_final_summary.csv", index=False)

    cohorts = [
        summary.loc[summary.collection.eq("natural100")].assign(cohort="natural100"),
        summary.loc[
            summary.collection.eq("designed_all")
            & summary.in_designed_primary100.fillna(False)
        ].assign(cohort="designed_primary100"),
        summary.loc[summary.collection.eq("designed_all")].assign(cohort="designed_all"),
    ]
    comparison = pd.concat(cohorts, ignore_index=True)
    fractions = (
        comparison.groupby("cohort")
        .agg(
            modules=("module_id", "nunique"),
            hurdler_usable_modules=("hurdler_usable", "sum"),
        )
        .reset_index()
    )
    fractions["hurdler_usable_fraction"] = (
        fractions.hurdler_usable_modules / fractions.modules
    )
    fractions.to_csv(args.output_dir / "module_hurdler_usable_fraction.csv", index=False)

    orderability_rows = []
    for cap in CAPS:
        part = summary.loc[summary.hurdler_usable].copy()
        part["fragment_limit_bp"] = cap
        part["max_orderable_module_copies"] = part[
            f"cap{cap}_max_orderable_module_copies"
        ]
        part["module_type"] = part.collection.map(
            {"natural100": "Natural", "designed_all": "Designed"}
        )
        part["orderability_class"] = "repeat_construct_orderable"
        part.loc[
            part.max_orderable_module_copies.eq(1), "orderability_class"
        ] = "single_module_only_not_repeat_construct"
        part.loc[
            part.max_orderable_module_copies.eq(0), "orderability_class"
        ] = "no_idt_orderable_construct_found"
        orderability_rows.append(part)
    orderability = pd.concat(orderability_rows, ignore_index=True)
    orderability.to_csv(
        args.output_dir / "module_copy_orderability_status.csv", index=False
    )
    # A zero is a search outcome, not a physical zero-module construct, and a
    # single module is not a repeat construct.  Retain both in the status table
    # but do not misrepresent either as a maximum-repeat scatter point.
    scatter = orderability.loc[
        orderability.max_orderable_module_copies.ge(2)
    ].copy()
    scatter.to_csv(args.output_dir / "module_length_copy_scatter_data.csv", index=False)

    sns.set_theme(style="whitegrid")
    figure, axes = plt.subplots(1, 2, figsize=(14, 5.8), sharey=False, facecolor="white")
    palette = {"Natural": "#4B2E83", "Designed": "#E57200"}
    for axis, cap in zip(axes, CAPS, strict=True):
        data = scatter.loc[scatter.fragment_limit_bp.eq(cap)]
        cap_status = orderability.loc[orderability.fragment_limit_bp.eq(cap)]
        class_counts = cap_status.orderability_class.value_counts()
        sns.scatterplot(
            data=data,
            x="unit_length_aa",
            y="max_orderable_module_copies",
            hue="module_type",
            palette=palette,
            alpha=0.78,
            s=48,
            edgecolor="white",
            linewidth=0.35,
            ax=axis,
        )
        axis.set_title(
            f"{cap:,} bp cap: {len(data)} orderable repeat constructs\n"
            "single-module only: "
            f"{int(class_counts.get('single_module_only_not_repeat_construct', 0))}; "
            "no IDT-orderable construct: "
            f"{int(class_counts.get('no_idt_orderable_construct_found', 0))}"
        )
        axis.set_xlabel("Selected middle repeat module length (AA)")
        axis.set_ylabel("Maximum IDT-orderable module copies (at least 2)")
        axis.set_ylim(bottom=1.5)
    sns.despine()
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(
            args.figure_dir / f"module_length_vs_max_orderable_copies.{suffix}",
            dpi=300,
            facecolor="white",
        )
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7.5, 5.2), facecolor="white")
    sns.barplot(
        data=fractions,
        x="cohort",
        y="hurdler_usable_fraction",
        color="#4B2E83",
        ax=axis,
    )
    axis.set_ylim(0, 1)
    axis.set_xlabel("")
    axis.set_ylabel("Fraction usable with HURDLER")
    axis.tick_params(axis="x", rotation=18)
    sns.despine()
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(
            args.figure_dir / f"module_hurdler_usable_fraction.{suffix}",
            dpi=300,
            facecolor="white",
        )
    plt.close(figure)

    top_lines = []
    for collection, data in summary.groupby("collection"):
        top = data.sort_values(
            ["cap3000_max_orderable_module_copies", "module_id"],
            ascending=[False, True],
        ).iloc[0]
        top_lines.append(
            f"- {collection}: {top.module_id}, "
            f"{int(top.cap3000_max_orderable_module_copies)} copies at 3,000 bp"
        )
    fraction_markdown = [
        "| cohort | modules | HURDLER usable | fraction |",
        "|---|---:|---:|---:|",
        *[
            f"| {row.cohort} | {int(row.modules)} | "
            f"{int(row.hurdler_usable_modules)} | "
            f"{float(row.hurdler_usable_fraction):.4f} |"
            for row in fractions.itertuples(index=False)
        ],
    ]
    middle_difference_counts = {
        str(collection): int(
            data.selected_module_sequence.ne(data.first_module_sequence).sum()
        )
        for collection, data in summary.groupby("collection")
    }
    orderability_lines = []
    for cap in CAPS:
        counts = (
            orderability.loc[orderability.fragment_limit_bp.eq(cap)]
            .orderability_class.value_counts()
        )
        orderability_lines.append(
            f"- {cap:,} bp: "
            f"{int(counts.get('repeat_construct_orderable', 0))} repeat constructs "
            f"(>=2 modules), "
            f"{int(counts.get('single_module_only_not_repeat_construct', 0))} "
            "single-module only, "
            f"{int(counts.get('no_idt_orderable_construct_found', 0))} with no "
            "IDT-orderable construct found"
        )
    markdown = [
        "# Final module summary",
        "",
        "Every analyzed AA sequence is the real middle module from its repeat region (ties choose the earlier central module). HURDLER usability means at least one successful result across the eight maintained plasmids. Maximum copies require an IDT zero-violation DNA and either the mathematical fragment ceiling or an explicit failed local+IDT evaluation for the next copy at 100 generations.",
        "",
        *fraction_markdown,
        "",
        "## Middle-copy audit",
        "",
        *[
            f"- {collection}: {count} selected middle sequences differ from the first repeat copy"
            for collection, count in middle_difference_counts.items()
        ],
        "",
        "## IDT-orderable repeat constructs",
        "",
        *orderability_lines,
        "",
        "The main scatter includes only constructs containing at least two modules. Zero-copy search outcomes and single-module-only outcomes remain in `module_copy_orderability_status.csv`; neither is presented as a repeat-protein maximum.",
        "",
        "## Largest 3,000 bp orderable constructs",
        "",
        *top_lines,
        "",
        "The complete AA sequences, final GA weights, exact IDT-passed DNA sequences, HURDLER schemes, and boundary evidence are in `module_final_summary.parquet` and `module_final_summary.csv`.",
    ]
    (args.output_dir / "module_final_summary.md").write_text("\n".join(markdown) + "\n")

    validation = {
        "module_rows": len(summary),
        "natural_rows": int(summary.collection.eq("natural100").sum()),
        "designed_rows": int(summary.collection.eq("designed_all").sum()),
        "unique_module_ids": int(summary.module_id.nunique()),
        "all_unit_lengths_match": bool(
            summary.unit_length_aa.eq(summary.unit_sequence.str.len()).all()
        ),
        "module_selection_policy": MODULE_SELECTION_POLICY,
        "ga_re_site_policy": GA_RE_SITE_POLICY,
        "all_caps_use_soft_nonselected_re_site_policy": bool(
            all(
                row[f"cap{cap}_ga_re_site_policy"] == GA_RE_SITE_POLICY
                for _, row in summary.iterrows()
                for cap in CAPS
            )
        ),
        "all_selected_modules_are_middle_policy": bool(
            summary.selected_module_policy.eq(MODULE_SELECTION_POLICY).all()
            and summary.module_selection_policy.eq(MODULE_SELECTION_POLICY).all()
        ),
        "all_unit_sequences_are_selected_module": bool(
            summary.unit_sequence.eq(summary.selected_module_sequence).all()
            and summary.unit_start.eq(summary.selected_module_start).all()
            and summary.unit_end.eq(summary.selected_module_end).all()
        ),
        "selected_middle_differs_from_first_by_collection": middle_difference_counts,
        "hurdler_usable_construct_caps": int(summary.hurdler_usable.sum() * len(CAPS)),
        "repeat_construct_scatter_rows": len(scatter),
        "scatter_minimum_copies": int(scatter.max_orderable_module_copies.min()),
        "orderability_status_rows": len(orderability),
        "hurdler_usable_caps_with_maximum_boundary_proof": int(
            sum(
                bool(row[f"cap{cap}_boundary_proven"])
                for _, row in summary.loc[summary.hurdler_usable].iterrows()
                for cap in CAPS
            )
        ),
        "idt_passed_dna_hash_mismatches": int(
            sum(
                dna_sha256(row[f"cap{cap}_idt_passed_dna_sequence"])
                != row[f"cap{cap}_idt_passed_dna_sha256"]
                for _, row in summary.iterrows()
                for cap in CAPS
            )
        ),
        "passed": bool(
            len(summary) == len(catalog)
            and summary.module_id.nunique() == len(catalog)
            and summary.unit_length_aa.eq(summary.unit_sequence.str.len()).all()
            and summary.selected_module_policy.eq(MODULE_SELECTION_POLICY).all()
            and summary.module_selection_policy.eq(MODULE_SELECTION_POLICY).all()
            and summary.unit_sequence.eq(summary.selected_module_sequence).all()
            and len(orderability) == int(summary.hurdler_usable.sum() * len(CAPS))
            and scatter.max_orderable_module_copies.ge(2).all()
            and all(
                row[f"cap{cap}_ga_re_site_policy"] == GA_RE_SITE_POLICY
                for _, row in summary.iterrows()
                for cap in CAPS
            )
            and all(
                bool(row[f"cap{cap}_boundary_proven"])
                for _, row in summary.loc[summary.hurdler_usable].iterrows()
                for cap in CAPS
            )
        ),
    }
    (args.output_dir / "module_final_summary_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )
    if not validation["passed"]:
        raise RuntimeError(json.dumps(validation, indent=2))
    print(json.dumps({**validation, "fractions": fractions.to_dict(orient="records")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
