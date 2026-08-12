import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    from pathlib import Path

    import marimo as mo
    import pandas as pd

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

    return (
        CompatibilityQuery,
        DESIGN_SCHEMA_VERSION_V2,
        DesignRequestV2,
        DesignSelection,
        IDTComplexityScorer,
        Path,
        configure_idt_credentials,
        configure_idt_credentials_from_bytes,
        configure_idt_credentials_from_values,
        design_construct_v2,
        design_query,
        json,
        mo,
        pd,
        write_design_outputs_v2,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # Annotation-aware HURDLER designer v2

    Protein RE-pair matching is performed before plasmid filtering. Seven
    physical annotated vectors provide eight profiles and four MCS cut
    schemes each. This local page produces design files only—never orders.
    """)
    return


@app.cell
def _(mo):
    input_mode = mo.ui.dropdown({"N-cap + module + C-cap": "split", "Complete exact protein": "full"}, value="N-cap + module + C-cap", label="Input mode")
    n_cap = mo.ui.text(value="M", label="N-cap AA")
    module = mo.ui.text_area(value="ACDEFGHIKLMNPQRSTVWY", label="Repeat module AA")
    copies = mo.ui.number(value=3, start=2, stop=10000, label="Repeat copies")
    c_cap = mo.ui.text(value="G", label="C-cap AA")
    full_protein = mo.ui.text_area(label="Complete protein AA/FASTA")
    repeat_start = mo.ui.number(value=1, start=1, label="Repeat start (1-based)")
    repeat_end = mo.ui.number(value=0, start=0, label="Repeat end (inclusive; 0 = infer first)")
    repeat_period = mo.ui.number(value=0, start=0, label="Repeat period (0 = infer first)")
    site_i = mo.ui.text(label="Site-I whitelist (comma-separated)")
    site_ii = mo.ui.text(label="Site-II whitelist (comma-separated)")
    site_iii = mo.ui.text(label="Site-III whitelist (comma-separated)")
    profiles = mo.ui.multiselect(options=["pGEX-4T-1", "pMAL-c5X", "pET-21a(+)", "pET-28a(+)", "pET-28a(+)_start_codon", "pCold_I", "pUC18", "pQE-3"], label="Optional plasmid profiles")
    allow_left = mo.ui.checkbox(False, label="Allow left cutter in HURDLER pair as final fallback")
    allow_right = mo.ui.checkbox(False, label="Allow right cutter in HURDLER pair as final fallback")
    query_button = mo.ui.run_button(label="1. Query all RE pairs and annotated vectors")
    mo.vstack([
        input_mode,
        mo.hstack([n_cap, copies, c_cap]),
        module,
        full_protein,
        mo.hstack([repeat_start, repeat_end, repeat_period]),
        mo.hstack([site_i, site_ii, site_iii]),
        profiles,
        mo.hstack([allow_left, allow_right]),
        query_button,
    ])
    return (
        allow_left,
        allow_right,
        c_cap,
        copies,
        full_protein,
        input_mode,
        module,
        n_cap,
        profiles,
        query_button,
        repeat_end,
        repeat_period,
        repeat_start,
        site_i,
        site_ii,
        site_iii,
    )


@app.cell
def _(
    CompatibilityQuery,
    DESIGN_SCHEMA_VERSION_V2,
    allow_left,
    allow_right,
    c_cap,
    copies,
    design_query,
    full_protein,
    input_mode,
    module,
    n_cap,
    profiles,
    query_button,
    repeat_end,
    repeat_period,
    repeat_start,
    site_i,
    site_ii,
    site_iii,
):
    def split_csv(value):
        return tuple(part.strip() for part in value.split(",") if part.strip())

    current_query = CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode=input_mode.value,
        sequence_id="marimo_design",
        n_cap=n_cap.value,
        repeat_module=module.value,
        repeat_copies=int(copies.value),
        c_cap=c_cap.value,
        full_protein_sequence=full_protein.value,
        repeat_region_start=int(repeat_start.value) if repeat_end.value and repeat_period.value else None,
        repeat_region_end=int(repeat_end.value) if repeat_end.value else None,
        repeat_period=int(repeat_period.value) if repeat_period.value else None,
        site_i_allowlist=split_csv(site_i.value),
        site_ii_allowlist=split_csv(site_ii.value),
        site_iii_allowlist=split_csv(site_iii.value),
        plasmid_allowlist=tuple(profiles.value),
        allow_left_cutter_in_hurdler_pair=allow_left.value,
        allow_right_cutter_in_hurdler_pair=allow_right.value,
    )
    query_result = design_query(current_query) if query_button.value else None
    return current_query, query_result


@app.cell
def _(mo, pd, query_result):
    mo.stop(query_result is None)
    route_table = pd.DataFrame(query_result.vector_routes)
    candidate_table = pd.DataFrame(query_result.protein_candidates)
    route_options = {
        f"#{row['rank']} {row['site_i_enzyme']}/{row['site_ii_enzyme']} → {row['profile_id']} {row['cut_scheme']}": index
        for index, row in enumerate(query_result.vector_routes)
    }
    route_choice = mo.ui.dropdown(route_options, value=next(iter(route_options), None), label="Confirmed route") if route_options else None
    mo.vstack([
        mo.md(f"**{query_result.status}** — {query_result.message}"),
        mo.ui.table(candidate_table.head(200)) if not candidate_table.empty else mo.md("No protein pair."),
        mo.ui.table(route_table.head(500)) if not route_table.empty else mo.md("No vector route."),
        route_choice if route_choice is not None else mo.md("Resolve the query before optimization."),
    ])
    return (route_choice,)


@app.cell
def _(json, mo, query_result, route_choice):
    mo.stop(query_result is None or not query_result.vector_routes or route_choice is None)
    validation = mo.ui.dropdown({"No optimization": "none", "Live IDT API score": "api", "IDT Bulk Input files": "batch"}, value="No optimization", label="Validation")
    credential_mode = mo.ui.dropdown({"External mode-600 env file": "path", "Manual OAuth": "manual", "Temporary env upload": "upload"}, value="External mode-600 env file", label="Credentials")
    credential_path = mo.ui.text(label="External credential path")
    client_id = mo.ui.text(kind="password", label="Client ID")
    client_secret = mo.ui.text(kind="password", label="Client secret")
    username = mo.ui.text(kind="password", label="Username")
    password = mo.ui.text(kind="password", label="Password")
    access_token = mo.ui.text(kind="password", label="Access token")
    upload = mo.ui.file(filetypes=[".env"], multiple=False, label="Temporary env upload")
    population = mo.ui.slider(4, 256, value=16, step=4, label="GA population")
    max_copies = mo.ui.number(value=20, start=2, stop=10000, label="Maximum repeat copies / resource limit")
    mutation = mo.ui.slider(0.001, 0.5, value=0.08, step=0.001, label="Mutation")
    crossover = mo.ui.slider(0.0, 1.0, value=0.75, step=0.01, label="Crossover")
    elite = mo.ui.slider(0.01, 0.5, value=0.15, step=0.01, label="Elite fraction")
    seed = mo.ui.number(value=42, label="Seed")
    auto_feedback = mo.ui.checkbox(True, label="Auto-adjust GA weights from IDT positive rules")
    weight_json = mo.ui.text_area(value=json.dumps({"selected_re_site_excess": 1e9, "repeated_re_site_excess": 1e4, "gc_window_violation": 1e9, "gc_window_soft_violation": 100, "repeated_8mer": 5, "repeated_13mer": 100, "repeated_14mer": 250, "hairpin_10mer_proxy": 25, "homopolymer_excess": 250, "terminal_repeat_proxy": 100, "negative_log_cai": 50}, sort_keys=True), label="GA score weights (JSON)")
    output_dir = mo.ui.text(value="output/marimo_vector_aware_design", label="Output directory")
    design_button = mo.ui.run_button(label="2. Optimize / export design files")
    mo.vstack([
        validation, credential_mode, credential_path,
        mo.hstack([client_id, client_secret]), mo.hstack([username, password]), access_token, upload,
        mo.hstack([population, max_copies, mutation]), mo.hstack([crossover, elite]), mo.hstack([seed, auto_feedback]),
        weight_json, output_dir, design_button,
    ])
    return (
        access_token,
        auto_feedback,
        client_id,
        client_secret,
        credential_mode,
        credential_path,
        crossover,
        design_button,
        elite,
        max_copies,
        mutation,
        output_dir,
        password,
        population,
        seed,
        upload,
        username,
        validation,
        weight_json,
    )


@app.cell
def _(
    DESIGN_SCHEMA_VERSION_V2,
    DesignRequestV2,
    DesignSelection,
    IDTComplexityScorer,
    Path,
    access_token,
    auto_feedback,
    client_id,
    client_secret,
    configure_idt_credentials,
    configure_idt_credentials_from_bytes,
    configure_idt_credentials_from_values,
    credential_mode,
    credential_path,
    crossover,
    current_query,
    design_button,
    design_construct_v2,
    elite,
    json,
    max_copies,
    mo,
    mutation,
    output_dir,
    password,
    population,
    query_result,
    route_choice,
    seed,
    upload,
    username,
    validation,
    weight_json,
    write_design_outputs_v2,
):
    mo.stop(not design_button.value)
    route = query_result.vector_routes[route_choice.value]
    if validation.value == "api":
        if credential_mode.value == "path":
            configure_idt_credentials(mode="path", path=credential_path.value, include_path_in_status=False)
        elif credential_mode.value == "manual":
            values = {"IDT_CLIENT_ID": client_id.value, "IDT_CLIENT_SECRET": client_secret.value, "IDT_USERNAME": username.value, "IDT_PASSWORD": password.value, "IDT_ACCESS_TOKEN": access_token.value}
            configure_idt_credentials_from_values(values)
            values.clear()
        else:
            if not upload.value:
                raise ValueError("Choose an env file")
            configure_idt_credentials_from_bytes(bytes(upload.value[0].contents))
    scorer = IDTComplexityScorer(Path(output_dir.value) / "idt_audit.jsonl") if validation.value == "api" else None
    request = DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=current_query,
        selection=DesignSelection(route["candidate_id"], route["profile_id"], route["scheme_id"], route["site_iii_options"][0]),
        validation_mode=validation.value,
        assembly_strategy="legacy_adaptive_max" if current_query.input_mode == "split" else "single_exact",
        max_repeat_copies=int(max_copies.value) if current_query.input_mode == "split" else None,
        population_size=int(population.value),
        mutation_rate=float(mutation.value),
        crossover_rate=float(crossover.value),
        elite_fraction=float(elite.value),
        seed=int(seed.value),
        generation_schedule=(10, 20, 40, 60, 80, 100),
        score_weights=json.loads(weight_json.value),
        auto_adjust_weights_from_idt=auto_feedback.value,
    )
    design_result = design_construct_v2(request, idt_scorer=scorer)
    exported = write_design_outputs_v2(design_result, output_dir.value)
    mo.vstack([mo.md(f"**{design_result.status}** — {design_result.message}"), mo.ui.table(design_result.primary_fragments), mo.ui.table(design_result.cloning_steps), mo.md(f"Files: `{exported}`")])
    return


if __name__ == "__main__":
    app.run()
