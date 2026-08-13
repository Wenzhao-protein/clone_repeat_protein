#!/usr/bin/env python3
"""Generate the output-free native-Forms exact-DNA HURDLER Colab."""

from __future__ import annotations

from pathlib import Path
import textwrap

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "notebooks" / "workflows" / "03_colab_exact_dna_hurdler_designer.ipynb"


RF00059 = (
    "ATCATCCACTAGGGGGGCCTTTAGAAAGGCTGAGATCAAAGTGTGCCTTTGAGACCCTTAGCACCTGATCT"
    "GGGTAATGCCAGCGTAGGGAAGTGGAGGAGCAGCACA"
)


def markdown(source: str):
    return nbf.v4.new_markdown_cell(textwrap.dedent(source).strip())


def code(source: str, *, cell_id: str, title: str, tags: list[str] | None = None):
    cell = nbf.v4.new_code_cell(textwrap.dedent(source).strip())
    cell.id = cell_id
    cell.metadata.update(
        {
            "id": cell_id,
            "cellView": "form",
            "colab": {},
            "tags": tags or ["colab-form"],
            "jupyter": {"source_hidden": True},
        }
    )
    cell.source = f'#@title {title} {{ display-mode: "form" }}\n' + cell.source
    cell.outputs = []
    cell.execution_count = None
    return cell


def main() -> int:
    book = nbf.v4.new_notebook()
    book.metadata.update(
        {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3.11"},
        }
    )
    book.cells = [
        markdown(
            """
            # Exact-DNA HURDLER designer

            Design an exact arbitrary DNA insert or a regulatory-element array through
            active/one-base-latent restriction sites. The final insert is never altered.

            Fill the visible forms, then choose **Runtime → Run all**. The default run is
            offline and performs no IDT request. After the route table appears, confirm one
            route and optionally choose Live IDT or Bulk Input export.

            The default is a four-copy array of the 108-bp Rfam RF00059 TPP riboswitch
            element. Every eligible active/one-base-latent RE, every maintained
            Site-III adapter enzyme, and all eight plasmid profiles start selected.
            The array is a derived
            cloning example, not a claim that four copies occur naturally. The final
            route is deliberately left unselected for manual confirmation.
            """
        ),
        code(
            """
            repository_url = "https://github.com/Wenzhao-protein/clone_repeat_protein.git" #@param {type:"string"}
            repository_ref = "agent/vector-aware-designer-v2" #@param {type:"string"}
            force_fresh_clone = True #@param {type:"boolean"}

            import importlib, os, shutil, subprocess, sys
            from pathlib import Path

            checkout = Path("/content/clone_repeat_protein") if Path("/content").is_dir() else Path.cwd()
            hosted_checkout = checkout != Path.cwd()
            if force_fresh_clone and checkout.exists() and hosted_checkout:
                shutil.rmtree(checkout)
            if not (checkout / "pyproject.toml").is_file():
                subprocess.run(["git", "clone", "--depth", "1", "--branch", repository_ref, repository_url, str(checkout)], check=True)
            elif hosted_checkout:
                subprocess.run(["git", "-C", str(checkout), "fetch", "origin", repository_ref, "--depth", "1"], check=True)
                subprocess.run(["git", "-C", str(checkout), "checkout", "--force", "-B", repository_ref, "FETCH_HEAD"], check=True)
                subprocess.run(["git", "-C", str(checkout), "clean", "-fd"], check=True)
            if str(checkout / "src") not in sys.path:
                sys.path.insert(0, str(checkout / "src"))
            try:
                import google.colab  # noqa: F401
                in_colab = True
            except ModuleNotFoundError:
                in_colab = False
            if in_colab:
                subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-e", f"{checkout}[notebooks]"], check=True)
            # A repeated Run all in one runtime must not reuse modules imported
            # from an older checkout that has just been replaced above.
            for module_name in list(sys.modules):
                if module_name == "hurdler" or module_name.startswith("hurdler."):
                    del sys.modules[module_name]
            importlib.invalidate_caches()
            try:
                import hurdler
            except ModuleNotFoundError as exc:
                raise RuntimeError(
                    'Local execution requires the prepared environment; run `python -m pip install -e ".[notebooks]"` inside it.'
                ) from exc
            print(f"HURDLER {hurdler.__version__} loaded from {Path(hurdler.__file__).resolve()}")
            """,
            cell_id="exact-dna-bootstrap",
            title="0. Install and load HURDLER",
            tags=["colab-native-form", "bootstrap"],
        ),
        code(
            f'''
            #@markdown **Array mode is the default. Only the selected mode is read.**
            input_mode = "Repeat unit × copies" #@param ["Repeat unit × copies", "Complete exact DNA / FASTA"]
            sequence_id = "RF00059_TPP_riboswitch_4copy" #@param {{type:"string"}}
            repeat_unit = "{RF00059}" #@param {{type:"string"}}
            optional_spacer = "" #@param {{type:"string"}}
            repeat_copies = 4 #@param {{type:"integer"}}
            complete_exact_dna_or_fasta = "" #@param {{type:"string", placeholder:"Paste one exact A/C/G/T sequence or one FASTA record"}}

            #@markdown The workflow strips case/whitespace only. U and ambiguity codes are rejected rather than silently converted.
            ''',
            cell_id="exact-dna-input-form",
            title="1. Exact DNA or repeat-array input",
            tags=["parameters", "colab-native-form"],
        ),
        code(
            """
            #@markdown Default search limits are suitable for the regulatory-array example. Budget exhaustion is reported as `search_incomplete`, never incompatible.
            use_advanced_search_settings = False #@param {type:"boolean"}
            max_purchase_bp = 3000 #@param {type:"integer"}
            max_search_states = 10000 #@param {type:"integer"}
            search_timeout_seconds = 600 #@param {type:"integer"}
            paths_per_state = 3 #@param {type:"integer"}
            maximum_complete_routes_per_group = 25 #@param {type:"integer"}
            allow_left_cutter_fallback = False #@param {type:"boolean"}
            allow_right_cutter_fallback = False #@param {type:"boolean"}
            route_confirmation_mode = "Automatically use top-ranked route" #@param ["Automatically use top-ranked route", "Select route manually after query"]
            auto_download_results_zip = True #@param {type:"boolean"}
            """,
            cell_id="exact-dna-search-form",
            title="2. Search and annotation-aware cutter policy",
            tags=["parameters", "colab-native-form"],
        ),
        code(
            """
            #@markdown Choose IDT behavior before the molecular query. No mode ever submits an order or changes the exact target DNA.
            idt_validation_choice = "No API — export IDT Bulk Input" #@param ["No API — export IDT Bulk Input", "Live IDT API — require every purchase fragment to score <10"]
            idt_credential_source = "Colab Secrets" #@param ["Colab Secrets", "Temporary idt.env upload"]
            #@markdown Live API tries alternate breakpoints/routes within the confirmed RE/plasmid/cut-scheme selection. If none pass, the output automatically contains IDT Bulk Input files.
            """,
            cell_id="exact-dna-idt-policy",
            title="3. Choose IDT validation policy",
            tags=["parameters", "colab-native-form"],
        ),
        code(
            """
            import hashlib, json, shutil, tempfile, time, traceback
            from dataclasses import asdict
            from pathlib import Path

            import ipywidgets as widgets
            import pandas as pd
            from IPython.display import Markdown, clear_output, display

            from hurdler.constants import PLASMIDS
            from hurdler.exact_dna_design import (
                EXACT_DNA_SCHEMA_VERSION, ExactDNAQuery, ExactDNASelection,
                confirm_best_exact_dna_route, confirm_exact_dna_route,
                load_exact_dna_enzyme_catalog,
                query_exact_dna, write_exact_dna_outputs,
            )
            from hurdler.idt import (
                IDTComplexityScorer, clear_idt_secret_environment,
                configure_idt_credentials_from_bytes,
                configure_idt_credentials_from_values,
            )
            geometries = load_exact_dna_enzyme_catalog()
            site_i_names = sorted(name for name, item in geometries.items() if item.site_i_eligible)
            site_ii_names = sorted(name for name, item in geometries.items() if item.site_ii_eligible)
            site_iii_names = sorted(name for name, item in geometries.items() if item.site_iii_eligible)
            print(f"Loaded {len(site_i_names)} Site-I, {len(site_ii_names)} Site-II, and {len(site_iii_names)} Site-III enzyme choices.")
            """,
            cell_id="exact-dna-engine",
            title="4. Load exact-DNA engine",
        ),
        code(
            """
            def checkbox_panel(names, heading, default_names=None):
                defaults = set(names if default_names is None else default_names)
                boxes = {name: widgets.Checkbox(value=name in defaults, description=name, indent=False, layout=widgets.Layout(width="145px")) for name in names}
                select_all = widgets.Button(description="Select all")
                select_none = widgets.Button(description="Select none")
                select_all.on_click(lambda _button: [setattr(box, "value", True) for box in boxes.values()])
                select_none.on_click(lambda _button: [setattr(box, "value", False) for box in boxes.values()])
                grid = widgets.GridBox(
                    list(boxes.values()),
                    layout=widgets.Layout(grid_template_columns="repeat(2, 145px)", grid_gap="2px 12px"),
                )
                return boxes, widgets.VBox([widgets.HTML(f"<b>{heading}</b>"), widgets.HBox([select_all, select_none]), grid])

            site_i_boxes, site_i_panel = checkbox_panel(site_i_names, "Site I enzymes — all selected")
            site_ii_boxes, site_ii_panel = checkbox_panel(site_ii_names, "Site II enzymes — all selected")
            site_iii_boxes, site_iii_panel = checkbox_panel(site_iii_names, "Site III adapter enzymes — all selected")
            plasmid_boxes, plasmid_panel = checkbox_panel(list(PLASMIDS), "Plasmid profiles")
            display(widgets.HTML("All eligible active/latent enzymes, all Site-III adapters, and all plasmid profiles are selected by default."))
            display(widgets.HBox([site_i_panel, site_ii_panel], layout=widgets.Layout(align_items="flex-start")))
            display(site_iii_panel)
            display(plasmid_panel)
            """,
            cell_id="exact-dna-re-plasmid-selection",
            title="5. Select individual RE enzymes and plasmids",
        ),
        code(
            """
            progress_bar = widgets.IntProgress(value=0, min=0, max=1, description="States")
            progress_label = widgets.HTML("Ready")
            query_button = widgets.Button(description="Run exact route query", button_style="primary", icon="search")
            query_output = widgets.Output()
            pair_dropdown = widgets.Dropdown(options=[("Choose RE pair", None)], description="RE pair", disabled=True, layout=widgets.Layout(width="47%"))
            plasmid_dropdown = widgets.Dropdown(options=[("Choose plasmid", None)], description="Plasmid", disabled=True, layout=widgets.Layout(width="47%"))
            scheme_dropdown = widgets.Dropdown(options=[("Choose cut scheme", None)], description="Cut scheme", disabled=True, layout=widgets.Layout(width="47%"))
            route_dropdown = widgets.Dropdown(options=[("Choose route", None)], description="Route", disabled=True, layout=widgets.Layout(width="47%"))
            confirm_button = widgets.Button(description="Confirm selected route", button_style="success", disabled=True)
            confirmation_output = widgets.Output()
            query_result = None
            query_fingerprint = ""
            confirmed_route_id = ""
            search_event_count = 0

            def selected(boxes):
                return tuple(name for name, box in boxes.items() if box.value)

            def current_query():
                chosen_i = selected(site_i_boxes)
                chosen_ii = selected(site_ii_boxes)
                chosen_iii = selected(site_iii_boxes)
                chosen_plasmids = selected(plasmid_boxes)
                if not chosen_i or not chosen_ii or not chosen_iii or not chosen_plasmids:
                    raise ValueError("Select at least one Site-I, Site-II, Site-III enzyme, and plasmid profile")
                advanced = bool(use_advanced_search_settings)
                return ExactDNAQuery(
                    schema_version=EXACT_DNA_SCHEMA_VERSION,
                    input_mode="array" if input_mode == "Repeat unit × copies" else "exact",
                    sequence_id=sequence_id,
                    repeat_unit=repeat_unit,
                    spacer=optional_spacer,
                    repeat_copies=int(repeat_copies),
                    exact_dna=complete_exact_dna_or_fasta,
                    site_i_allowlist=chosen_i,
                    site_ii_allowlist=chosen_ii,
                    site_iii_allowlist=chosen_iii,
                    plasmid_allowlist=chosen_plasmids,
                    allow_left_cutter_in_hurdler_pair=bool(allow_left_cutter_fallback),
                    allow_right_cutter_in_hurdler_pair=bool(allow_right_cutter_fallback),
                    max_purchase_bp=int(max_purchase_bp),
                    max_states=int(max_search_states if advanced else 10000),
                    timeout_seconds=int(search_timeout_seconds if advanced else 600),
                    paths_per_state=int(paths_per_state if advanced else 3),
                    max_complete_routes=int(maximum_complete_routes_per_group if advanced else 25),
                )

            def fingerprint(query):
                return hashlib.sha256(json.dumps(asdict(query), sort_keys=True).encode()).hexdigest()

            def invalidate_confirmation():
                global confirmed_route_id
                confirmed_route_id = ""
                confirm_button.disabled = not bool(route_dropdown.value)
                if "validate_button" in globals():
                    validate_button.disabled = True

            def routes_for_current_selection(*, through="route"):
                if query_result is None:
                    return []
                rows = list(query_result.route_candidates)
                pair = pair_dropdown.value
                if pair:
                    rows = [row for row in rows if any(
                        (item["site_i_enzyme"], item["site_ii_enzyme"]) == tuple(pair)
                        for item in row["pairs"]
                    )]
                if through in {"plasmid", "scheme", "route"} and plasmid_dropdown.value:
                    rows = [row for row in rows if row["profile_id"] == plasmid_dropdown.value]
                if through in {"scheme", "route"} and scheme_dropdown.value:
                    rows = [row for row in rows if row["scheme_id"] == scheme_dropdown.value]
                return rows

            def refresh_plasmids(_change=None):
                invalidate_confirmation()
                pair = pair_dropdown.value
                choices = sorted({row["profile_id"] for row in routes_for_current_selection(through="pair")}) if pair else []
                plasmid_dropdown.options = [("Choose plasmid", None), *[(item, item) for item in choices]]
                plasmid_dropdown.value = None
                plasmid_dropdown.disabled = not bool(choices)

            def refresh_schemes(_change=None):
                invalidate_confirmation()
                choices = sorted({
                    (row["scheme_id"], row["cut_scheme"])
                    for row in routes_for_current_selection(through="plasmid")
                }) if plasmid_dropdown.value else []
                scheme_dropdown.options = [("Choose cut scheme", None), *[(f"{label} · {identifier}", identifier) for identifier, label in choices]]
                scheme_dropdown.value = None
                scheme_dropdown.disabled = not bool(choices)

            def refresh_routes(_change=None):
                invalidate_confirmation()
                rows = routes_for_current_selection(through="scheme") if scheme_dropdown.value else []
                options = [("Choose route", None)]
                for rank, row in enumerate(rows, 1):
                    options.append((
                        f"{rank}: {row['transition_count']} transitions · {row['hurdler_step_count']} cycles · {row['unique_purchase_count']} unique purchases",
                        row["route_id"],
                    ))
                route_dropdown.options = options
                route_dropdown.value = None
                route_dropdown.disabled = not bool(rows)

            def route_changed(_change=None):
                invalidate_confirmation()

            pair_dropdown.observe(refresh_plasmids, names="value")
            plasmid_dropdown.observe(refresh_schemes, names="value")
            scheme_dropdown.observe(refresh_routes, names="value")
            route_dropdown.observe(route_changed, names="value")
            for selector_box in [*site_i_boxes.values(), *site_ii_boxes.values(), *site_iii_boxes.values(), *plasmid_boxes.values()]:
                selector_box.observe(lambda _change: invalidate_confirmation(), names="value")

            def on_progress(event):
                global search_event_count
                details = dict(event.details)
                if event.status == "state_completed":
                    search_event_count += 1
                    progress_bar.value = min(progress_bar.max, search_event_count)
                progress_label.value = (
                    f"<b>{event.stage}</b> · {event.status} · {event.message} "
                    f"· state={details.get('state_length_bp', event.copies or '—')} "
                    f"· edges={details.get('edge_count', '—')}"
                )

            def run_query(_button=None):
                global query_result, query_fingerprint, confirmed_route_id, search_event_count
                query_button.disabled = True
                pair_dropdown.disabled = True
                plasmid_dropdown.disabled = True
                scheme_dropdown.disabled = True
                route_dropdown.disabled = True
                confirm_button.disabled = True
                confirmed_route_id = ""
                with query_output:
                    clear_output(wait=True)
                    try:
                        query = current_query()
                        search_event_count = 0
                        progress_bar.max = max(1, int(query.repeat_copies) - 1) if query.input_mode == "array" else int(query.max_states)
                        progress_bar.value = 0
                        progress_bar.bar_style = "info"
                        display(progress_bar, progress_label)
                        started = time.monotonic()
                        query_result = query_exact_dna(query, progress_callback=on_progress)
                        query_fingerprint = fingerprint(query)
                        progress_bar.bar_style = "success" if query_result.status == "hurdler_compatible_molecular" else "warning"
                        display(Markdown(f"**Status:** `{query_result.status}` — {query_result.message}"))
                        display(Markdown(
                            f"Target: **{query_result.target_length_bp:,} bp** · "
                            f"active hits: **{sum(row['state']=='active' for row in query_result.restriction_hits)}** · "
                            f"latent hits: **{sum(row['state']=='latent' for row in query_result.restriction_hits)}** · "
                            f"routes: **{len(query_result.route_candidates)}** · elapsed: **{time.monotonic()-started:.1f}s**"
                        ))
                        if query_result.restriction_hits:
                            hit_columns = [
                                "canonical_enzyme", "state", "orientation", "start", "end",
                                "observed", "site", "mismatch_position", "target_base", "active_base",
                            ]
                            display(Markdown("**Active and one-base latent sites in the exact target**"))
                            display(pd.DataFrame(query_result.restriction_hits)[hit_columns])
                        if query_result.pair_candidates:
                            display(pd.DataFrame(query_result.pair_candidates))
                        if query_result.route_candidates:
                            table = pd.DataFrame(query_result.route_candidates).drop(columns=["seed", "pairs", "silencing_decisions"], errors="ignore")
                            rows_per_page = 25
                            page_count = max(1, (len(table) + rows_per_page - 1) // rows_per_page)
                            route_page = widgets.BoundedIntText(value=1, min=1, max=page_count, description="Route page")
                            route_page_output = widgets.Output()
                            def render_route_page(_change=None):
                                with route_page_output:
                                    clear_output(wait=True)
                                    first = (int(route_page.value) - 1) * rows_per_page
                                    display(table.iloc[first:first + rows_per_page])
                                    display(Markdown(f"Page **{route_page.value}/{page_count}** · all **{len(table)}** grouped routes remain selectable below."))
                            route_page.observe(render_route_page, names="value")
                            display(route_page, route_page_output)
                            render_route_page()
                            pairs = sorted({
                                (item["site_i_enzyme"], item["site_ii_enzyme"])
                                for row in query_result.route_candidates for item in row["pairs"]
                            })
                            pair_dropdown.options = [("Choose RE pair", None), *[(f"{left} / {right}", (left, right)) for left, right in pairs]]
                            pair_dropdown.value = None
                            pair_dropdown.disabled = False
                    except Exception as exc:
                        progress_bar.bar_style = "danger"
                        display(Markdown(f"**Query failed safely:** `{type(exc).__name__}: {exc}`"))
                query_button.disabled = False

            def confirm_route(_button=None):
                global confirmed_route_id
                with confirmation_output:
                    clear_output(wait=True)
                    try:
                        query = current_query()
                        if query_result is None or fingerprint(query) != query_fingerprint:
                            raise RuntimeError("Inputs or enzyme/plasmid selections changed; re-run the query")
                        confirmed_route_id = str(route_dropdown.value)
                        if not route_dropdown.value:
                            raise RuntimeError("Select an RE pair, plasmid, cut scheme, and exact route")
                        route = next(row for row in query_result.route_candidates if row["route_id"] == confirmed_route_id)
                        preview = confirm_exact_dna_route(query_result, ExactDNASelection(
                            confirmed_route_id, "none",
                            plasmid_profile=str(plasmid_dropdown.value),
                            cut_scheme_id=str(scheme_dropdown.value),
                            site_i_enzyme=str(pair_dropdown.value[0]),
                            site_ii_enzyme=str(pair_dropdown.value[1]),
                        ))
                        display(Markdown(
                            f"**Confirmed:** `{confirmed_route_id}` · "
                            f"{route['profile_id']} · {route['cut_scheme']}. "
                            "The validation/export panel below is now enabled."
                        ))
                        display(Markdown("**Verified active/latent transitions**"))
                        display(pd.DataFrame(preview.latent_transitions))
                        display(Markdown("**Exact purchase topology (IDT not called)**"))
                        display(pd.DataFrame(preview.purchase_fragments).drop(
                            columns=["purchase_sequence", "secondary_purchase_sequence", "primer_forward_5to3", "primer_reverse_5to3"], errors="ignore"
                        ))
                        if "validate_button" in globals():
                            validate_button.disabled = False
                    except Exception as exc:
                        confirmed_route_id = ""
                        if "validate_button" in globals():
                            validate_button.disabled = True
                        display(Markdown(f"**Confirmation failed:** `{type(exc).__name__}: {exc}`"))

            query_button.on_click(run_query)
            confirm_button.on_click(confirm_route)
            display(widgets.VBox([
                query_button, query_output,
                widgets.HBox([pair_dropdown, plasmid_dropdown]),
                widgets.HBox([scheme_dropdown, route_dropdown]),
                confirm_button, confirmation_output,
            ]))

            def select_and_confirm_top_route():
                if query_result is None or not query_result.route_candidates:
                    return False
                top = query_result.route_candidates[0]
                pair = top["pairs"][0]
                pair_dropdown.value = (pair["site_i_enzyme"], pair["site_ii_enzyme"])
                plasmid_dropdown.value = top["profile_id"]
                scheme_dropdown.value = top["scheme_id"]
                route_dropdown.value = top["route_id"]
                confirm_route()
                return bool(confirmed_route_id)

            # Default Runtime → Run all completes query and route confirmation.
            # Manual mode preserves the cascading selectors for expert review.
            run_query()
            if route_confirmation_mode == "Automatically use top-ranked route":
                select_and_confirm_top_route()
            """,
            cell_id="exact-dna-query-and-confirm",
            title="6. Search, inspect, and confirm an exact route",
        ),
        code(
            """
            validation_mode_value = "api" if idt_validation_choice.startswith("Live IDT") else "batch"
            credential_source_value = "secrets" if idt_credential_source == "Colab Secrets" else "upload"
            credential_upload = widgets.FileUpload(accept=".env,text/plain", multiple=False, description="Upload temporary idt.env")
            default_output_directory = (
                "/content/exact_dna_hurdler_design"
                if in_colab else
                str(Path(tempfile.gettempdir()) / "exact_dna_hurdler_design")
            )
            output_directory = widgets.Text(value=default_output_directory, description="Output", layout=widgets.Layout(width="98%"))
            validate_button = widgets.Button(
                description="Validate / export", button_style="primary",
                disabled=not bool(confirmed_route_id),
            )
            download_button = widgets.Button(description="Download ZIP", icon="download", disabled=True)
            audit_download_button = widgets.Button(description="Download technical audit ZIP", icon="download", disabled=True)
            validation_progress = widgets.HTML("Confirm a route above first.")
            validation_output = widgets.Output()
            output_zip = None
            audit_zip = None

            def secret_value(reader, name):
                try:
                    return str(reader.get(name) or "").strip()
                except Exception:
                    return ""

            def configure_secrets():
                from google.colab import userdata
                token = secret_value(userdata, "IDT_ACCESS_TOKEN")
                if token:
                    return configure_idt_credentials_from_values({"IDT_ACCESS_TOKEN": token}, auth_method="access_token")
                values = {name: secret_value(userdata, name) for name in ("IDT_CLIENT_ID", "IDT_CLIENT_SECRET", "IDT_USERNAME", "IDT_PASSWORD")}
                return configure_idt_credentials_from_values(values, auth_method="password")

            def uploaded_bytes():
                value = credential_upload.value
                if not value:
                    raise RuntimeError("Upload an external idt.env file or use Colab Secrets")
                item = next(iter(value.values())) if isinstance(value, dict) else value[0]
                return bytes(item["content"] if isinstance(item, dict) else item.content)

            def clear_upload_control():
                global credential_upload
                credential_upload = widgets.FileUpload(accept=".env,text/plain", multiple=False, description="Upload temporary idt.env")
                credential_row.children = (credential_upload,)

            def configure_credentials():
                if credential_source_value == "secrets":
                    return configure_secrets()
                try:
                    return configure_idt_credentials_from_bytes(uploaded_bytes())
                finally:
                    clear_upload_control()

            def validation_event(event):
                validation_progress.value = f"<b>{event.stage}</b> · {event.status} · {event.message}"

            def run_validation(_button=None):
                global output_zip, audit_zip
                validate_button.disabled = True
                download_button.disabled = True
                with validation_output:
                    clear_output(wait=True)
                    try:
                        query = current_query()
                        if query_result is None or fingerprint(query) != query_fingerprint or not confirmed_route_id:
                            raise RuntimeError("The route is missing or stale; re-run and confirm the query")
                        destination = Path(output_directory.value).expanduser()
                        if destination.exists():
                            shutil.rmtree(destination)
                        destination.mkdir(parents=True)
                        scorer = None
                        with tempfile.TemporaryDirectory(prefix="hurdler-idt-audit-") as temporary:
                            if validation_mode_value == "api":
                                configure_credentials()
                                scorer = IDTComplexityScorer(Path(temporary) / "raw.jsonl")
                            selected = ExactDNASelection(
                                confirmed_route_id, validation_mode_value,
                                plasmid_profile=str(plasmid_dropdown.value),
                                cut_scheme_id=str(scheme_dropdown.value),
                                site_i_enzyme=str(pair_dropdown.value[0]),
                                site_ii_enzyme=str(pair_dropdown.value[1]),
                            )
                            confirmer = confirm_best_exact_dna_route if validation_mode_value == "api" else confirm_exact_dna_route
                            result = confirmer(
                                query_result, selected, idt_scorer=scorer, progress_callback=validation_event
                            )
                        files = write_exact_dna_outputs(result, destination)
                        audit_zip = Path(files.pop("technical_audit_zip"))
                        output_zip = Path(shutil.make_archive(str(destination.resolve()), "zip", root_dir=destination.resolve()))
                        download_button.disabled = False
                        audit_download_button.disabled = not audit_zip.is_file()
                        display(Markdown(f"**Result:** `{result.status}` — {result.message}"))
                        display(Markdown(f"**Cloning steps:** {len(result.cloning_steps)}"))
                        display(pd.DataFrame(result.cloning_steps))
                        display(Markdown(
                            f"Generated **{len(files)}** user-facing files. Technical hashes and route audit are in a separate optional ZIP. No order was submitted."
                        ))
                    except Exception as exc:
                        validation_progress.value = f"<b>failed</b> · {type(exc).__name__}"
                        display(Markdown(f"**Validation/export failed safely:** `{type(exc).__name__}: {exc}`"))
                    finally:
                        clear_idt_secret_environment()
                        validate_button.disabled = not bool(confirmed_route_id)

            def download_zip(_button=None):
                if output_zip is None or not output_zip.is_file():
                    raise FileNotFoundError("Run validation/export first")
                from google.colab import files
                files.download(str(output_zip))

            def download_audit_zip(_button=None):
                if audit_zip is None or not audit_zip.is_file():
                    raise FileNotFoundError("Run validation/export first")
                from google.colab import files
                files.download(str(audit_zip))

            validate_button.on_click(run_validation)
            download_button.on_click(download_zip)
            audit_download_button.on_click(download_audit_zip)
            credential_row = widgets.HBox([credential_upload])
            display(widgets.VBox([
                widgets.HTML("<b>IDT is score-only. It never optimizes DNA and never submits an order.</b>"),
                widgets.HTML(f"Policy selected above: <b>{idt_validation_choice}</b>"),
                credential_row if validation_mode_value == "api" and credential_source_value == "upload" else widgets.HTML("Credentials will be read only when Live API validation starts."),
                output_directory, validation_progress,
                widgets.HBox([validate_button, download_button]), audit_download_button, validation_output,
            ]))

            # Default Run all produces the complete offline Bulk package. Live
            # API mode likewise runs automatically after credentials are read.
            if route_confirmation_mode == "Automatically use top-ranked route" and confirmed_route_id:
                run_validation()
                if auto_download_results_zip and in_colab and output_zip is not None:
                    download_zip()
            """,
            cell_id="exact-dna-idt-export",
            title="7. Independent validation, IDT scoring, and export",
        ),
    ]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(book, OUTPUT)
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
