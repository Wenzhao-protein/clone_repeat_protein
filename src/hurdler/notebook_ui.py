"""ipywidgets presentation layer for the interactive HURDLER designer."""

from __future__ import annotations

import getpass
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .constants import PLASMIDS
from .design import (
    DesignRequest,
    DesignResult,
    analyze_repeat_sequence,
    boundary_confirmation_token,
    bundled_index_dir,
    confirm_repeat_boundary,
    design_construct,
    enumerate_design_candidates,
    parse_protein_input,
    role_enzyme_options,
    write_design_outputs,
)
from .idt import (
    IDTComplexityScorer,
    clear_idt_secret_environment,
    configure_idt_credentials,
    configure_idt_credentials_from_bytes,
    configure_idt_credentials_from_values,
)
from .index import PatternIndex


RPB1_CTD_FASTA = """>S_cerevisiae_Rpb1_CTD_26_repeats_plus_C_terminal_tip
FSPTSPTYSPTSPAYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPT
SPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPAY
SPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPNYSPTSPSYSPTSPGYSPGSPAYSPKQ
DEQKHNENENSR"""


class PassingMockIDTScorer:
    """Explicit test-only scorer for notebook headless smoke execution."""

    def score(self, name: str, sequence: str) -> dict[str, Any]:
        sequence_sha = hashlib.sha256(sequence.encode()).hexdigest()
        response_sha = hashlib.sha256(f"mock|{name}|{sequence_sha}".encode()).hexdigest()
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
            "idt_response_sha256": response_sha,
            "mock_scorer": True,
        }


def run_headless_designer_smoke(
    output_dir: str | Path | None = None,
) -> dict[str, DesignResult]:
    """Run the user-supplied Rpb1 CTD golden case and an orderable control."""
    index = PatternIndex.load(bundled_index_dir())
    rpb1_analysis = analyze_repeat_sequence(RPB1_CTD_FASTA)
    if rpb1_analysis.proposed_period != 7:
        raise AssertionError("Rpb1 CTD golden input did not recover the 7-AA heptad period")
    rpb1_request = DesignRequest(
        sequence_id=rpb1_analysis.sequence_id,
        full_protein_sequence=rpb1_analysis.full_protein_sequence,
        target_repeat_copies=26,
        plasmid="pGEX-4T-1",
        confirmed_repeat_start=1,
        confirmed_repeat_end=182,
        confirmed_period=7,
        confirmation_token=boundary_confirmation_token(
            rpb1_analysis.full_protein_sequence, 1, 182, 7
        ),
        optimize=False,
    )
    rpb1_result = design_construct(rpb1_request, index=index)
    if rpb1_result.status != "hurdler_incompatible":
        raise AssertionError(
            "Rpb1 CTD YSPTSPS middle heptad golden result unexpectedly changed"
        )
    if any(
        enumerate_design_candidates("YSPTSPS", maintained_plasmid, index)
        for maintained_plasmid in PLASMIDS
    ):
        raise AssertionError("Rpb1 CTD middle heptad unexpectedly passed another maintained plasmid")

    # A separate positive control is required to exercise GA, fragment-level
    # IDT policy, exact assembly, and export even though the biological golden
    # example is correctly reported incompatible by the frozen index.
    module = "HIKLMNPQRST"
    protein = "M" + module * 4 + "K"
    start = 2
    end = 1 + len(module) * 4
    request = DesignRequest(
        sequence_id="headless_smoke",
        full_protein_sequence=protein,
        target_repeat_copies=4,
        plasmid="pGEX-4T-1",
        confirmed_repeat_start=start,
        confirmed_repeat_end=end,
        confirmed_period=len(module),
        confirmation_token=boundary_confirmation_token(protein, start, end, len(module)),
        optimize=True,
        population_size=4,
    )
    result = design_construct(request, index=index, idt_scorer=PassingMockIDTScorer())
    if output_dir is not None:
        write_design_outputs(rpb1_result, Path(output_dir) / "rpb1_ctd_example")
        write_design_outputs(result, output_dir)
    if not result.orderable:
        raise AssertionError(f"Headless designer smoke did not complete: {result.status}")
    return {"rpb1_ctd": rpb1_result, "compatible_control": result}


@dataclass
class DesignerApp:
    widget: Any
    index: PatternIndex
    last_result: DesignResult | None = None
    confirmation_token: str | None = None
    credentials_ready: bool = False


def build_designer_ui(*, index_dir: str | Path | None = None) -> DesignerApp:
    """Build the interactive UI; scientific work remains in :mod:`hurdler.design`."""
    try:
        import ipywidgets as widgets
        from IPython.display import Markdown, clear_output, display
    except ImportError as exc:  # pragma: no cover - exercised by notebook environment
        raise RuntimeError("Install hurdler[notebooks] to use the interactive designer") from exc

    index = PatternIndex.load(index_dir or bundled_index_dir())
    enzyme_options = role_enzyme_options(index)

    sequence = widgets.Textarea(
        description="Protein",
        value=RPB1_CTD_FASTA,
        placeholder="Paste one raw AA sequence or one FASTA record",
        layout=widgets.Layout(width="100%", height="160px"),
    )
    target_copies = widgets.BoundedIntText(value=26, min=2, max=10000, description="Target copies")
    plasmid = widgets.Dropdown(options=PLASMIDS, value=PLASMIDS[0], description="Plasmid")
    infer_button = widgets.Button(description="1 · Infer repeat boundary", button_style="primary")
    infer_output = widgets.Output()

    start = widgets.BoundedIntText(value=1, min=1, max=1_000_000, description="Start (1-based)")
    end = widgets.BoundedIntText(value=1, min=1, max=1_000_000, description="End (inclusive)")
    period = widgets.BoundedIntText(value=1, min=1, max=10000, description="Period (AA)")
    confirm_button = widgets.Button(description="2 · Confirm these boundaries", button_style="warning")
    confirmation_status = widgets.HTML("<b>Not confirmed.</b> Sequence changes invalidate confirmation.")

    site_i = widgets.SelectMultiple(
        options=enzyme_options["site_i"], description="Site I", rows=7,
        layout=widgets.Layout(width="32%"),
    )
    site_ii = widgets.SelectMultiple(
        options=enzyme_options["site_ii"], description="Site II", rows=7,
        layout=widgets.Layout(width="32%"),
    )
    site_iii = widgets.SelectMultiple(
        options=enzyme_options["site_iii"], description="Site III", rows=7,
        layout=widgets.Layout(width="32%"),
    )
    query_button = widgets.Button(description="3 · Enumerate HURDLER routes", button_style="info")
    query_output = widgets.Output()
    selected_candidate = widgets.Dropdown(options=[], description="Route", layout=widgets.Layout(width="100%"))

    optimize = widgets.Checkbox(value=False, description="Run codon optimization and live IDT scoring")
    credential_mode = widgets.Dropdown(
        options=[
            ("Manual OAuth fields", "manual_password"),
            ("Manual access token", "manual_token"),
            ("Repo-external mode-600 env path (hidden prompt)", "path"),
            ("Temporary env upload", "upload"),
        ],
        value="path",
        description="Credentials",
    )
    client_id = widgets.Password(description="Client ID")
    client_secret = widgets.Password(description="Client secret")
    username = widgets.Password(description="Username")
    password = widgets.Password(description="Password")
    access_token = widgets.Password(description="Access token")
    upload = widgets.FileUpload(accept=".env,text/plain", multiple=False, description="Upload env")
    configure_button = widgets.Button(description="Configure IDT credentials")
    credential_status = widgets.HTML("Credentials are not configured.")
    output_dir = widgets.Text(value="output/interactive_design", description="Output", layout=widgets.Layout(width="100%"))
    run_button = widgets.Button(description="4 · Optimize, score, and export", button_style="success")
    run_output = widgets.Output()

    app = DesignerApp(widget=None, index=index)

    def invalidate(_change: Any = None) -> None:
        app.confirmation_token = None
        app.last_result = None
        confirmation_status.value = "<b>Not confirmed.</b> Sequence changes invalidate confirmation."
        selected_candidate.options = []

    sequence.observe(invalidate, names="value")
    for boundary_widget in (start, end, period):
        boundary_widget.observe(invalidate, names="value")

    def on_infer(_button: Any) -> None:
        invalidate()
        with infer_output:
            clear_output(wait=True)
            try:
                analysis = analyze_repeat_sequence(sequence.value)
            except Exception as exc:
                display(Markdown(f"**Inference failed:** `{type(exc).__name__}: {exc}`"))
                return
            if analysis.proposed_start is not None:
                start.value = int(analysis.proposed_start)
                end.value = int(analysis.proposed_end)
                period.value = int(analysis.proposed_period)
            display(Markdown(
                "Candidate periods are suggestions only. Review the complete region and confirm or edit the coordinates."
            ))
            frame = pd.DataFrame([candidate.to_dict() for candidate in analysis.candidates])
            display(frame.head(25) if not frame.empty else Markdown("No confident period; enter boundaries manually."))

    infer_button.on_click(on_infer)

    def on_confirm(_button: Any) -> None:
        with infer_output:
            try:
                _, normalized = parse_protein_input(sequence.value)
                token = boundary_confirmation_token(normalized, start.value, end.value, period.value)
                # The core validator rejects partial units and out-of-range coordinates now.
                confirmed = confirm_repeat_boundary(
                    normalized, start=start.value, end=end.value, period=period.value,
                    expected_token=token,
                )
            except Exception as exc:
                app.confirmation_token = None
                confirmation_status.value = f"<b>Confirmation failed:</b> {type(exc).__name__}: {exc}"
                return
        app.confirmation_token = token
        confirmation_status.value = (
            f"<b>Confirmed.</b> Middle unit {confirmed.middle_unit_index}/{confirmed.repeat_count}: "
            f"<code>{confirmed.middle_module}</code> ({confirmed.period} AA)."
        )
        with infer_output:
            display(
                pd.DataFrame(
                    {
                        "unit_index": range(1, confirmed.repeat_count + 1),
                        "sequence": confirmed.unit_sequences,
                        "selected_middle": [
                            index == confirmed.middle_unit_index
                            for index in range(1, confirmed.repeat_count + 1)
                        ],
                    }
                )
            )
            display(
                pd.DataFrame(
                    {
                        "module_position_1based": range(1, confirmed.period + 1),
                        "consensus_residue": list(confirmed.consensus_module),
                        "conservation": confirmed.position_conservation,
                        "fixed_at_80_percent": [
                            position in confirmed.fixed_positions_1based
                            for position in range(1, confirmed.period + 1)
                        ],
                    }
                )
            )

    confirm_button.on_click(on_confirm)

    def request_for_current_state(*, run_optimization: bool) -> DesignRequest:
        if app.confirmation_token is None:
            raise RuntimeError("Confirm the boundary after the most recent sequence edit")
        _, normalized = parse_protein_input(sequence.value)
        return DesignRequest(
            full_protein_sequence=normalized,
            target_repeat_copies=target_copies.value,
            plasmid=plasmid.value,
            confirmed_repeat_start=start.value,
            confirmed_repeat_end=end.value,
            confirmed_period=period.value,
            confirmation_token=app.confirmation_token,
            site_i_allowlist=tuple(site_i.value),
            site_ii_allowlist=tuple(site_ii.value),
            site_iii_allowlist=tuple(site_iii.value),
            selected_candidate_id=selected_candidate.value or None,
            optimize=run_optimization,
        )

    def on_query(_button: Any) -> None:
        with query_output:
            clear_output(wait=True)
            try:
                request = request_for_current_state(run_optimization=False)
                result = design_construct(request, index=app.index)
            except Exception as exc:
                display(Markdown(f"**HURDLER query failed:** `{type(exc).__name__}: {exc}`"))
                return
            app.last_result = result
            selected_candidate.options = [
                ("Auto · try all allowed routes in deterministic rank order", ""),
                *[
                    (
                        f"#{row['rank']} {row['site_i_enzyme']} / {row['site_ii_enzyme']} / "
                        f"{row['site_iii_enzyme']} · {row['direction']} · positions "
                        f"{int(row['site_i_position']) + 1}/{int(row['site_ii_position']) + 1}",
                        row["candidate_id"],
                    )
                    for row in result.candidates
                ],
            ]
            display(Markdown(f"**{len(result.candidates):,} allowed routes.** Status: `{result.status}`"))
            columns = [
                "rank", "candidate_id", "direction", "site_i_position", "site_ii_position",
                "site_i_enzyme", "site_ii_enzyme", "site_iii_enzyme",
                "site_i_recognition_site", "site_ii_recognition_site",
                "site_iii_recognition_site", "site_i_ovhg", "site_ii_ovhg", "site_iii_ovhg", "orthogonality",
                "site_i_9mer_bp", "site_ii_9mer_bp_original", "site_ii_9mer_bp_mutated",
                "minimum_target_repeat_copies_for_locked_windows",
                "requested_target_geometry_supported",
            ]
            display(pd.DataFrame(result.candidates)[columns].head(100) if result.candidates else Markdown(result.message))

    query_button.on_click(on_query)

    secret_widgets = (client_id, client_secret, username, password, access_token)

    def clear_secret_controls() -> None:
        for widget in secret_widgets:
            widget.value = ""
        try:
            upload.value = ()
        except (AttributeError, TypeError):
            pass

    def on_configure(_button: Any) -> None:
        status: dict[str, object] | None = None
        try:
            if credential_mode.value == "path":
                private_path = getpass.getpass("Repo-external IDT env path: ")
                status = configure_idt_credentials(
                    mode="path",
                    path=private_path,
                    headless=False,
                    include_path_in_status=False,
                )
                private_path = ""
            elif credential_mode.value == "upload":
                if not upload.value:
                    raise RuntimeError("Choose one env file first")
                entry = upload.value[0]
                payload = bytes(entry["content"])
                status = configure_idt_credentials_from_bytes(payload)
                payload = b""
            elif credential_mode.value == "manual_token":
                status = configure_idt_credentials_from_values(
                    {"IDT_ACCESS_TOKEN": access_token.value}, auth_method="access_token"
                )
            else:
                status = configure_idt_credentials_from_values(
                    {
                        "IDT_CLIENT_ID": client_id.value,
                        "IDT_CLIENT_SECRET": client_secret.value,
                        "IDT_USERNAME": username.value,
                        "IDT_PASSWORD": password.value,
                    },
                    auth_method="password",
                )
            app.credentials_ready = True
            credential_status.value = (
                f"<b>Configured in memory.</b> Mode: {status['credential_mode']}; "
                f"authentication: {status['auth_method']}."
            )
        except Exception:
            app.credentials_ready = False
            credential_status.value = (
                "<b>Credential configuration failed.</b> Verify the selected format, owner-only file permission, "
                "and that the file is outside this repository. No path or secret was retained."
            )
        finally:
            clear_secret_controls()

    configure_button.on_click(on_configure)

    def on_run(_button: Any) -> None:
        with run_output:
            clear_output(wait=True)
            if not optimize.value:
                display(Markdown(
                    "Optimization is disabled. Use **Enumerate HURDLER routes** for the topology draft; no orderable DNA will be written."
                ))
                return
            if not app.credentials_ready:
                display(Markdown("**Configure IDT credentials first.**"))
                return
            try:
                request = request_for_current_state(run_optimization=True)
                destination = Path(output_dir.value).expanduser()
                scorer = IDTComplexityScorer(destination / "idt_audit.jsonl")
                result = design_construct(request, index=app.index, idt_scorer=scorer)
                files = write_design_outputs(result, destination)
                app.last_result = result
                display(Markdown(f"**Status:** `{result.status}` — {result.message}"))
                display(pd.DataFrame(result.purchase_fragments))
                display(pd.DataFrame(result.cloning_steps))
                display(Markdown("Exported files:\n" + "\n".join(f"- `{path}`" for path in files.values())))
            except Exception as exc:
                display(Markdown(f"**Design failed:** `{type(exc).__name__}: {exc}`"))
            finally:
                clear_idt_secret_environment()
                app.credentials_ready = False
                credential_status.value = "Credentials cleared from the process; configure again for another live run."
                clear_secret_controls()

    run_button.on_click(on_run)

    credential_box = widgets.Accordion(
        children=[
            widgets.VBox(
                [
                    widgets.HTML(
                        "Credentials are used only for IDT scoring. Uploads and password fields are cleared immediately; "
                        "the path mode uses a hidden prompt and does not store the path in widget state."
                    ),
                    credential_mode,
                    widgets.HBox([client_id, client_secret]),
                    widgets.HBox([username, password]),
                    access_token,
                    upload,
                    configure_button,
                    credential_status,
                ]
            )
        ]
    )
    credential_box.set_title(0, "Optional optimization · private IDT credentials")
    app.widget = widgets.VBox(
        [
            widgets.HTML(
                "<h2>Interactive HURDLER construct designer</h2>"
                "<p>Fast sequence-only inference is a suggestion. You must confirm or edit the repeat boundary. "
                "DSSP, Foldseek, and structure prediction are intentionally not run here.</p>"
            ),
            sequence,
            widgets.HBox([target_copies, plasmid]),
            infer_button,
            infer_output,
            widgets.HBox([start, end, period]),
            confirm_button,
            confirmation_status,
            widgets.HTML("<h3>Restriction-enzyme whitelists</h3><p>Leave a role empty to allow every supported enzyme.</p>"),
            widgets.HBox([site_i, site_ii, site_iii]),
            query_button,
            query_output,
            selected_candidate,
            optimize,
            credential_box,
            output_dir,
            run_button,
            run_output,
        ]
    )
    return app
