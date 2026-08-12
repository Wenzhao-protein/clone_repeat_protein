#!/usr/bin/env python3
"""Validate canonical PDF/PNG pairs and render a visual contact sheet."""

from __future__ import annotations

import argparse
import base64
import hashlib
import html
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from PIL import Image, ImageOps, ImageDraw
import plotly.io as pio


FIGURE_SOURCES = {
    "success_rate_1_60": "notebooks/tasks/02_success_rate_1_60.ipynb",
    "success_rate_1_60_scan_3x": "notebooks/tasks/02_success_rate_1_60.ipynb",
    "module_benchmark": "notebooks/tasks/03_repeat_module_benchmark.ipynb",
    "adaptive_maximum_search": "notebooks/tasks/03_repeat_module_benchmark.ipynb",
    "module_length_vs_max_orderable_copies": "studies/hurdler_validation/scripts/build_module_summary.py",
    "module_hurdler_usable_fraction": "studies/hurdler_validation/scripts/build_module_summary.py",
    "source_vs_primitive_module_length": "notebooks/tasks/05_module_boundary_inference.ipynb",
    "module_harmonic_ratio": "notebooks/tasks/05_module_boundary_inference.ipynb",
    "module_fixed_fraction": "notebooks/tasks/05_module_boundary_inference.ipynb",
    "secondary_structure_coverage": "notebooks/tasks/05_module_boundary_inference.ipynb",
    "sequence_vs_secondary_structure_evidence": "notebooks/tasks/05_module_boundary_inference.ipynb",
    "module_boundary_examples": "notebooks/tasks/05_module_boundary_inference.ipynb",
    "module_compatibility_by_length": "notebooks/tasks/06_module_compatibility_by_length.ipynb",
    "module_length_vs_maximum_copies": "notebooks/tasks/07_adaptive_selected_pair_capacity.ipynb",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_source(value: object) -> str:
    return str(value).split("clone_repeat_protein/", 1)[-1]


def extract_legacy_figures(repo: Path, destination: Path) -> tuple[list[Path], dict[str, str]]:
    """Materialize every static-capable figure from passed legacy notebooks."""
    paths: list[Path] = []
    sources: dict[str, str] = {}
    legacy_root = repo / "studies" / "hurdler_validation" / "step05_reproducibility" / "legacy_notebooks"
    destination.mkdir(parents=True, exist_ok=True)
    for status_path in sorted(legacy_root.glob("*/status.json")):
        status = json.loads(status_path.read_text())
        if status.get("status") != "passed":
            continue
        notebook = Path(str(status.get("output_path", "")))
        if not notebook.is_file():
            candidates = sorted(status_path.parent.glob("*_executed.ipynb"))
            if not candidates:
                continue
            notebook = candidates[0]
        book = json.loads(notebook.read_text())
        for cell_index, cell in enumerate(book.get("cells", [])):
            for output_index, output in enumerate(cell.get("outputs", [])):
                data = output.get("data", {})
                figure_id = f"legacy_{status_path.parent.name}_c{cell_index:03d}_o{output_index:02d}"
                png = destination / f"{figure_id}.png"
                pdf = destination / f"{figure_id}.pdf"
                if "image/png" in data:
                    encoded = data["image/png"]
                    if isinstance(encoded, list):
                        encoded = "".join(encoded)
                    with Image.open(io.BytesIO(base64.b64decode(encoded))) as image:
                        rgb = Image.new("RGB", image.size, "white")
                        if "A" in image.getbands():
                            rgb.paste(image.convert("RGBA"), mask=image.convert("RGBA").getchannel("A"))
                        else:
                            rgb.paste(image.convert("RGB"))
                        rgb.save(png, dpi=(300, 300))
                        rgb.save(pdf, "PDF", resolution=300)
                elif "application/vnd.plotly.v1+json" in data:
                    figure = pio.from_json(json.dumps(data["application/vnd.plotly.v1+json"]))
                    figure.write_image(png, width=1200, height=750, scale=2)
                    figure.write_image(pdf, width=1200, height=750)
                else:
                    continue
                paths.append(png)
                sources[figure_id] = normalized_source(status.get("source", ""))
    return paths, sources


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--required-figure",
        action="append",
        help="Require only the named figure(s); default requires every canonical figure.",
    )
    args = parser.parse_args()
    repo = args.repo.resolve()
    roots = [
        repo / "studies" / "hurdler_validation" / "step02_success_landscape" / "figures",
        repo / "studies" / "hurdler_validation" / "step03_module_corpus" / "figures",
        repo / "studies" / "hurdler_validation" / "step04_module_optimization" / "figures",
    ]
    if args.required_figure:
        extracted, extracted_sources = [], {}
    else:
        extracted, extracted_sources = extract_legacy_figures(
            repo, args.output_dir / "extracted_legacy"
        )
    rows = []
    thumbnails: list[tuple[str, Image.Image]] = []
    canonical_by_stem: dict[str, Path] = {}
    for root in roots:
        for png in sorted(root.rglob("*.png")):
            previous = canonical_by_stem.get(png.stem)
            # The versioned middle-module analysis supersedes the top-level
            # first-unit figure with the same stable figure ID.
            if (
                previous is None
                or "periodic_v4" in png.parts
                or len(png.parts) < len(previous.parts)
            ):
                canonical_by_stem[png.stem] = png
    if args.required_figure:
        requested = set(args.required_figure)
        figure_paths = [
            canonical_by_stem[key]
            for key in sorted(canonical_by_stem)
            if key in requested
        ]
    else:
        figure_paths = [canonical_by_stem[key] for key in sorted(canonical_by_stem)] + extracted
    for png in figure_paths:
        pdf = png.with_suffix(".pdf")
        if not pdf.is_file() or not pdf.stat().st_size:
            raise RuntimeError(f"Missing non-empty PDF pair for {png}")
        with Image.open(png) as image:
            image.verify()
        with Image.open(png) as image:
            width, height = image.size
            thumbnail = ImageOps.contain(image.convert("RGB"), (700, 430))
            thumbnails.append((png.stem, thumbnail.copy()))
        source = FIGURE_SOURCES[png.stem] if png.stem in FIGURE_SOURCES else extracted_sources[png.stem]
        rows.append(
            {
                "figure_id": png.stem,
                "source_notebook": source,
                "source_sha256": sha256(repo / source),
                "png": str(png),
                "pdf": str(pdf),
                "png_sha256": sha256(png),
                "pdf_sha256": sha256(pdf),
                "width_px": width,
                "height_px": height,
                "png_bytes": png.stat().st_size,
                "pdf_bytes": pdf.stat().st_size,
                "generated_at": datetime.fromtimestamp(png.stat().st_mtime, timezone.utc).isoformat(),
                "style": (
                    "original-repo continuous lines and circular markers"
                    if png.stem in {"success_rate_1_60", "success_rate_1_60_scan_3x"}
                    else "white-background UW/IPD palette"
                ),
                "analysis_version": (
                    "historical-notebook-success-v1-three-copy-scan"
                    if png.stem == "success_rate_1_60_scan_3x"
                    else "historical-notebook-success-v1"
                    if png.stem == "success_rate_1_60"
                    else "expanded-middle-repeatsdb-foldseek-v1"
                    if "expanded-middle-repeatsdb-foldseek-v1" in png.parts
                    else "periodic_v4_middle_unit"
                    if "periodic_v4" in png.parts
                    else "legacy_or_unversioned"
                ),
                "status": "passed",
            }
        )
    expected = set(args.required_figure or FIGURE_SOURCES)
    observed = {row["figure_id"] for row in rows}
    if not expected.issubset(observed):
        raise RuntimeError(f"Figure set mismatch: missing={expected-observed}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.output_dir / "figure_manifest.csv", index=False)
    frame.to_json(args.output_dir / "figure_manifest.json", orient="records", indent=2)

    sheet_width = 800
    panel_width = 400
    panel_height = 330
    sheet = Image.new("RGB", (sheet_width, panel_height * ((len(thumbnails) + 1) // 2)), "white")
    draw = ImageDraw.Draw(sheet)
    for index, (label, thumbnail) in enumerate(thumbnails):
        column = index % 2
        row = index // 2
        left = column * panel_width
        top = row * panel_height
        thumbnail = ImageOps.contain(thumbnail, (360, 270))
        sheet.paste(thumbnail, (left + (panel_width - thumbnail.width) // 2, top + 45))
        draw.text((left + 15, top + 12), label, fill="#4B2E83")
    sheet_path = args.output_dir / "contact_sheet.png"
    sheet.save(sheet_path, dpi=(150, 150))

    html_rows = "\n".join(
        f"<tr><td>{html.escape(row['figure_id'])}</td><td>{html.escape(row['source_notebook'])}</td>"
        f"<td>{row['width_px']} x {row['height_px']}</td><td>{row['status']}</td></tr>"
        for row in rows
    )
    report = f"""<!doctype html><meta charset='utf-8'><title>HURDLER figure report</title>
<style>body{{font-family:Arial,sans-serif;margin:2rem;color:#222}}table{{border-collapse:collapse}}
td,th{{border:1px solid #bbb;padding:.5rem}}h1{{color:#4B2E83}}img{{max-width:760px}}</style>
<h1>HURDLER figure validation</h1><p>All canonical figures have non-empty PDF and PNG pairs.</p>
<table><thead><tr><th>Figure</th><th>Source</th><th>Pixels</th><th>Status</th></tr></thead>
<tbody>{html_rows}</tbody></table><h2>Contact sheet</h2><img src='contact_sheet.png'>"""
    (args.output_dir / "figure_report.html").write_text(report)
    summary = {"figure_count": len(rows), "passed": True, "contact_sheet_sha256": sha256(sheet_path)}
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
