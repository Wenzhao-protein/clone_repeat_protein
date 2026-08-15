"""Backend for V2 notebook 12: Colab-compatible agarose gel analysis."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any, Mapping

import matplotlib
import numpy as np
import pandas as pd
from PIL import Image
from scipy.signal import find_peaks, peak_widths

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "12_agarose",
    "Agarose gel analysis",
    "Read PNG/TIFF/SCN gels, quantify lanes and bands, and report missing historical SCN files explicitly.",
    default_request={"lane_count": 8},
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _source(request: Mapping[str, Any]) -> Path:
    default = repo_root() / "agarose_gel_analysis/data/wenzhao_runs/2023-04-16_dArmRP_triangle_step1_BamHI_EcoRI_2.png"
    return Path(str(request.get("image", default)))


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    source = _source(request)
    if not source.is_file():
        status = "blocked_missing_input" if source.suffix.lower() == ".scn" else "failed"
        return {"status": status, "missing": str(source.absolute())}
    if int(request.get("lane_count", 8)) < 1:
        raise ValueError("lane_count must be positive")
    return {"status": "passed", "source": str(source.absolute())}


def _read_scn(path: Path) -> np.ndarray:
    module_path = repo_root() / "agarose_gel_analysis/src/scn_reader.py"
    spec = importlib.util.spec_from_file_location("hurdler_legacy_scn_reader", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the retained SCN reader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    image, _dimensions, error = module.read_scn_file(path)
    if error:
        raise ValueError(error)
    return np.asarray(image)


def _gray(path: Path) -> np.ndarray:
    image = _read_scn(path) if path.suffix.lower() == ".scn" else np.asarray(Image.open(path))
    if image.ndim == 3:
        image = image[..., :3].astype(float).mean(axis=2)
    image = image.astype(float)
    span = float(image.max() - image.min())
    return (image - image.min()) / span if span else np.zeros_like(image)


def _analyze(image: np.ndarray, lane_count: int, invert: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    signal_image = 1.0 - image if invert else image
    height, width = signal_image.shape
    edges = np.linspace(0, width, lane_count + 1, dtype=int)
    bands: list[dict[str, Any]] = []
    lanes: list[dict[str, Any]] = []
    for lane in range(lane_count):
        start, end = int(edges[lane]), int(edges[lane + 1])
        profile = signal_image[:, start:end].mean(axis=1)
        peaks, properties = find_peaks(profile, prominence=max(0.03, float(profile.std()) * 0.5), distance=max(3, height // 50))
        widths = peak_widths(profile, peaks, rel_height=0.5)[0] if len(peaks) else []
        for ordinal, (peak, width_value) in enumerate(zip(peaks, widths, strict=True), 1):
            bands.append({
                "lane": lane + 1, "band": ordinal, "relative_peak_position": peak / max(1, height - 1),
                "relative_band_width": float(width_value) / height, "peak_intensity": float(profile[peak]),
                "prominence": float(properties["prominences"][ordinal - 1]),
            })
        lanes.append({
            "lane": lane + 1, "x_start": start, "x_end": end, "band_count": len(peaks),
            "mean_intensity": float(profile.mean()), "maximum_intensity": float(profile.max()),
        })
    return pd.DataFrame(bands), pd.DataFrame(lanes)


def _plot(image: np.ndarray, lanes: pd.DataFrame, bands: pd.DataFrame, output: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.imshow(image, cmap="gray", aspect="auto")
    for row in lanes.itertuples(index=False):
        ax.axvline(row.x_start, color="#E69F00", linewidth=0.6)
        ax.text((row.x_start + row.x_end) / 2, 10, str(row.lane), color="white", ha="center", va="top", fontsize=7)
    height = image.shape[0]
    for row in bands.itertuples(index=False):
        lane = lanes.iloc[int(row.lane) - 1]
        ax.plot((lane.x_start + lane.x_end) / 2, row.relative_peak_position * (height - 1), "o", color="#4B2E83", markersize=3)
    ax.set_title("Agarose gel lane and band analysis")
    ax.axis("off")
    fig.tight_layout()
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix in ("png", "pdf", "svg"):
        path = output / f"agarose_gel_annotated.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    info = preflight(context, request)
    if info["status"] != "passed":
        report = context.directory("reports") / "agarose_blocked_missing_input.txt"
        report.write_text(f"Required historical input is missing: {info['missing']}\n")
        return result_from_paths(context, backend_id=SPEC.notebook_id, request=request, paths=[report], metrics=info,
                                 status=info["status"], limitations=["No gel measurements were fabricated."])
    image = _gray(Path(info["source"]))
    bands, lanes = _analyze(image, int(request.get("lane_count", 8)), bool(request.get("invert", True)))
    band_paths = write_frame(bands, context.directory("tables") / "agarose_bands")
    lane_paths = write_frame(lanes, context.directory("tables") / "agarose_lanes")
    figures = _plot(image, lanes, bands, context.directory("figures"))
    return result_from_paths(context, backend_id=SPEC.notebook_id, request=request,
                             paths=[*band_paths, *lane_paths, *figures],
                             metrics={"lane_count": len(lanes), "band_count": len(bands), "image_height": image.shape[0], "image_width": image.shape[1]})


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
