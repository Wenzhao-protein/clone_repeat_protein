"""Backend for V2 notebook 13: SEC MAT parsing and quantitative plots."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any, Mapping

import matplotlib
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.signal import find_peaks

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ..notebook_workspace import NotebookContext, NotebookResult, ProgressCallback
from .common import BackendSpec, repo_root, result_from_paths, write_frame


SPEC = BackendSpec(
    "13_sec",
    "SEC analysis",
    "Read retained or uploaded MAT files, normalize traces and calculate peak/area summaries.",
)


def get_spec() -> dict[str, Any]:
    return SPEC.to_dict()


def _source(request: Mapping[str, Any]) -> Path:
    return Path(str(request.get("mat_file", repo_root() / "SEC/input/20250530_sec_result.mat")))


def preflight(context: NotebookContext, request: Mapping[str, Any]) -> dict[str, Any]:
    source = _source(request)
    if not source.is_file():
        raise FileNotFoundError(source)
    return {"status": "passed", "mat_file": str(source.absolute())}


def _legacy_module():
    path = repo_root() / "SEC/src/sec_utils.py"
    spec = importlib.util.spec_from_file_location("hurdler_retained_sec_utils", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load retained SEC utilities")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _tidy(source: Path) -> pd.DataFrame:
    module = _legacy_module()
    frame = module.read_mat(source)
    frame = module.filter_by_sample_name(frame)
    rows = []
    for row in frame.itertuples(index=False):
        time = np.asarray(row.time, dtype=float).ravel()
        signal = np.asarray(row.signal, dtype=float).ravel()
        for x, y in zip(time, signal, strict=False):
            rows.append({"datetime": row.datetime, "abs_wl": row.abs_wl, "sample_name": row.sample_name, "time": x, "signal": y})
    return pd.DataFrame(rows)


def run(
    context: NotebookContext,
    request: Mapping[str, Any],
    progress_callback: ProgressCallback | None = None,
) -> NotebookResult:
    context.prepare()
    info = preflight(context, request)
    tidy = _tidy(Path(info["mat_file"]))
    if tidy.empty:
        raise ValueError("SEC input contained no usable traces")
    tidy["normalized_signal"] = tidy.groupby(["sample_name", "abs_wl"], sort=False).signal.transform(
        lambda values: (values - values.min()) / (values.max() - values.min()) if values.max() > values.min() else 0.0
    )
    summaries = []
    for (sample, wavelength), rows in tidy.groupby(["sample_name", "abs_wl"], sort=True):
        rows = rows.sort_values("time")
        peaks, _ = find_peaks(rows.normalized_signal.to_numpy(), prominence=0.05)
        summaries.append({
            "sample_name": sample, "abs_wl": wavelength, "point_count": len(rows),
            "peak_count": len(peaks), "primary_peak_time": float(rows.iloc[peaks[np.argmax(rows.iloc[peaks].normalized_signal)]].time) if len(peaks) else np.nan,
            "normalized_area": float(trapezoid(rows.normalized_signal, rows.time)),
        })
    summary = pd.DataFrame(summaries)
    tidy_paths = write_frame(tidy, context.directory("tables") / "sec_tidy")
    summary_paths = write_frame(summary, context.directory("tables") / "sec_peak_summary")
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for (sample, wavelength), rows in tidy.groupby(["sample_name", "abs_wl"], sort=True):
        rows = rows.sort_values("time")
        ax.plot(rows.time, rows.normalized_signal, linewidth=1.0, alpha=0.8, label=f"{sample} | {wavelength}")
    ax.set_xlabel("Elution time")
    ax.set_ylabel("Normalized signal")
    ax.set_title("SEC traces")
    if tidy.groupby(["sample_name", "abs_wl"]).ngroups <= 16:
        ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    figure_paths = []
    for suffix in ("png", "pdf", "svg"):
        path = context.directory("figures") / f"sec_normalized_traces.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white")
        figure_paths.append(path)
    plt.close(fig)
    return result_from_paths(context, backend_id=SPEC.notebook_id, request=request,
                             paths=[*tidy_paths, *summary_paths, *figure_paths],
                             metrics={"trace_count": int(tidy.groupby(["sample_name", "abs_wl"]).ngroups), "point_count": len(tidy), "peak_count": int(summary.peak_count.sum())})


def write_outputs(context: NotebookContext, result: NotebookResult) -> dict[str, Any]:
    return result.to_dict()
