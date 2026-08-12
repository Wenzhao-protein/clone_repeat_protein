"""Structured, credential-free progress events for interactive design runs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Mapping


@dataclass(frozen=True)
class DesignProgressEvent:
    """One stable progress update emitted by GA, IDT, and RDL orchestration."""

    stage: str
    status: str
    message: str = ""
    fragment_kind: str = ""
    copies: int | None = None
    phase: str = ""
    generations: int | None = None
    generation: int | None = None
    ga_score: float | None = None
    selected_pair_re_site_excess: int | None = None
    elapsed_seconds: float | None = None
    details: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


ProgressCallback = Callable[[DesignProgressEvent], None]


def emit_progress(
    callback: ProgressCallback | None,
    *,
    stage: str,
    status: str,
    **values: Any,
) -> None:
    """Emit an event only when a caller requested progress reporting."""
    if callback is not None:
        callback(DesignProgressEvent(stage=stage, status=status, **values))
