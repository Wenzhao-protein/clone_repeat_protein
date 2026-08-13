"""Structured, credential-free progress events for interactive design runs."""

from __future__ import annotations

import threading
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Mapping


class DesignRunStopped(RuntimeError):
    """Raised at a cooperative safe point after a user requests Stop."""


class DesignRunControl:
    """Thread-safe cooperative pause/resume/stop control for long design runs.

    A pause never mutates the GA state: the worker waits at the next explicit
    safe point and resumes from the same population and RNG state.  Stop is
    deliberately cooperative as well, so an in-flight HTTP request or one
    process-pool fitness batch is allowed to finish before this class raises.
    """

    def __init__(self) -> None:
        self._condition = threading.Condition()
        self._paused = False
        self._stopped = False

    @property
    def paused(self) -> bool:
        with self._condition:
            return self._paused and not self._stopped

    @property
    def stopped(self) -> bool:
        with self._condition:
            return self._stopped

    def pause(self) -> None:
        with self._condition:
            if not self._stopped:
                self._paused = True

    def resume(self) -> None:
        with self._condition:
            self._paused = False
            self._condition.notify_all()

    def stop(self) -> None:
        with self._condition:
            self._stopped = True
            self._paused = False
            self._condition.notify_all()

    def safe_point(self) -> None:
        """Wait while paused and raise once Stop has been requested."""
        with self._condition:
            while self._paused and not self._stopped:
                self._condition.wait(timeout=0.25)
            if self._stopped:
                raise DesignRunStopped("GA run stopped by the user at a safe point")


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
    feedback_round: int | None = None
    max_feedback_rounds: int | None = None
    population_size: int | None = None
    mutation_rate: float | None = None
    crossover_rate: float | None = None
    idt_score: float | None = None
    idt_positive_rules: tuple[str, ...] = ()
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
