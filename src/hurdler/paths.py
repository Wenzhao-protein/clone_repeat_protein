"""Repository and exact-scratch path handling."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProjectPaths:
    root: Path

    @classmethod
    def discover(cls, start: str | Path | None = None) -> "ProjectPaths":
        current = Path(start or __file__).absolute()
        if current.is_file():
            current = current.parent
        for candidate in (current, *current.parents):
            if (candidate / "pyproject.toml").exists() and (candidate / "src" / "hurdler").exists():
                return cls(candidate)
        raise FileNotFoundError("Could not locate the HURDLER repository root")

    @property
    def output(self) -> Path:
        return self.root / "output"

    @property
    def reference_output(self) -> Path:
        return self.root / "data" / "reference_output"

    @property
    def scratch_root(self) -> Path:
        absolute = self.root.absolute()
        try:
            relative = absolute.relative_to("/home")
        except ValueError as exc:
            raise ValueError(f"Project is not below /home: {absolute}") from exc
        return Path("/net/scratch") / relative
