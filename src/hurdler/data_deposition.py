"""Deterministic file-tree rendering for the SI data deposition."""

from __future__ import annotations

from pathlib import Path


HIDDEN_METADATA_NOTE = (
    "# Hidden CLC-native metadata files (.acl, .clcinfo, .orderlist2) are present "
    "in the repository but omitted from this listing."
)
REPOSITORY_ONLY_NOTE = (
    "# CLC project files (.clc), navigation README.md files, VALIDATION.md, and "
    "this generated listing are repository-only and omitted from the SI tree."
)
EXCLUDED_FILENAMES = frozenset(
    {"README.md", "VALIDATION.md", "directory_structure.txt"}
)


def include_in_si_tree(path: Path) -> bool:
    """Return whether *path* belongs in the human-readable SI tree."""
    if path.name.startswith("."):
        return False
    if path.is_file() and path.suffix.lower() == ".clc":
        return False
    if path.is_file() and path.name in EXCLUDED_FILENAMES:
        return False
    return True


def _tree_lines(directory: Path, prefix: str = "") -> list[str]:
    items = [item for item in directory.iterdir() if include_in_si_tree(item)]
    items.sort(key=lambda item: (not item.is_dir(), item.name.casefold(), item.name))
    lines: list[str] = []
    for index, item in enumerate(items):
        last = index == len(items) - 1
        branch = "└── " if last else "├── "
        continuation = "    " if last else "│   "
        lines.append(f"{prefix}{branch}{item.name}{'/' if item.is_dir() else ''}")
        if item.is_dir():
            lines.extend(_tree_lines(item, prefix + continuation))
    return lines


def build_directory_structure(directory: str | Path) -> str:
    """Render the stable SI listing, including explicit exclusion notes."""
    root = Path(directory).resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    lines = [
        HIDDEN_METADATA_NOTE,
        REPOSITORY_ONLY_NOTE,
        f"{root.name}/",
        *_tree_lines(root),
    ]
    return "\n".join(lines) + "\n"


def write_directory_structure(directory: str | Path) -> Path:
    root = Path(directory).resolve()
    destination = root / "directory_structure.txt"
    destination.write_text(build_directory_structure(root), encoding="utf-8")
    return destination
