"""JSON document I/O for the regression harness."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from support.errors import DocumentError


def load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise DocumentError(f"cannot read {label} {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise DocumentError(f"invalid JSON in {label} {path}: {exc.msg}") from exc
    if not isinstance(document, dict):
        raise DocumentError(f"{label} must contain a JSON object")
    return document


def write_json_direct(path: Path, document: dict[str, Any]) -> None:
    """Write JSON directly when the surrounding operation owns publication."""
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def write_json_atomic(path: Path, document: dict[str, Any], label: str) -> None:
    path = path.expanduser().resolve()
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)
    except OSError as exc:
        raise DocumentError(f"cannot write {label} {path}: {exc}") from exc
