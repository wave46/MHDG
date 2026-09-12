"""JSON reading, schema validation and direct/atomic record writing."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .support import BundleError, DocumentError

def load_json(path: Path, label: str) -> dict[str, Any]:
    """Read one JSON object and reject arrays or scalar documents."""
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




def load_validated_json(
    path: Path,
    schema_path: Path,
    label: str,
    *,
    definition: str | None = None,
) -> dict[str, Any]:
    """Load a JSON object and validate it against a regression schema."""
    try:
        from jsonschema import Draft202012Validator, FormatChecker
        from jsonschema.exceptions import SchemaError
    except ImportError as exc:  # pragma: no cover - depends on the local environment
        raise SystemExit(
            "jsonschema is required; install regression_tests/requirements.txt"
        ) from exc
    document = load_json(path, label)
    schema = load_json(schema_path, f"schema {schema_path.name}")
    if definition is not None:
        schema = {"$defs": schema["$defs"], "$ref": f"#/$defs/{definition}"}
    try:
        Draft202012Validator.check_schema(schema)
    except SchemaError as exc:
        raise BundleError(
            f"invalid regression schema {schema_path}: {exc.message}"
        ) from exc

    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    errors = sorted(
        validator.iter_errors(document),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if not errors:
        return document
    error = errors[0]
    location = ".".join(str(part) for part in error.absolute_path)
    where = f".{location}" if location else ""
    raise BundleError(f"{label}{where}: {error.message}")
