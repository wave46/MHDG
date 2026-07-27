"""Load JSON documents validated against regression schemas."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from support.documents import load_json
from support.errors import BundleError

try:
    from jsonschema import Draft202012Validator, FormatChecker
    from jsonschema.exceptions import SchemaError
except ImportError as exc:  # pragma: no cover - depends on the local environment
    raise SystemExit(
        "jsonschema is required; install regression_tests/requirements.txt"
    ) from exc


def load_validated_json(
    path: Path,
    schema_path: Path,
    label: str,
) -> dict[str, Any]:
    """Load a JSON object and validate it against a regression schema."""
    document = load_json(path, label)
    _validate_json(document, schema_path, label)
    return document


def _validate_json(document: Any, schema_path: Path, label: str) -> None:
    schema = load_json(schema_path, f"schema {schema_path.name}")
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
        return
    error = errors[0]
    location = ".".join(str(part) for part in error.absolute_path)
    where = f".{location}" if location else ""
    raise BundleError(f"{label}{where}: {error.message}")
