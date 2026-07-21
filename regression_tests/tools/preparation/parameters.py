"""Render selected values into copied MHDG parameter files."""

from __future__ import annotations

import math
import re
from pathlib import Path

from support.errors import BundleError


ASSIGNMENT_RE = re.compile(
    r"^(?P<prefix>\s*(?P<key>[A-Za-z][A-Za-z0-9_]*)\s*=\s*).*$"
)


def render_parameter_file(
    source: Path,
    destination: Path,
    replacements: dict[str, Path | str],
    parameter_overrides: dict[str, bool | float | int | str] | None = None,
) -> None:
    """Render selected path and parameter assignments in a copied input file."""
    lines = _read_parameter_lines(source)
    values = _replacement_values(replacements, parameter_overrides)
    rendered, counts = _render_assignments(lines, values)
    _require_single_assignment(counts)
    _write_parameter_file(destination, rendered)


def _read_parameter_lines(source: Path) -> list[str]:
    try:
        return source.read_text(encoding="utf-8").splitlines(keepends=True)
    except OSError as exc:
        raise BundleError(f"cannot read parameter file {source}: {exc}") from exc


def _replacement_values(
    replacements: dict[str, Path | str],
    parameter_overrides: dict[str, bool | float | int | str] | None,
) -> dict[str, str]:
    values = {}
    for key, raw_value in replacements.items():
        value = str(raw_value)
        if "'" in value:
            raise BundleError(f"cannot render a path containing a quote: {value}")
        values[key.lower()] = f"'{value}'"

    for key, value in (parameter_overrides or {}).items():
        normalized = key.lower()
        if normalized in values:
            raise BundleError(f"duplicate parameter replacement: {key}")
        values[normalized] = _format_parameter_value(key, value)
    return values


def _format_parameter_value(key: str, value: bool | float | int | str) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise BundleError(f"parameter override must be finite: {key}")
        return repr(value)
    if isinstance(value, str):
        if "'" in value or "\n" in value or "\r" in value:
            raise BundleError(f"parameter override contains invalid text: {key}")
        return f"'{value}'"
    raise BundleError(f"unsupported parameter override: {key}")


def _render_assignments(
    lines: list[str],
    values: dict[str, str],
) -> tuple[list[str], dict[str, int]]:
    counts = dict.fromkeys(values, 0)
    rendered = []
    for line in lines:
        body = line.rstrip("\r\n")
        ending = line[len(body) :]
        code, marker, comment = body.partition("!")
        match = ASSIGNMENT_RE.match(code)
        key = match.group("key").lower() if match else ""
        if key not in values:
            rendered.append(line)
            continue

        suffix = f" !{comment}" if marker else ""
        rendered.append(f"{match.group('prefix')}{values[key]}{suffix}{ending}")
        counts[key] += 1
    return rendered, counts


def _require_single_assignment(counts: dict[str, int]) -> None:
    invalid = [key for key, count in counts.items() if count != 1]
    if invalid:
        details = ", ".join(f"{key} ({counts[key]} matches)" for key in invalid)
        raise BundleError(f"parameter assignments must appear once: {details}")


def _write_parameter_file(destination: Path, rendered: list[str]) -> None:
    try:
        destination.write_text("".join(rendered), encoding="utf-8")
    except OSError as exc:
        raise BundleError(f"cannot write parameter file {destination}: {exc}") from exc
