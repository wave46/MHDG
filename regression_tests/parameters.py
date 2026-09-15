"""Render selected values into copied MHDG parameter files."""

from __future__ import annotations

import math
import re
from pathlib import Path

from .support import BundleError


ASSIGNMENT_RE = re.compile(
    r"^(?P<prefix>\s*(?P<key>[A-Za-z][A-Za-z0-9_]*)\s*=\s*).*$"
)
def render_parameter_file(
    source: Path,
    destination: Path,
    replacements: dict[str, Path | str],
    parameter_overrides: dict[str, bool | float | int | str] | None = None,
    parameter_namelists: dict[str, str] | None = None,
) -> None:
    """Render selected path and parameter assignments in a copied input file."""
    lines = _read_parameter_lines(source)
    values = _replacement_values(replacements, parameter_overrides)
    rendered, counts = _render_assignments(lines, values)
    rendered = _insert_missing_assignments(rendered, values, counts, parameter_namelists or {})
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
        code, marker, comment = body, "", ""
        for token in re.finditer(r"'[^']*'|\"[^\"]*\"|!", body):
            if token.group() == "!":
                code, marker, comment = body[:token.start()], "!", body[token.end():]
                break
        match = ASSIGNMENT_RE.match(code)
        key = match.group("key").lower() if match else ""
        if key not in values:
            rendered.append(line)
            continue

        scalar_code = re.sub(r"'[^']*'|\"[^\"]*\"", "''", code)
        if re.search(r",\s*[A-Za-z][A-Za-z0-9_]*(?:\([^)]*\))?\s*=", scalar_code) or scalar_code.rstrip().endswith("&"):
            raise BundleError(f"parameter replacement requires one complete assignment per line: {key}")

        suffix = f" !{comment}" if marker else ""
        rendered.append(f"{match.group('prefix')}{values[key]}{suffix}{ending}")
        counts[key] += 1
    return rendered, counts


def _require_single_assignment(counts: dict[str, int]) -> None:
    invalid = [key for key, count in counts.items() if count != 1]
    if invalid:
        details = ", ".join(f"{key} ({counts[key]} matches)" for key in invalid)
        raise BundleError(f"parameter assignments must appear once: {details}")


def _insert_missing_assignments(
    lines: list[str],
    values: dict[str, str],
    counts: dict[str, int],
    namelists: dict[str, str],
) -> list[str]:
    """Insert absent scalars only into explicitly declared namelists."""
    rendered = list(lines)
    for spelling, namelist in namelists.items():
        key = spelling.lower()
        if key not in values or counts[key] != 0:
            continue
        header = re.compile(rf"^\s*&{re.escape(namelist)}\s*(?:!.*)?$", re.IGNORECASE)
        matches = [
            index
            for index, line in enumerate(rendered)
            if header.match(line.rstrip("\r\n"))
        ]
        if len(matches) != 1:
            continue
        index = matches[0] + 1
        ending = "\r\n" if rendered[matches[0]].endswith("\r\n") else "\n"
        rendered.insert(index, f"    {spelling} = {values[key]}{ending}")
        counts[key] = 1
    return rendered


def _write_parameter_file(destination: Path, rendered: list[str]) -> None:
    try:
        destination.write_text("".join(rendered), encoding="utf-8")
    except OSError as exc:
        raise BundleError(f"cannot write parameter file {destination}: {exc}") from exc


def read_selected_input_values(path, names):
    """Read selected scalar/list assignments, returning case-insensitive lowercase keys."""
    names = {name.lower() for name in names}
    values = {}
    for line in _read_parameter_lines(path):
        match = ASSIGNMENT_RE.match(line)
        key = match.group("key").lower() if match else ""
        if key not in names:
            continue
        if key in values:
            raise BundleError(f"parameter assignments must appear once: {key}")
        raw = line[len(match.group("prefix")):].split("!", 1)[0].strip().rstrip(",")
        tokens = re.findall(r"'[^']*'|\"[^\"]*\"|[^,\s]+", raw)
        parsed = []
        for token in tokens:
            if token.lower() in (".true.", ".false."):
                parsed.append(token.lower() == ".true.")
            elif token[:1] in ("'", '"'):
                parsed.append(token[1:-1])
            else:
                try:
                    parsed.append(float(token.replace("D", "E").replace("d", "e")))
                except ValueError as exc:
                    raise BundleError(f"unsupported parameter assignment in {path.name}: {key}") from exc
        if not parsed:
            raise BundleError(f"empty parameter assignment: {key}")
        values[key] = parsed[0] if len(parsed) == 1 else parsed
    return values
