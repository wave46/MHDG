"""Read regression settings and resolve their external bundle root."""

from __future__ import annotations

import re
from pathlib import Path

from support.errors import BundleError


SETTING_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")


def read_settings(settings_path: Path) -> dict[str, str]:
    """Parse KEY=VALUE settings without executing shell code."""
    try:
        lines = settings_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise BundleError(f"cannot read settings file {settings_path}: {exc}") from exc

    settings: dict[str, str] = {}
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise BundleError(
                f"{settings_path}:{line_number}: expected a KEY=VALUE setting"
            )
        key, value = (part.strip() for part in line.split("=", 1))
        if not SETTING_RE.fullmatch(key):
            raise BundleError(f"{settings_path}:{line_number}: invalid setting name")
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
            value = value[1:-1]
        settings[key] = value
    return settings


def bundle_root_from_settings(settings: dict[str, str]) -> Path:
    """Resolve the required bundle root from parsed settings."""
    if settings.get("MHDG_REGRESSION_SETTINGS_VERSION") != "2":
        raise BundleError("settings must define MHDG_REGRESSION_SETTINGS_VERSION=2")
    data_root = settings.get("MHDG_REGRESSION_DATA_ROOT")
    if not data_root:
        raise BundleError("settings must define MHDG_REGRESSION_DATA_ROOT")

    root = Path(data_root).expanduser()
    if not root.is_absolute():
        raise BundleError("MHDG_REGRESSION_DATA_ROOT must be an absolute path")
    try:
        root = root.resolve(strict=True)
    except FileNotFoundError as exc:
        raise BundleError(f"bundle root does not exist: {root}") from exc
    if not root.is_dir():
        raise BundleError(f"bundle root is not a directory: {root}")
    return root
