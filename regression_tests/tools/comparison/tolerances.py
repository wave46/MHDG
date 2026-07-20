"""Load comparison tolerance profiles declared by regression workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from support.documents import load_json
from support.errors import ComparisonError


def load_fixed_tolerances(
    path: Path,
    workflow: dict[str, Any],
    layout_id: str,
    tolerance_profile_override: str | None,
) -> tuple[str, dict[str, Any]]:
    """Select and validate the tolerance profile for a fixed-mesh run."""
    document = load_json(path, "tolerance definitions")
    profile_id = tolerance_profile_override or workflow.get("tolerance_profile")
    if (
        tolerance_profile_override is None
        and layout_id != workflow.get("default_layout")
    ):
        profile_id = workflow.get("cross_layout_tolerance_profile", profile_id)

    profiles = document.get("profiles", {})
    if not profile_id or profile_id not in profiles:
        raise ComparisonError(f"unknown tolerance profile: {profile_id}")
    profile = profiles[profile_id]

    required = {
        "newton_error_max",
        "mesh_coordinate_atol",
        "relative_l2_max",
        "normalized_linf_max",
    }
    missing = sorted(required - profile.keys())
    if missing:
        raise ComparisonError(f"tolerance profile is missing: {', '.join(missing)}")
    return profile_id, profile
