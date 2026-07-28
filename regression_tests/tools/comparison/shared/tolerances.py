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
    profile_id = tolerance_profile_override or workflow.get("tolerance_profile")
    if (
        tolerance_profile_override is None
        and layout_id != workflow.get("default_layout")
    ):
        profile_id = workflow.get("cross_layout_tolerance_profile", profile_id)

    if not profile_id:
        raise ComparisonError(f"unknown tolerance profile: {profile_id}")
    profile = _load_profile(path, profile_id)

    required = {
        "newton_error_max",
        "mesh_connectivity",
        "mesh_coordinate_atol",
        "relative_l2_max",
        "normalized_linf_max",
    }
    missing = sorted(required - profile.keys())
    if missing:
        raise ComparisonError(f"tolerance profile is missing: {', '.join(missing)}")
    if profile["mesh_connectivity"] not in ("exact", "numbering_invariant"):
        raise ComparisonError(
            f"invalid mesh connectivity mode: {profile['mesh_connectivity']}"
        )
    return profile_id, profile


def load_adaptive_tolerances(
    path: Path,
    workflow: dict[str, Any],
    tolerance_profile_override: str | None = None,
) -> tuple[str, dict[str, Any]]:
    """Select and validate the tolerance profile for an adaptive run."""
    profile_id = tolerance_profile_override or workflow.get("tolerance_profile")
    profile = _load_profile(path, profile_id)
    required = {
        "newton_error_max",
        "samples_per_element",
        "minimum_point_coverage",
        "solution",
        "gradient",
    }
    if not isinstance(profile, dict) or not required <= profile.keys():
        raise ComparisonError(f"invalid adaptive tolerance profile: {profile_id}")
    for dataset in ("solution", "gradient"):
        limits = profile[dataset]
        if not isinstance(limits, dict) or not {
            "relative_l2_max",
            "normalized_linf_max",
        } <= limits.keys():
            raise ComparisonError(
                f"adaptive tolerance profile has invalid {dataset} limits"
            )
    return profile_id, profile


def _load_profile(path: Path, profile_id: str) -> dict[str, Any]:
    document = load_json(path, "tolerance definitions")
    if document.get("schema_version") != 2:
        raise ComparisonError("tolerance definitions require schema_version 2")
    profiles = document.get("profiles")
    if not isinstance(profiles, dict) or profile_id not in profiles:
        raise ComparisonError(f"unknown tolerance profile: {profile_id}")
    declaration = profiles[profile_id]
    defaults = document.get("defaults")
    if not isinstance(declaration, dict) or not isinstance(defaults, dict):
        raise ComparisonError(f"invalid tolerance profile: {profile_id}")

    profile = {**defaults, **declaration}
    if "solution" not in declaration:
        fixed_defaults = document.get("fixed_defaults")
        if not isinstance(fixed_defaults, dict):
            raise ComparisonError("fixed tolerance defaults are missing")
        profile = {**defaults, **fixed_defaults, **declaration}
    return profile
