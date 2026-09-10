"""Shared output selection, norms, convergence and tolerance loading."""

from __future__ import annotations

import re
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Any, Literal

import numpy as np

from support.documents import load_json
from support.errors import ComparisonError
from support.paths import recorded_file, require_file


TIME_SAVE_RE = re.compile(r"_\d{4}\.h5$")


OUTPUT_RE = re.compile(r"Output written to file\s+(.+\.h5)\s*$", re.MULTILINE)


def select_candidate(
    run_directory: Path,
    metadata: dict[str, Any],
    override: Path | None = None,
) -> Path:
    """Select the final HDF5 result recorded for a completed run."""
    if override is not None:
        return resolve_run_file(run_directory, override, "", "candidate")

    recorded_outputs = metadata.get("hdf5_outputs")
    if not isinstance(recorded_outputs, list) or not recorded_outputs:
        raise ComparisonError("run metadata contains no HDF5 output")
    declared = [
        recorded_file(value, "HDF5 output", run_directory)
        for value in recorded_outputs
    ]

    logged_output = _last_logged_output(run_directory, declared)
    if logged_output is not None:
        return logged_output

    final_outputs = [path for path in declared if not TIME_SAVE_RE.search(path.name)]
    if len(final_outputs) == 1:
        return final_outputs[0]
    if len(declared) == 1:
        return declared[0]
    raise ComparisonError("cannot select one final HDF5 output; use --candidate")


def resolve_run_file(
    run_directory: Path,
    override: Path | None,
    default: str,
    label: str,
) -> Path:
    """Resolve an override or run-relative default file."""
    path = override if override is not None else Path(default)
    path = path.expanduser()
    if not path.is_absolute():
        path = run_directory / path
    return require_file(path, label)


def _last_logged_output(
    run_directory: Path,
    declared: list[Path],
) -> Path | None:
    stdout_path = run_directory / "stdout.log"
    if not stdout_path.is_file():
        return None

    text = stdout_path.read_text(encoding="utf-8", errors="replace")
    for match in reversed(OUTPUT_RE.findall(text)):
        path = Path(match.strip()).expanduser()
        path = (
            path.resolve()
            if path.is_absolute()
            else (run_directory / path).resolve()
        )
        if path in declared:
            return path
    return None


@dataclass(frozen=True)
class ErrorNorms:
    compatible: bool
    finite: bool
    relative_l2: float | None
    normalized_linf: float | None

    @property
    def available(self) -> bool:
        return (
            self.compatible
            and self.finite
            and self.relative_l2 is not None
            and self.normalized_linf is not None
        )


def calculate_error_norms(
    reference: np.ndarray, candidate: np.ndarray
) -> ErrorNorms:
    """Calculate error norms when both arrays are compatible and finite."""
    first = np.asarray(reference, dtype=float)
    second = np.asarray(candidate, dtype=float)
    compatible = first.shape == second.shape
    finite = bool(np.isfinite(first).all() and np.isfinite(second).all())
    if not compatible or not finite:
        return ErrorNorms(compatible, finite, None, None)
    if first.size == 0:
        return ErrorNorms(True, True, None, None)

    difference = second - first
    tiny = np.finfo(float).tiny
    relative_l2 = float(
        np.linalg.norm(difference.ravel())
        / max(float(np.linalg.norm(first.ravel())), tiny)
    )
    normalized_linf = float(
        np.max(np.abs(difference)) / max(float(np.max(np.abs(first))), tiny)
    )
    return ErrorNorms(True, True, relative_l2, normalized_linf)


NEWTON_CONVERGENCE_FAILURE = "final Newton error exceeds tolerance"


NEWTON_FINITE_FAILURE = "final Newton error is missing or non-finite"


_ERROR_RE = re.compile(r"^\s*Error:\s*([-+0-9.eE]+)\s*$", re.MULTILINE)


NewtonCheck = Literal["bounded", "finite_only"]


@dataclass(frozen=True)
class NewtonConvergence:
    """The final Newton error and its acceptance threshold."""

    final_error: float | None
    maximum: float | None

    @property
    def passed(self) -> bool:
        return self.failure is None

    @property
    def failure(self) -> str | None:
        """Explain why the recorded error does not satisfy this check."""
        if self.final_error is None or not isfinite(self.final_error):
            return NEWTON_FINITE_FAILURE
        if self.maximum is not None and self.final_error > self.maximum:
            return NEWTON_CONVERGENCE_FAILURE
        return None

    def as_report(self) -> dict[str, bool | float | None]:
        """Return the stable convergence section used in comparison reports."""
        return {
            "passed": self.passed,
            "final_newton_error": self.final_error,
            "maximum": self.maximum,
        }


def read_newton_convergence(
    log_path: Path,
    maximum: float | None,
) -> NewtonConvergence:
    """Read the last Newton error and optionally enforce an upper bound."""
    try:
        solver_output = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ComparisonError(f"cannot read solver log {log_path}: {exc}") from exc

    errors = _ERROR_RE.findall(solver_output)
    final_error = float(errors[-1]) if errors else None
    return NewtonConvergence(final_error=final_error, maximum=maximum)


def effective_newton_maximum(
    configured_maximum: float | None,
    check: NewtonCheck,
) -> float | None:
    """Apply the stage's declared Newton-convergence policy."""
    if check == "bounded":
        return configured_maximum
    if check == "finite_only":
        return None
    raise ComparisonError(f"unsupported Newton check: {check}")


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
        "mesh_coordinate_atol",
        "relative_l2_max",
        "normalized_linf_max",
    }
    missing = sorted(required - profile.keys())
    if missing:
        raise ComparisonError(f"tolerance profile is missing: {', '.join(missing)}")
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
