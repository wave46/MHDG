"""Load suite declarations and resolve their summary directories."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition
from bundle.schemas import load_validated_json
from catalogs.layouts import load_layouts
from support.errors import BundleError


def load_suite_definition(
    suite_id: str,
    suites_path: Path,
    layouts_path: Path,
    case_directory: Path,
) -> dict[str, Any]:
    """Load one suite and validate its workflow and layout references."""
    schema_path = suites_path.parent / "schemas" / "suites.schema.json"
    document = load_validated_json(suites_path, schema_path, "suite definitions")
    declaration = document["suites"].get(suite_id)
    if declaration is None:
        available = ", ".join(sorted(document["suites"]))
        raise BundleError(f"unknown suite {suite_id}; available: {available}")
    suite = {
        "description": declaration["description"],
        "case_id": declaration["case"],
        "workflow_ids": declaration["workflows"],
        "layouts": declaration.get("layouts") or _paired_layouts(declaration),
    }
    if "layout_comparisons" in declaration:
        suite["layout_comparisons"] = declaration["layout_comparisons"]
        suite["tolerance_profile"] = declaration["tolerance_profile"]
        for option in ("layout_comparison_policy", "reference_comparisons"):
            if option in declaration:
                suite[option] = declaration[option]

    layouts = load_layouts(layouts_path)
    unknown_layouts = [
        layout for layout in suite["layouts"] if layout not in layouts
    ]
    if unknown_layouts:
        raise BundleError(
            f"suite {suite_id} has unknown layouts: {', '.join(unknown_layouts)}"
        )

    case = load_case_definition(suite["case_id"], case_directory)
    unknown_workflows = [
        workflow_id
        for workflow_id in suite["workflow_ids"]
        if workflow_id not in case["workflows"]
    ]
    if unknown_workflows:
        raise BundleError(
            f"suite {suite_id} has unknown workflows: "
            f"{', '.join(unknown_workflows)}"
        )
    return suite


def _paired_layouts(declaration: dict[str, Any]) -> list[str]:
    pairs = declaration["layout_comparisons"]
    for pair in pairs:
        if pair["baseline"] == pair["candidate"]:
            raise BundleError("a layout cannot be compared with itself")
    return list(
        dict.fromkeys(
            layout
            for pair in pairs
            for layout in (pair["baseline"], pair["candidate"])
        )
    )


def require_bundle_class(
    bundle_root: Path,
    case_directory: Path,
    required: str,
) -> None:
    """Require the source bundle to declare the requested publication class."""
    schema = case_directory.parent / "schemas" / "bundle-manifest.schema.json"
    manifest = load_validated_json(
        bundle_root / "manifest.json",
        schema,
        "bundle manifest",
    )
    actual = manifest.get("bundle_class", "unspecified")
    if actual != required:
        raise BundleError(
            f"golden-check requires bundle_class={required}; found {actual}"
        )


def suite_directory(
    settings: dict[str, str],
    suite_id: str,
    run_id: str,
    resume: bool,
) -> Path:
    """Resolve and create or validate the suite summary directory."""
    configured = settings.get("MHDG_REGRESSION_RUN_ROOT")
    if not configured:
        raise BundleError("settings must define MHDG_REGRESSION_RUN_ROOT")
    root = Path(configured).expanduser()
    if not root.is_absolute():
        raise BundleError("MHDG_REGRESSION_RUN_ROOT must be an absolute path")
    path = root.resolve() / "suites" / suite_id / run_id
    if resume:
        if not path.is_dir():
            raise BundleError(f"suite summary directory does not exist: {path}")
        return path
    try:
        path.mkdir(parents=True)
    except FileExistsError as exc:
        raise BundleError(f"suite summary directory already exists: {path}") from exc
    except OSError as exc:
        raise BundleError(
            f"cannot create suite summary directory {path}: {exc}"
        ) from exc
    return path
