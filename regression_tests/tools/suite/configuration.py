"""Load suite declarations and resolve their summary directories."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition
from bundle.schemas import load_validated_json
from catalogs.layouts import layout_pairs, load_layouts
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
    layouts = load_layouts(layouts_path)
    defaults = document["defaults"]
    relations = declaration.get("relations", [])
    selected = declaration.get(
        "layouts", "all" if relations else [defaults["layout"]]
    )
    selected = list(layouts) if selected == "all" else selected
    unknown_layouts = [layout for layout in selected if layout not in layouts]
    if unknown_layouts:
        raise BundleError(
            f"suite {suite_id} has unknown layouts: {', '.join(unknown_layouts)}"
        )

    pairs = layout_pairs({name: layouts[name] for name in selected}, relations)
    suite = {
        "description": declaration["description"],
        "case_id": declaration.get("case", defaults["case"]),
        "workflow_ids": declaration["workflows"],
        "layouts": list(dict.fromkeys(
            layout for pair in pairs for layout in pair.values()
        )) if pairs else selected,
        "reference_comparisons": declaration.get("reference_comparisons", not pairs),
    }
    if pairs:
        suite["layout_comparisons"] = pairs
        suite["tolerance_profile"] = declaration["tolerance_profile"]
        if "layout_comparison_policy" in declaration:
            suite["layout_comparison_policy"] = declaration["layout_comparison_policy"]

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
