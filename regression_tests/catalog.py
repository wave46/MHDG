"""Resolve case/workflow inheritance and generate execution layouts and relations.

Callers share a plain path-to-definition dictionary for one operation. Only
repository declarations belong in it; mutable bundles/builds/runs are read at
their validation boundaries. Omitting it starts a fresh catalog snapshot.
"""

from __future__ import annotations

import re
from copy import deepcopy
from itertools import combinations
from pathlib import Path
from typing import Any

from .documents import load_validated_json
from .support import BundleError
from .config import DEFAULT_MODEL, build_key

MEDIA_TYPES = {
    ".h5": "application/x-hdf5",
    ".json": "application/json",
    ".msh": "application/x-gmsh",
    ".geo": "text/plain",
    ".nml": "text/plain",
    ".txt": "text/plain",
}


def normalize_case_definition(
    case_id: str,
    declaration: dict[str, Any],
    shared: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Derive machine-oriented fields from one readable case declaration."""
    shared = shared or {}
    reference = declaration["reference"]
    return {
        "schema_version": declaration["schema_version"],
        "case_id": case_id,
        "description": declaration["description"],
        "reference_branch": reference["branch"],
        "reference_revision": reference["revision"],
        "bundle_files": _bundle_files(case_id, declaration.get("files", {})),
        "workflows": _workflows(
            declaration["workflows"],
            shared.get("workflows", {}),
            shared.get("sequences", {}),
            {**shared.get("parameter_namelists", {}), **declaration.get("parameter_namelists", {})},
        ),
    }


def _bundle_files(
    case_id: str,
    files: dict[str, dict[str, str]],
) -> dict[str, dict[str, Any]]:
    normalized = {}
    for optional, section in (
        (False, files.get("required", {})),
        (True, files.get("optional", {})),
    ):
        for role, filename in section.items():
            if role in normalized:
                raise BundleError(f"case file role is declared twice: {role}")
            suffix = Path(filename).suffix.lower()
            try:
                media_type = MEDIA_TYPES[suffix]
            except KeyError as exc:
                raise BundleError(
                    f"cannot infer media type for case file: {filename}"
                ) from exc
            normalized[role] = {
                "filename": filename,
                "artifact_id": f"{case_id}_{role}",
                "media_type": media_type,
                "optional": optional,
            }
    return normalized


def _workflows(
    declarations: dict[str, dict[str, Any]],
    shared: dict[str, dict[str, Any]],
    sequences: dict[str, list[dict[str, Any]]],
    parameter_namelists: dict[str, str],
) -> dict[str, dict[str, Any]]:
    # Overlay case choices before resolving parents, so adaptive/feature
    # workflows inherit this case's base inputs and stage adjustments.
    available = {
        **shared,
        **{
            name: _extend_workflow(shared.get(name, {}), changes)
            for name, changes in declarations.items()
        },
    }
    resolved: dict[str, dict[str, Any]] = {}
    active: set[str] = set()

    def resolve(workflow_id: str) -> dict[str, Any]:
        if workflow_id in resolved:
            return resolved[workflow_id]
        if workflow_id in active:
            raise BundleError(f"workflow inheritance cycle at {workflow_id}")
        active.add(workflow_id)
        declaration = available[workflow_id]
        parent_id = declaration.get("extends")
        if parent_id:
            if parent_id not in available:
                raise BundleError(
                    f"workflow {workflow_id} extends unknown workflow {parent_id}"
                )
            changes = deepcopy(
                {
                    key: value
                    for key, value in declaration.items()
                    if key != "extends"
                }
            )
            merged = _extend_workflow(resolve(parent_id), changes)
        else:
            merged = deepcopy(declaration)
        active.remove(workflow_id)
        resolved[workflow_id] = merged
        return merged

    return {
        workflow_id: _expand_workflow(
            _extend_workflow({"parameter_namelists": parameter_namelists}, resolve(workflow_id)),
            sequences,
        )
        for workflow_id in declarations
    }


def _extend_workflow(
    base: dict[str, Any],
    changes: dict[str, Any],
) -> dict[str, Any]:
    """Apply one derived declaration while retaining base parameter values."""
    extended = deepcopy(base)
    extended.update(changes)
    for key in ("parameter_overrides", "parameter_namelists"):
        if key in changes:
            extended[key] = {**base.get(key, {}), **changes[key]}
    if "stage_overrides" in changes:
        previous = base.get("stage_overrides", {})
        extended["stage_overrides"] = {
            **deepcopy(previous),
            **{
                name: _extend_workflow(previous.get(name, {}), stage)
                for name, stage in changes["stage_overrides"].items()
            },
        }
    return extended


def _expand_workflow(
    declaration: dict[str, Any],
    sequences: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    """Expand stage sequences and defaults without renaming declaration fields."""
    missing = [key for key in ("type", "description") if key not in declaration]
    if missing:
        raise BundleError("resolved workflow is missing: " + ", ".join(missing))
    workflow = deepcopy(declaration)
    workflow.pop("extends", None)
    workflow.setdefault("model", DEFAULT_MODEL)
    workflow.setdefault("inputs", [])
    stages = []
    for entry in workflow.pop("stages", []):
        if isinstance(entry, str):
            if entry not in sequences:
                raise BundleError(f"unknown stage sequence: {entry}")
            stages.extend(sequences[entry])
        else:
            stages.append(entry)
    stage_changes = workflow.pop("stage_overrides", {})
    stage_ids = {stage["id"] for stage in stages}
    unknown = sorted(stage_changes.keys() - stage_ids)
    if unknown:
        raise BundleError("stage_overrides contains unknown stages: " + ", ".join(unknown))
    adaptive_stages = set(workflow.pop("adaptive_stages", []))
    if stages:
        unknown = sorted(adaptive_stages - stage_ids)
        if unknown:
            raise BundleError("adaptive_stages contains unknown stages: " + ", ".join(unknown))
        workflow["stages"] = []
        for index, declaration in enumerate(stages):
            stage = _extend_workflow(declaration, stage_changes.get(declaration["id"], {}))
            stage["restart_from"] = "analytical" if index == 0 else "previous_stage"
            stage.setdefault("newton_check", "bounded")
            if stage["id"] in adaptive_stages:
                stage.setdefault("parameter_overrides", {})["rest_adapt"] = True
            workflow["stages"].append(stage)
    return workflow


def load_case_definition(case_id: str, case_directory: Path, *, catalog=None) -> dict[str, Any]:
    """Load and validate one tracked case definition by its identifier."""
    if not re.fullmatch(r"[a-z][a-z0-9_]*", case_id):
        raise BundleError(f"invalid case identifier: {case_id}")

    catalog = {} if catalog is None else catalog
    path = case_directory / f"{case_id}.json"
    if path in catalog:
        return catalog[path]
    schema_path = case_directory.parent / "schemas" / "case.schema.json"
    declaration = load_validated_json(
        path,
        schema_path,
        f"case definition {path.name}",
    )
    shared_path = case_directory.parent / "workflows.json"
    if shared_path not in catalog:
        catalog[shared_path] = load_validated_json(
            shared_path, schema_path, "shared workflows", definition="workflow_catalog",
        )
    shared = catalog[shared_path]
    case = normalize_case_definition(case_id, declaration, shared)
    _validate_case_workflows(case)
    catalog[path] = case
    return case


def required_case_roles(
    case: dict[str, Any],
    workflow_id: str | None = None,
    *, references: bool = True,
) -> set[str]:
    """Return base package roles or the roles required by one workflow."""
    if workflow_id is None:
        return {
            role
            for role, spec in case.get("bundle_files", {}).items()
            if not spec.get("optional", False)
        }
    try:
        workflow = case["workflows"][workflow_id]
    except KeyError as exc:
        raise BundleError(
            f"case {case['case_id']} has no workflow {workflow_id}"
        ) from exc
    return workflow_required_roles(workflow, references=references)


def workflow_required_roles(workflow: dict[str, Any], *, references: bool = True) -> set[str]:
    """Resolve execution inputs, optionally including comparison references."""
    roles = set(workflow.get("inputs", []))
    for name in (
        "mesh",
        "parameters",
        "transport",
        "restart",
        "impurity_configuration",
    ):
        if workflow.get(name):
            roles.add(workflow[name])
    if references and workflow.get("reference"):
        roles.add(workflow["reference"])
    for stage in workflow.get("stages", []):
        roles.add(stage["parameters"])
        if stage.get("transport"):
            roles.add(stage["transport"])
    return roles


def _validate_case_workflows(case: dict[str, Any]) -> None:
    declared_roles = set(case.get("bundle_files", {}))
    for workflow_id, workflow in case["workflows"].items():
        used_roles = workflow_required_roles(workflow) | set(
            workflow.get("outputs", [])
        )
        unknown = (
            sorted(used_roles - declared_roles)
            if declared_roles
            else []
        )
        if unknown:
            raise BundleError(
                f"workflow {workflow_id} uses undefined artifact roles: "
                + ", ".join(unknown)
            )
        if workflow["type"] not in {
            "staged_fixed_mesh",
            "staged_adaptive_mesh",
        }:
            continue

        missing = [
            name
            for name in ("mesh", "stages")
            if not workflow.get(name)
        ]
        if missing:
            raise BundleError(
                f"workflow {workflow_id} is missing: {', '.join(missing)}"
            )
        stages = workflow["stages"]
        stage_ids = [stage["id"] for stage in stages]
        if len(stage_ids) != len(set(stage_ids)):
            raise BundleError(f"workflow {workflow_id} has duplicate stage IDs")


SERIAL_LAYOUT = re.compile(r"serial_omp([1-9][0-9]*)")


MPI_LAYOUT = re.compile(r"mpi([1-9][0-9]*)_omp([1-9][0-9]*)")


def load_layouts(path: Path, *, catalog=None) -> dict[str, dict[str, Any]]:
    """Load every declared layout and derive its execution parameters."""
    catalog = {} if catalog is None else catalog
    if path in catalog:
        return catalog[path]
    schema = path.parent / "schemas/layouts.schema.json"
    document = load_validated_json(path, schema, "layout definitions")
    catalog[path] = {layout_id: _layout(layout_id) for layout_id in document["layouts"]}
    return catalog[path]


def load_layout(layout_id: str, path: Path, *, catalog=None) -> dict[str, Any]:
    """Return one layout or report the available identifiers."""
    layouts = load_layouts(path, catalog=catalog)
    try:
        return layouts[layout_id]
    except KeyError as exc:
        available = ", ".join(layouts)
        raise BundleError(
            f"unknown layout {layout_id}; available: {available}"
        ) from exc


def layout_pairs(
    layouts: dict[str, dict[str, Any]], relations: list[str],
) -> list[dict[str, str]]:
    """Generate requested comparisons, sharing each unordered pair once."""
    pairs = []
    seen = set()
    for relation in relations:
        if relation == "all_pairs":
            matches = list(combinations(layouts, 2))
        elif relation in {"openmp", "mpi", "hybrid"}:
            matches = []
            for baseline, first in layouts.items():
                if first["omp_threads"] != 1:
                    continue
                for candidate, second in layouts.items():
                    if relation == "openmp":
                        match = (
                            first["execution"] == second["execution"] == "serial"
                            and second["omp_threads"] > 1
                        )
                    elif relation == "mpi":
                        match = (
                            first["execution"] == "serial"
                            and second["execution"] == "mpi"
                            and second["omp_threads"] == 1
                        )
                    else:
                        match = (
                            first["execution"] == second["execution"] == "mpi"
                            and first["mpi_ranks"] == second["mpi_ranks"]
                            and second["omp_threads"] > 1
                        )
                    if match:
                        matches.append((baseline, candidate))
        else:
            raise BundleError(f"unknown layout relation: {relation}")
        if not matches:
            raise BundleError(f"layout relation {relation} has no matching layouts")
        for baseline, candidate in matches:
            identity = frozenset((baseline, candidate))
            if identity not in seen:
                pairs.append({"baseline": baseline, "candidate": candidate})
                seen.add(identity)
    return pairs


def _layout(layout_id: str) -> dict[str, Any]:
    serial = SERIAL_LAYOUT.fullmatch(layout_id)
    if serial:
        threads = int(serial.group(1))
        return {
            "description": f"Serial solver with {threads} OpenMP thread(s)",
            "execution": "serial",
            "mpi_ranks": 1,
            "omp_threads": threads,
        }

    mpi = MPI_LAYOUT.fullmatch(layout_id)
    if mpi:
        ranks, threads = (int(value) for value in mpi.groups())
        return {
            "description": (
                f"Parallel solver with {ranks} MPI rank(s) and "
                f"{threads} OpenMP thread(s)"
            ),
            "execution": "mpi",
            "mpi_ranks": ranks,
            "omp_threads": threads,
        }
    raise BundleError(f"invalid layout identifier: {layout_id}")


def load_suite_definition(
    suite_id: str,
    suites_path: Path,
    layouts_path: Path,
    case_directory: Path,
    *, case_id: str | None = None, catalog=None,
) -> dict[str, Any]:
    """Load one suite and validate its workflow and layout references."""
    catalog = {} if catalog is None else catalog
    document = _suite_document(suites_path, catalog)
    declaration = document["suites"].get(suite_id)
    if declaration is None:
        available = ", ".join(sorted(document["suites"]))
        raise BundleError(f"unknown suite {suite_id}; available: {available}")
    layouts = load_layouts(layouts_path, catalog=catalog)
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
        "case_id": case_id or declaration.get("case", defaults["case"]),
        "diagnostics": declaration.get("diagnostics", "off"),
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

    case = load_case_definition(suite["case_id"], case_directory, catalog=catalog)
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


def _suite_document(path, catalog):
    if path not in catalog:
        catalog[path] = load_validated_json(path, path.parent / "schemas/suites.schema.json", "suite definitions")
    return catalog[path]


def load_selection(name, suites_path, layouts_path, case_directory, *, case_id=None, catalog=None):
    """Expand a profile into unique case/suite selections; focused suites use the same loader."""
    catalog = {} if catalog is None else catalog
    document = _suite_document(suites_path, catalog)
    profiles = document.get("profiles", {})
    if profiles.keys() & document["suites"].keys():
        raise BundleError("profile and suite names must be distinct")
    is_profile = name in profiles
    if is_profile and case_id:
        raise BundleError("--case applies to a focused suite; profiles declare their cases")

    def expand(profile, ancestors=()):
        if profile in ancestors:
            raise BundleError(f"cyclic profile inclusion: {' -> '.join((*ancestors, profile))}")
        if profile not in profiles:
            raise BundleError(f"unknown included profile: {profile}")
        for parent in profiles[profile].get("include", []):
            yield from expand(parent, (*ancestors, profile))
        yield from profiles[profile]["checks"]

    entries = expand(name) if is_profile else [{"suite": name, "case": case_id}]
    selected = {}
    for entry in entries:
        suite = load_suite_definition(
            entry["suite"], suites_path, layouts_path, case_directory, case_id=entry.get("case"), catalog=catalog,
        )
        selected[(entry["suite"], suite["case_id"])] = {"suite_id": entry["suite"], **suite}
    if not selected:
        raise BundleError(f"selection {name} contains no checks")
    return is_profile, list(selected.values())


def required_builds(case, workflows, layouts):
    """Derive the unique model/execution combinations for selected cells."""
    for name in workflows:
        if name not in case["workflows"]:
            raise BundleError(f"case {case['case_id']} has no workflow {name}")
    requirements = {(case["workflows"][name]["model"], layout["execution"])
                    for name in workflows for layout in layouts}
    for requirement in requirements:
        build_key(*requirement)
    return requirements


def selection_builds(checks, case_directory, layouts, *, catalog=None):
    catalog = {} if catalog is None else catalog
    return set().union(*(required_builds(load_case_definition(check["case_id"], case_directory, catalog=catalog),
                                        check["workflow_ids"], [layouts[name] for name in check["layouts"]])
                         for check in checks))
