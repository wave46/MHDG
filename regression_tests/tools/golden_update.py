#!/usr/bin/env python3
"""Run or resume an ordered golden-reference update."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from build.configuration import parse_build_jobs
from build.workflow import build_solver
from bundle.promotion import promote_bundle, promote_mapped_bundle
from bundle.schemas import load_validated_json
from bundle.settings import bundle_root_from_settings, read_settings
from bundle.validation import validate_bundle_root
from suite.configuration import require_bundle_class
from suite.pairs import compare_layout_pairs
from suite.runner import run_suite
from suite.verification import verify_suite
from support.documents import load_json, write_json_atomic
from support.errors import BundleError, HarnessError
from support.files import file_identity
from support.identifiers import IDENTIFIER_RE
from support.paths import require_file
from support.time import utc_now


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = REGRESSION_ROOT.parent
STATE_FILE = "campaign.json"


def main(argv: list[str] | None = None) -> int:
    parser = _argument_parser()
    args = parser.parse_args(argv)
    try:
        if args.action == "status":
            state = load_json(args.workspace / STATE_FILE, "golden campaign state")
        else:
            state = update_campaign(args)
    except HarnessError as exc:
        print(f"golden update failed: {exc}", file=sys.stderr)
        return 1
    _print_status(state)
    return 0


def update_campaign(args: argparse.Namespace) -> dict[str, Any]:
    """Create or continue one campaign until its next acceptance gate."""
    if not IDENTIFIER_RE.fullmatch(args.run_id):
        raise BundleError(f"invalid golden run identifier: {args.run_id}")
    if not args.bundle_version.strip():
        raise BundleError("golden bundle version must not be empty")

    declaration = _load_declaration(args.campaigns, args.case_id)
    source_settings = require_file(args.settings, "source settings")
    settings = read_settings(source_settings)
    source_bundle = bundle_root_from_settings(settings)
    validate_bundle_root(source_bundle, args.cases)
    require_bundle_class(source_bundle, args.cases, "golden")
    source_manifest = load_json(source_bundle / "manifest.json", "bundle manifest")
    if source_manifest["case_id"] != args.case_id:
        raise BundleError("campaign and source bundle use different cases")

    selected = _selected_components(args.only, declaration)
    partial = args.only is not None
    workspace = _workspace(args.workspace, settings, args.run_id)
    inputs = {
        "case_id": args.case_id,
        "run_id": args.run_id,
        "source_settings": _file_record(source_settings),
        "source_bundle": str(source_bundle),
        "source_manifest": _file_record(source_bundle / "manifest.json"),
        "campaign_catalog": _file_record(args.campaigns),
        "output": str(args.output.expanduser().resolve()),
        "bundle_version": args.bundle_version,
        "build_jobs": args.build_jobs,
        "only": sorted(selected) if partial else None,
    }
    state = _load_or_create(
        workspace,
        inputs,
        declaration,
        args.output,
        selected,
        partial,
    )
    if args.accept:
        _accept(state, args.accept)
    return _advance(state, args)


def _load_or_create(
    workspace: Path,
    inputs: dict[str, Any],
    declaration: dict[str, Any],
    output: Path,
    selected: set[str],
    partial: bool,
) -> dict[str, Any]:
    state_path = workspace / STATE_FILE
    if state_path.is_file():
        state = load_json(state_path, "golden campaign state")
        if state.get("inputs") != inputs:
            raise BundleError("golden campaign inputs changed")
        return state
    if workspace.exists():
        raise BundleError(f"campaign workspace already exists: {workspace}")
    if output.expanduser().resolve().exists():
        raise BundleError(f"golden output already exists: {output}")

    workspace.mkdir(parents=True)
    state = {
        "schema_version": 1,
        "workspace": str(workspace),
        "status": "ready",
        "inputs": inputs,
        "build": {"status": "pending"},
        "active_bundle": inputs["source_bundle"],
        "active_bundle_class": "golden",
        "active_settings": None,
        "warnings": _selection_warnings(declaration, selected, partial),
        "stages": [
            {
                **stage,
                "status": (
                    "pending"
                    if not partial
                    or stage["kind"] == "verification"
                    or stage.get("component") in selected
                    else "skipped"
                ),
            }
            for stage in declaration["stages"]
        ],
    }
    _save(state)
    return state


def _advance(state: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    if state["build"]["status"] != "completed":
        _run_build(state, args)

    for stage in state["stages"]:
        if stage["status"] == "awaiting_acceptance":
            state["status"] = "awaiting_acceptance"
            _save(state)
            return state
        if stage["status"] in {"completed", "skipped"}:
            continue
        if stage["status"] == "failed":
            raise BundleError(f"golden stage failed: {stage['id']}")
        _run_stage(state, stage, args)
        if stage["status"] == "awaiting_acceptance":
            return state

    state["status"] = "stages_completed"
    _save(state)
    return state


def _run_build(state: dict[str, Any], args: argparse.Namespace) -> None:
    state["build"]["status"] = "running"
    state["status"] = "running"
    _save(state)
    result = build_solver(args.settings, args.repository_root, args.build_jobs)
    state["build"] = {
        "status": "completed",
        "directory": str(result.path),
        "settings": _file_record(result.settings_path),
        "metadata": str(result.metadata_path),
    }
    state["active_settings"] = state["build"]["settings"]
    _save(state)


def _run_stage(
    state: dict[str, Any],
    stage: dict[str, Any],
    args: argparse.Namespace,
) -> None:
    resume = stage["status"] == "running"
    stage["status"] = "running"
    state["status"] = "running"
    _save(state)

    settings_record = state["active_settings"]
    settings_path = Path(settings_record["path"])
    if _file_record(settings_path) != settings_record:
        raise BundleError("campaign settings changed")
    if stage["kind"] == "verification":
        _record_verification_candidate(state)
    summary_path, summary = run_suite(
        settings_path,
        stage["suite"],
        args.cases,
        args.layouts,
        args.suites,
        args.tolerances,
        f"{state['inputs']['run_id']}-{stage['id']}",
        state["active_bundle_class"],
        compare=stage["kind"] not in {"matrix", "reference"},
        resume=resume,
    )
    passed = summary["status"] == "passed"
    stage["summary"] = str(summary_path)
    if stage["kind"] == "matrix" and passed:
        passed = _record_matrix_reports(state, stage, summary_path, summary, args)
    elif (
        stage["kind"] == "reference"
        and stage["acceptance_required"]
        and passed
    ):
        old_golden_path, _ = verify_suite(
            summary_path,
            args.cases,
            args.tolerances,
        )
        stage["old_golden_report"] = str(old_golden_path)
    if not passed:
        stage["status"] = "failed"
        state["status"] = "failed"
        _save(state)
        raise BundleError(f"golden stage failed: {stage['id']}")

    if stage["kind"] in {"matrix", "reference"}:
        _compose_stage(state, stage, summary_path, settings_path, args)

    stage["status"] = (
        "awaiting_acceptance"
        if stage["acceptance_required"]
        else "completed"
    )
    if stage["status"] == "completed" and stage.get("candidate"):
        _activate_candidate(state, stage)
    state["status"] = stage["status"]
    _save(state)


def _record_matrix_reports(
    state: dict[str, Any],
    stage: dict[str, Any],
    summary_path: Path,
    summary: dict[str, Any],
    args: argparse.Namespace,
) -> bool:
    comparisons = compare_layout_pairs(summary, args.cases, args.tolerances)
    report = {
        "schema_version": 1,
        "status": (
            "passed"
            if all(item["status"] == "passed" for item in comparisons)
            else "failed"
        ),
        "comparisons": comparisons,
    }
    report_path = Path(state["workspace"]) / "reports" / f"{stage['id']}.json"
    write_json_atomic(report_path, report, "matrix layout-pair report")
    stage["layout_pair_report"] = str(report_path)
    if report["status"] != "passed":
        return False
    old_golden_path, _ = verify_suite(
        summary_path,
        args.cases,
        args.tolerances,
        include_layout_pairs=False,
    )
    stage["old_golden_report"] = str(old_golden_path)
    return True


def _accept(state: dict[str, Any], stage_id: str) -> None:
    matching = [stage for stage in state["stages"] if stage["id"] == stage_id]
    if len(matching) != 1 or matching[0]["status"] != "awaiting_acceptance":
        raise BundleError(f"golden stage is not awaiting acceptance: {stage_id}")
    matching[0]["status"] = "completed"
    matching[0]["accepted_utc"] = utc_now()
    if matching[0].get("candidate"):
        _activate_candidate(state, matching[0])
    state["status"] = "ready"
    _save(state)


def _load_declaration(path: Path, case_id: str) -> dict[str, Any]:
    document = load_validated_json(
        path,
        path.parent / "schemas/golden-campaigns.schema.json",
        "golden campaign catalog",
    )
    declaration = document["campaigns"].get(case_id)
    if declaration is None:
        raise BundleError(f"no golden campaign is declared for {case_id}")
    if declaration["case"] != case_id:
        raise BundleError(f"golden campaign {case_id} declares another case")
    ids = [stage["id"] for stage in declaration["stages"]]
    if len(ids) != len(set(ids)):
        raise BundleError("golden campaign contains duplicate stage identifiers")
    if any(
        stage["kind"] == "reference" and not stage.get("role_mappings")
        for stage in declaration["stages"]
    ):
        raise BundleError("reference stages require role_mappings")
    return declaration


def _compose_stage(
    state: dict[str, Any],
    stage: dict[str, Any],
    summary_path: Path,
    settings_path: Path,
    args: argparse.Namespace,
) -> None:
    candidate = Path(state["workspace"]) / "candidates" / stage["id"]
    stage["candidate"] = str(candidate)
    _save(state)
    if candidate.exists():
        validate_bundle_root(candidate, args.cases)
    elif stage["kind"] == "matrix":
        promote_bundle(
            settings_path,
            summary_path,
            candidate,
            f"{state['inputs']['bundle_version']}-{stage['id']}",
            args.cases,
            bundle_class="candidate",
            matrix_warm_roles=("warm_restart",),
        )
    else:
        promote_mapped_bundle(
            settings_path,
            summary_path,
            stage["role_mappings"],
            candidate,
            f"{state['inputs']['bundle_version']}-{stage['id']}",
            args.cases,
        )

    generated = Path(state["workspace"]) / "settings" / f"{stage['id']}.env"
    _write_bundle_settings(settings_path, generated, candidate)
    stage["candidate_settings"] = _file_record(generated)
    _save(state)


def _activate_candidate(state: dict[str, Any], stage: dict[str, Any]) -> None:
    state["active_bundle"] = stage["candidate"]
    state["active_bundle_class"] = "candidate"
    state["active_settings"] = stage["candidate_settings"]


def _record_verification_candidate(state: dict[str, Any]) -> None:
    root = Path(state["active_bundle"])
    candidate = {
        "root": str(root),
        "manifest": _file_record(root / "manifest.json"),
    }
    recorded = state.get("verification_candidate")
    if recorded is not None and recorded != candidate:
        raise BundleError("verification candidate manifest changed")
    state["verification_candidate"] = candidate
    _save(state)


def _write_bundle_settings(template: Path, output: Path, bundle: Path) -> None:
    settings = read_settings(template)
    settings["MHDG_REGRESSION_DATA_ROOT"] = str(bundle.resolve())
    contents = "".join(
        f"{key}={value}\n" for key, value in sorted(settings.items())
    )
    if output.exists():
        if output.read_text(encoding="utf-8") != contents:
            raise BundleError(f"generated campaign settings changed: {output}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".env.tmp")
    try:
        temporary.write_text(contents, encoding="utf-8")
        temporary.replace(output)
    except OSError as exc:
        raise BundleError(f"cannot write campaign settings {output}: {exc}") from exc


def _workspace(
    configured: Path | None,
    settings: dict[str, str],
    run_id: str,
) -> Path:
    if configured is not None:
        path = configured.expanduser()
    else:
        run_root = settings.get("MHDG_REGRESSION_RUN_ROOT")
        if not run_root:
            raise BundleError("settings must define MHDG_REGRESSION_RUN_ROOT")
        path = Path(run_root).expanduser() / "golden_campaigns" / run_id
    if not path.is_absolute():
        raise BundleError("golden campaign workspace must be absolute")
    return path.resolve()


def _selected_components(
    requested: list[str] | None,
    declaration: dict[str, Any],
) -> set[str]:
    available = {
        stage["component"]
        for stage in declaration["stages"]
        if "component" in stage
    }
    if requested is None:
        return available
    selected = set(requested)
    unknown = sorted(selected - available)
    if unknown:
        raise BundleError(
            "unknown golden update component: " + ", ".join(unknown)
        )
    return selected


def _selection_warnings(
    declaration: dict[str, Any],
    selected: set[str],
    partial: bool,
) -> list[str]:
    if not partial:
        return []
    warnings = []
    for group in declaration.get("update_together", []):
        chosen = selected.intersection(group)
        if chosen and chosen != set(group):
            warnings.append(
                "components normally updated together were split: "
                + ", ".join(group)
            )
    return warnings


def _file_record(path: Path) -> dict[str, Any]:
    path = path.expanduser().resolve()
    return {"path": str(path), **file_identity(path)}


def _save(state: dict[str, Any]) -> None:
    write_json_atomic(
        Path(state["workspace"]) / STATE_FILE,
        state,
        "golden campaign state",
    )


def _print_status(state: dict[str, Any]) -> None:
    print(f"campaign: {state['inputs']['run_id']} ({state['status']})")
    print(f"workspace: {state['workspace']}")
    print(f"build: {state['build']['status']}")
    for stage in state["stages"]:
        print(f"{stage['id']}: {stage['status']}")
    for warning in state.get("warnings", []):
        print(f"warning: {warning}")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    update = actions.add_parser("update", help="start or continue an update")
    update.add_argument("case_id", metavar="CASE")
    update.add_argument("--settings", required=True, type=Path)
    update.add_argument("--run-id", required=True)
    update.add_argument("--output", required=True, type=Path)
    update.add_argument("--bundle-version", required=True, metavar="VERSION")
    update.add_argument("--workspace", type=Path)
    update.add_argument("--accept", metavar="STAGE")
    update.add_argument("--only", action="append", metavar="COMPONENT")
    update.add_argument("--build-jobs", type=parse_build_jobs, metavar="N")
    update.add_argument("--repository-root", type=Path, default=REPOSITORY_ROOT)
    update.add_argument(
        "--campaigns",
        type=Path,
        default=REGRESSION_ROOT / "golden_campaigns.json",
    )
    update.add_argument("--cases", type=Path, default=REGRESSION_ROOT / "cases")
    update.add_argument(
        "--layouts",
        type=Path,
        default=REGRESSION_ROOT / "layouts.json",
    )
    update.add_argument("--suites", type=Path, default=REGRESSION_ROOT / "suites.json")
    update.add_argument(
        "--tolerances",
        type=Path,
        default=REGRESSION_ROOT / "tolerances.json",
    )
    status = actions.add_parser("status", help="show persisted campaign state")
    status.add_argument("workspace", type=Path)
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
