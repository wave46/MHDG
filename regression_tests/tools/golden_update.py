#!/usr/bin/env python3
"""Run or resume an ordered golden-reference update."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from build.configuration import parse_build_jobs
from build.workflow import build_solver
from check_balance_diagnostics import check_suite as check_balance_diagnostics
from bundle.promotion import (
    promote_bundle,
    promote_mapped_bundle,
    publish_campaign_bundle,
)
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
ATTEMPT_RECORD_FIELDS = (
    "summary",
    "layout_pair_report",
    "old_golden_report",
    "balance_diagnostics_report",
)
STAGE_RESULT_FIELDS = ATTEMPT_RECORD_FIELDS + (
    "candidate",
    "candidate_settings",
    "accepted_utc",
    "failed_utc",
)
DECLARATION_STAGE_FIELDS = (
    "id",
    "kind",
    "suite",
    "component",
    "acceptance_required",
    "checker",
    "role_mappings",
)


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
    """Create or continue one campaign until completion or final review."""
    if not IDENTIFIER_RE.fullmatch(args.run_id):
        raise BundleError(f"invalid golden run identifier: {args.run_id}")
    if not args.bundle_version.strip():
        raise BundleError("golden bundle version must not be empty")

    declaration = _load_declaration(args.campaigns, args.case_id)
    source_settings = require_file(args.settings, "source settings")
    settings = read_settings(source_settings)
    source_bundle = bundle_root_from_settings(settings)
    validate_bundle_root(source_bundle, args.cases)
    source_bundle_class = "candidate" if args.bootstrap_candidate else "golden"
    require_bundle_class(source_bundle, args.cases, source_bundle_class)
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
        "source_bundle_class": source_bundle_class,
        "bootstrap_candidate": args.bootstrap_candidate,
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
        allow_catalog_update=args.retry_from is not None,
    )
    if args.retry_from:
        _retry_from_stage(
            state,
            declaration,
            args.retry_from,
            inputs["campaign_catalog"],
        )
    if args.retry_failed:
        _retry_failed_stage(state)
    if args.accept:
        _accept_campaign(state, args.accept)
    return _advance(state, args)


def _retry_failed_stage(state: dict[str, Any]) -> None:
    """Archive one failed stage attempt and make the stage runnable again."""
    failed = [stage for stage in state["stages"] if stage["status"] == "failed"]
    if state["status"] != "failed" or len(failed) != 1:
        raise BundleError("golden campaign has no single failed stage to retry")

    stage = failed[0]
    retry_count = stage.get("retry_count", 0) + 1
    attempt = {
        "attempt": retry_count - 1,
        "archived_utc": utc_now(),
    }
    for name in ATTEMPT_RECORD_FIELDS:
        if name in stage:
            attempt[name] = stage.pop(name)
    if "failed_utc" in stage:
        attempt["failed_utc"] = stage.pop("failed_utc")
    stage.setdefault("failed_attempts", []).append(attempt)
    stage["retry_count"] = retry_count
    stage["status"] = "pending"
    state["status"] = "ready"
    _save(state)


def _retry_from_stage(
    state: dict[str, Any],
    declaration: dict[str, Any],
    stage_id: str,
    campaign_catalog: dict[str, Any],
) -> None:
    """Rewind a failed campaign to a corrected declaration stage."""
    if state["status"] != "failed":
        raise BundleError("golden campaign is not failed")
    stages = state["stages"]
    declared = declaration["stages"]
    if [stage["id"] for stage in stages] != [stage["id"] for stage in declared]:
        raise BundleError("golden campaign stage order changed")
    indices = [index for index, stage in enumerate(stages) if stage["id"] == stage_id]
    if not indices:
        raise BundleError(f"unknown golden retry stage: {stage_id}")
    start = indices[0]
    if not any(stage["status"] == "failed" for stage in stages[start:]):
        raise BundleError("golden retry stage is after the recorded failure")

    _restore_preceding_candidate(state, start)
    for stage in stages[start:]:
        if stage["status"] == "skipped":
            continue
        if any(name in stage for name in ATTEMPT_RECORD_FIELDS):
            _archive_stage_attempt(stage)
            stage["retry_count"] = stage.get("retry_count", 0) + 1
        for name in STAGE_RESULT_FIELDS:
            stage.pop(name, None)
        stage["status"] = "pending"
    for stage, current in zip(stages, declared):
        _update_stage_declaration(stage, current)

    previous_catalog = state["inputs"]["campaign_catalog"]
    state.setdefault("campaign_amendments", []).append(
        {
            "amended_utc": utc_now(),
            "retry_from": stage_id,
            "previous_campaign_catalog": previous_catalog,
            "campaign_catalog": campaign_catalog,
        }
    )
    state["inputs"]["campaign_catalog"] = campaign_catalog
    state["acceptance"] = {
        "status": "pending",
        "required_stages": [
            stage["id"]
            for stage in stages
            if stage["status"] != "skipped" and stage["acceptance_required"]
        ],
    }
    state["publication"] = {"status": "pending"}
    state.pop("verification_candidate", None)
    state["status"] = "ready"
    _save(state)


def _restore_preceding_candidate(state: dict[str, Any], start: int) -> None:
    state["active_bundle"] = state["inputs"]["source_bundle"]
    state["active_bundle_class"] = state["inputs"]["source_bundle_class"]
    state["active_settings"] = state["build"].get("settings")
    for stage in reversed(state["stages"][:start]):
        if stage.get("candidate") and stage["status"] == "completed":
            _activate_candidate(state, stage)
            return


def _update_stage_declaration(
    stage: dict[str, Any], declaration: dict[str, Any]
) -> None:
    for name in DECLARATION_STAGE_FIELDS:
        stage.pop(name, None)
        if name in declaration:
            stage[name] = declaration[name]


def _archive_stage_attempt(stage: dict[str, Any]) -> None:
    attempt = {
        "attempt": stage.get("retry_count", 0),
        "status": stage["status"],
        "archived_utc": utc_now(),
        "declaration": {
            name: stage[name]
            for name in DECLARATION_STAGE_FIELDS
            if name in stage
        },
    }
    for name in STAGE_RESULT_FIELDS:
        if name in stage:
            attempt[name] = stage[name]
    stage.setdefault("archived_attempts", []).append(attempt)


def _load_or_create(
    workspace: Path,
    inputs: dict[str, Any],
    declaration: dict[str, Any],
    output: Path,
    selected: set[str],
    partial: bool,
    allow_catalog_update: bool = False,
) -> dict[str, Any]:
    state_path = workspace / STATE_FILE
    if state_path.is_file():
        state = load_json(state_path, "golden campaign state")
        if state.get("inputs") != inputs:
            catalog_only = _only_campaign_catalog_changed(
                state.get("inputs"), inputs
            )
            if catalog_only and allow_catalog_update:
                return state
            if catalog_only and _declaration_matches_state(state, declaration):
                _refresh_campaign_catalog(
                    state,
                    inputs["campaign_catalog"],
                )
            else:
                raise BundleError("golden campaign inputs changed")
        return state
    if workspace.exists():
        raise BundleError(f"campaign workspace already exists: {workspace}")
    if output.expanduser().resolve().exists():
        raise BundleError(f"golden output already exists: {output}")

    workspace.mkdir(parents=True)
    stages = [
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
    ]
    review_stages = [
        stage["id"]
        for stage in stages
        if stage["status"] != "skipped" and stage["acceptance_required"]
    ]
    state = {
        "schema_version": 1,
        "workspace": str(workspace),
        "status": "ready",
        "inputs": inputs,
        "build": {"status": "pending"},
        "active_bundle": inputs["source_bundle"],
        "active_bundle_class": inputs["source_bundle_class"],
        "active_settings": None,
        "warnings": _selection_warnings(declaration, selected, partial),
        "acceptance": {
            "status": "pending" if review_stages else "not_required",
            "required_stages": review_stages,
        },
        "publication": {"status": "pending"},
        "stages": stages,
    }
    _save(state)
    return state


def _only_campaign_catalog_changed(
    recorded: dict[str, Any] | None, current: dict[str, Any]
) -> bool:
    if recorded is None:
        return False
    old_catalog = recorded.get("campaign_catalog", {})
    new_catalog = current.get("campaign_catalog", {})
    if old_catalog.get("path") != new_catalog.get("path"):
        return False
    old_inputs = {**recorded, "campaign_catalog": new_catalog}
    return old_inputs == current


def _declaration_matches_state(
    state: dict[str, Any], declaration: dict[str, Any]
) -> bool:
    recorded_stages = [
        {
            name: stage[name]
            for name in DECLARATION_STAGE_FIELDS
            if name in stage
        }
        for stage in state.get("stages", [])
    ]
    current_stages = [
        {
            name: stage[name]
            for name in DECLARATION_STAGE_FIELDS
            if name in stage
        }
        for stage in declaration["stages"]
    ]
    return recorded_stages == current_stages


def _refresh_campaign_catalog(
    state: dict[str, Any], campaign_catalog: dict[str, Any]
) -> None:
    """Record an unrelated shared-catalog update without rewinding stages."""
    previous_catalog = state["inputs"]["campaign_catalog"]
    state.setdefault("campaign_catalog_refreshes", []).append(
        {
            "refreshed_utc": utc_now(),
            "case_id": state["inputs"]["case_id"],
            "reason": "selected case declaration unchanged",
            "previous_campaign_catalog": previous_catalog,
            "campaign_catalog": campaign_catalog,
        }
    )
    state["inputs"]["campaign_catalog"] = campaign_catalog
    _save(state)


def _advance(state: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    if state["build"]["status"] != "completed":
        _run_build(state, args)

    for stage in state["stages"]:
        if stage["status"] in {"completed", "skipped"}:
            continue
        if stage["status"] == "failed":
            raise BundleError(f"golden stage failed: {stage['id']}")
        _run_stage(state, stage, args)

    _record_verification_candidate(state)
    if state["acceptance"]["status"] == "pending":
        state["status"] = "awaiting_acceptance"
        _save(state)
        return state
    return _publish(state, args)


def _run_build(state: dict[str, Any], args: argparse.Namespace) -> None:
    state["build"]["status"] = "running"
    state["status"] = "running"
    _save(state)
    result = build_solver(args.settings, args.repository_root, args.build_jobs)
    state["build"] = {
        "status": "completed",
        "directory": str(result.path),
        "settings": _file_record(result.settings_path),
        "metadata": _file_record(result.metadata_path),
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
        _stage_run_id(state, stage),
        state["active_bundle_class"],
        compare=stage["kind"] not in {"matrix", "reference"},
        resume=resume,
    )
    passed = summary["status"] == "passed"
    stage["summary"] = _file_record(summary_path)
    if passed and stage.get("checker") == "balance_diagnostics":
        report = check_balance_diagnostics(summary_path)
        report_path = (
            Path(state["workspace"])
            / "reports"
            / f"{_stage_attempt_name(stage)}-balance-diagnostics.json"
        )
        write_json_atomic(
            report_path,
            report,
            "balance diagnostics report",
        )
        stage["balance_diagnostics_report"] = _file_record(report_path)
        passed = report["status"] == "passed"
    if stage["kind"] == "matrix" and passed:
        passed = _record_matrix_reports(state, stage, summary_path, summary, args)
    elif (
        stage["kind"] == "reference"
        and stage["acceptance_required"]
        and passed
    ):
        old_golden_path, old_golden = verify_suite(
            summary_path,
            args.cases,
            args.tolerances,
        )
        stage["old_golden_report"] = _file_record(old_golden_path)
        passed = _reference_runs_converged(old_golden)
    if not passed:
        stage["status"] = "failed"
        stage["failed_utc"] = utc_now()
        state["status"] = "failed"
        _save(state)
        raise BundleError(f"golden stage failed: {stage['id']}")

    if stage["kind"] in {"matrix", "reference"}:
        _compose_stage(state, stage, summary_path, settings_path, args)

    stage["status"] = "completed"
    if stage.get("candidate"):
        _activate_candidate(state, stage)
    state["status"] = stage["status"]
    _save(state)


def _reference_runs_converged(verification: dict[str, Any]) -> bool:
    """Accept changed reference fields only when every producer converged."""
    results = verification.get("results")
    return bool(results) and all(
        result.get("convergence_status") == "passed" for result in results
    )


def _stage_run_id(state: dict[str, Any], stage: dict[str, Any]) -> str:
    run_id = f"{state['inputs']['run_id']}-{_stage_attempt_name(stage)}"

    return run_id


def _stage_attempt_name(stage: dict[str, Any]) -> str:
    name = stage["id"]
    retry_count = stage.get("retry_count", 0)
    if retry_count:
        name += f"-retry-{retry_count}"
    return name


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
    report_path = (
        Path(state["workspace"])
        / "reports"
        / f"{_stage_attempt_name(stage)}.json"
    )
    write_json_atomic(report_path, report, "matrix layout-pair report")
    stage["layout_pair_report"] = _file_record(report_path)
    if report["status"] != "passed":
        return False
    old_golden_path, _ = verify_suite(
        summary_path,
        args.cases,
        args.tolerances,
        include_layout_pairs=False,
    )
    stage["old_golden_report"] = _file_record(old_golden_path)
    return True


def _accept_campaign(state: dict[str, Any], acceptance: str) -> None:
    if acceptance != "campaign":
        raise BundleError("golden acceptance target must be 'campaign'")
    if (
        state["status"] != "awaiting_acceptance"
        or state["acceptance"]["status"] != "pending"
    ):
        raise BundleError("golden campaign is not awaiting acceptance")
    _record_verification_candidate(state)
    accepted_utc = utc_now()
    state["acceptance"].update(
        {"status": "accepted", "accepted_utc": accepted_utc}
    )
    required = set(state["acceptance"]["required_stages"])
    for stage in state["stages"]:
        if stage["id"] in required:
            stage["accepted_utc"] = accepted_utc
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
    attempt_name = _stage_attempt_name(stage)
    candidate = Path(state["workspace"]) / "candidates" / attempt_name
    stage["candidate"] = str(candidate)
    _save(state)
    if candidate.exists():
        validate_bundle_root(candidate, args.cases)
    elif stage["kind"] == "matrix":
        promote_bundle(
            settings_path,
            summary_path,
            candidate,
            f"{state['inputs']['bundle_version']}-{attempt_name}",
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
            f"{state['inputs']['bundle_version']}-{attempt_name}",
            args.cases,
        )

    generated = Path(state["workspace"]) / "settings" / f"{attempt_name}.env"
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


def _publish(state: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    _record_verification_candidate(state)
    output = Path(state["inputs"]["output"])
    publication = state["publication"]
    if publication["status"] == "completed":
        if _file_record(output / "manifest.json") != publication["manifest"]:
            raise BundleError("published golden manifest changed")
        _validate_published(state, output, args.cases)
        state["status"] = "published"
        return state

    if output.exists():
        if publication["status"] != "publishing":
            raise BundleError(f"golden output already exists: {output}")
        _validate_published(state, output, args.cases)
    else:
        publication.update({"status": "publishing", "started_utc": utc_now()})
        state["status"] = "publishing"
        _save(state)
        publish_campaign_bundle(
            Path(state["active_bundle"]),
            output,
            state["inputs"]["bundle_version"],
            args.cases,
            _campaign_provenance_files(state),
        )

    publication.update(
        {
            "status": "completed",
            "finished_utc": utc_now(),
            "manifest": _file_record(output / "manifest.json"),
        }
    )
    state["status"] = "published"
    _save(state)
    return state


def _validate_published(
    state: dict[str, Any],
    output: Path,
    case_directory: Path,
) -> None:
    validate_bundle_root(output, case_directory)
    manifest = load_json(output / "manifest.json", "published bundle manifest")
    candidate = load_json(
        Path(state["active_bundle"]) / "manifest.json",
        "candidate bundle manifest",
    )
    expected = {
        "bundle_id": candidate["bundle_id"],
        "bundle_version": state["inputs"]["bundle_version"],
        "bundle_class": "golden",
        "case_id": state["inputs"]["case_id"],
    }
    if any(manifest.get(name) != value for name, value in expected.items()):
        raise BundleError("published golden bundle does not match campaign")
    records = [
        artifact
        for artifact in manifest["artifacts"].values()
        if artifact["path"] == "provenance/golden_campaign/campaign.json"
    ]
    if len(records) != 1:
        raise BundleError("published golden bundle has no campaign state")
    recorded = load_json(output / records[0]["path"], "published campaign state")
    if (
        recorded.get("workspace") != state["workspace"]
        or recorded.get("verification_candidate")
        != state["verification_candidate"]
    ):
        raise BundleError("published golden bundle belongs to another campaign")


def _campaign_provenance_files(state: dict[str, Any]) -> list[tuple[str, Path]]:
    workspace = Path(state["workspace"])
    build_record = state["build"]["metadata"]
    build_metadata = Path(build_record["path"])
    if _file_record(build_metadata) != build_record:
        raise BundleError("campaign build metadata changed")
    files = [
        ("campaign.json", workspace / STATE_FILE),
        ("declaration.json", Path(state["inputs"]["campaign_catalog"]["path"])),
        ("build/build_metadata.json", build_metadata),
    ]
    for stage in state["stages"]:
        prefix = f"stages/{stage['id']}"
        if stage.get("summary") is not None:
            _append_attempt_provenance(
                files, prefix, stage, build_record, stage["id"]
            )
        for attempt in stage.get("failed_attempts", []):
            attempt_number = attempt["attempt"] + 1
            _append_attempt_provenance(
                files,
                f"{prefix}/failed_attempts/{attempt_number:03d}",
                attempt,
                build_record,
                f"{stage['id']} failed attempt {attempt_number}",
            )
        for attempt in stage.get("archived_attempts", []):
            attempt_number = attempt["attempt"] + 1
            _append_attempt_provenance(
                files,
                f"{prefix}/archived_attempts/{attempt_number:03d}",
                attempt,
                build_record,
                f"{stage['id']} archived attempt {attempt_number}",
            )
    return files


def _append_attempt_provenance(
    files: list[tuple[str, Path]],
    prefix: str,
    attempt: dict[str, Any],
    build_record: dict[str, Any],
    label: str,
) -> None:
    summary_path = _recorded_path(attempt["summary"], "suite summary")
    documents = [("suite_summary.json", summary_path)]
    documents.extend(
        (f"{name}.json", _recorded_path(attempt[name], name))
        for name in (
            "layout_pair_report",
            "old_golden_report",
            "balance_diagnostics_report",
        )
        if name in attempt
    )
    report_paths = set()
    for name, path in documents:
        files.append((f"{prefix}/{name}", path))
        report_paths.update(_comparison_report_paths(load_json(path, name)))

    summary = load_json(summary_path, "suite summary")
    if summary.get("execution_inputs", {}).get("build_manifest") != build_record:
        raise BundleError(f"golden stage used another build: {label}")
    for number, result in enumerate(summary.get("results", []), start=1):
        run = Path(result["run_directory"])
        run_prefix = f"{prefix}/runs/{number:03d}"
        for name in ("run_plan.json", "run_metadata.json"):
            files.append((f"{run_prefix}/{name}", run / name))
            for nested in sorted(run.glob(f"stages/*/{name}")):
                files.append(
                    (
                        f"{run_prefix}/stages/{nested.parent.name}/{name}",
                        nested,
                    )
                )
    for number, path in enumerate(sorted(report_paths), start=1):
        files.append((f"{prefix}/comparisons/{number:03d}.json", path))


def _comparison_report_paths(document: Any) -> set[Path]:
    if isinstance(document, dict):
        paths = {
            Path(value)
            for name, value in document.items()
            if name == "comparison_report" and isinstance(value, str)
        }
        for value in document.values():
            paths.update(_comparison_report_paths(value))
        return paths
    if isinstance(document, list):
        paths = set()
        for value in document:
            paths.update(_comparison_report_paths(value))
        return paths
    return set()


def _recorded_path(record: dict[str, Any], label: str) -> Path:
    path = Path(record["path"])
    if _file_record(path) != record:
        raise BundleError(f"campaign {label} changed")
    return path


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
    acceptance = state.get("acceptance")
    if acceptance is not None:
        print(f"campaign acceptance: {acceptance['status']}")
    for warning in state.get("warnings", []):
        print(f"warning: {warning}")
    if state.get("publication", {}).get("status") == "completed":
        print(f"output: {state['inputs']['output']}")


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
    recovery = update.add_mutually_exclusive_group()
    recovery.add_argument("--accept", choices=("campaign",))
    recovery.add_argument(
        "--retry-failed",
        action="store_true",
        help="archive and rerun the campaign's failed stage",
    )
    recovery.add_argument(
        "--retry-from",
        metavar="STAGE",
        help="archive and rerun a failed campaign from a corrected stage",
    )
    update.add_argument(
        "--bootstrap-candidate",
        action="store_true",
        help="allow a candidate source for the first promotion of a new case",
    )
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
