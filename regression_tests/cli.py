"""Public commands; repository catalogs are discovered beside this module."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import reporting

ROOT = Path(__file__).resolve().parent


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="python -m regression_tests", description="MHDG scientific regression checks.",
        allow_abbrev=False,
    )
    commands = result.add_subparsers(dest="command", required=True)

    def command(name, help, *, parent=commands, settings=False):
        child = parent.add_parser(name, help=help, description=help, allow_abbrev=False)
        if settings:
            child.add_argument(
                "--settings", type=Path,
                help="machine JSON (default: settings.local.json or MHDG_REGRESSION_SETTINGS)",
            )
        return child

    def selection(child):
        child.add_argument("--bundle", type=Path, help="external bundle; overrides the selected case default")
        child.add_argument("--build-manifest", type=Path, help="build_metadata.json; overrides the default build")

    check = command("check", "Run a suite and its scientific comparisons.", settings=True)
    selection(check)
    check.add_argument("suite", nargs="?", default="warm", metavar="SUITE")
    check.add_argument("--run-id")
    check.add_argument("--resume", action="store_true", help="resume an existing suite run")
    check.add_argument("--run-only", action="store_true", help="record runs without comparisons")
    check.add_argument("--allow-candidate", action="store_true", help="also accept candidate bundles")
    check.add_argument("--build", action="store_true", help="build executables before checking")
    check.add_argument("--build-jobs", type=_positive_integer, metavar="N")
    check.add_argument("--diagnostics", choices=("off", "summary", "equations", "detailed"))

    for name, help in (
        ("run", "Prepare and execute one workflow; use check for comparisons."),
        ("prepare", "Debug: prepare one workflow without execution."),
    ):
        run = command(name, help, settings=True)
        selection(run)
        run.add_argument("case", metavar="CASE")
        run.add_argument("workflow", metavar="WORKFLOW")
        run.add_argument("--layout", default="mpi4_omp4", metavar="LAYOUT")
        run.add_argument("--run-id")

    listing = command("list", "Show available scientific selections without machine settings.")
    lists = listing.add_subparsers(dest="listing", required=True)
    for name in ("cases", "suites", "layouts", "workflows"):
        item = command(name, f"List {name}.", parent=lists)
        if name == "workflows":
            item.add_argument("case", nargs="?", metavar="CASE")

    build = command("build", "Build regression executables.", settings=True)
    build.add_argument("--jobs", type=_positive_integer, metavar="N")

    doctor = command("doctor", "Check selected catalogs, data, build and machine prerequisites.", settings=True)
    selection(doctor)
    doctor.add_argument("suite", nargs="?", default="warm", metavar="SUITE")

    compare = command("compare", "Debug: compare a saved run or suite summary.")
    compare.add_argument("path", type=Path, metavar="RUN_OR_SUITE_SUMMARY")
    compare.add_argument("--suite", action="store_true", help="compare a saved suite summary")
    compare.add_argument("--candidate", type=Path)
    compare.add_argument("--reference", type=Path)
    compare.add_argument("--tolerance-profile")
    compare.add_argument("--report", type=Path)

    bundle = command("bundle", "Create or validate external regression data.")
    bundles = bundle.add_subparsers(dest="action", required=True)
    create = command("create", "Copy prepared inputs into a new candidate bundle.", parent=bundles)
    create.add_argument("--case", required=True, metavar="CASE")
    create.add_argument("--source", required=True, type=Path)
    create.add_argument("--output", required=True, type=Path)
    create.add_argument("--bundle-version", default="1.0.0", metavar="VERSION")
    validate = command("validate", "Validate a bundle, including artifact checksums.", parent=bundles)
    validate.add_argument("path", type=Path, metavar="BUNDLE")
    return result


def _positive_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        number = 0
    if number < 1:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return number


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv[1:] if argv is None else argv
    command_parser = parser()
    if not arguments:
        command_parser.print_help()
        return 0
    args = command_parser.parse_args(arguments)
    # Temporary bridge to existing internals, removed as their owning modules
    # migrate in steps 3-6. Help remains usable without scientific dependencies.
    tools_path = str(ROOT / "tools")
    if tools_path not in sys.path:
        sys.path.insert(0, tools_path)
    from support.errors import HarnessError

    try:
        return _dispatch(args)
    except (HarnessError, OSError) as exc:
        return reporting.error(str(exc))
    except KeyboardInterrupt:
        reporting.error("interrupted")
        return 130


def _dispatch(args: argparse.Namespace) -> int:
    if args.command == "list":
        return _list(args)
    if args.command == "check":
        return _check(args)
    if args.command == "doctor":
        from .doctor import diagnose

        return diagnose(args)
    if args.command in {"run", "prepare"}:
        from catalogs.layouts import load_layout
        from .prepare import prepare_run
        from .execute import execute_prepared

        settings = _settings(args, args.case)
        from .config import runtime_settings

        runtime_settings(settings, [load_layout(args.layout, ROOT / "layouts.json")])
        prepared = prepare_run(
            settings, args.case, args.workflow, args.layout,
            ROOT / "cases", ROOT / "layouts.json", args.run_id,
        )
        reporting.prepared(prepared)
        if args.command == "prepare":
            return 0
        result = execute_prepared(prepared, settings)
        return reporting.run_result(result)
    if args.command == "build":
        from .build import build_solver
        from .config import settings

        result = build_solver(settings(args.settings, use_build=False), ROOT.parent, args.jobs)
        reporting.status("build", "completed", result.path)
        print(f"select with --build-manifest {result.metadata_path}")
        return 0
    if args.command == "compare":
        return _compare(args)
    if args.command == "bundle":
        from bundle.creation import create_bundle
        from bundle.validation import validate_bundle_root

        if args.action == "create":
            summary = create_bundle(
                args.case, args.source, args.output, ROOT / "cases", args.bundle_version,
            )
            reporting.status("bundle", "created", args.output)
        else:
            summary = validate_bundle_root(args.path, ROOT / "cases")
            reporting.status("bundle", "valid", args.path)
        reporting.bundle(summary)
        return 0
    raise AssertionError(f"unhandled command: {args.command}")


def _check(args: argparse.Namespace) -> int:
    from .build import build_solver
    from suite.runner import run_suite
    from suite.reporting import print_run_summary
    from support.errors import BundleError
    from suite.configuration import load_suite_definition
    from catalogs.layouts import load_layouts
    from .config import build_settings, runtime_settings

    if args.build and args.resume:
        raise BundleError("--build cannot be used with --resume; reuse the original build settings")
    if args.build_jobs is not None and not args.build:
        raise BundleError("--build-jobs requires --build")
    if args.build and args.build_manifest:
        raise BundleError("select --build or --build-manifest, not both")
    suite = load_suite_definition(
        args.suite, ROOT / "suites.json", ROOT / "layouts.json", ROOT / "cases",
    )
    settings = _settings(args, suite["case_id"], use_build=not args.build)
    layouts = load_layouts(ROOT / "layouts.json")
    if args.build:
        variants = {"serial" if layouts[name]["execution"] == "serial" else "parallel"
                    for name in suite["layouts"]}
        build = build_solver(settings, ROOT.parent, args.build_jobs, variants=variants)
        settings.update(build_settings(build.metadata_path))
        reporting.status("build", "completed", build.path)
        print(f"select with --build-manifest {build.metadata_path}")
    runtime_settings(settings, [layouts[name] for name in suite["layouts"]])
    path, summary = run_suite(
        settings, args.suite, ROOT / "cases", ROOT / "layouts.json",
        ROOT / "suites.json", ROOT / "tolerances.json", args.run_id,
        required_bundle_class=None if args.allow_candidate else "golden",
        compare=not args.run_only, resume=args.resume,
        parameter_overrides=(
            {"balance_diagnostics_mode": args.diagnostics} if args.diagnostics else None
        ),
    )
    print_run_summary(summary, path)
    return 0 if summary["status"] == "passed" else 1


def _settings(args, case, *, use_build=True):
    from .config import settings
    from support.errors import BundleError

    values = settings(
        args.settings, case=case, bundle=args.bundle,
        build_manifest=args.build_manifest, use_build=use_build,
    )
    if not values.get("MHDG_REGRESSION_DATA_ROOT"):
        raise BundleError(f"no bundle selected for {case}; pass --bundle DIR or set defaults.bundles")
    print(f"case: {case}")
    print(f"bundle: {values['MHDG_REGRESSION_DATA_ROOT']}")
    selected_build = values.get("MHDG_BUILD_MANIFEST", "prebuilt executable settings")
    print(f"build: {selected_build if use_build else 'new build requested'}")
    if use_build and values.get("MHDG_SOLVER_REVISION"):
        print(f"solver revision: {values['MHDG_SOLVER_REVISION']}")
    return values


def _compare(args: argparse.Namespace) -> int:
    from support.errors import BundleError

    if args.suite:
        from suite.verification import verify_suite
        from suite.reporting import print_verification_summary

        if any((args.candidate, args.reference, args.tolerance_profile, args.report)):
            raise BundleError("--suite cannot be combined with individual comparison overrides")
        path, summary = verify_suite(args.path, ROOT / "cases", ROOT / "tolerances.json")
        print_verification_summary(summary, path)
        return 0 if summary["status"] == "passed" else 1
    from comparison.workflow import compare_completed_run

    policy, path, report = compare_completed_run(
        args.path, ROOT / "cases", ROOT / "tolerances.json",
        candidate_override=args.candidate, reference_override=args.reference,
        tolerance_profile_override=args.tolerance_profile, report_override=args.report,
    )
    reporting.comparison(policy, path, report)
    return 0 if report["status"] == "passed" else 1


def _list(args: argparse.Namespace) -> int:
    from bundle.cases import load_case_definition
    from bundle.schemas import load_validated_json
    from catalogs.layouts import load_layouts

    if args.listing == "layouts":
        entries = load_layouts(ROOT / "layouts.json")
    elif args.listing == "suites":
        entries = load_validated_json(
            ROOT / "suites.json", ROOT / "schemas/suites.schema.json", "suites",
        )["suites"]
    else:
        names = [args.case] if getattr(args, "case", None) else [
            path.stem for path in sorted((ROOT / "cases").glob("*.json"))
        ]
        entries = {}
        for name in names:
            case = load_case_definition(name, ROOT / "cases")
            if args.listing == "cases":
                entries[name] = case
            else:
                for workflow, definition in case["workflows"].items():
                    label = workflow if args.case else f"{name}/{workflow}"
                    entries[label] = definition
    for name, definition in entries.items():
        print(f"{name}: {definition['description']}")
    return 0
