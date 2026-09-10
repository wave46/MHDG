"""Read-only checks of the selected regression setup; never build or run MHDG."""

import importlib
import os
from pathlib import Path

from . import config
from support.errors import HarnessError, BundleError


def diagnose(args):
    failures = []

    def check(label, action):
        try:
            detail = action()
        except (HarnessError, OSError, ImportError) as exc:
            failures.append(label)
            print(f"FAIL {label}: {exc}")
            return None
        print(f"PASS {label}: {detail}")
        return detail

    for dependency in ("numpy", "h5py", "jsonschema"):
        check(f"Python {dependency}", lambda name=dependency: importlib.import_module(name).__name__)
    if "Python jsonschema" in failures:
        print("Install regression_tests/requirements.txt with this Python interpreter.")
        return 1

    from bundle.validation import validate_bundle_root
    from bundle.settings import bundle_root_from_settings
    from bundle.cases import load_case_definition
    from catalogs.layouts import load_layouts
    from suite.configuration import load_suite_definition

    def catalogs():
        for path in (config.ROOT / "cases").glob("*.json"):
            load_case_definition(path.stem, config.ROOT / "cases")
        return load_suite_definition(
            args.suite, config.ROOT / "suites.json", config.ROOT / "layouts.json",
            config.ROOT / "cases",
        )

    # Keep detailed objects out of the terminal while sharing the actual loader.
    try:
        suite = catalogs()
        layouts = load_layouts(config.ROOT / "layouts.json")
        values = config.settings(
            args.settings, case=suite["case_id"], bundle=args.bundle,
            build_manifest=args.build_manifest,
        )
    except (HarnessError, OSError) as exc:
        print(f"FAIL configuration: {exc}")
        return 1
    print(f"PASS catalogs: {args.suite} / {suite['case_id']}")
    case = load_case_definition(suite["case_id"], config.ROOT / "cases")
    if any(case["workflows"][name].get("comparison_policy") == "mesh_independent"
           for name in suite["workflow_ids"]):
        check("Python hdg_postprocess", lambda: importlib.import_module("hdg_postprocess").__name__)

    def bundle():
        if not values.get("MHDG_REGRESSION_DATA_ROOT"):
            raise BundleError("no bundle selected; pass --bundle DIR or set defaults.bundles")
        root = bundle_root_from_settings(values)
        summary = validate_bundle_root(root, config.ROOT / "cases")
        if summary.case_id != suite["case_id"]:
            raise BundleError(f"selected suite needs {suite['case_id']}, bundle contains {summary.case_id}")
        return f"{root} ({summary.verified_artifact_count} artifacts verified)"

    def runtime():
        config.runtime_settings(values, [layouts[name] for name in suite["layouts"]])
        return values.get("MHDG_BUILD_MANIFEST", "legacy executable settings")

    check("bundle", bundle)
    check("build/runtime", runtime)
    check("run root", lambda: writable_root(Path(values["MHDG_REGRESSION_RUN_ROOT"])))
    print(f"doctor: {'failed' if failures else 'passed'}")
    return 1 if failures else 0


def writable_root(path):
    path = path.expanduser()
    if not path.is_absolute():
        raise BundleError(f"scratch root must be absolute: {path}")
    parent = path
    while not parent.exists() and parent != parent.parent:
        parent = parent.parent
    if not parent.is_dir() or not os.access(parent, os.W_OK | os.X_OK):
        raise BundleError(f"scratch root is not writable: {path}")
    return str(path)
