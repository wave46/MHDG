"""Read-only checks of the selected regression setup; never build or run MHDG."""

import importlib
import os
from pathlib import Path

from . import config
from .support import HarnessError, BundleError


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

    from .bundles import validate_bundle_root
    from .config import bundle_root_from_settings
    from .catalog import load_case_definition, load_selection
    from .catalog import load_layouts
    try:
        _, checks = load_selection(
            args.suite, config.ROOT / "suites.json", config.ROOT / "layouts.json",
            config.ROOT / "cases", case_id=args.case,
        )
        layouts = load_layouts(config.ROOT / "layouts.json")
        cases = list(dict.fromkeys(item["case_id"] for item in checks))
        if args.bundle and len(cases) != 1:
            raise BundleError("--bundle requires a single-case selection; configure defaults.bundles")
    except (HarnessError, OSError) as exc:
        print(f"FAIL configuration: {exc}")
        return 1
    print(f"PASS catalogs: {args.suite} / {', '.join(cases)}")
    needs_interpolation = False
    for case_id in cases:
        case = load_case_definition(case_id, config.ROOT / "cases")
        selected = [item for item in checks if item["case_id"] == case_id]
        needs_interpolation |= any(case["workflows"][name].get("comparison_policy") == "mesh_independent"
                                   for item in selected for name in item["workflow_ids"])
        try:
            values = config.settings(args.settings, case=case_id, bundle=args.bundle,
                                     build_manifest=args.build_manifest)
        except (HarnessError, OSError) as exc:
            failures.append(case_id)
            print(f"FAIL {case_id} settings: {exc}")
            continue

        def bundle():
            root = bundle_root_from_settings(values)
            workflows = {name for item in selected for name in item["workflow_ids"]}
            summary = validate_bundle_root(root, config.ROOT / "cases", workflows=workflows)
            if summary.case_id != case_id:
                raise BundleError(f"selected check needs {case_id}, bundle contains {summary.case_id}")
            return f"{root} ({summary.verified_artifact_count} artifacts verified)"

        def runtime():
            names = {name for item in selected for name in item["layouts"]}
            config.runtime_settings(values, [layouts[name] for name in names])
            return values["MHDG_BUILD_MANIFEST"]

        check(f"{case_id} bundle", bundle)
        check(f"{case_id} build/runtime", runtime)
        check("run root", lambda: writable_root(Path(values["MHDG_REGRESSION_RUN_ROOT"])))
    if needs_interpolation:
        check("Python hdg_postprocess", lambda: importlib.import_module("hdg_postprocess").__name__)
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
