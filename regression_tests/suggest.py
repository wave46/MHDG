"""Advisory file-to-check mapping; never execute checks or infer solver semantics."""

from fnmatch import fnmatchcase
from pathlib import Path
import subprocess

from .catalog import load_selection
from .support import BundleError

ROOT = Path(__file__).resolve().parent
# First matching rule owns a file. Keep the mapping explicit and small; unknown
# files need human assessment, especially models outside the current case tree.
RULES = (
    (("regression_tests/tests/*", "regression_tests/clean.py", "regression_tests/suggest.py",
      "regression_tests/reporting.py", "regression_tests/doctor.py", "regression_tests/build.py"),
     ("pytest",), "Harness tests, tooling or reporting; verify their Python behavior"),
    (("regression_tests/diagnostics.py",), ("pytest", "routine-extended", "sources", "neutralgamma"),
     "Harness diagnostic checks on real parallel/source outputs, plus Python behavior"),
    (("regression_tests/*",), ("pytest", "full"),
     "Shared harness workflows, comparisons or scientific inputs; verify Python behavior and both real cases"),
    (("lib/Makefile*", "lib/Make.inc/*"), ("routine-extended",), "serial/MPI builds and parallel execution"),
    (("src/Utils/Diagnostics/*",), ("routine-extended", "sources", "neutralgamma", "diagnostic-tests"),
     "diagnostic output and totals on existing parallel/source checks; output-format checker tests"),
    (("src/Models/neutral_flux_limiter.f90",), ("neutrals", "neutral_parallel"),
     "neutral limiter/projection goldens and parallel consistency"),
    (("src/Models/NGammaTiTe/impurity_radiation/*",), ("warm", "impurities"),
     "tungsten warm reference plus disabled/nitrogen/mixed-impurity references"),
    (("src/Models/NGammaTiTe/transport_1d/*",), ("routine-extended",),
     "transport_1d warm reconvergence and short fixed/adaptive parallel checks"),
    (("src/HDG/initialization.f90", "src/HDG/preprocess.f90", "src/InOut/read_input.f90",
      "src/MHDG.f90", "src/Convergence.f90", "src/Definitions/*", "src/Models/adimensionalization.f90",
      "src/Models/magnetic_*.f90", "src/Models/NGammaTiTe/physics.f90", "src/Models/NGammaTiTe/analytical.f90"),
     ("full",), "shared physics, initialization or convergence: both topologies, cold chains and feature coverage"),
    (("src/Adaptivity/*",), ("routine-extended",),
     "short adaptive execution, field and generated-mesh agreement across layouts"),
    (("src/MPI_OMP/*", "src/HDG/*", "src/LinearAlgebra/*", "src/InOut/*"), ("routine-extended",),
     "shared assembly, ownership or I/O: warm and fixed/adaptive parallel consistency"),
)


def _git(repository, *arguments):
    result = subprocess.run(["git", *arguments], cwd=repository, capture_output=True, text=True)
    if result.returncode:
        raise BundleError(f"cannot inspect Git changes: {result.stderr.strip()}")
    return result.stdout


def changed_paths(repository, base=None):
    """Include both sides of renames and each local layer, even if edits cancel."""
    paths = set()
    if base is not None:
        revision = _git(repository, "rev-parse", "--verify", "--end-of-options", f"{base}^{{commit}}").strip()
        paths.update(_git(repository, "diff", "--name-only", "-z", "--no-renames", f"{revision}...HEAD", "--").split("\0"))
    for arguments in (("diff", "--name-only", "-z", "--no-renames", "--cached", "--"),
                      ("diff", "--name-only", "-z", "--no-renames", "--"),
                      ("ls-files", "--others", "--exclude-standard", "-z")):
        paths.update(_git(repository, *arguments).split("\0"))
    return sorted(paths - {""})


def recommendations(paths, *, catalog_root=ROOT):
    matched, unmapped, documentation = {}, [], []
    paths = sorted(set(paths))
    needs_build = any(path.startswith(("src/", "lib/")) and Path(path).suffix.lower() not in {".md", ".rst"}
                      for path in paths)
    for path in paths:
        if Path(path).suffix.lower() in {".md", ".rst"}:
            documentation.append(path)
            continue
        for patterns, checks, reason in RULES:
            if any(fnmatchcase(path, pattern) for pattern in patterns):
                for check in checks:
                    matched.setdefault(check, {})[path] = reason
                break
        else:
            unmapped.append(path)
    # Catalog expansion owns profile membership. Suppress a smaller suggestion
    # only when the broader one includes all of its actual case/suite selections.
    coverage, catalog = {}, {}
    for name in matched:
        if name in {"pytest", "diagnostic-tests"}:
            continue
        _, selections = load_selection(name, catalog_root / "suites.json", catalog_root / "layouts.json",
                                       catalog_root / "cases", catalog=catalog)
        coverage[name] = {(item["case_id"], item["suite_id"]) for item in selections}
    for name in sorted(coverage, key=lambda item: len(coverage[item])):
        covering = next((other for other in matched if other in coverage and coverage[name] < coverage[other]), None)
        if covering:
            matched[covering].update(matched.pop(name))
    if "pytest" in matched and "diagnostic-tests" in matched:
        matched["pytest"].update(matched.pop("diagnostic-tests"))
    checks = []
    for name, reasons in sorted(matched.items(), key=lambda item: (item[0] not in {"pytest", "diagnostic-tests"}, item[0])):
        if name in {"pytest", "diagnostic-tests"}:
            target = "regression_tests/tests" if name == "pytest" else "regression_tests/tests/test_balance_diagnostics.py"
            command = ["python", "-m", "pytest", "-q", target]
        else:
            command = ["python", "-m", "regression_tests", "check", name, *(["--build"] if needs_build else [])]
        checks.append({"command": command, "reasons": reasons})
    return {"checks": checks, "unmapped": unmapped, "documentation": documentation}


def suggest(repository, base=None):
    paths = changed_paths(repository, base)
    return {"base": base, "paths": paths, **recommendations(paths)}
