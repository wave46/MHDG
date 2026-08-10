#!/usr/bin/env python3
"""Check relocated neutral wall sources in a completed regression suite."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

from support.documents import load_json, write_json_atomic
from support.errors import HarnessError
from support.paths import require_directory, require_file


CONSERVATION_TOLERANCE = 1.0e-12
LAYOUT_TOLERANCE = 5.0e-8
MARKER = "NEUTRAL_WALL_SOURCE_CONSERVATION"
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
RECORD = re.compile(
    rf"^\s*{MARKER}\s+(?P<source>puff|pump)\s+"
    rf"wall=\s*(?P<wall>{NUMBER})\s+volume=\s*(?P<volume>{NUMBER})\s*$"
)
SOURCES = ("puff", "pump")


def check_run(run_directory: Path) -> dict:
    """Check every wall/volume pair in every completed stage log."""
    run_directory = require_directory(run_directory, "neutral wall-source run")
    logs = _stdout_logs(run_directory)
    failures = []
    counts = dict.fromkeys(SOURCES, 0)
    final = {}
    if not logs:
        failures.append("run metadata contains no completed stage logs")

    for path in logs:
        seen = set()
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            for line_number, line in enumerate(stream, start=1):
                if MARKER not in line:
                    continue
                match = RECORD.match(line)
                if match is None:
                    failures.append(f"{path}:{line_number}: malformed record")
                    continue
                source = match.group("source")
                wall = _number(match.group("wall"))
                volume = _number(match.group("volume"))
                difference = _relative_difference(wall, volume)
                if not (math.isfinite(wall) and math.isfinite(volume)):
                    failures.append(f"{path}:{line_number}: non-finite {source} total")
                elif difference > CONSERVATION_TOLERANCE:
                    failures.append(
                        f"{path}:{line_number}: {source} wall/volume relative "
                        f"difference {difference:.6e} exceeds "
                        f"{CONSERVATION_TOLERANCE:.6e}"
                    )
                seen.add(source)
                counts[source] += 1
                final[source] = {"wall": wall, "volume": volume}
        failures.extend(
            f"{path}: missing {source} conservation records"
            for source in SOURCES
            if source not in seen
        )

    missing = [source for source in SOURCES if source not in final]
    if missing:
        failures.append("missing final source totals: " + ", ".join(missing))
    for source in SOURCES:
        if source not in final:
            continue
        nonpositive = [
            f"{measure}={final[source][measure]:.16e}"
            for measure in ("wall", "volume")
            if math.isfinite(final[source][measure]) and final[source][measure] <= 0.0
        ]
        if nonpositive:
            failures.append(
                f"final {source} totals must be positive: " + ", ".join(nonpositive)
            )
    return {
        "run_directory": str(run_directory),
        "record_counts": counts,
        "final_totals": final,
        "status": "passed" if not failures else "failed",
        "failures": failures,
    }


def check_suite(summary_path: Path) -> dict:
    """Check all suite runs and all declared cross-layout pairs."""
    summary_path = require_file(summary_path, "suite summary")
    summary = load_json(summary_path, "suite summary")
    checked = {}
    run_reports = []
    failures = []

    for result in summary.get("results", []):
        key = (result.get("workflow_id"), result.get("layout_id"))
        if result.get("run_status") == "completed":
            report = check_run(Path(result["run_directory"]))
            checked[key] = report
        else:
            report = {
                "status": "failed",
                "failures": [f"solver run did not complete: {key[0]}/{key[1]}"],
            }
        run_reports.append({"workflow_id": key[0], "layout_id": key[1], **report})
        failures.extend(report["failures"])

    pair_reports = []
    for workflow in summary.get("workflow_ids", []):
        for pair in summary.get("layout_comparisons", []):
            first = pair["baseline"]
            second = pair["candidate"]
            if (workflow, first) in checked and (workflow, second) in checked:
                report = _compare_totals(checked[(workflow, first)], checked[(workflow, second)])
            else:
                report = {
                    "status": "failed",
                    "failures": [f"missing completed layout pair: {workflow}/{first}->{second}"],
                }
            pair_reports.append(
                {
                    "workflow_id": workflow,
                    "baseline_layout_id": first,
                    "candidate_layout_id": second,
                    **report,
                }
            )
            failures.extend(report["failures"])

    return {
        "suite_summary": str(summary_path),
        "conservation_tolerance": CONSERVATION_TOLERANCE,
        "layout_tolerance": LAYOUT_TOLERANCE,
        "runs": run_reports,
        "comparisons": pair_reports,
        "status": "passed" if not failures else "failed",
        "failures": failures,
    }


def _compare_totals(baseline: dict, candidate: dict) -> dict:
    failures = [
        *(f"baseline: {item}" for item in baseline["failures"]),
        *(f"candidate: {item}" for item in candidate["failures"]),
    ]
    differences = {}
    if not failures:
        for source in SOURCES:
            differences[source] = {}
            for measure in ("wall", "volume"):
                difference = _relative_difference(
                    baseline["final_totals"][source][measure],
                    candidate["final_totals"][source][measure],
                )
                differences[source][measure] = difference
                if difference > LAYOUT_TOLERANCE:
                    failures.append(
                        f"final {source} {measure} relative difference "
                        f"{difference:.6e} exceeds {LAYOUT_TOLERANCE:.6e}"
                    )
    return {
        "relative_differences": differences,
        "status": "passed" if not failures else "failed",
        "failures": failures,
    }


def _stdout_logs(run_directory: Path) -> list[Path]:
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    if not metadata.get("stages"):
        return [require_file(run_directory / "stdout.log", "solver stdout")]
    return [
        require_file(Path(stage["run_directory"]) / "stdout.log", "stage stdout")
        for stage in metadata["stages"]
        if stage.get("status") == "completed"
    ]


def _number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _relative_difference(first: float, second: float) -> float:
    scale = max(abs(first), abs(second))
    return abs(first - second) / scale if scale else 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_summary", type=Path)
    args = parser.parse_args()
    try:
        report = check_suite(args.suite_summary)
        path = args.suite_summary.parent / "neutral_wall_source_check.json"
        write_json_atomic(path, report, "neutral wall-source report")
    except HarnessError as exc:
        parser.error(str(exc))
    print(f"neutral wall-source check: {report['status']}")
    print(f"report: {path}")
    for failure in report["failures"]:
        print(f"failure: {failure}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
