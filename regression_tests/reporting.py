"""Shared terminal conventions; scientific details retain their existing reports."""

import shlex
import sys


def error(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 1


def status(subject, outcome, path) -> None:
    print(f"{subject} {outcome}: {path}", flush=True)


def prepared(run) -> None:
    from preparation.models import PreparedRun

    status("run", "prepared", run.path)
    if isinstance(run, PreparedRun):
        print(f"command: {shlex.join(run.command)}", flush=True)
    else:
        print(f"stages: {len(run.stages)}", flush=True)


def run_result(result) -> int:
    print(f"status: {result.status}")
    print(f"solver exit code: {result.exit_code}")
    print(f"runtime: {result.duration_seconds:.3f} s")
    print(f"HDF5 outputs: {len(result.hdf5_outputs)}")
    return 0 if result.status == "completed" else 1


def bundle(summary) -> None:
    print(f"case: {summary.case_id}; version: {summary.bundle_version}")
    print(
        f"verified {summary.verified_artifact_count}/{summary.artifact_count} artifacts "
        f"({summary.verified_bytes} bytes)"
    )
    for warning in summary.warnings:
        print(f"warning: {warning}", file=sys.stderr)


def comparison(policy, path, report) -> None:
    from comparison.shared.reporting import print_adaptive_summary, print_fixed_summary

    print(f"comparison policy: {policy}")
    status("comparison", report["status"], path)
    if policy in {"fixed_hdf5", "mesh_independent"}:
        files = report if policy == "fixed_hdf5" else report["files"]
        print(f"candidate: {files['candidate']}")
        print(f"reference: {files['reference']}")
        print(f"tolerance profile: {report['tolerance_profile']['id']}")
        if policy == "fixed_hdf5":
            print_fixed_summary(report)
        else:
            print_adaptive_summary(report)
    else:
        print(f"stages checked: {report['checked_stage_count']}/{report['total_stage_count']}")
        for failure in report["failures"]:
            print(f"FAIL: {failure}")
