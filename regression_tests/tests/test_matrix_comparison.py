from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from bundle.cases import load_case_definition  # noqa: E402
from comparison.workflow import compare_completed_run  # noqa: E402


class ReferenceMatrixComparisonTests(unittest.TestCase):
    @patch("comparison.matrix.compare_fixed_run")
    def test_fixed_comparison_stops_at_first_divergent_stage(self, compare) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run, references = self._fixture(root, "cold_fixed")

            def result(*args, **kwargs):
                stage_id = args[0].run_directory.name.split("_", 1)[1]
                failures = (
                    ["solution/u exceeds tolerance"]
                    if stage_id == "diffusion_reduction"
                    else []
                )
                return kwargs["report_path"], {
                    "status": "failed" if failures else "passed",
                    "failures": failures,
                }

            compare.side_effect = result
            policy, path, report = compare_completed_run(
                run,
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "tolerances.json",
            )

        self.assertEqual(policy, "reference_matrix")
        self.assertEqual(path.name, "matrix_comparison.json")
        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["checked_stage_count"], 2)
        self.assertEqual(report["first_failed_stage"], "diffusion_reduction")
        self.assertEqual(compare.call_count, 2)
        self.assertEqual(
            compare.call_args.kwargs["overrides"].reference,
            references["diffusion_reduction"],
        )
        self.assertEqual(
            compare.call_args.kwargs["overrides"].tolerance_profile,
            "fixed_stage_reference",
        )

    @patch("comparison.matrix.compare_adaptive_run")
    def test_adaptive_comparison_checks_every_matching_stage(self, compare) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run, _ = self._fixture(root, "cold_adaptive")
            compare.side_effect = lambda *args, **kwargs: (
                kwargs["report_path"],
                {"status": "passed", "failures": []},
            )

            policy, _, report = compare_completed_run(
                run,
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "tolerances.json",
            )

        self.assertEqual(policy, "reference_matrix")
        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["checked_stage_count"], 7)
        self.assertIsNone(report["first_failed_stage"])
        self.assertEqual(compare.call_count, 7)
        self.assertEqual(
            compare.call_args.kwargs["overrides"].tolerance_profile,
            "adaptive_reference",
        )

    @staticmethod
    def _fixture(root: Path, workflow_id: str) -> tuple[Path, dict[str, Path]]:
        case = load_case_definition("legacy_case", REGRESSION_ROOT / "cases")
        stage_ids = [
            stage["stage_id"] for stage in case["workflows"][workflow_id]["stages"]
        ]
        bundle = root / "golden"
        references_dir = bundle / "references"
        references_dir.mkdir(parents=True)
        references = {}
        artifacts = {}
        entries = []
        for stage_id in stage_ids:
            reference = references_dir / f"{stage_id}.h5"
            reference.write_text(f"golden {stage_id}\n", encoding="utf-8")
            artifact_id = f"golden_{stage_id}"
            references[stage_id] = reference.resolve()
            artifacts[artifact_id] = ReferenceMatrixComparisonTests._artifact(
                reference.relative_to(bundle), "application/x-hdf5"
            )
            entries.append(
                {
                    "workflow_id": workflow_id,
                    "layout_id": "mpi4_omp4",
                    "stage_id": stage_id,
                    "artifact_id": artifact_id,
                }
            )

        index = references_dir / "index.json"
        index.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "created_utc": "2026-07-20T00:00:00Z",
                    "case_id": "legacy_case",
                    "suite_id": "cold_matrix",
                    "suite_run_id": "golden-1",
                    "source_bundle": {
                        "bundle_id": "legacy_bundle",
                        "bundle_version": "candidate-1",
                    },
                    "tracked_reference": {
                        "branch": "develop",
                        "revision": "0" * 40,
                    },
                    "references": entries,
                }
            ),
            encoding="utf-8",
        )
        artifacts["golden_index"] = ReferenceMatrixComparisonTests._artifact(
            index.relative_to(bundle), "application/json"
        )
        (bundle / "manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "bundle_id": "golden_bundle",
                    "bundle_version": "golden-1",
                    "bundle_class": "golden",
                    "created_utc": "2026-07-20T00:00:00Z",
                    "case_id": "legacy_case",
                    "roles": {"reference_matrix": "golden_index"},
                    "artifacts": artifacts,
                }
            ),
            encoding="utf-8",
        )

        run = root / "run"
        stage_records = []
        plan_stages = []
        for number, stage_id in enumerate(stage_ids, start=1):
            stage_dir = run / "stages" / f"{number:02d}_{stage_id}"
            stage_dir.mkdir(parents=True)
            candidate = stage_dir / "candidate.h5"
            candidate.write_text(f"candidate {stage_id}\n", encoding="utf-8")
            (stage_dir / "run_metadata.json").write_text(
                json.dumps({"status": "completed"}),
                encoding="utf-8",
            )
            stage_records.append(
                {
                    "stage_id": stage_id,
                    "run_directory": str(stage_dir),
                    "selected_hdf5": str(candidate),
                    "status": "completed",
                }
            )
            plan_stages.append({"stage_id": stage_id})
        (run / "run_plan.json").write_text(
            json.dumps(
                {
                    "case_id": "legacy_case",
                    "workflow_id": workflow_id,
                    "layout_id": "mpi4_omp4",
                    "bundle": {
                        "root": str(bundle),
                        "bundle_id": "golden_bundle",
                        "bundle_version": "golden-1",
                    },
                    "stages": plan_stages,
                }
            ),
            encoding="utf-8",
        )
        (run / "run_metadata.json").write_text(
            json.dumps({"status": "completed", "stages": stage_records}),
            encoding="utf-8",
        )
        return run, references

    @staticmethod
    def _artifact(path: Path, media_type: str) -> dict[str, object]:
        return {
            "path": path.as_posix(),
            "sha256": "0" * 64,
            "size_bytes": 1,
            "media_type": media_type,
        }
