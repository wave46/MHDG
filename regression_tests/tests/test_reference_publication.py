from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from bundle.cases import load_case_definition  # noqa: E402
from bundle.creation import create_bundle  # noqa: E402
from bundle.promotion import promote_bundle  # noqa: E402
from bundle.validation import validate_bundle_root  # noqa: E402
from support.errors import BundleError  # noqa: E402
from support.files import file_identity  # noqa: E402
from tests.fixtures.case_data import write_case_source  # noqa: E402
from tests.fixtures.harness import run_command  # noqa: E402


class ReferencePublicationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        source = write_case_source(self.root / "source")
        (source / "reference_mpi4_omp4.h5").write_text(
            "old reference\n",
            encoding="utf-8",
        )
        self.candidate = self.root / "candidate"
        create_bundle(
            "legacy_case",
            source,
            self.candidate,
            REGRESSION_ROOT / "cases",
        )
        self.settings = self.root / "settings.env"
        self.settings.write_text(
            "MHDG_REGRESSION_SETTINGS_VERSION=2\n"
            f"MHDG_REGRESSION_DATA_ROOT={self.candidate}\n",
            encoding="utf-8",
        )
        self.summary = self._write_canonical_summary()
        self.output = self.root / "golden"

    def test_publication_creates_a_valid_golden_without_mutating_source(self) -> None:
        source_manifest = _load_json(self.candidate / "manifest.json")
        source_reference = self.candidate / source_manifest["artifacts"][
            "legacy_case_warm_reference"
        ]["path"]

        completed = run_command(
            "bundle",
            "promote",
            str(self.summary),
            "--settings",
            str(self.settings),
            "--output",
            str(self.output),
            "--bundle-version",
            "golden-1",
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        validate_bundle_root(self.output, REGRESSION_ROOT / "cases")
        manifest = _load_json(self.output / "manifest.json")
        reference = self.output / manifest["artifacts"][
            "legacy_case_warm_reference"
        ]["path"]
        self.assertEqual(manifest["bundle_class"], "golden")
        self.assertEqual(manifest["bundle_version"], "golden-1")
        self.assertEqual(reference.read_text(encoding="utf-8"), "new golden\n")
        self.assertTrue(
            (self.output / "provenance/golden_reference/comparison.json").is_file()
        )
        self.assertEqual(source_reference.read_text(encoding="utf-8"), "old reference\n")
        self.assertEqual(
            _load_json(self.candidate / "manifest.json")["bundle_class"],
            "candidate",
        )

    def test_matrix_publication_collects_every_stage(self) -> None:
        summary = self._write_matrix_summary()

        promote_bundle(
            self.settings,
            summary,
            self.output,
            "matrix-golden-1",
            REGRESSION_ROOT / "cases",
        )

        validate_bundle_root(self.output, REGRESSION_ROOT / "cases")
        manifest = _load_json(self.output / "manifest.json")
        index = self.output / manifest["artifacts"]["golden_matrix_index"]["path"]
        references = _load_json(index)["references"]
        self.assertEqual(len(references), 7)
        self.assertEqual(len({entry["artifact_id"] for entry in references}), 7)
        for entry in references:
            self.assertTrue(
                (self.output / manifest["artifacts"][entry["artifact_id"]]["path"]).is_file()
            )
        for role in ("warm_restart", "warm_reference"):
            artifact = manifest["artifacts"][manifest["roles"][role]]
            self.assertEqual(
                (self.output / artifact["path"]).read_text(encoding="utf-8"),
                "continuation_05\n",
            )

    def test_publication_applies_multiple_summaries_in_order(self) -> None:
        matrix_summary = self._write_matrix_summary()

        completed = run_command(
            "bundle",
            "promote",
            str(matrix_summary),
            str(self.summary),
            "--settings",
            str(self.settings),
            "--output",
            str(self.output),
            "--bundle-version",
            "composed-golden-1",
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        validate_bundle_root(self.output, REGRESSION_ROOT / "cases")
        manifest = _load_json(self.output / "manifest.json")
        warm_reference = manifest["artifacts"][manifest["roles"]["warm_reference"]]
        warm_restart = manifest["artifacts"][manifest["roles"]["warm_restart"]]
        self.assertEqual(
            (self.output / warm_reference["path"]).read_text(encoding="utf-8"),
            "new golden\n",
        )
        self.assertEqual(
            (self.output / warm_restart["path"]).read_text(encoding="utf-8"),
            "continuation_05\n",
        )
        self.assertTrue(
            (self.output / "provenance/golden_matrix/suite_summary.json").is_file()
        )
        self.assertTrue(
            (self.output / "provenance/golden_reference/suite_summary.json").is_file()
        )

    def test_failed_suite_is_not_publishable(self) -> None:
        summary = _load_json(self.summary)
        summary["status"] = "failed"
        _write_json(self.summary, summary)

        with self.assertRaisesRegex(BundleError, "only a passing"):
            promote_bundle(
                self.settings,
                self.summary,
                self.output,
                "golden-1",
                REGRESSION_ROOT / "cases",
            )

    def _write_canonical_summary(self) -> Path:
        run = self.root / "canonical_run"
        (run / "outputs").mkdir(parents=True)
        solution = run / "outputs/result.h5"
        solution.write_text("new golden\n", encoding="utf-8")
        manifest = _load_json(self.candidate / "manifest.json")
        reference = self.candidate / manifest["artifacts"][
            "legacy_case_warm_reference"
        ]["path"]
        identity = {
            "case_id": "legacy_case",
            "workflow_id": "warm",
            "layout_id": "mpi4_omp4",
        }
        plan = {
            **identity,
            "bundle": {
                "root": str(self.candidate),
                "bundle_id": manifest["bundle_id"],
                "bundle_version": manifest["bundle_version"],
            },
        }
        metadata = {
            "status": "completed",
            "hdf5_outputs": ["outputs/result.h5"],
            "solver": {"revision": "abc123"},
            "executable": file_identity(solution),
        }
        comparison = {
            **identity,
            "status": "passed",
            "candidate": str(solution),
            "reference": str(reference),
            "files": {
                "candidate": file_identity(solution),
                "reference": file_identity(reference),
            },
        }
        evidence = {
            "run_plan.json": plan,
            "run_metadata.json": metadata,
            "comparison.json": comparison,
        }
        for filename, document in evidence.items():
            _write_json(run / filename, document)
        (run / "stdout.log").write_text("completed\n", encoding="utf-8")
        (run / "stderr.log").write_text("", encoding="utf-8")

        path = self.root / "suite_summary.json"
        _write_json(
            path,
            {
                "schema_version": 2,
                "status": "passed",
                "suite_id": "warm",
                "run_id": "suite-test",
                "case_id": "legacy_case",
                "workflow_ids": ["warm"],
                "layout_ids": ["mpi4_omp4"],
                "comparison_mode": "immediate",
                "results": [
                    {
                        "layout_id": "mpi4_omp4",
                        "status": "passed",
                        "run_status": "completed",
                        "comparison_status": "passed",
                        "run_directory": str(run),
                    }
                ],
            },
        )
        return path

    def _write_matrix_summary(self) -> Path:
        workflow_id = "cold_fixed"
        layout_id = "mpi4_omp4"
        case = load_case_definition("legacy_case", REGRESSION_ROOT / "cases")
        stage_ids = [
            stage["stage_id"] for stage in case["workflows"][workflow_id]["stages"]
        ]
        manifest = _load_json(self.candidate / "manifest.json")
        run = self.root / "matrix_run"
        plan_stages = []
        metadata_stages = []
        for number, stage_id in enumerate(stage_ids, start=1):
            stage = run / "stages" / f"{number:02d}_{stage_id}"
            (stage / "outputs").mkdir(parents=True)
            solution = stage / "outputs/result.h5"
            solution.write_text(f"{stage_id}\n", encoding="utf-8")
            plan_stages.append({"stage_id": stage_id})
            metadata_stages.append(
                {
                    "stage_id": stage_id,
                    "run_directory": str(stage),
                    "selected_hdf5": str(solution),
                    "status": "completed",
                }
            )
        _write_json(
            run / "run_plan.json",
            {
                "case_id": "legacy_case",
                "workflow_id": workflow_id,
                "layout_id": layout_id,
                "bundle": {
                    "root": str(self.candidate),
                    "bundle_id": manifest["bundle_id"],
                    "bundle_version": manifest["bundle_version"],
                },
                "stages": plan_stages,
            },
        )
        _write_json(
            run / "run_metadata.json",
            {"status": "completed", "stages": metadata_stages},
        )
        path = self.root / "matrix_summary.json"
        _write_json(
            path,
            {
                "schema_version": 2,
                "status": "passed",
                "suite_id": "cold_matrix",
                "run_id": "matrix-test",
                "case_id": "legacy_case",
                "workflow_ids": [workflow_id],
                "layout_ids": [layout_id],
                "comparison_mode": "deferred",
                "results": [
                    {
                        "workflow_id": workflow_id,
                        "layout_id": layout_id,
                        "status": "passed",
                        "run_status": "completed",
                        "comparison_status": "not_run",
                        "run_directory": str(run),
                    }
                ],
            },
        )
        return path


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, document: dict) -> None:
    path.write_text(json.dumps(document), encoding="utf-8")
