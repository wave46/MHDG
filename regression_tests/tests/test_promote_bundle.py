from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from check_bundle import BundleError, validate_bundle_root  # noqa: E402
from create_bundle import create_bundle  # noqa: E402
from promote_bundle import promote_bundle  # noqa: E402


REQUIRED_FILES = (
    "mesh.msh",
    "geometry.geo",
    "equilibrium.h5",
    "current_density.h5",
    "param.txt",
    "transport_model.nml",
    "restart.h5",
    "reference_mpi4_omp4.h5",
)


class GoldenBundlePromotionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.source_bundle = self._create_bundle()
        self.settings = self.root / "settings.env"
        self.settings.write_text(
            "MHDG_REGRESSION_SETTINGS_VERSION=1\n"
            f"MHDG_REGRESSION_DATA_ROOT={self.source_bundle}\n",
            encoding="utf-8",
        )
        self.summary = self._create_summary()
        self.output = self.root / "golden"

    def test_command_creates_complete_valid_golden_bundle(self) -> None:
        completed = subprocess.run(
            [
                str(REGRESSION_ROOT / "regression.sh"),
                "--settings",
                str(self.settings),
                "bundle",
                "promote",
                str(self.summary),
                "--output",
                str(self.output),
                "--bundle-version",
                "golden-1",
            ],
            check=False,
            capture_output=True,
            env={**os.environ, "PYTHON": sys.executable},
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        validate_bundle_root(self.output, REGRESSION_ROOT / "cases")
        manifest = _load_json(self.output / "manifest.json")
        reference = self.output / manifest["artifacts"][
            "legacy_warm_reference_mpi4_omp4"
        ]["path"]
        self.assertEqual(manifest["bundle_class"], "golden")
        self.assertEqual(manifest["bundle_version"], "golden-1")
        self.assertEqual(reference.read_text(encoding="utf-8"), "new golden\n")
        self.assertTrue(
            (self.output / "case_data/legacy_fixed/restart.h5").is_file()
        )
        self.assertTrue(
            (self.output / "provenance/golden_reference/run_metadata.json").is_file()
        )
        self.assertIn("golden bundle created:", completed.stdout)

    def test_promotion_leaves_source_bundle_unchanged(self) -> None:
        source_manifest = _load_json(self.source_bundle / "manifest.json")
        source_reference = self.source_bundle / source_manifest["artifacts"][
            "legacy_warm_reference_mpi4_omp4"
        ]["path"]

        promote_bundle(
            self.settings,
            self.summary,
            self.output,
            "golden-1",
            REGRESSION_ROOT / "cases",
        )

        self.assertEqual(source_reference.read_text(encoding="utf-8"), "old reference\n")
        self.assertEqual(
            _load_json(self.source_bundle / "manifest.json")["bundle_class"],
            "candidate",
        )

    def test_matrix_promotion_collects_every_stage_without_manual_roles(self) -> None:
        summary = self._create_matrix_summary()

        promote_bundle(
            self.settings,
            summary,
            self.output,
            "matrix-golden-1",
            REGRESSION_ROOT / "cases",
        )

        validate_bundle_root(self.output, REGRESSION_ROOT / "cases")
        manifest = _load_json(self.output / "manifest.json")
        roles = manifest["case_data"]["legacy_fixed"]["roles"]
        self.assertEqual(roles["reference_matrix"], "golden_matrix_index")
        index_path = self.output / manifest["artifacts"][
            "golden_matrix_index"
        ]["path"]
        references = _load_json(index_path)["references"]
        self.assertEqual(len(references), 28)
        self.assertEqual(
            len({entry["artifact_id"] for entry in references}), 28
        )
        for entry in references:
            artifact = manifest["artifacts"][entry["artifact_id"]]
            self.assertTrue((self.output / artifact["path"]).is_file())
        reference = self.output / manifest["artifacts"][
            "legacy_warm_reference_mpi4_omp4"
        ]["path"]
        self.assertEqual(reference.read_text(encoding="utf-8"), "old reference\n")

    def test_failed_suite_is_rejected(self) -> None:
        summary = _load_json(self.summary)
        summary["status"] = "failed"
        self.summary.write_text(json.dumps(summary), encoding="utf-8")

        with self.assertRaisesRegex(BundleError, "only a passing"):
            promote_bundle(
                self.settings,
                self.summary,
                self.output,
                "golden-1",
                REGRESSION_ROOT / "cases",
            )

    def test_existing_output_is_not_replaced(self) -> None:
        self.output.mkdir()
        marker = self.output / "keep"
        marker.write_text("keep\n", encoding="utf-8")

        with self.assertRaisesRegex(BundleError, "output already exists"):
            promote_bundle(
                self.settings,
                self.summary,
                self.output,
                "golden-1",
                REGRESSION_ROOT / "cases",
            )
        self.assertEqual(marker.read_text(encoding="utf-8"), "keep\n")

    def _create_bundle(self) -> Path:
        prepared = self.root / "prepared"
        prepared.mkdir()
        for filename in REQUIRED_FILES:
            contents = "old reference\n" if filename.startswith("reference") else filename
            (prepared / filename).write_text(contents, encoding="utf-8")
        bundle = self.root / "candidate"
        create_bundle("legacy_fixed", prepared, bundle, REGRESSION_ROOT / "cases")
        return bundle

    def _create_summary(self) -> Path:
        result = self._create_canonical_run()
        summary = {
            "schema_version": 1,
            "status": "passed",
            "suite_id": "warm",
            "run_id": "suite-test",
            "case_id": "legacy_fixed",
            "workflow_id": "warm",
            "results": [result],
        }
        path = self.root / "suite_summary.json"
        path.write_text(json.dumps(summary), encoding="utf-8")
        return path

    def _create_canonical_run(self) -> dict[str, object]:
        run = self.root / "run"
        (run / "outputs").mkdir(parents=True)
        solution = run / "outputs/result.h5"
        solution.write_text("new golden\n", encoding="utf-8")
        manifest = _load_json(self.source_bundle / "manifest.json")
        reference = self.source_bundle / manifest["artifacts"][
            "legacy_warm_reference_mpi4_omp4"
        ]["path"]
        plan = {
            "case_id": "legacy_fixed",
            "workflow_id": "warm",
            "layout_id": "mpi4_omp4",
            "bundle": {
                "root": str(self.source_bundle),
                "bundle_id": manifest["bundle_id"],
                "bundle_version": manifest["bundle_version"],
            },
        }
        metadata = {
            "status": "completed",
            "hdf5_outputs": ["outputs/result.h5"],
            "solver": {"revision": "abc123", "build_description": "test"},
            "executable": _identity(solution),
        }
        comparison = {
            "status": "passed",
            "case_id": "legacy_fixed",
            "workflow_id": "warm",
            "layout_id": "mpi4_omp4",
            "candidate": str(solution),
            "reference": str(reference),
            "files": {
                "candidate": _identity(solution),
                "reference": _identity(reference),
            },
        }
        files = {
            "run_plan.json": json.dumps(plan),
            "run_metadata.json": json.dumps(metadata),
            "comparison.json": json.dumps(comparison),
            "stdout.log": "Output written to file outputs/result.h5\n",
            "stderr.log": "",
        }
        for name, contents in files.items():
            (run / name).write_text(contents, encoding="utf-8")
        return {
            "layout_id": "mpi4_omp4",
            "status": "passed",
            "run_status": "completed",
            "comparison_status": "passed",
            "run_directory": str(run),
        }

    def _create_matrix_summary(self) -> Path:
        manifest = _load_json(self.source_bundle / "manifest.json")
        case = _load_json(REGRESSION_ROOT / "cases/legacy_fixed.json")
        workflows = ["cold_fixed", "cold_adaptive"]
        layouts = ["serial_omp1", "mpi2_omp4"]
        results = []

        for workflow_id in workflows:
            stage_ids = [
                stage["stage_id"]
                for stage in case["workflows"][workflow_id]["stages"]
            ]
            for layout_id in layouts:
                run = self.root / f"run_{workflow_id}_{layout_id}"
                run.mkdir()
                plan_stages = []
                metadata_stages = []
                for number, stage_id in enumerate(stage_ids, start=1):
                    stage_dir = run / "stages" / f"{number:02d}_{stage_id}"
                    outputs = stage_dir / "outputs"
                    outputs.mkdir(parents=True)
                    solution = outputs / "result.h5"
                    solution.write_text(
                        f"{workflow_id} {layout_id} {stage_id}\n",
                        encoding="utf-8",
                    )
                    plan_stages.append({"stage_id": stage_id})
                    metadata_stages.append(
                        {
                            "stage_id": stage_id,
                            "run_directory": str(stage_dir),
                            "selected_hdf5": str(solution),
                            "status": "completed",
                        }
                    )

                plan = {
                    "case_id": "legacy_fixed",
                    "workflow_id": workflow_id,
                    "layout_id": layout_id,
                    "bundle": {
                        "root": str(self.source_bundle),
                        "bundle_id": manifest["bundle_id"],
                        "bundle_version": manifest["bundle_version"],
                    },
                    "stages": plan_stages,
                }
                metadata = {"status": "completed", "stages": metadata_stages}
                (run / "run_plan.json").write_text(
                    json.dumps(plan), encoding="utf-8"
                )
                (run / "run_metadata.json").write_text(
                    json.dumps(metadata), encoding="utf-8"
                )
                results.append(
                    {
                        "workflow_id": workflow_id,
                        "layout_id": layout_id,
                        "status": "passed",
                        "run_status": "completed",
                        "comparison_status": "not_run",
                        "run_directory": str(run),
                    }
                )

        summary = {
            "schema_version": 1,
            "status": "passed",
            "suite_id": "cold_matrix",
            "run_id": "matrix-test",
            "case_id": "legacy_fixed",
            "workflow_ids": workflows,
            "layout_ids": layouts,
            "comparison_mode": "deferred",
            "results": results,
        }
        path = self.root / "matrix_summary.json"
        path.write_text(json.dumps(summary), encoding="utf-8")
        return path


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _identity(path: Path) -> dict[str, object]:
    contents = path.read_bytes()
    return {
        "size_bytes": len(contents),
        "sha256": hashlib.sha256(contents).hexdigest(),
    }


if __name__ == "__main__":
    unittest.main()
