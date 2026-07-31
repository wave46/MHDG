from __future__ import annotations

import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

import golden_update  # noqa: E402
from support.errors import BundleError  # noqa: E402
from tests.fixtures.harness import create_harness, run_command  # noqa: E402


class GoldenUpdateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        harness_root = self.root / "harness"
        harness_root.mkdir()
        self.harness = create_harness(harness_root)
        self.harness.set_bundle_class("golden")
        self.workspace = self.root / "campaign"
        self.output = self.root / "new-golden"
        self.build = SimpleNamespace(
            path=self.root / "build",
            settings_path=self.harness.settings,
            metadata_path=self.root / "build/metadata.json",
        )
        self.build.metadata_path.parent.mkdir()
        self.build.metadata_path.write_text("{}\n", encoding="utf-8")
        self.build_solver = self._patch("build_solver", return_value=self.build)
        self.run_suite = self._patch("run_suite", side_effect=self._passing_suite)
        self.compare_pairs = self._patch(
            "compare_layout_pairs",
            return_value=[{"status": "passed"}],
        )
        old_report = self.root / "old-golden.json"
        old_report.write_text("{}\n", encoding="utf-8")
        self.verify_suite = self._patch(
            "verify_suite",
            return_value=(old_report, {"status": "failed"}),
        )
        self.promote_bundle = self._patch(
            "promote_bundle",
            side_effect=lambda *args, **kwargs: self._create_candidate(args[2]),
        )
        self.promote_mapped_bundle = self._patch(
            "promote_mapped_bundle",
            side_effect=lambda *args, **kwargs: self._create_candidate(args[3]),
        )
        self.publish_campaign_bundle = self._patch(
            "publish_campaign_bundle",
            side_effect=lambda *args, **kwargs: self._create_candidate(args[1]),
        )
        self.validate_published = self._patch("_validate_published")

    def test_update_stops_at_gate_then_continues_in_order(self) -> None:
        suites = []

        def run_suite(*args, **kwargs):
            suite_id = args[1]
            suites.append((suite_id, kwargs["compare"], kwargs["resume"]))
            return self._passing_suite(*args, **kwargs)

        self.run_suite.side_effect = run_suite
        self.assertEqual(self._run(), 0)
        self.assertEqual(suites, [("cold_matrix", False, False)])
        state = self._state()
        self.assertEqual(state["status"], "awaiting_acceptance")
        self.assertEqual(state["stages"][0]["status"], "awaiting_acceptance")
        self.assertEqual(self._run("--accept", "cold_matrix"), 0)
        self.assertEqual(suites[-1][0], "warm")
        self.assertEqual(
            self._state()["stages"][1]["status"],
            "awaiting_acceptance",
        )
        self.assertEqual(self._run("--accept", "warm_reference"), 0)
        self.assertEqual(
            [suite for suite, _, _ in suites[-2:]],
            ["impurity_references", "impurity_references"],
        )
        self.assertEqual(
            self._state()["stages"][3]["status"],
            "awaiting_acceptance",
        )
        self.assertEqual(self._run("--accept", "impurity_references"), 0)

        self.assertEqual(
            [suite for suite, _, _ in suites],
            [
                "cold_matrix",
                "warm",
                "impurity_references",
                "impurity_references",
                "initialization_smoke",
                "race",
                "warm",
                "impurity_mixture",
            ],
        )
        self.assertEqual(self.promote_bundle.call_count, 1)
        self.assertEqual(self.promote_mapped_bundle.call_count, 3)
        calls = self.run_suite.call_args_list
        self.assertEqual(
            [call.args[7] for call in calls],
            ["golden", *("candidate" for _ in range(7))],
        )
        self.assertEqual(
            [Path(call.args[0]).name for call in calls[1:]],
            [
                "cold_matrix.env",
                "warm_reference.env",
                "impurity_restarts.env",
                "impurity_references.env",
                "impurity_references.env",
                "impurity_references.env",
                "impurity_references.env",
            ],
        )
        self.assertEqual(
            self.promote_bundle.call_args.kwargs["matrix_warm_roles"],
            ("warm_restart",),
        )
        state = self._state()
        self.assertEqual(state["status"], "published")
        self.assertEqual(
            state["verification_candidate"]["root"],
            state["active_bundle"],
        )
        self.assertIsNotNone(state["stages"][0].get("accepted_utc"))
        self.assertTrue(self.output.exists())
        self.assertEqual(self.publish_campaign_bundle.call_count, 1)
        published_files = dict(self.publish_campaign_bundle.call_args.args[4])
        self.assertIn("campaign.json", published_files)
        self.assertIn("build/build_metadata.json", published_files)
        self.assertIn(
            "stages/verify_impurity_mixture/suite_summary.json",
            published_files,
        )
        self.assertEqual(self._run(), 0)
        self.assertEqual(self.publish_campaign_bundle.call_count, 1)
        self.validate_published.assert_called_once()
        completed = run_command("golden", "status", str(self.workspace))
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("cold_matrix: completed", completed.stdout)

    def test_partial_warm_update_warns_and_preserves_other_components(self) -> None:
        self.assertEqual(self._run("--only", "warm"), 0)
        state = self._state()
        self.assertEqual(state["stages"][0]["status"], "skipped")
        self.assertEqual(state["stages"][1]["status"], "awaiting_acceptance")
        self.assertEqual(len(state["warnings"]), 1)

        self.assertEqual(
            self._run("--only", "warm", "--accept", "warm_reference"),
            0,
        )
        self.assertEqual(
            [call.args[1] for call in self.run_suite.call_args_list],
            ["warm", "warm", "impurity_mixture"],
        )
        state = self._state()
        self.assertEqual(state["status"], "published")
        self.assertEqual(state["stages"][2]["status"], "skipped")
        self.assertEqual(state["stages"][3]["status"], "skipped")
        self.assertEqual(self.promote_bundle.call_count, 0)
        self.assertEqual(self.promote_mapped_bundle.call_count, 1)

    def test_interrupted_suite_uses_existing_resume(self) -> None:
        calls = []

        def run_suite(*args, **kwargs):
            calls.append(kwargs["resume"])
            if len(calls) == 1:
                raise BundleError("synthetic interruption")
            return self._passing_suite(*args, **kwargs)

        self.run_suite.side_effect = run_suite
        with redirect_stderr(StringIO()):
            self.assertEqual(self._run(), 1)
            self.assertEqual(self._run(), 0)

        self.assertEqual(calls, [False, True])
        self.assertEqual(self._state()["stages"][0]["status"], "awaiting_acceptance")

    def test_changed_source_settings_reject_resume(self) -> None:
        self.assertEqual(self._run(), 0)
        with self.harness.settings.open("a", encoding="utf-8") as stream:
            stream.write("# changed\n")
        with redirect_stderr(StringIO()) as errors:
            self.assertEqual(self._run(), 1)

        self.assertIn("campaign inputs changed", errors.getvalue())

    def _arguments(self, *extra: str) -> list[str]:
        return [
            "update",
            "legacy_case",
            "--settings",
            str(self.harness.settings),
            "--run-id",
            "golden-test",
            "--workspace",
            str(self.workspace),
            "--output",
            str(self.output),
            "--bundle-version",
            "golden-test-1",
            *extra,
        ]

    def _run(self, *extra: str) -> int:
        with redirect_stdout(StringIO()):
            return golden_update.main(self._arguments(*extra))

    def _passing_suite(self, *args, **kwargs):
        path = self.root / f"{args[1]}-summary.json"
        summary = {
            "status": "passed",
            "execution_inputs": {
                "build_manifest": golden_update._file_record(
                    self.build.metadata_path,
                )
            },
            "results": [],
        }
        path.write_text(json.dumps(summary) + "\n", encoding="utf-8")
        return path, summary

    def _create_candidate(self, path: Path) -> None:
        path.mkdir(parents=True)
        (path / "manifest.json").write_text("{}\n", encoding="utf-8")

    def _patch(self, name: str, **kwargs):
        patcher = patch.object(golden_update, name, **kwargs)
        self.addCleanup(patcher.stop)
        return patcher.start()

    def _state(self) -> dict:
        return json.loads(
            (self.workspace / golden_update.STATE_FILE).read_text(encoding="utf-8")
        )


if __name__ == "__main__":
    unittest.main()
