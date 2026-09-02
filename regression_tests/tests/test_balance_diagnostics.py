from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from check_balance_diagnostics import check_suite  # noqa: E402


class BalanceDiagnosticsCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)

    def test_checks_schema_identities_and_terminal_history(self) -> None:
        summary, solution = self._write_suite()
        report = check_suite(summary)
        self.assertEqual(report["status"], "passed")
        history = report["runs"][0]["stages"][0]["terminal_history"]
        self.assertEqual(len(history), 2)
        self.assertEqual(history[-1]["time_iteration"], 2)
        self.assertEqual(history[-1]["newton_iteration"], 3)

        with h5py.File(solution, "r+") as handle:
            handle["diagnostics/discrete/total_n/residual"][()] = 8.5
        failed = check_suite(summary)
        self.assertEqual(failed["status"], "failed")
        self.assertTrue(
            any("discrete/total_n/residual" in item for item in failed["failures"])
        )

    def test_checks_warm_run_output(self) -> None:
        summary, solution = self._write_suite()
        run = self.root / "run"
        (run / "stdout.log").write_text(
            self._terminal_block(1, 2) + self._terminal_block(2, 3),
            encoding="utf-8",
        )
        (run / "run_metadata.json").write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "status": "completed",
                    "hdf5_outputs": [str(solution.relative_to(run))],
                }
            ),
            encoding="utf-8",
        )

        report = check_suite(summary)

        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["runs"][0]["stages"][0]["stage_id"], "warm")

    def test_rejects_nonfinite_detailed_component(self) -> None:
        summary, solution = self._write_suite()
        with h5py.File(solution, "r+") as handle:
            handle[
                "diagnostics/physical/n/exchange/charge_exchange_rate"
            ][()] = float("nan")
        report = check_suite(summary)
        self.assertEqual(report["status"], "failed")
        self.assertTrue(any("charge_exchange" in item for item in report["failures"]))

    def test_requires_neutralgamma_components_for_neutralgamma_model(self) -> None:
        summary, solution = self._write_suite()
        with h5py.File(solution, "r+") as handle:
            names = "simulation_parameters/physics/conservative_variable_names"
            del handle[names]
            handle.create_dataset(
                names, data=[b"n", b"Gamma", b"Ei", b"Ee", b"nn", b"Gamman"]
            )
        report = check_suite(summary)
        self.assertEqual(report["status"], "failed")
        self.assertTrue(
            any("neutral_gamma_convection" in item for item in report["failures"])
        )

    def test_checks_relocated_puff_and_pump_placement(self) -> None:
        summary, solution = self._write_suite()
        with h5py.File(solution, "r+") as handle:
            handle[
                "simulation_parameters/switches/neutral_wall_sources_in_elements"
            ][()] = 1
            source = "diagnostics/bc/n_n/source_components"
            handle[f"{source}/puff_source"][()] = 0.0
            handle[f"{source}/pump_sink"][()] = 0.0
            handle["diagnostics/bc/n_n/imposed_source_inward"][()] = 3.8
            handle["diagnostics/bc/n_n/physical_flux_inward"][()] = 3.7
            handle[
                "diagnostics/bc/n_n/physical_flux_components_inward/limited_diffusion"
            ][()] = 3.6
        terminal = solution.parent / "stdout.log"
        text = terminal.read_text(encoding="utf-8")
        for old, new in {
            "puff source 5.000E-01": "puff source 0.000E+00",
            "pump sink (subtracted) 1.000E-01": (
                "pump sink (subtracted) 0.000E+00"
            ),
            "imposed source inward 4.200E+00": "imposed source inward 3.800E+00",
            "limited diffusion 4.000E+00": "limited diffusion 3.600E+00",
            "physical flux inward 4.100E+00": "physical flux inward 3.700E+00",
        }.items():
            text = text.replace(old, new)
        terminal.write_text(text, encoding="utf-8")

        report = check_suite(summary)
        self.assertEqual(report["status"], "passed")

        with h5py.File(solution, "r+") as handle:
            handle[
                "diagnostics/physical/n_n/volume_components/puff_source"
            ][()] = 0.75
        failed = check_suite(summary)
        self.assertEqual(failed["status"], "failed")
        self.assertTrue(any("puff_source" in item for item in failed["failures"]))

    def _write_suite(self) -> tuple[Path, Path]:
        run = self.root / "run"
        stage = run / "stages" / "time_init"
        stage.mkdir(parents=True)
        solution = stage / "solution.h5"
        self._write_solution(solution)
        (stage / "stdout.log").write_text(
            self._terminal_block(1, 2) + self._terminal_block(2, 3),
            encoding="utf-8",
        )
        (run / "run_metadata.json").write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "stages": [
                        {
                            "stage_id": "time_init",
                            "status": "completed",
                            "run_directory": str(stage),
                            "selected_hdf5": str(solution),
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        summary = self.root / "suite_summary.json"
        summary.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "results": [
                        {
                            "workflow_id": "cold_fixed_balance_diagnostics",
                            "layout_id": "mpi4_omp4",
                            "run_status": "completed",
                            "run_directory": str(run),
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        return summary, solution

    def _write_solution(self, path: Path) -> None:
        physical = {
            "n": (10.0, 5.0, 2.0, -3.0, 6.0),
            "n_n": (2.0, 1.0, -1.5, -0.5, 3.0),
            "total_n": (12.0, 6.0, 0.5, -3.5, 9.0),
        }
        physical_fields = (
            "content",
            "temporal",
            "volume",
            "boundary_physical_inward",
            "physical_imbalance",
        )
        discrete = {
            "n": (1.0, 5.0),
            "n_n": (0.5, 2.5),
            "total_n": (1.5, 7.5),
        }
        components = {
            "n/volume_components": {
                "ionization": 3.0,
                "recombination": -1.0,
                "prescribed_source": 0.0,
            },
            "n/boundary_components_inward": {
                "parallel": -4.0,
                "diffusion": 1.0,
                "pinch": 0.0,
            },
            "n_n/volume_components": {
                "ionization": -3.0,
                "recombination": 1.0,
                "prescribed_source": 0.0,
                "puff_source": 1.0,
                "pump_source": -0.5,
            },
            "n_n/boundary_components_inward": {
                "limited_diffusion": -1.0,
                "limited_pressure": 0.5,
            },
        }
        bc = {
            "n": {
                "diffusion_inward": -1.0,
                "hdg_tau_inward": 1.0,
                "residual": 0.0,
            },
            "n_n": {
                "imposed_source_inward": 4.2,
                "physical_flux_inward": 4.1,
                "hdg_tau_inward": 0.1,
                "residual": 0.0,
            },
            "n_n/source_components": {
                "recycling_parallel_inward": 4.0,
                "recycling_diffusion_inward": -0.2,
                "recycling_pinch_inward": 0.0,
                "puff_source": 0.5,
                "pump_sink": 0.1,
            },
            "n_n/physical_flux_components_inward": {
                "limited_diffusion": 4.0,
                "limited_pressure": 0.1,
            },
        }

        with h5py.File(path, "w") as handle:
            handle.create_dataset(
                "simulation_parameters/switches/balance_diagnostics_mode",
                data="detailed",
            )
            handle.create_dataset(
                "simulation_parameters/switches/neutral_wall_sources_in_elements",
                data=0,
            )
            handle.create_dataset("simulation_parameters/physics/puff", data=1.0)
            handle.create_dataset(
                "simulation_parameters/physics/conservative_variable_names",
                data=[b"n", b"Gamma", b"Ei", b"Ee", b"nn"],
            )
            for equation, fields in physical.items():
                prefix = f"diagnostics/physical/{equation}"
                handle.create_dataset(f"{prefix}/content_units", data="particles")
                handle.create_dataset(f"{prefix}/rate_units", data="particles/s")
                for name, value in zip(physical_fields, fields):
                    handle.create_dataset(f"{prefix}/{name}", data=value)
            for equation, fields in discrete.items():
                prefix = f"diagnostics/discrete/{equation}"
                handle.create_dataset(f"{prefix}/units", data="particles/s")
                for name, value in zip(("hdg_tau_inward", "residual"), fields):
                    handle.create_dataset(f"{prefix}/{name}", data=value)
            for group, fields in components.items():
                prefix = f"diagnostics/physical/{group}"
                handle.create_dataset(f"{prefix}/units", data="particles/s")
                for name, value in fields.items():
                    handle.create_dataset(f"{prefix}/{name}", data=value)
            handle.create_dataset(
                "diagnostics/physical/n/exchange/units", data="particles/s"
            )
            handle.create_dataset(
                "diagnostics/physical/n/exchange/charge_exchange_rate", data=7.0
            )
            for group, fields in bc.items():
                prefix = f"diagnostics/bc/{group}"
                handle.create_dataset(f"{prefix}/units", data="particles/s")
                for name, value in fields.items():
                    handle.create_dataset(f"{prefix}/{name}", data=value)

    @staticmethod
    def _terminal_block(time_iteration: int, newton_iteration: int) -> str:
        return f"""Time iteration = {time_iteration}
NR iteration: {newton_iteration}
Balance diagnostics (detailed)
  Physical n
      content [particles] 1.000E+01
      temporal [particles/s] 5.000E+00
      volume 2.000E+00
      boundary physical inward -3.000E+00
      physical imbalance 6.000E+00
    Discrete
      HDG tau inward 1.000E+00
      residual 5.000E+00
  Physical n_n
      content [particles] 2.000E+00
      temporal [particles/s] 1.000E+00
      volume -1.500E+00
      boundary physical inward -5.000E-01
      physical imbalance 3.000E+00
    Discrete
      HDG tau inward 5.000E-01
      residual 2.500E+00
  Physical total_n (derived)
      content [particles] 1.200E+01
      temporal [particles/s] 6.000E+00
      volume 5.000E-01
      boundary physical inward -3.500E+00
      physical imbalance 9.000E+00
    Discrete
      HDG tau inward 1.500E+00
      residual 7.500E+00
  Independent boundary-condition checks [particles/s]
    n: diffusion + HDG tau = 0
      diffusion inward -1.000E+00
      HDG tau inward 1.000E+00
      residual 0.000E+00
    n_n: imposed source - physical flux - HDG tau = 0
      Imposed-source components: recycling + puff - pump
        recycling parallel inward 4.000E+00
        recycling diffusion inward -2.000E-01
        recycling pinch inward 0.000E+00
        puff source 5.000E-01
        pump sink (subtracted) 1.000E-01
      imposed source inward 4.200E+00
      Physical-flux components inward
        limited diffusion 4.000E+00
        limited pressure 1.000E-01
      physical flux inward 4.100E+00
      HDG tau inward 1.000E-01
      residual 0.000E+00
"""


if __name__ == "__main__":
    unittest.main()
