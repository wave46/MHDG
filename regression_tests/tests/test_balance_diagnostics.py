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

    def test_checks_equation_contract_and_terminal_history(self) -> None:
        summary, _ = self._write_suite()

        report = check_suite(summary)

        self.assertEqual(report["status"], "passed")
        history = report["runs"][0]["stages"][0]["terminal_history"]
        self.assertEqual(len(history), 2)
        self.assertEqual(history[-1]["time_iteration"], 2)
        self.assertEqual(history[-1]["newton_iteration"], 3)

    def test_rejects_broken_derived_residual(self) -> None:
        summary, solution = self._write_suite()
        with h5py.File(solution, "r+") as handle:
            handle["diagnostics/equations/total_n/discrete/residual"][()] = 8.5

        report = check_suite(summary)

        self.assertEqual(report["status"], "failed")
        self.assertTrue(any("total_n/discrete/residual" in failure for failure in report["failures"]))

    def test_rejects_broken_internal_energy_exchange(self) -> None:
        summary, solution = self._write_suite()
        with h5py.File(solution, "r+") as handle:
            handle[
                "diagnostics/equations/nEe/physical/volume_components/temperature_exchange"
            ][()] = 0.85
            handle[
                "diagnostics/equations/nEe/physical/volume_components/prescribed_source"
            ][()] = 3.8

        report = check_suite(summary)

        self.assertEqual(report["status"], "failed")
        self.assertTrue(any("temperature_exchange" in failure for failure in report["failures"]))

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
                            "workflow_id": "balance_diagnostics_warm",
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
            "n": (10.0, 5.0, 2.0, -3.0, -6.0),
            "nu": (4.0, 1.0, -2.0, 0.5, -2.5),
            "nEi": (20.0, 6.0, 3.0, -2.0, -5.0),
            "nEe": (30.0, 7.0, 4.0, -1.0, -4.0),
            "n_n": (2.0, 1.0, -1.5, -0.5, -3.0),
            "total_n": (12.0, 6.0, 0.5, -3.5, -9.0),
            "total_E": (50.0, 13.0, 7.0, -3.0, -9.0),
        }
        discrete = {
            "n": (-3.0, 1.0, -2.0, -5.0),
            "nu": (0.5, 0.5, 1.0, -2.0),
            "nEi": (-4.0, 1.0, -3.0, -6.0),
            "nEe": (-3.0, 0.25, -2.75, -5.75),
            "n_n": (-0.5, 0.5, 0.0, -2.5),
            "total_n": (-3.5, 1.5, -2.0, -7.5),
            "total_E": (-7.0, 1.25, -5.75, -11.75),
        }
        units = {
            "n": ("particles", "particles/s"),
            "nu": ("kg m s^-1", "N"),
            "nEi": ("J", "W"),
            "nEe": ("J", "W"),
            "n_n": ("particles", "particles/s"),
            "total_n": ("particles", "particles/s"),
            "total_E": ("J", "W"),
        }
        volume = {
            "n": {"ionization": 3.0, "recombination": -1.0, "prescribed_source": 0.0},
            "nu": {
                "ionization": 0.2,
                "recombination": -0.4,
                "charge_exchange": -0.6,
                "pressure_divergence": -1.2,
                "prescribed_source": 0.0,
            },
            "nEi": {
                "ionization": 2.0,
                "recombination": -1.0,
                "charge_exchange": -0.5,
                "parallel_electric_work": -0.25,
                "temperature_exchange": -0.75,
                "prescribed_source": 3.5,
            },
            "nEe": {
                "ionization": -2.0,
                "recombination": 0.5,
                "radiation": -0.4,
                "ohmic": 1.0,
                "parallel_electric_work": 0.25,
                "temperature_exchange": 0.75,
                "prescribed_source": 3.9,
            },
            "n_n": {
                "ionization": -3.0,
                "recombination": 1.0,
                "prescribed_source": 0.0,
                "puff": 1.0,
                "pump": -0.5,
            },
        }
        physical_boundary = {
            "n": {"parallel_convection": -3.0, "pinch": 0.0},
            "nEi": {"sheath": -2.0, "pinch": 0.0},
            "nEe": {"sheath": -1.0, "pinch": 0.0},
            "n_n": {
                "recycling_parallel": -1.0,
                "recycling_diffusion": 0.0,
                "recycling_pinch": 0.0,
                "puff": 1.0,
                "pump": -0.5,
            },
        }
        equation_boundary = {
            "n": {"parallel_convection": -4.0, "diffusion": 1.0, "pinch": 0.0},
            "nu": {"convection": 0.3, "diffusion": 0.1, "pinch": 0.1},
            "nEi": {
                "convection": -2.0,
                "diffusion": -0.5,
                "parallel_conduction": -1.5,
                "pinch": 0.0,
            },
            "nEe": {
                "convection": -0.5,
                "diffusion": -0.5,
                "parallel_conduction": -2.0,
                "pinch": 0.0,
            },
            "n_n": {"limited_diffusion": -1.0, "limited_pressure": 0.5},
        }

        with h5py.File(path, "w") as handle:
            self._put(handle, "simulation_parameters/switches/balance_diagnostics_mode", "detailed")
            self._put(handle, "simulation_parameters/switches/neutral_wall_sources_in_elements", 0)
            self._put(handle, "simulation_parameters/physics/puff", 1.0)
            self._put(
                handle,
                "simulation_parameters/physics/conservative_variable_names",
                [b"n", b"Gamma", b"Ei", b"Ee", b"nn"],
            )
            for equation, fields in physical.items():
                prefix = f"diagnostics/equations/{equation}"
                content, temporal, source, boundary, imbalance = fields
                content_units, rate_units = units[equation]
                for name, value in {
                    "content": content,
                    "content_units": content_units,
                    "rate_units": rate_units,
                    "physical/units": rate_units,
                    "physical/temporal": temporal,
                    "physical/volume": source,
                    "physical/boundary_inward": boundary,
                    "physical/imbalance": imbalance,
                    "discrete/units": rate_units,
                }.items():
                    self._put(handle, f"{prefix}/{name}", value)
                for name, value in zip(
                    (
                        "equation_boundary_inward",
                        "tau_stabilization_inward",
                        "numerical_boundary_inward",
                        "residual",
                    ),
                    discrete[equation],
                ):
                    self._put(handle, f"{prefix}/discrete/{name}", value)
            for equation, components in volume.items():
                self._write_components(handle, equation, "physical/volume_components", components)
            for equation, components in physical_boundary.items():
                self._write_components(handle, equation, "physical/boundary_components_inward", components)
            for equation, components in equation_boundary.items():
                self._write_components(handle, equation, "discrete/boundary_components_inward", components)
            self._put(handle, "diagnostics/equations/n/physical/exchange/charge_exchange_rate", 7.0)
            self._put(handle, "diagnostics/equations/n/bc/units", "particles/s")
            self._put(handle, "diagnostics/equations/n/bc/diffusion_inward", -1.0)
            self._put(
                handle,
                "diagnostics/equations/n/bc/tau_stabilization_inward",
                1.0,
            )
            self._put(handle, "diagnostics/equations/n/bc/residual", 0.0)
            plasma_bc = {
                "nu": {
                    "units": "N",
                    "perpendicular_diffusion_inward": 0.2,
                    "split_diffusion_inward": 0.1,
                    "tau_stabilization_inward": -0.3,
                    "residual": 0.0,
                },
                "nEi": {
                    "units": "W",
                    "perpendicular_diffusion_inward": 0.4,
                    "split_diffusion_inward": 0.1,
                    "parallel_conduction_inward": -0.2,
                    "sheath_minus_bulk_inward": -0.1,
                    "tau_stabilization_inward": -0.2,
                    "residual": 0.0,
                },
                "nEe": {
                    "units": "W",
                    "perpendicular_diffusion_inward": 0.3,
                    "split_diffusion_inward": 0.1,
                    "parallel_conduction_inward": -0.15,
                    "sheath_minus_bulk_inward": -0.05,
                    "tau_stabilization_inward": -0.2,
                    "residual": 0.0,
                },
                "total_E": {
                    "units": "W",
                    "perpendicular_diffusion_inward": 0.7,
                    "split_diffusion_inward": 0.2,
                    "parallel_conduction_inward": -0.35,
                    "sheath_minus_bulk_inward": -0.15,
                    "tau_stabilization_inward": -0.4,
                    "residual": 0.0,
                },
            }
            for equation, components in plasma_bc.items():
                for name, value in components.items():
                    self._put(handle, f"diagnostics/equations/{equation}/bc/{name}", value)
            neutral_bc = {
                "units": "particles/s",
                "imposed_source_inward": 4.2,
                "physical_flux_inward": 4.1,
                "tau_stabilization_inward": 0.1,
                "residual": 0.0,
                "source_components/recycling_parallel_inward": 4.0,
                "source_components/recycling_diffusion_inward": -0.2,
                "source_components/recycling_pinch_inward": 0.0,
                "source_components/puff": 0.5,
                "source_components/pump": 0.1,
                "physical_flux_components_inward/limited_diffusion": 4.0,
                "physical_flux_components_inward/limited_pressure": 0.1,
            }
            for name, value in neutral_bc.items():
                self._put(handle, f"diagnostics/equations/n_n/bc/{name}", value)

    @staticmethod
    def _put(handle: h5py.File, path: str, value: object) -> None:
        parent, name = path.rsplit("/", 1)
        handle.require_group(parent).create_dataset(name, data=value)

    def _write_components(
        self, handle: h5py.File, equation: str, group: str, values: dict[str, float]
    ) -> None:
        for name, value in values.items():
            self._put(handle, f"diagnostics/equations/{equation}/{group}/{name}", value)

    @staticmethod
    def _terminal_block(time_iteration: int, newton_iteration: int) -> str:
        return f"""Time iteration = {time_iteration}
NR iteration: {newton_iteration}
Balance diagnostics (detailed)
  Content
    particles [particles]  n   1.0000E+01  n_n 2.0000E+00  n+n_n 1.2000E+01
    momentum [kg m s^-1]  nu 4.0000E+00
    plasma energy [J]  nEi 2.0000E+01  nEe 3.0000E+01  nEi+nEe 5.0000E+01
  Physical balances
    volume + physical boundary inward - temporal = imbalance
    n [particles/s]
                5.0000E+00  2.0000E+00 -3.0000E+00 -6.0000E+00
    nu [N]
                1.0000E+00 -2.0000E+00  5.0000E-01 -2.5000E+00
    nEi [W]
                6.0000E+00  3.0000E+00 -2.0000E+00 -5.0000E+00
    nEe [W]
                7.0000E+00  4.0000E+00 -1.0000E+00 -4.0000E+00
    n_n [particles/s]
                1.0000E+00 -1.5000E+00 -5.0000E-01 -3.0000E+00
    n+n_n [particles/s] (derived)
                6.0000E+00  5.0000E-01 -3.5000E+00 -9.0000E+00
    nEi+nEe [W] (derived)
                1.3000E+01  7.0000E+00 -3.0000E+00 -9.0000E+00
  Discrete equations
    volume + numerical boundary - temporal = residual
    equation boundary + tau stabilization = numerical boundary
    n [particles/s]
                5.0000E+00  2.0000E+00 -3.0000E+00  1.0000E+00 -2.0000E+00 -5.0000E+00
    nu [N]
                1.0000E+00 -2.0000E+00  5.0000E-01  5.0000E-01  1.0000E+00 -2.0000E+00
    nEi [W]
                6.0000E+00  3.0000E+00 -4.0000E+00  1.0000E+00 -3.0000E+00 -6.0000E+00
    nEe [W]
                7.0000E+00  4.0000E+00 -3.0000E+00  2.5000E-01 -2.7500E+00 -5.7500E+00
    n_n [particles/s]
                1.0000E+00 -1.5000E+00 -5.0000E-01  5.0000E-01  0.0000E+00 -2.5000E+00
    n+n_n [particles/s] (derived)
                6.0000E+00  5.0000E-01 -3.5000E+00  1.5000E+00 -2.0000E+00 -7.5000E+00
    nEi+nEe [W] (derived)
                1.3000E+01  7.0000E+00 -7.0000E+00  1.2500E+00 -5.7500E+00 -1.1750E+01
"""


if __name__ == "__main__":
    unittest.main()
