"""Write compact grouped or legacy-flat synthetic MHDG solutions."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np


def write_solution(
    path: Path,
    *,
    grouped: bool = True,
    solution_offset: float = 0.0,
) -> None:
    mesh_values = {
        "X": np.array([[1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 0.0, 1.0]]),
        "T": np.array([[1, 2], [2, 3]], dtype=np.int32),
        "Tlin": np.array([[1, 2], [2, 3]], dtype=np.int32),
        "Tb": np.array([[1, 2], [3, 4]], dtype=np.int32),
        "Nelems": np.array([2], dtype=np.int32),
        "Nnodesperelem": np.array([2], dtype=np.int32),
    }
    solution_values = {
        "u": np.arange(1.0, 9.0) + solution_offset,
        "q": np.arange(1.0, 17.0) + solution_offset,
        "u_tilde": np.arange(1.0, 13.0) + solution_offset,
    }

    with h5py.File(path, "w") as handle:
        mesh = handle.create_group("mesh") if grouped else handle
        solution = handle.create_group("solution") if grouped else handle
        for name, values in mesh_values.items():
            mesh.create_dataset(name, data=values)
        for name, values in solution_values.items():
            solution.create_dataset(name, data=values)

        if grouped:
            parameters = handle.create_group("simulation_parameters")
            parameters.create_dataset("Neq", data=np.array([2], dtype=np.int32))
            physics = parameters.create_group("physics")
            physics.create_dataset(
                "conservative_variable_names",
                data=np.array([b"rho", b"Gamma"]),
            )
        else:
            handle.create_dataset("Neq", data=np.array([2], dtype=np.int32))
            handle.create_dataset(
                "conservative_variable_names",
                data=np.array([b"rho", b"Gamma"]),
            )

        transport = handle.create_group("transport_1d")
        transport.create_group("coefficients").create_dataset(
            "d_fs",
            data=np.array([1.0, 2.0]),
        )
        transport.create_group("profiles").create_dataset(
            "rho_grid",
            data=np.array([0.0, 1.0]),
        )
