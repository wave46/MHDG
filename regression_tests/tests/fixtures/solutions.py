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
    model: bool = False,
) -> None:
    mesh_values = {
        "X": np.array([[0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0]]),
        "T": np.array([[1, 2], [2, 4], [3, 3]], dtype=np.int32),
        "Tlin": np.array([[1, 2], [2, 4], [3, 3]], dtype=np.int32),
        "Tb": np.array([[1, 2, 4, 3], [2, 4, 3, 1]], dtype=np.int32),
        "intfaces": np.array([[1], [2], [2], [3], [2]], dtype=np.int32),
        "extfaces": np.array(
            [[1, 2, 2, 1], [1, 1, 2, 3]], dtype=np.int32
        ),
        "boundaryFlag": np.array([1, 2, 3, 4], dtype=np.int32),
        "Ndim": np.array([2], dtype=np.int32),
        "elemType": np.array([0], dtype=np.int32),
        "F": np.array([[2, 3], [1, 4], [5, 1]], dtype=np.int32),
        "Nnodes": np.array([4], dtype=np.int32),
        "Nelems": np.array([2], dtype=np.int32),
        "Nnodesperelem": np.array([3], dtype=np.int32),
        "Nnodesperface": np.array([2], dtype=np.int32),
        "Nintfaces": np.array([1], dtype=np.int32),
        "Nextfaces": np.array([4], dtype=np.int32),
        "Nfaces": np.array([5], dtype=np.int32),
    }
    names = [b"rho", b"Gamma", b"nEi", b"nEe", b"rhon"] if model else [b"rho", b"Gamma"]
    if model == "NGammaTiTeNeutralGamma":
        names.append(b"Gamman")
    count = len(names)
    solution_values = {
        "u": np.arange(1.0, 6 * count + 1.) + solution_offset,
        "q": np.arange(1.0, 12 * count + 1.) + solution_offset,
        "u_tilde": np.arange(1.0, 10 * count + 1.) + solution_offset,
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
            parameters.create_dataset("Neq", data=np.array([count], dtype=np.int32))
            physics = parameters.create_group("physics")
            physics.create_dataset(
                "conservative_variable_names",
                data=np.array(names),
            )
        else:
            handle.create_dataset("Neq", data=np.array([count], dtype=np.int32))
            handle.create_dataset(
                "conservative_variable_names",
                data=np.array(names),
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


def write_solver_output(model="NGammaTiTeNeutral", build_id="fixture-build", revision="fixture-revision", dirty=False, *, path="outputs/result.h5"):
    """Known output for a tiny executable; never infer values from its run plan."""
    write_solution(Path(path), model=model)
    with h5py.File(path, "r+") as h:
        for key, value in {
            "simulation_parameters/model": "N-Gamma-Ti-Te-Neutral" + ("Gamma" if model.endswith("Gamma") else ""),
            "simulation_parameters/Ndim": 2, "simulation_parameters/switches/balance_diagnostics_mode": "off",
            "simulation_parameters/switches/impurity_radiation": int(model != "NGammaTiTeNeutralGamma"),
            "simulation_parameters/numerics/NeutralP_lambda": 0.,
            "simulation_parameters/physics/impurity_names": "W",
            "simulation_parameters/physics/impurity_concentrations": 1e-4,
            "magnetic/jtor_source": "bicubic_psi", "provenance/build_id": build_id,
            "provenance/git_commit": revision, "provenance/git_dirty": int(dirty),
        }.items():
            h[key] = value
