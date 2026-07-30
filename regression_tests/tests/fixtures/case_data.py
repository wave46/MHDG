"""Create the smallest prepared directory accepted by the legacy case."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path


PARAMETERS = """&INPUT_LST
    transport_model_path = '/old/transport_model.nml'
    field_path = '/old/equilibrium.h5'
    jtor_path = '/old/current_density.h5'
    save_folder = '/old/output/'
/
&SWITCH_LST
    steady = .false.
    saveNR = .true.
    impurity_radiation = .true.
/
&PHYS_LST
    impurity_name = 'W'
    impurity_concentration = 1e-4
/
&NUMER_LST
    nrp = 40
/
&ADAPT_LST
    adaptivity = .true.
    time_adapt = .true.
    NR_adapt = .true.
    div_adapt = .true.
    rest_adapt = .true.
    osc_adapt = .true.
    geometry_path = '/old/geometry.geo'
/
&TIME_LST
    nts = 10
/
"""

COLD_PARAMETER_FILES = (
    "param_cold_fixed_time_init.txt",
    "param_cold_fixed_diffusion_reduction.txt",
    *(f"param_cold_fixed_continuation_{index:02d}.txt" for index in range(1, 6)),
)
COLD_TRANSPORT_FILES = (
    "transport_cold_fixed_initial.nml",
    *(
        f"transport_cold_fixed_continuation_{index:02d}.nml"
        for index in range(1, 6)
    ),
)
BASE_FILES = (
    "mesh.msh",
    "mesh_adaptive_initial.msh",
    "geometry.geo",
    "equilibrium.h5",
    "current_density.h5",
    "transport_model.nml",
    "restart.h5",
    "restart_impurity_off.h5",
    "reference_impurity_off_mpi4_omp4.h5",
    "restart_impurity_n.h5",
    "reference_impurity_n_mpi4_omp4.h5",
)


def write_case_source(
    directory: Path,
    reference_writer: Callable[[Path], None] | None = None,
) -> Path:
    """Write generic warm and cold artifacts for ``legacy_case``."""
    directory.mkdir()
    for filename in BASE_FILES:
        (directory / filename).write_text(
            f"synthetic {filename}\n",
            encoding="utf-8",
        )
    (directory / "param.txt").write_text(PARAMETERS, encoding="utf-8")
    for filename in COLD_PARAMETER_FILES:
        (directory / filename).write_text(PARAMETERS, encoding="utf-8")
    for filename in COLD_TRANSPORT_FILES:
        (directory / filename).write_text(
            f"synthetic {filename}\n",
            encoding="utf-8",
        )

    reference = directory / "reference_mpi4_omp4.h5"
    if reference_writer is None:
        reference.write_text("synthetic reference\n", encoding="utf-8")
    else:
        reference_writer(reference)
    return directory
