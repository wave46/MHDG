"""Create the smallest prepared directory accepted by the legacy case."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path


PARAMETERS = """&INPUT_LST
    transport_model_path = '/old/transport_model.nml'
    impurity_model_path = '/old/impurity_model.nml'
    field_path = '/old/equilibrium.h5'
    jtor_path = '/old/current_density.h5'
    compute_from_flux = .false.
    save_folder = '/old/output/'
/
&SWITCH_LST
    steady = .false.
    saveNR = .true.
    impurity_radiation = .true.
/
&PHYS_LST
/
&NUMER_LST
    nrp = 40
    tau(6) = 1.0
    tau(7) = 1.0
    neutralp_lambda = 0.0
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
    "impurity_model_w.nml",
    "impurity_model_n.nml",
    "impurity_model_nw.nml",
    "restart.h5",
    "restart_neutral_sources_in_elements.h5",
    "reference_neutral_sources_in_elements_mpi4_omp4.h5",
    "reference_neutral_pressure_mpi4_omp4.h5",
    "reference_neutral_perpendicular_mpi4_omp4.h5",
    "reference_neutral_limiter_fixed_mpi4_omp4.h5",
    "reference_neutral_limiter_ti_mpi4_omp4.h5",
    "restart_impurity_off.h5",
    "reference_impurity_off_mpi4_omp4.h5",
    "restart_impurity_n.h5",
    "reference_impurity_n_mpi4_omp4.h5",
    "restart_impurity_nw.h5",
    "reference_impurity_nw_mpi4_omp4.h5",
)

IMPURITY_CONFIGURATIONS = {
    "impurity_model_w.nml": ("'W'", "1.0d-4"),
    "impurity_model_n.nml": ("'N'", "1.0d-2"),
    "impurity_model_nw.nml": ("'N', 'W'", "1.0d-2, 1.0d-4"),
}


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
    for filename, (names, concentrations) in IMPURITY_CONFIGURATIONS.items():
        count = names.count("'") // 2
        (directory / filename).write_text(
            "&IMPURITY_RADIATION_LST\n"
            f"  n_impurities = {count}\n"
            f"  impurity_names = {names}\n"
            f"  impurity_concentrations = {concentrations}\n"
            "/\n",
            encoding="utf-8",
        )

    reference = directory / "reference_mpi4_omp4.h5"
    if reference_writer is None:
        reference.write_text("synthetic reference\n", encoding="utf-8")
    else:
        reference_writer(reference)
    return directory
