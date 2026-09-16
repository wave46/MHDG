"""Create the smallest prepared directory accepted by the legacy case."""

from __future__ import annotations

import json
import shutil
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
&UTILS_LST
/
&NUMER_LST
    nrp = 40
    tNR = 2e-4
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

FILES = ("mesh.msh", "mesh_adaptive_initial.msh", "geometry.geo", "equilibrium.h5",
         "current_density.h5", "transport_model.nml", "transport_cold_fixed_initial.nml")


def write_case_source(directory):
    """Only warm and initial-stage inputs; no inventory of production features."""
    from .solutions import write_solver_output
    directory.mkdir()
    for name in FILES:
        (directory / name).write_text(f"fixture {name}\n")
    for name in ("param.txt", "param_cold_fixed_time_init.txt"):
        (directory / name).write_text(PARAMETERS)
    (directory / "impurity_model_w.nml").write_text(
        "&IMPURITY_RADIATION_LST\n impurity_names = 'W'\n impurity_concentrations = 1.0d-4\n/\n")
    write_solver_output(path=directory / "reference_mpi4_omp4.h5")
    shutil.copy2(directory / "reference_mpi4_omp4.h5", directory / "restart.h5")
    return directory


def write_catalog(root):
    """An explicit small mechanics catalog, independent of production cold recipes."""
    repository = Path(__file__).resolve().parents[2]
    shutil.copytree(repository / "schemas", root / "schemas")
    for name in ("layouts.json", "tolerances.json"):
        shutil.copy2(repository / name, root / name)
    (root / "cases").mkdir()
    required = dict(zip(("mesh", "coarse_mesh", "geometry", "equilibrium_magnetic_field",
                        "equilibrium_current_density", "transport_configuration", "initial_transport"), FILES))
    required.update(warm_parameters="param.txt", impurity_configuration="impurity_model_w.nml")
    optional = {"warm_restart": "restart.h5", "warm_reference": "reference_mpi4_omp4.h5",
                "initial_parameters": "param_cold_fixed_time_init.txt"}
    warm = {"description": "Warm fixture", "type": "warm_same_state", "inputs": list(required),
            "impurity_configuration": "impurity_configuration", "layout": "mpi4_omp4",
            "restart": "warm_restart", "reference": "warm_reference", "outputs": ["warm_reference"],
            "parameter_overrides": {"compute_from_flux": True},
            "comparison": {"method": "fixed_hdf5", "profile": "fixed_same_layout", "cross_layout_profile": "fixed_cross_layout"}}
    stages = [{"id": name, "parameters": "initial_parameters", "transport": "initial_transport",
               "newton_check": "finite_only" if name == "initial" else "bounded"}
              for name in ("initial", "continued", "final")]
    workflows = {
        "warm": warm,
        "cold_fixed": {"description": "Fixed fixture", "type": "staged_fixed_mesh",
                       "inputs": list(required), "mesh": "mesh", "layout": "mpi4_omp4",
                       "impurity_configuration": "impurity_configuration", "outputs": ["warm_restart"],
                       "reference": "warm_reference", "stages": stages,
                       "parameter_overrides": {"compute_from_flux": True, "rest_adapt": False},
                       "comparison": {"method": "fixed_hdf5", "profile": "cold_fixed_reference", "stage_profile": "fixed_stage_reference"}},
        "cold_adaptive": {"extends": "cold_fixed", "type": "staged_adaptive_mesh", "mesh": "coarse_mesh",
                          "adaptive_stages": ["initial", "continued"], "outputs": [],
                          "comparison": {"method": "mesh_independent", "profile": "adaptive_reference",
                                         "stage_profile": "adaptive_reference", "direct_stage_profile": "fixed_stage_reference"}},
        "cold_step_fixed": {"extends": "cold_fixed", "layout": "serial_omp1", "stages": stages[:1],
                            "outputs": [], "comparison": {"method": "fixed_hdf5", "profile": "race_step", "stage_profile": "race_step"}},
        "cold_step_adaptive": {"extends": "cold_step_fixed", "type": "staged_adaptive_mesh", "adaptive_stages": ["initial"]},
        "cold_step_neutralgamma": {"extends": "cold_step_fixed", "model": "NGammaTiTeNeutralGamma",
                                   "parameter_overrides": {"impurity_radiation": False}},
    }
    documents = {"workflows.json": {"schema_version": 2, "workflows": workflows,
                                    "parameter_namelists": {"balance_diagnostics_mode": "utils_lst"}},
                 "suites.json": {"schema_version": 2, "defaults": {"case": "diverted_case", "layout": "mpi4_omp4"},
                    "suites": {
                        "warm": {"description": "Warm", "workflows": ["warm"]},
                        "warm_parallelism": {"description": "Resume cells", "workflows": ["warm"], "layouts": "all"},
                        "parallel": {"description": "One pair", "workflows": ["cold_step_fixed"],
                                     "layouts": ["serial_omp1", "mpi4_omp4"], "relations": ["all_pairs"], "tolerance_profile": "race_step"},
                        "initialization": {"description": "File validity", "workflows": ["cold_step_fixed"], "reference_comparisons": False},
                        "bootstrap": {"description": "Stage convergence", "workflows": ["cold_adaptive"], "reference_comparisons": False},
                        "neutralgamma": {"description": "Model selection", "case": "legacy_case", "workflows": ["cold_step_neutralgamma"], "layouts": ["serial_omp1"], "reference_comparisons": False}},
                    "profiles": {"routine-extended": {"description": "Fixture profile", "checks": [{"suite": "warm"}, {"suite": "parallel"}]},
                                 "full": {"description": "Mixed models", "include": ["routine-extended"], "checks": [{"suite": "neutralgamma"}]}}}}
    for case in ("legacy_case", "diverted_case"):
        documents[f"cases/{case}.json"] = {"schema_version": 2, "description": "Mechanics fixture",
            "reference": {"branch": "test", "revision": "a" * 40}, "files": {"required": required, "optional": optional},
            "workflows": {name: {} for name in workflows}}
    for name, document in documents.items():
        (root / name).write_text(json.dumps(document))
    return root
