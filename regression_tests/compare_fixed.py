"""Direct comparison of mesh, conservative fields, transport and magnetic data."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from .support import ComparisonError
from .compare_common import calculate_error_norms


def compare_hdf5_files(
    reference_path: Path,
    candidate_path: Path,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
    """Compare grouped solutions on the same fixed mesh."""
    failures: list[str] = []
    report: dict[str, Any] = {"failures": failures}

    try:
        with h5py.File(reference_path, "r") as reference, h5py.File(
            candidate_path, "r"
        ) as candidate:
            report["mesh"] = compare_mesh(
                reference,
                candidate,
                tolerances,
                failures,
            )
            if not report["mesh"]["passed"]:
                report["solution"] = {
                    "passed": False,
                    "reason": "mesh comparison failed",
                }
            else:
                report["solution"] = compare_solution(
                    reference,
                    candidate,
                    tolerances,
                    failures,
                )
            report["transport_1d"] = compare_transport(
                reference,
                candidate,
                tolerances,
                failures,
            )
            report["magnetic"] = compare_magnetic(
                reference,
                candidate,
                tolerances,
                failures,
            )
    except (OSError, KeyError) as exc:
        raise ComparisonError(f"cannot read HDF5 comparison data: {exc}") from exc
    except ComparisonError as exc:
        failures.append(str(exc))

    report["status"] = "passed" if not failures else "failed"
    return report


def validate_solution_file(path: Path) -> list[str]:
    """Check one solution's mesh, field storage and finite optional numeric data."""
    failures = []

    def finite(values, label):
        if values.size == 0 or not _finite(values):
            failures.append(f"{label} is empty or contains non-finite values")

    try:
        with h5py.File(path, "r") as handle:
            read_discrete_mesh(handle)
            equations = _equation_count(handle)
            for name in ("u", "q", "u_tilde"):
                finite(read_solution_field(handle, name, equations), f"solution/{name}")
            transport = optional_group(handle, "transport_1d")
            for name, dataset in transport_datasets(transport):
                finite(np.asarray(dataset), f"transport_1d/{name}")
            magnetic = optional_group(handle, "magnetic")
            if magnetic is not None:
                for name, dataset in magnetic.items():
                    if isinstance(dataset, h5py.Dataset) and dataset.dtype.kind not in "biuOSU":
                        finite(np.asarray(dataset), f"magnetic/{name}")
    except (OSError, KeyError, TypeError, ValueError) as exc:
        failures.append(str(exc))
    return failures


def required_array(handle: h5py.File, group: str, name: str) -> np.ndarray:
    """Read a required dataset from grouped storage."""
    array = optional_array(handle, group, name)
    if array is None:
        raise ComparisonError(f"required dataset is missing: {group}/{name}")
    return array


def required_scalar(handle: h5py.File, group: str, name: str) -> int | float:
    """Read the scalar value stored in a required MHDG dataset."""
    return required_array(handle, group, name).reshape(-1)[0].item()


def optional_array(
    handle: h5py.File,
    group: str,
    name: str,
) -> np.ndarray | None:
    """Read an optional dataset from grouped storage."""
    path = f"{group}/{name}"
    return np.asarray(handle[path]) if path in handle and isinstance(handle[path], h5py.Dataset) else None


def optional_group(handle: h5py.File, path: str) -> h5py.Group | None:
    """Return an optional HDF5 group."""
    return handle[path] if path in handle and isinstance(handle[path], h5py.Group) else None


def read_discrete_mesh(handle: h5py.File) -> dict[str, Any]:
    """Read a 2D triangular mesh after checking counts, connectivity and geometry."""
    path = handle.filename
    values = {}
    for name in ("Ndim", "elemType", "Nnodes", "Nelems", "Nnodesperelem",
                 "Nnodesperface", "Nfaces", "Nintfaces", "Nextfaces"):
        value = required_array(handle, "mesh", name)
        if value.size != 1 or value.dtype.kind not in "iu" or value.item() < 0:
            raise ComparisonError(f"invalid mesh/{name} in {path}")
        values[name] = int(value.item())
    dim, nodes, elements = values["Ndim"], values["Nnodes"], values["Nelems"]
    npe, npf = values["Nnodesperelem"], values["Nnodesperface"]
    faces, internal, external = values["Nfaces"], values["Nintfaces"], values["Nextfaces"]
    if dim != 2 or values["elemType"] != 0:
        raise ComparisonError(f"mesh comparison currently requires 2D triangles: {path}")
    if min(nodes, elements, faces) < 1 or npf < 2 or npe != npf * (npf + 1) // 2:
        raise ComparisonError(f"invalid mesh counts or polynomial order in {path}")
    if faces != internal + external or 3 * elements != 2 * internal + external:
        raise ComparisonError(f"inconsistent mesh face counts in {path}")
    for name, shape in {
        "X": (2, nodes), "T": (npe, elements), "Tlin": (3, elements),
        "Tb": (npf, external), "F": (3, elements),
        "intfaces": (5, internal), "extfaces": (2, external),
    }.items():
        array = required_array(handle, "mesh", name)
        if array.shape != shape or not np.isfinite(array).all():
            raise ComparisonError(f"invalid mesh/{name} shape or non-finite values in {path}")
        if name != "X" and array.dtype.kind not in "iu":
            raise ComparisonError(f"mesh/{name} must contain integer indices in {path}")
        values[name] = array
    for name, maximum in (("T", nodes), ("Tlin", nodes), ("Tb", nodes), ("F", faces)):
        if ((values[name] < 1) | (values[name] > maximum)).any():
            raise ComparisonError(f"mesh/{name} contains invalid indices in {path}")
    if (np.diff(np.sort(values["T"], axis=0), axis=0) == 0).any():
        raise ComparisonError(f"mesh/T contains repeated element nodes in {path}")
    if not (values["Tlin"][:, None, :] == values["T"][None, :, :]).any(axis=1).all():
        raise ComparisonError(f"mesh/Tlin vertices do not belong to their elements in {path}")
    incidence = np.bincount(values["F"].ravel().astype(np.intp), minlength=faces + 1)[1:]
    if np.count_nonzero(incidence == 1) != external or np.count_nonzero(incidence == 2) != internal:
        raise ComparisonError(f"mesh/F has inconsistent face incidence in {path}")
    for name, row, maximum in (("intfaces", [0, 2], elements), ("intfaces", [1, 3], 3),
                               ("extfaces", [0], elements), ("extfaces", [1], 3)):
        indices = values[name][row]
        if ((indices < 1) | (indices > maximum)).any():
            raise ComparisonError(f"mesh/{name} contains invalid indices in {path}")
    vertices = values["X"][:, values["Tlin"] - 1]
    first, second = vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]
    area = first[0] * second[1] - first[1] * second[0]
    if not np.isfinite(area).all() or (area == 0).any():
        raise ComparisonError(f"mesh/Tlin contains degenerate triangles in {path}")
    return values


def numeric_comparison(
    reference: np.ndarray,
    candidate: np.ndarray,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
    """Build the standard fixed-mesh numeric tolerance report."""
    norms = calculate_error_norms(reference, candidate)
    if not norms.available:
        return {
            "passed": False,
            "finite": norms.finite,
            "reference_shape": list(reference.shape),
            "candidate_shape": list(candidate.shape),
            "relative_l2": None,
            "normalized_linf": None,
        }

    relative_l2 = norms.relative_l2
    normalized_linf = norms.normalized_linf
    assert relative_l2 is not None and normalized_linf is not None
    passed = (
        relative_l2 <= tolerances["relative_l2_max"]
        and normalized_linf <= tolerances["normalized_linf_max"]
    )
    return {
        "passed": passed,
        "finite": norms.finite,
        "reference_shape": list(reference.shape),
        "candidate_shape": list(candidate.shape),
        "relative_l2": relative_l2,
        "relative_l2_max": tolerances["relative_l2_max"],
        "normalized_linf": normalized_linf,
        "normalized_linf_max": tolerances["normalized_linf_max"],
    }


def compare_mesh(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare mesh coordinates and connectivity in storage order."""
    connectivity = {}
    for name in ("T", "Tlin", "Tb", "F", "intfaces", "extfaces", "Ndim", "elemType",
                 "Nnodes", "Nelems", "Nnodesperelem", "Nnodesperface", "Nfaces", "Nintfaces", "Nextfaces"):
        details = _compare_connectivity(reference, candidate, name)
        if details is None:
            continue
        connectivity[name] = details

    coordinates = _compare_coordinates(reference, candidate, tolerances)
    report = {
        "mode": "exact",
        "connectivity": connectivity,
        "coordinates": coordinates,
    }
    report["passed"] = coordinates["passed"] and all(
        details["passed"] for details in connectivity.values()
    )
    _record_mesh_failures(report, failures)
    return report


def _record_mesh_failures(
    report: dict[str, Any], failures: list[str]
) -> None:
    coordinates = report.get("coordinates", {})
    if not coordinates.get("passed"):
        failures.append(
            f"mesh/X: {coordinates.get('reason', 'coordinates differ')}"
        )
        return
    connectivity = report.get("connectivity", {})
    for name, details in connectivity.items():
        if not details.get("passed"):
            reason = details.get("reason", "connectivity differs")
            failures.append(f"mesh/{name}: {reason}")


def _compare_connectivity(
    reference: h5py.File,
    candidate: h5py.File,
    name: str,
) -> dict[str, Any] | None:
    first = optional_array(reference, "mesh", name)
    second = optional_array(candidate, "mesh", name)
    if first is None and second is None and name == "Tb":
        return None
    if first is None or second is None:
        return {"passed": False, "reason": "dataset missing"}

    return {
        "passed": first.shape == second.shape and np.array_equal(first, second),
        "reference_shape": list(first.shape),
        "candidate_shape": list(second.shape),
    }


def _compare_coordinates(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
    first = required_array(reference, "mesh", "X")
    second = required_array(candidate, "mesh", "X")
    finite = _finite(first) and _finite(second)
    same_shape = first.shape == second.shape
    maximum_error = (
        float(np.max(np.abs(second - first))) if finite and same_shape else None
    )
    passed = (
        finite
        and same_shape
        and maximum_error is not None
        and maximum_error <= tolerances["mesh_coordinate_atol"]
    )
    return {
        "passed": passed,
        "finite": finite,
        "reference_shape": list(first.shape),
        "candidate_shape": list(second.shape),
        "maximum_absolute_error": maximum_error,
        "absolute_tolerance": tolerances["mesh_coordinate_atol"],
    }


def _finite(array: np.ndarray) -> bool:
    return bool(np.isfinite(np.asarray(array, dtype=float)).all())


def compare_solution(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare each conservative field equation by equation."""
    equation_count = _matching_equation_count(reference, candidate)
    names, names_match = _equation_names(reference, candidate, equation_count)
    if not names_match:
        failures.append("conservative variable names differ")

    datasets = {
        name: _compare_dataset(
            reference,
            candidate,
            name,
            equation_count,
            names,
            tolerances,
            failures,
        )
        for name in ("u", "q", "u_tilde")
    }
    return {
        "passed": all(dataset["passed"] for dataset in datasets.values())
        and names_match,
        "equation_count": equation_count,
        "equation_names": names,
        "datasets": datasets,
    }


def _compare_dataset(
    reference: h5py.File,
    candidate: h5py.File,
    dataset_name: str,
    equation_count: int,
    equation_names: list[str],
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    reference_values = read_solution_field(reference, dataset_name, equation_count)
    candidate_values = read_solution_field(candidate, dataset_name, equation_count)
    report: dict[str, Any] = {
        "reference_size": int(reference_values.size),
        "candidate_size": int(candidate_values.size),
        "equations": {},
    }

    reference_by_equation = np.moveaxis(reference_values, 2, -1).reshape(
        -1, equation_count
    )
    candidate_by_equation = np.moveaxis(candidate_values, 2, -1).reshape(
        -1, equation_count
    )
    for index, name in enumerate(equation_names):
        metrics = numeric_comparison(
            reference_by_equation[:, index],
            candidate_by_equation[:, index],
            tolerances,
        )
        report["equations"][name] = metrics
        if not metrics["passed"]:
            failures.append(f"solution/{dataset_name}/{name} exceeds tolerance")
    report["passed"] = all(
        result["passed"] for result in report["equations"].values()
    )
    return report


def read_solution_field(
    handle: h5py.File,
    name: str,
    equation_count: int,
) -> np.ndarray:
    """Read a conservative field with the size and storage order required by its mesh."""
    values = required_array(handle, "solution", name)
    element_count = int(required_scalar(handle, "mesh", "Nelems"))
    nodes_per_element = int(
        required_scalar(handle, "mesh", "Nnodesperelem")
    )
    if name == "u":
        shape = (element_count, nodes_per_element, equation_count)
    elif name == "q":
        dimension = int(required_scalar(handle, "mesh", "Ndim"))
        shape = (element_count, nodes_per_element, equation_count, dimension)
    else:
        face_count = int(required_scalar(handle, "mesh", "Nfaces"))
        nodes_per_face = int(
            required_scalar(handle, "mesh", "Nnodesperface")
        )
        shape = (face_count, nodes_per_face, equation_count)
    if equation_count < 1 or values.size != np.prod(shape):
        raise ComparisonError(f"solution/{name} size is inconsistent with the mesh in {handle.filename}")
    return values.reshape(shape)


def _matching_equation_count(reference: h5py.File, candidate: h5py.File) -> int:
    reference_count = _equation_count(reference)
    candidate_count = _equation_count(candidate)
    if reference_count != candidate_count:
        raise ComparisonError(
            f"equation counts differ: {reference_count} != {candidate_count}"
        )
    return reference_count


def _equation_names(
    reference: h5py.File,
    candidate: h5py.File,
    count: int,
) -> tuple[list[str], bool]:
    reference_names = _valid_equation_names(reference, count)
    candidate_names = _valid_equation_names(candidate, count)
    names = reference_names or candidate_names or [
        f"equation_{index + 1}" for index in range(count)
    ]
    names_match = not (
        reference_names and candidate_names and reference_names != candidate_names
    )
    return names, names_match


def _equation_count(handle: h5py.File) -> int:
    if "simulation_parameters/Neq" in handle:
        return int(np.asarray(handle["simulation_parameters/Neq"]).reshape(-1)[0])

    names = _optional_names(handle)
    if names:
        return len(names)

    values = required_array(handle, "solution", "u").size
    elements = int(required_scalar(handle, "mesh", "Nelems"))
    nodes = int(required_scalar(handle, "mesh", "Nnodesperelem"))
    if elements <= 0 or nodes <= 0 or values % (elements * nodes):
        raise ComparisonError("cannot infer the number of equations")
    return values // (elements * nodes)


def _valid_equation_names(handle: h5py.File, count: int) -> list[str]:
    names = _optional_names(handle)
    return names if len(names) == count else []


def _optional_names(handle: h5py.File) -> list[str]:
    path = "simulation_parameters/physics/conservative_variable_names"
    if path not in handle:
        return []
    return [
        value.decode(errors="replace").strip(" \x00")
        if isinstance(value, bytes)
        else str(value).strip()
        for value in np.asarray(handle[path]).reshape(-1)
    ]


def compare_transport(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare reference transport datasets when they are present."""
    reference_group = optional_group(reference, "transport_1d")
    candidate_group = optional_group(candidate, "transport_1d")
    if reference_group is None:
        return {"present": False, "passed": True, "datasets": {}}
    if candidate_group is None:
        failures.append("candidate is missing transport_1d data")
        return {"present": True, "passed": False, "datasets": {}}

    datasets = {}
    for label, reference_dataset in transport_datasets(reference_group):
        candidate_dataset = candidate_group.get(label)
        if not isinstance(candidate_dataset, h5py.Dataset):
            datasets[label] = {"passed": False, "reason": "dataset missing"}
            failures.append(f"transport_1d/{label} is missing")
            continue
        metrics = numeric_comparison(np.asarray(reference_dataset), np.asarray(candidate_dataset), tolerances)
        datasets[label] = metrics
        if not metrics["passed"]:
            failures.append(f"transport_1d/{label} exceeds tolerance")
    return {
        "present": True,
        "passed": all(result["passed"] for result in datasets.values()),
        "datasets": datasets,
    }


def transport_datasets(group):
    """Iterate the optional coefficient/profile datasets used by both checks."""
    for section in ("coefficients", "profiles"):
        if group is None or section not in group:
            continue
        if not isinstance(group[section], h5py.Group):
            raise ComparisonError(f"transport_1d/{section} must be a group")
        for name, dataset in group[section].items():
            if isinstance(dataset, h5py.Dataset):
                yield f"{section}/{name}", dataset


def compare_magnetic(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare every direct dataset in the reference magnetic group."""
    reference_group = optional_group(reference, "magnetic")
    candidate_group = optional_group(candidate, "magnetic")
    if reference_group is None:
        return {"present": False, "passed": True, "datasets": {}}
    if candidate_group is None:
        failures.append("candidate is missing magnetic data")
        return {"present": True, "passed": False, "datasets": {}}

    datasets: dict[str, dict[str, Any]] = {}
    for name, reference_dataset in reference_group.items():
        if not isinstance(reference_dataset, h5py.Dataset):
            continue
        candidate_dataset = candidate_group.get(name)
        if not isinstance(candidate_dataset, h5py.Dataset):
            datasets[name] = {"passed": False, "reason": "dataset missing"}
            failures.append(f"magnetic/{name} is missing")
            continue

        reference_array = np.asarray(reference_dataset)
        candidate_array = np.asarray(candidate_dataset)
        if reference_array.dtype.kind in "biuOSU":
            metrics = _exact_comparison(reference_array, candidate_array)
        else:
            metrics = numeric_comparison(
                reference_array,
                candidate_array,
                tolerances,
            )
        datasets[name] = metrics
        if not metrics["passed"]:
            failures.append(f"magnetic/{name} differs")

    return {
        "present": True,
        "passed": all(result["passed"] for result in datasets.values()),
        "datasets": datasets,
    }


def _exact_comparison(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict[str, Any]:
    """Compare categorical and integer magnetic metadata exactly."""
    return {
        "passed": np.array_equal(reference, candidate),
        "mode": "exact",
        "reference_shape": list(reference.shape),
        "candidate_shape": list(candidate.shape),
    }


def check_output_contract(path, directory, model, provenance, overrides):
    """Check the requested model, build provenance and explicit feature switches."""
    from .config import MODELS
    from .parameters import read_selected_input_values

    names = {"impurity_radiation", "transport_1d", "neutral_wall_sources_in_elements",
             "neutral_perpendicular_diffusion", "neutralp_lambda", "neutral_flux_limiter_mode",
             "neutral_flux_limiter_tn_source", "neutral_flux_limiter_tn_ev", "compute_from_flux"}
    overrides = {**read_selected_input_values(directory / "param.txt", names),
                 **{key.lower(): value for key, value in overrides.items()}}
    label, equations = MODELS[model]
    expected = {"simulation_parameters/model": label, "simulation_parameters/Ndim": 2,
                "simulation_parameters/Neq": len(equations),
                "simulation_parameters/physics/conservative_variable_names": list(equations),
                **{f"provenance/{key}": value for key, value in provenance.items()}}
    for name in ("impurity_radiation", "neutral_wall_sources_in_elements", "neutral_perpendicular_diffusion"):
        if name in overrides:
            expected[f"simulation_parameters/switches/{name}"] = overrides[name]
    for name in ("neutral_flux_limiter_mode", "neutral_flux_limiter_tn_source"):
        if name in overrides:
            expected[f"simulation_parameters/physics/{name}"] = overrides[name]
    if "neutralp_lambda" in overrides:
        expected["simulation_parameters/numerics/NeutralP_lambda"] = overrides["neutralp_lambda"]
    if "compute_from_flux" in overrides:
        expected["magnetic/jtor_source"] = "bicubic_psi" if overrides["compute_from_flux"] else "stored_hdf5"
    if overrides.get("impurity_radiation") and (directory / "inputs/impurity_model.nml").is_file():
        impurity = read_selected_input_values(directory / "inputs/impurity_model.nml", {"impurity_names", "impurity_concentrations"})
        expected.update({f"simulation_parameters/physics/{name}": value for name, value in impurity.items()})
    failures = []
    try:
        with h5py.File(path) as handle:
            if overrides.get("transport_1d") and "transport_1d" not in handle:
                failures.append("missing expected output: transport_1d")
            if overrides.get("neutral_flux_limiter_tn_source") == "fixed" and "neutral_flux_limiter_tn_ev" in overrides:
                scale = float(np.asarray(handle["simulation_parameters/adimensionalization/temperature_scale"]).item())
                expected["simulation_parameters/physics/neutral_flux_limiter_tn"] = overrides["neutral_flux_limiter_tn_ev"] / scale
            for name, wanted in expected.items():
                if name not in handle:
                    failures.append(f"missing expected output: {name}")
                    continue
                actual = [value.decode().strip(" \x00") if isinstance(value, bytes) else value.item() if isinstance(value, np.generic) else value
                          for value in np.asarray(handle[name]).reshape(-1)]
                wanted = wanted if isinstance(wanted, list) else [wanted]
                matches = actual == wanted if any(isinstance(value, str) for value in wanted) else (
                    len(actual) == len(wanted) and np.allclose(actual, wanted, rtol=1e-12, atol=0))
                if not matches:
                    failures.append(f"{name}: expected {wanted}, got {actual}")
    except (OSError, KeyError, ValueError, TypeError, ZeroDivisionError) as exc:
        failures.append(f"cannot verify output contract: {exc}")
    return {"status": "failed" if failures else "passed", "failures": failures}
