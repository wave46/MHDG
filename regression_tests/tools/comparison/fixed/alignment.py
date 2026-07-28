"""Align MHDG meshes that differ only in global numbering."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any

import h5py
import numpy as np

from comparison.fixed.data import mesh_rows, required_array, required_scalar
from support.errors import ComparisonError


@dataclass(frozen=True)
class MeshAlignment:
    """Candidate element and face indices expressed in reference ordering."""

    candidate_element_to_reference: np.ndarray
    candidate_face_to_reference: np.ndarray
    reverse_candidate_faces: np.ndarray


@dataclass(frozen=True)
class _Mesh:
    coordinates: np.ndarray
    element_nodes: np.ndarray
    element_vertices: np.ndarray
    boundary_face_nodes: np.ndarray
    boundary_flags: np.ndarray
    interior_faces: np.ndarray
    nodes_per_face: int


def compare_numbering_invariant(
    reference: h5py.File,
    candidate: h5py.File,
    coordinate_tolerance: float,
) -> tuple[dict[str, Any], MeshAlignment | None]:
    """Require the same high-order mesh after global-tag renumbering."""
    reference_mesh = _read_mesh(reference)
    candidate_mesh = _read_mesh(candidate)

    try:
        candidate_node_to_reference, maximum_error = _match_nodes(
            reference_mesh.coordinates,
            candidate_mesh.coordinates,
            coordinate_tolerance,
        )
    except ComparisonError as exc:
        return _report(coordinate_tolerance, None, str(exc)), None

    try:
        candidate_element_to_reference = _match_elements(
            reference_mesh, candidate_mesh, candidate_node_to_reference
        )
        _check_boundary_faces(
            reference_mesh, candidate_mesh, candidate_node_to_reference
        )
        candidate_face_to_reference, reverse_candidate_faces = _match_faces(
            reference_mesh, candidate_mesh, candidate_node_to_reference
        )
    except ComparisonError as exc:
        return _report(coordinate_tolerance, maximum_error, str(exc)), None

    report = _report(
        coordinate_tolerance,
        maximum_error,
        None,
        element_count=candidate_element_to_reference.size,
        face_count=candidate_face_to_reference.size,
        reversed_face_count=int(np.count_nonzero(reverse_candidate_faces)),
    )
    alignment = MeshAlignment(
        candidate_element_to_reference,
        candidate_face_to_reference,
        reverse_candidate_faces,
    )
    return report, alignment


def _read_mesh(handle: h5py.File) -> _Mesh:
    return _Mesh(
        coordinates=mesh_rows(handle, "X", float),
        element_nodes=mesh_rows(handle, "T", np.int64),
        element_vertices=mesh_rows(handle, "Tlin", np.int64),
        boundary_face_nodes=mesh_rows(handle, "Tb", np.int64),
        boundary_flags=required_array(
            handle, "mesh", "boundaryFlag"
        ).reshape(-1),
        interior_faces=mesh_rows(handle, "intfaces", np.int64),
        nodes_per_face=int(required_scalar(handle, "mesh", "Nnodesperface")),
    )


def _match_nodes(
    reference: np.ndarray,
    candidate: np.ndarray,
    tolerance: float,
) -> tuple[np.ndarray, float]:
    if reference.shape[0] != candidate.shape[0]:
        raise ComparisonError("node counts differ")
    if not np.isfinite(reference).all() or not np.isfinite(candidate).all():
        raise ComparisonError("mesh coordinates are not finite")

    reference_cells: dict[tuple[int, ...], list[int]] = {}
    for index, point in enumerate(reference):
        reference_cells.setdefault(_cell(point, tolerance), []).append(index)

    candidate_node_to_reference = np.empty(candidate.shape[0], dtype=np.int64)
    reference_used = np.zeros(reference.shape[0], dtype=bool)
    maximum_error = 0.0
    neighbor_offsets = tuple(product((-1, 0, 1), repeat=reference.shape[1]))

    for candidate_index, point in enumerate(candidate):
        cell = _cell(point, tolerance)
        nearby = (
            reference_index
            for offset in neighbor_offsets
            for reference_index in reference_cells.get(
                tuple(value + shift for value, shift in zip(cell, offset)), ()
            )
        )
        matches = []
        for reference_index in nearby:
            error = float(np.max(np.abs(point - reference[reference_index])))
            if error <= tolerance:
                matches.append((reference_index, error))
        if len(matches) != 1:
            raise ComparisonError(
                f"candidate node {candidate_index + 1} has {len(matches)} matches"
            )

        reference_index, error = matches[0]
        if reference_used[reference_index]:
            raise ComparisonError(
                f"reference node {reference_index + 1} is matched twice"
            )
        candidate_node_to_reference[candidate_index] = reference_index
        reference_used[reference_index] = True
        maximum_error = max(maximum_error, error)

    return candidate_node_to_reference, maximum_error


def _cell(point: np.ndarray, width: float) -> tuple[int, ...]:
    return tuple(np.floor(point / width).astype(np.int64))


def _match_elements(
    reference: _Mesh,
    candidate: _Mesh,
    candidate_node_to_reference: np.ndarray,
) -> np.ndarray:
    candidate_element_to_reference = _match_ordered_entities(
        reference.element_nodes,
        candidate.element_nodes,
        candidate_node_to_reference,
        "elements",
    )
    candidate_vertices = _use_reference_node_ids(
        candidate.element_vertices, candidate_node_to_reference
    )
    if not np.array_equal(
        candidate_vertices,
        reference.element_vertices[candidate_element_to_reference],
    ):
        raise ComparisonError("element vertex ordering differs")
    return candidate_element_to_reference


def _check_boundary_faces(
    reference: _Mesh,
    candidate: _Mesh,
    candidate_node_to_reference: np.ndarray,
) -> None:
    candidate_boundary_to_reference = _match_ordered_entities(
        reference.boundary_face_nodes,
        candidate.boundary_face_nodes,
        candidate_node_to_reference,
        "boundary faces",
    )
    if not np.array_equal(
        candidate.boundary_flags,
        reference.boundary_flags[candidate_boundary_to_reference],
    ):
        raise ComparisonError("boundary flags differ")


def _match_faces(
    reference: _Mesh,
    candidate: _Mesh,
    candidate_node_to_reference: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    reference_nodes = _ordered_face_nodes(reference)
    candidate_nodes = _use_reference_node_ids(
        _ordered_face_nodes(candidate), candidate_node_to_reference
    )
    candidate_face_to_reference = _match_rows_by_node_set(
        reference_nodes, candidate_nodes, "faces"
    )
    matched_reference = reference_nodes[candidate_face_to_reference]

    same_direction = np.all(candidate_nodes == matched_reference, axis=1)
    reverse_direction = np.all(
        candidate_nodes[:, ::-1] == matched_reference, axis=1
    )
    if not np.all(same_direction | reverse_direction):
        raise ComparisonError("face node ordering is neither equal nor reversed")
    return candidate_face_to_reference, ~same_direction


def _ordered_face_nodes(mesh: _Mesh) -> np.ndarray:
    interior_count = mesh.interior_faces.shape[0]
    ordered = np.empty(
        (
            interior_count + mesh.boundary_face_nodes.shape[0],
            mesh.nodes_per_face,
        ),
        dtype=np.int64,
    )
    local_faces = _triangle_face_nodes(mesh.nodes_per_face)
    for face_index, face in enumerate(mesh.interior_faces):
        element = face[0] - 1
        local_face = face[1] - 1
        ordered[face_index] = mesh.element_nodes[
            element, local_faces[local_face]
        ]
    ordered[interior_count:] = mesh.boundary_face_nodes
    return ordered


def _triangle_face_nodes(nodes_per_face: int) -> tuple[np.ndarray, ...]:
    edge_node_count = nodes_per_face - 2
    faces = []
    for face in range(3):
        first_edge_node = 3 + face * edge_node_count
        edge_nodes = range(first_edge_node, first_edge_node + edge_node_count)
        faces.append(
            np.asarray([face, *edge_nodes, (face + 1) % 3], dtype=np.int64)
        )
    return tuple(faces)


def _use_reference_node_ids(
    candidate_nodes: np.ndarray,
    candidate_node_to_reference: np.ndarray,
) -> np.ndarray:
    return candidate_node_to_reference[candidate_nodes - 1] + 1


def _match_ordered_entities(
    reference_nodes: np.ndarray,
    candidate_nodes: np.ndarray,
    candidate_node_to_reference: np.ndarray,
    label: str,
) -> np.ndarray:
    candidate_nodes = _use_reference_node_ids(
        candidate_nodes, candidate_node_to_reference
    )
    candidate_to_reference = _match_rows_by_node_set(
        reference_nodes, candidate_nodes, label
    )
    if not np.array_equal(
        candidate_nodes, reference_nodes[candidate_to_reference]
    ):
        raise ComparisonError(f"{label} have different local node ordering")
    return candidate_to_reference


def _match_rows_by_node_set(
    reference: np.ndarray,
    candidate: np.ndarray,
    label: str,
) -> np.ndarray:
    if reference.shape != candidate.shape:
        raise ComparisonError(f"{label} have different shapes")

    reference_rows = {
        tuple(sorted(row)): index for index, row in enumerate(reference)
    }
    if len(reference_rows) != reference.shape[0]:
        raise ComparisonError(f"reference {label} are not unique")

    row_map = np.asarray(
        [reference_rows.get(tuple(sorted(row)), -1) for row in candidate],
        dtype=np.int64,
    )
    if np.any(row_map < 0) or np.unique(row_map).size != row_map.size:
        raise ComparisonError(f"{label} differ after node renumbering")
    return row_map


def _report(
    tolerance: float,
    maximum_error: float | None,
    failure: str | None,
    **counts: int,
) -> dict[str, Any]:
    coordinates_passed = maximum_error is not None
    topology_passed = coordinates_passed and failure is None
    coordinates = {
        "passed": coordinates_passed,
        "maximum_absolute_error": maximum_error,
        "absolute_tolerance": tolerance,
    }
    if not coordinates_passed and failure:
        coordinates["reason"] = failure
    return {
        "mode": "numbering_invariant",
        "passed": topology_passed,
        "coordinates": coordinates,
        "connectivity": {
            "topology": {
                "passed": topology_passed,
                **({"reason": failure} if failure else counts),
            }
        },
    }
