"""Typed golden-reference records."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from support.errors import BundleError


REFERENCE_MATRIX_ROLE = "reference_matrix"
REFERENCE_MATRIX_ID = "golden_matrix_index"


@dataclass(frozen=True)
class StageReference:
    stage_id: str
    solution: Path


@dataclass(frozen=True)
class MatrixRun:
    workflow_id: str
    layout_id: str
    directory: Path
    stages: tuple[StageReference, ...]


@dataclass(frozen=True)
class ReferenceMatrix:
    bundle_root: Path
    bundle_id: str
    bundle_version: str
    references: dict[tuple[str, str, str], Path]

    def reference_for(
        self,
        workflow_id: str,
        layout_id: str,
        stage_id: str,
    ) -> Path:
        key = (workflow_id, layout_id, stage_id)
        try:
            return self.references[key]
        except KeyError as exc:
            cell = "/".join(key)
            raise BundleError(
                f"golden matrix has no reference for {cell}"
            ) from exc
