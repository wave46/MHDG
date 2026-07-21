"""Typed results produced by bundle operations."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ValidationSummary:
    bundle_id: str
    bundle_version: str
    artifact_count: int
    verified_artifact_count: int
    verified_bytes: int
    case_data: list[str]
    checked_cases: list[str]
    warnings: list[str]
