"""Shared relative L2 and normalized Linf calculations."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ErrorNorms:
    compatible: bool
    finite: bool
    relative_l2: float | None
    normalized_linf: float | None

    @property
    def available(self) -> bool:
        return (
            self.compatible
            and self.finite
            and self.relative_l2 is not None
            and self.normalized_linf is not None
        )


def calculate_error_norms(
    reference: np.ndarray, candidate: np.ndarray
) -> ErrorNorms:
    """Calculate error norms when both arrays are compatible and finite."""
    first = np.asarray(reference, dtype=float)
    second = np.asarray(candidate, dtype=float)
    compatible = first.shape == second.shape
    finite = bool(np.isfinite(first).all() and np.isfinite(second).all())
    if not compatible or not finite:
        return ErrorNorms(compatible, finite, None, None)
    if first.size == 0:
        return ErrorNorms(True, True, None, None)

    difference = second - first
    tiny = np.finfo(float).tiny
    relative_l2 = float(
        np.linalg.norm(difference.ravel())
        / max(float(np.linalg.norm(first.ravel())), tiny)
    )
    normalized_linf = float(
        np.max(np.abs(difference)) / max(float(np.max(np.abs(first))), tiny)
    )
    return ErrorNorms(True, True, relative_l2, normalized_linf)


def maximum_metric(
    metrics: Iterable[Mapping[str, float | None]], key: str
) -> float | None:
    """Return the largest available metric value."""
    values = (metric.get(key) for metric in metrics)
    return max((value for value in values if value is not None), default=None)


def format_metric(value: float | None) -> str:
    """Format an available metric or a readable missing-value marker."""
    return f"{value:.3e}" if value is not None else "n/a"
