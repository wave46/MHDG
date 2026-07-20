"""Shared identifier rules for regression-harness records."""

from __future__ import annotations

import re


IDENTIFIER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
