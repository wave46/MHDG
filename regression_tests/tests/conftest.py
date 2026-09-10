"""Import the package and remaining legacy helpers during the migration."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT.parent), str(ROOT), str(ROOT / "tools")]
