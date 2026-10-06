"""Check signs and separation using the production boundary accumulator."""

from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def check(handle, path, expected):
    actual = np.asarray(handle[path][()]).item()
    assert np.isclose(actual, expected, rtol=1e-13, atol=1e-13), (path, actual, expected)


executable = Path(sys.argv[1]).resolve()
cases = [(mode, wall, 0) for mode in (1, 2, 3) for wall in (0., 3.)]
cases += [(3, wall, 1) for wall in (0., 3.)] + [(0, 3., 1)]
for mode, wall, relocated in cases:
    with tempfile.TemporaryDirectory(prefix="neutral_wall_diagnostics_") as directory:
        result = subprocess.run([str(executable), str(mode), str(wall), str(relocated)],
                                cwd=directory, capture_output=True, text=True, check=True)
        with h5py.File(Path(directory) / "diagnostics.h5") as handle:
            if mode == 0:
                assert "diagnostics" not in handle
                assert "Balance diagnostics" not in result.stdout
                continue
            if mode == 1:
                root = "diagnostics/summary"
                check(handle, root + "/external_sources/particles/puff_in", 14.)
                check(handle, root + "/external_sources/particles/pump_out", 12.)
                check(handle, root + "/external_sources/particles/neutral_wall_absorption", 2*wall)
                check(handle, root + "/balances/total_n/physical_imbalance", -3.-2*wall)
            else:
                root = "diagnostics/equations"
                check(handle, root + "/n_n/physical/imbalance", 7.-2*wall)
                check(handle, root + "/total_n/physical/imbalance", -3.-2*wall)
                check(handle, root + "/n_n/bc/residual", 0.)
                check(handle, root + "/n_n/discrete/residual", 2. if relocated else 0.)
                if mode == 3:
                    bc = root + "/n_n/bc/source_components/"
                    physical = root + "/n_n/physical/boundary_components_inward/"
                    check(handle, bc + "neutral_wall_absorption", 2*wall)
                    check(handle, physical + "neutral_wall_absorption", -2*wall)
                    check(handle, bc + "recycling_parallel_inward", 5.)
                    check(handle, bc + "puff", 0. if relocated else 14.)
                    check(handle, bc + "pump", 0. if relocated else 12.)
                    if relocated:
                        volume = root + "/n_n/physical/volume_components/"
                        check(handle, volume + "puff", 14.)
                        check(handle, volume + "pump", -12.)
                    assert "neutral_wall_absorption" not in handle[root + "/n_n/discrete/boundary_components_inward"]
                    assert "neutral_wall_absorption" not in handle[root + "/n_n/physical/volume_components"]
            if mode in (1, 3):
                assert "neutral wall absorption" in result.stdout
print(f"Neutral wall diagnostic checks passed ({len(cases)} cases)")
