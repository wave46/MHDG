"""Run `make check-adaptivity-estimator` from lib/ for a 2D Ti/Te build.

The single-rank MPI driver checks affine reconstruction for 1/2/5/6 equations,
planar/axisymmetric geometry, length-scale invariance of physical-field errors,
the deliberate 10 cm cap, and unchanged solver gradients. It uses the generic
node table and no external case data, mesh generation or plasma solve.
"""
from pathlib import Path
import os
import subprocess
import sys
import tempfile


def main():
    executable = Path(sys.argv[1]).resolve()
    repository = Path(__file__).resolve().parents[1]
    environment = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    with tempfile.TemporaryDirectory(prefix="mhdg-estimator-") as directory:
        scratch = Path(directory)
        (scratch / "positionFeketeNodesTri2D.h5").symlink_to(
            repository / "test/positionFeketeNodesTri2D.h5"
        )
        subprocess.run([str(executable)], cwd=scratch, env=environment,
                       check=True, timeout=60)


if __name__ == "__main__":
    main()
