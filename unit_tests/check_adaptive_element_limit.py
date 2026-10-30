"""Check real candidate caps and failure propagation on a generic square."""
from pathlib import Path
import os
import shlex
import subprocess
import sys
import tempfile


GEOMETRY = """Point(1) = {1, 0, 0, 1};
Point(2) = {1.1, 0, 0, 1};
Point(3) = {1.1, 0.1, 0, 1};
Point(4) = {1, 0.1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
"""


def main():
    executable = Path(sys.argv[1]).resolve()
    repository = Path(__file__).resolve().parents[1]
    environment = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    launchers = [[]]
    if sys.argv[2] == "parall":
        launchers.append([*shlex.split(os.environ.get("MPIEXEC", "mpirun.openmpi")), "-n", "2"])
    with tempfile.TemporaryDirectory(prefix="mhdg-element-limit-") as directory:
        root = Path(directory)
        for index, launcher in enumerate(launchers):
            for mode in ["generate", "reject", "periodic"]:
                scratch = root / f"{index}-{mode}"
                scratch.mkdir()
                (scratch / "res").mkdir()
                (scratch / "square.geo").write_text(GEOMETRY)
                (scratch / "positionFeketeNodesTri2D.h5").symlink_to(
                    repository / "test/positionFeketeNodesTri2D.h5"
                )
                if mode != "generate":
                    (scratch / "res/temp.msh").write_text("preserve previous candidate\n")
                run = subprocess.run([*launcher, str(executable), mode], cwd=scratch,
                                     env=environment, text=True, capture_output=True, timeout=60)
                output = run.stdout + run.stderr
                if mode == "periodic":
                    passed = run.returncode != 0 and "Periodic mesh adaptation failed:" in output
                else:
                    passed = run.returncode == 0 and "ADAPTIVE_ELEMENT_LIMIT_PASS" in output
                if not passed:
                    raise SystemExit(f"FAIL: launcher={launcher} mode={mode}\n{output}")
                if mode == "generate":
                    assert (scratch / "res/cap_mesh_n1.msh").is_file()
                    assert (scratch / "res/cap_mesh_n2.msh").is_file()
                    assert not (scratch / "res/cap_mesh_n3.msh").exists()
                    assert (scratch / "res/temp.msh").read_text() == "preserve previous candidate\n"
                else:
                    assert (scratch / "res/temp.msh").read_text() == "preserve previous candidate\n"
                    assert sorted(path.name for path in (scratch / "res").iterdir()) == ["temp.msh"]
        print("Adaptive element cap: unlimited/exact/oversized candidates, retained mesh/fields, "
              "Newton failure propagation and periodic stop PASS")


if __name__ == "__main__":
    main()
