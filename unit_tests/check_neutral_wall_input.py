"""Check the production input loader without running a solver simulation."""

from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main():
    executable = Path(sys.argv[1]).resolve()
    repository = Path(__file__).resolve().parents[1]
    template = (repository / "test/param_initial.txt").read_text()
    valid = [
        ("omitted", "", (1.0, 1.0)),
        ("explicit_inactive", "Re_n=1, Re_n_pump=1", (1.0, 1.0)),
        ("inherit", "Re_n=0.99", (0.99, 0.99)),
        ("pump_override", "Re_n=0.99, Re_n_pump=0.95", (0.99, 0.95)),
        ("pump_only", "Re_n_pump=0.99", (1.0, 0.99)),
        ("wall_zero", "Re_n=0, Re_n_pump=1", (0.0, 1.0)),
        ("pump_zero", "Re_n=1, Re_n_pump=0", (1.0, 0.0)),
    ]
    invalid = [
        (f"{key}_{label}", f"{key}={value}",
         f"{key} must be finite and in [0,1]")
        for key in ("Re_n", "Re_n_pump")
        for label, value in (("negative", "-0.01"), ("too_large", "1.01"),
                             ("nan", "NaN"), ("infinity", "Infinity"))
    ]
    with tempfile.TemporaryDirectory(prefix="mhdg-neutral-wall-input-") as scratch:
        for name, entries, expected in valid + invalid:
            directory = Path(scratch) / name
            directory.mkdir()
            text = template.replace("&PHYS_LST\n", f"&PHYS_LST\n  {entries}\n")
            (directory / "param.txt").write_text(text)
            run = subprocess.run([str(executable)], cwd=directory,
                                 capture_output=True, text=True, timeout=15)
            output = run.stdout + run.stderr
            if isinstance(expected, tuple):
                match = re.search(r"^NEUTRAL_WALL_ALBEDOS\s+(\S+)\s+(\S+)",
                                  output, re.M)
                passed = run.returncode == 0 and match is not None
                if passed:
                    values = tuple(float(value) for value in match.groups())
                    passed = values == expected
            else:
                passed = (run.returncode != 0 and expected in output
                          and "NEUTRAL_WALL_ALBEDOS" not in output)
            if not passed:
                raise SystemExit(f"FAIL: {name}\n{output}")
    print(f"neutral wall input checks: {len(valid) + len(invalid)} PASS")


if __name__ == "__main__":
    main()
