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
    # Omission and inheritance cases must not inherit explicit template albedos.
    template = re.sub(r"(?mi)^\s*Re_n(?:_pump)?\s*=[^\n]*\n", "", template)
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
    invalid.append(("removed_iter_core", "bcflags(8)=6",
                    "Boundary condition 6 (ITER core) has been removed"))
    cases = [("PHYS_LST", *case) for case in valid + invalid]
    cases.append(("NUMER_LST", "legacy_bohm_zero", "bohmtypebc=0", (1.0, 1.0)))
    cases.extend(("NUMER_LST", f"removed_bohm_{value}", f"bohmtypebc={value}",
                  "Alternate Bohm boundary conditions have been removed")
                 for value in (1, 2))
    cases.extend([
        ("SWITCH_LST", "legacy_savetau_false", "saveTau=.false.", (1.0, 1.0)),
        ("SWITCH_LST", "removed_savetau_true", "saveTau=.true.",
         "saveTau output has been removed"),
    ])
    cases.extend(("SWITCH_LST", f"retired_iter_{value}", f"testcase={value}",
                  "ITER testcases 80-89 have been removed")
                 for value in range(80, 90))
    cases.extend(("SWITCH_LST", f"supported_{value}", f"testcase={value}", (1.0, 1.0))
                 for value in (52, 60))
    cases.extend([
        ("SWITCH_LST", "legacy_diffmin", "diffmin=0.1", (1.0, 1.0)),
        ("PHYS_LST", "legacy_core_sources", "part_source=1, ener_source=20e6", (1.0, 1.0)),
    ])
    with tempfile.TemporaryDirectory(prefix="mhdg-neutral-wall-input-") as scratch:
        for namelist, name, entries, expected in cases:
            directory = Path(scratch) / name
            directory.mkdir()
            text = re.sub(rf"(?ms)(&{namelist}\b.*?)(^\s*/\s*$)",
                          lambda match: f"{match[1]}  {entries}\n{match[2]}",
                          template, count=1)
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
    print(f"neutral wall input checks: {len(cases)} PASS")


if __name__ == "__main__":
    main()
