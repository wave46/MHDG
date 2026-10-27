"""Check the production input loader without running a solver simulation."""

from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main():
    executable = Path(sys.argv[1]).resolve()
    build_kind = subprocess.check_output(
        [str(executable), "--build-kind"], text=True, timeout=15).strip()
    if build_kind not in ("k", "non-k"):
        raise SystemExit(f"unexpected input-check build kind: {build_kind!r}")
    k_equation = build_kind == "k"
    repository = Path(__file__).resolve().parents[1]
    template = (repository / "test/param_initial.txt").read_text()
    # Omission and inheritance cases must not inherit explicit template albedos.
    template = re.sub(r"(?mi)^\s*Re_n(?:_pump)?\s*=[^\n]*\n", "", template)
    template = re.sub(r"(?mi)^\s*recycle_plasma_diffusion\s*=[^\n]*\n", "", template)
    if k_equation:
        template = template.replace("&PHYS_LST", "&PHYS_LST\n"
            "  diff_k_min=10, diff_k_max=10, k_max=1.e10", 1)
    valid = [
        ("omitted", "", (1.0, 1.0)),
        ("diffusion_false", "recycle_plasma_diffusion=.false.", (1.0, 1.0)),
        ("diffusion_true", "recycle_plasma_diffusion=.true.", (1.0, 1.0)),
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
    valid.append(("bohm_unsuppressed", "bohm_energy_thresh=0, bohmth=0", (1.0, 1.0)))
    invalid.extend((f"bohm_energy_{name}", f"bohm_energy_thresh={value}",
                    "bohm_energy_thresh must be finite and nonnegative")
                   for name, value in (("negative", "-1"), ("nan", "NaN"), ("infinite", "Infinity")))
    invalid.extend((f"bohm_angle_{name}", f"bohmth={value}", "bohmth must be finite")
                   for name, value in (("nan", "NaN"), ("infinite", "Infinity")))
    cases = [("PHYS_LST", *case) for case in valid + invalid]
    cases.extend([
        ("NUMER_LST", "stabilization_prescribed", "stab=1, tau=0.11,0.22,0.33,0.44,0.55", (1.0, 1.0)),
        ("NUMER_LST", "stabilization_automatic", "stab=5", (1.0, 1.0)),
    ])
    cases.append(("NUMER_LST", "legacy_bohm_zero", "bohmtypebc=0", (1.0, 1.0)))
    cases.extend(("NUMER_LST", f"removed_bohm_{value}", f"bohmtypebc={value}",
                  "Alternate Bohm boundary conditions have been removed")
                 for value in (1, 2))
    cases.extend([
        ("SWITCH_LST", "legacy_savetau_false", "saveTau=.false.", (1.0, 1.0)),
        ("SWITCH_LST", "removed_savetau_true", "saveTau=.true.",
         "saveTau output has been removed"),
        ("SWITCH_LST", "legacy_limrho_zero", "limrho=0", (1.0, 1.0)),
        ("NUMER_LST", "legacy_limrho_coefficients", "minrho=1.e-12, so_coe=0.01, df_coe=100", (1.0, 1.0)),
    ])
    cases.extend(("SWITCH_LST", f"removed_limrho_{value}", f"limrho={value}",
                  "limrho has been removed")
                 for value in (-1, 1, 2, 3))
    cases.extend(("SWITCH_LST", f"retired_iter_{value}", f"testcase={value}",
                  "ITER testcases 80-89 have been removed")
                 for value in range(80, 90))
    cases.extend(("SWITCH_LST", f"supported_{value}", f"testcase={value}", (1.0, 1.0))
                 for value in (52, 60))
    cases.extend([
        ("SWITCH_LST", "legacy_diffmin", "diffmin=0.1", (1.0, 1.0)),
        ("PHYS_LST", "legacy_core_sources", "part_source=1, ener_source=20e6", (1.0, 1.0)),
    ])
    cases.extend(("SWITCH_LST", f"k_range_{value}", f"testcase={value}",
                  "k-equation builds require testcase in 50-69"
                  if k_equation and not 50 <= value <= 69 else (1.0, 1.0))
                 for value in (2, 49, 50, 69, 70))
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
                    flag = re.search(r"^RECYCLE_PLASMA_DIFFUSION\s+([TF])", output, re.M)
                    passed = passed and flag is not None and (flag[1] == "T") == (name == "diffusion_true")
            else:
                passed = (run.returncode != 0 and expected in output
                          and "NEUTRAL_WALL_ALBEDOS" not in output)
            if passed and name == "stabilization_prescribed":
                selected = re.search(r"^STABILIZATION_TYPE\s+(\d+)", output, re.M)
                prescribed = re.search(r"^STABILIZATION_TAU\s+(.+)$", output, re.M)
                passed = selected is not None and selected[1] == "1" and prescribed is not None
                if passed:
                    values = tuple(float(value) for value in prescribed[1].split())
                    passed = all(abs(actual-expected) < 1.e-14
                                 for actual, expected in zip(values, (.11,.22,.33,.44,.55)))
            if not passed:
                raise SystemExit(f"FAIL: {name}\n{output}")
    print(f"neutral wall input checks: {len(cases)} PASS")


if __name__ == "__main__":
    main()
