"""Check the production Bohm blend and face assembly against independent residuals.

Requires gfortran and NumPy; no solver build or external case data is needed.
"""
import argparse
from contextlib import contextmanager
import ctypes
from dataclasses import dataclass
from itertools import product
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MACH_WIDTH = .20
NORMAL = np.array([.6, .8])


@dataclass(frozen=True)
class Model:
    isothermal: bool = False
    neutral_gamma: bool = False
    optimized: bool = False

    @property
    def neq(self):
        return 6 if self.neutral_gamma else 5

    @property
    def flags(self):
        optimization = ["-Ofast", "-march=native"] if self.optimized else ["-O0", "-g", "-fcheck=all"]
        flags = ["-cpp", "-DNGAMMA", "-fdefault-real-8", "-fdefault-double-8",
                 "-ffree-line-length-none", *optimization]
        if not self.isothermal:
            flags.append("-DTEMPERATURE")
        if self.neutral_gamma:
            flags += ["-DNEUTRAL", "-DNEUTRALGAMMA"]
        return flags


@dataclass(frozen=True)
class Fixture:
    label: str
    trace: np.ndarray
    magnetic: np.ndarray
    angle: float
    energy_threshold: float
    derivatives: bool = True


@dataclass
class Results:
    states: int = 0
    columns: int = 0
    residual: float = 0.
    derivative: float = 0.
    target: float = 0.


def routine(source, name):
    pattern = rf"^\s*SUBROUTINE {name}\(.*?^\s*END\s*SUBROUTINE {name}\b[^\n]*"
    return re.search(pattern, source, re.I | re.M | re.S)[0]


def conditional_block(source):
    depth, result = 0, []
    for line in source.splitlines():
        code = line.split("!", 1)[0].strip()
        if re.match(r"IF\s*\(.*\)\s*THEN\b", code, re.I):
            depth += 1
        result.append(line)
        if re.match(r"END\s*IF\b", code, re.I):
            depth -= 1
            if depth == 0:
                return "\n".join(result)
    raise ValueError("Unclosed production assembly block")


def production_wrapper(model):
    boundary_path = ROOT / "src/HDG/hdg_BC.f90"
    boundary = subprocess.check_output(
        ["gfortran", *model.flags, "-E", "-P", str(boundary_path)], text=True)
    assembly = routine(boundary, "assembly_bohm_bc")
    physics = (ROOT / "src/Models/NGammaTiTe/physics.f90").read_text()
    helper_names = ["softplus", "softplus_deriv", "compute_cs", "compute_dcs_du",
                    "compute_W2", "compute_dW2_dU", "compute_W3", "compute_dW3_dU",
                    "compute_W4", "compute_dW4_dU"]
    diagnostics = ROOT / "src/Utils/Diagnostics"
    flux = re.search(
        r"  MODULE FUNCTION perpendicular_diffusive_flux\(.*?  END FUNCTION perpendicular_diffusive_flux",
        (diagnostics / "balance_diagnostics_common.f90").read_text(), re.S)[0]
    blocks = {
        "HELPERS": "\n".join(routine(physics, name) for name in helper_names),
        "SELECTOR": "!" + routine(boundary, "set_Bohm_bc").split(
            "! Bohm selector:", 1)[1].split("! Shape functions", 1)[0],
        "STABILIZATION": assembly.split("! Begin diagonal Bohm stabilization", 1)[1].split("! End diagonal Bohm stabilization", 1)[0],
        "DIFFUSION": conditional_block(assembly.split("! Perpendicular diffusion", 1)[1]),
        "DIAGNOSTIC_TYPE": re.search(
            r"  TYPE, PUBLIC :: balance_plasma_bc_type.*?  END TYPE balance_plasma_bc_type",
            (diagnostics / "balance_diagnostics.f90").read_text(), re.S)[0],
        "DIAGNOSTIC_HELPERS": flux.replace("MODULE FUNCTION", "FUNCTION", 1) + "\n" + routine(
            (diagnostics / "balance_diagnostics_bc.f90").read_text(), "accumulate_momentum_closure"),
    }
    if not model.isothermal:
        expression = re.search(r"momentum_target=(.*?), &", boundary_path.read_text(), re.S)[1]
        blocks["DIAGNOSTIC_TARGET"] = "target=" + expression.replace("&", "").replace("\n", "")
    wrapper = (ROOT / "unit_tests/fortran/test_bohm_momentum.f90").read_text()
    for name, block in blocks.items():
        wrapper = wrapper.replace("! PRODUCTION_" + name, block)
    return wrapper


class BohmProbe:
    def __init__(self, library, model):
        self.model = model
        self.library = library
        neq = model.neq
        self.gradient = np.asfortranarray(np.linspace(-.25, .4, 2*neq).reshape(2, neq))
        self.iso = np.diag(np.linspace(.3, .8, neq))
        self.ani = np.diag(np.linspace(.06, .12, neq))
        array = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
        library.probe.argtypes = [array]*8 + [ctypes.c_double]*3 + [ctypes.c_int] + [array]*7
        library.probe.restype = None

    def evaluate(self, fixture, element, tau, tangent):
        neq = self.model.neq
        ju, jh = [np.zeros((neq, neq), order="F") for _ in range(2)]
        jq = np.zeros((neq, 2*neq), order="F")
        rhs, weights, target, diagnostic = np.zeros(neq), np.zeros(3), np.zeros(1), np.zeros(3)
        arrays = (element, fixture.trace, self.gradient, self.iso, self.ani,
                  fixture.magnetic, NORMAL, tau)
        lag_velocity = fixture.trace[1]/fixture.trace[0]
        self.library.probe(*[np.asfortranarray(a) for a in arrays], fixture.angle,
                           fixture.energy_threshold, lag_velocity, tangent,
                           ju, jh, jq, rhs, weights, target, diagnostic)
        return ju, jh, jq, rhs, weights, target[0], diagnostic


@contextmanager
def compiled_probe(model):
    with tempfile.TemporaryDirectory(prefix="bohm_momentum_") as directory:
        tmp = Path(directory)
        (tmp / "extracted.f90").write_text(production_wrapper(model))
        sources = [ROOT / "src/Definitions/prec_const.f90",
                   ROOT / "src/Utils/smooth_transition.f90",
                   ROOT / "src/Models/bohm_boundary.f90", tmp / "extracted.f90"]
        subprocess.run(["gfortran", *model.flags, "-shared", "-fPIC", "-J", str(tmp),
                        "-I", str(tmp), *map(str, sources), "-o", str(tmp / "probe.so")],
                       cwd=tmp, check=True)
        yield BohmProbe(ctypes.CDLL(str(tmp / "probe.so")), model)


def smoothstep(z):
    t = np.clip(z, 0., 1.)
    return t*t*(3.-2.*t)


def reference_target(model, fixture, element, trace):
    bn = fixture.magnetic[0]*NORMAL[0]+fixture.magnetic[1]*NORMAL[1]
    sign = 1. if bn >= 0 else -1.
    beta = max(trace[2], 1.e-20)/(max(trace[2], 1.e-20) + fixture.energy_threshold)
    if model.isothermal:
        original = max(fixture.trace[2], 1.e-20)
        beta = original/(original + fixture.energy_threshold)  # Existing frozen suppression.
        active = float(abs(trace[1]/trace[0]) <= 1./beta or abs(bn) <= fixture.angle)
        angle = 1.
        if abs(bn) <= fixture.angle:
            angle = abs(bn)/fixture.angle if fixture.angle > 0 else 0.
        weights = np.array([0., active, 1.-active])
        target = active*sign*beta*angle*element[0] + (1.-active)*element[1]
        return weights, target
    cs = np.sqrt(2./3./trace[0]*(trace[2]+trace[3]-.5*trace[1]**2/trace[0]))
    angle = smoothstep(abs(bn)/fixture.angle) if fixture.angle > 0 else float(bn != 0)
    active = 1.-smoothstep((sign*trace[1]/(trace[0]*cs)-1./beta)/MACH_WIDTH)
    weights = np.array([active*(1.-beta*angle), active*beta*angle, 1.-active])
    target = weights[1]*sign*cs*trace[0] + weights[2]*element[1]
    return weights, target


def reference_residual(probe, fixture, element, trace, gradient, tau, tangent):
    model = probe.model
    weights, momentum = reference_target(model, fixture, element, trace)
    target = element.copy()
    target[1] = momentum
    if model.neutral_gamma:
        target[5] = -(weights[0]+weights[1])*(fixture.trace[1]/fixture.trace[0])*element[4]
    flux = perpendicular_flux(probe, fixture, trace, gradient)
    if tangent:
        flux = NORMAL @ gradient
    if model.neutral_gamma:
        flux[4] = 0.  # Common assembly skips neutral density.
    return tau @ (target-trace) - flux


def perpendicular_flux(probe, fixture, trace, gradient):
    b = fixture.magnetic
    bn = b @ NORMAL
    perpendicular = NORMAL-bn*b
    diffusion = probe.iso.diagonal()
    flux = diffusion*(NORMAL @ gradient)-probe.ani.diagonal()*bn*(b @ gradient)
    flux[1] += (diffusion[0]-diffusion[1])*trace[1]/trace[0]*(perpendicular @ gradient[:, 0])
    if probe.model.isothermal:
        return flux
    rn, ru, ri, re = diffusion[:4]
    w3 = np.zeros(probe.model.neq)
    w3[0] = (rn-ri)*trace[2]/trace[0]+(ri-ru)*(trace[1]/trace[0])**2
    w3[1] = -(ri-ru)*trace[1]/trace[0]
    flux[2] += perpendicular @ gradient @ w3
    flux[3] += (rn-re)*trace[3]/trace[0]*(perpendicular @ gradient[:, 0])
    return flux


def mach_state(neq, sign, mach, density=.8, cs=1.3):
    trace = np.full(neq, .1)
    trace[:4] = [density, sign*density*cs*mach,
                 density*cs**2*(1.1+.5*mach**2), .4*density*cs**2]
    return trace


def transition_fixtures(model):
    for sign, angle_fraction, beta in product([-1., 1.], [0., .25, .5, .75, 1.], [1., .5, .01]):
        offsets = [-.2, .25, .75, 1.2] if beta >= .5 else [-.2]
        for offset in offsets:
            mach = 1./beta + MACH_WIDTH*offset if beta >= .5 else .4
            trace = mach_state(model.neq, sign, mach)
            label = f"sign={sign:g}, angle={angle_fraction:g}, beta={beta:g}, band={offset:g}"
            yield Fixture(label, trace, sign*angle_fraction*.1*NORMAL, .1, trace[2]*(1.-beta)/beta)


def special_fixtures(model):
    tail = [.1]*(model.neq-4)
    if model.isothermal:
        trace = np.array([.8, .2, 1.6, 1., *tail])
        yield Fixture('zero angular threshold at exact tangency', trace, np.zeros(2), 0., .2)
    for sign in [-1., 1.]:
        trace = np.array([.8, -sign*.3, 2., 1., *tail])
        yield Fixture('inward flow without angular threshold', trace, sign*.03*NORMAL, 0., .2)
        for beta in [1., .5]:
            trace = mach_state(model.neq, sign, -.035)
            yield Fixture('grazing inward flow', trace, sign*.01865*NORMAL, .1, trace[2]*(1.-beta)/beta)
    yield Fixture('ion-energy floor', np.array([.8, 0., 1.e-22, 1., *tail]), .05*NORMAL, .1, 1.)
    yield Fixture('low plasma density', np.array([.8e-8, .2e-8, 1.6e-8, 1.e-8, *tail]), .05*NORMAL, .1, 1., False)
    for mach in [1., 1.+MACH_WIDTH]:
        trace = mach_state(model.neq, 1., mach, density=1., cs=1.)
        yield Fixture('Mach support endpoint', trace, .1*NORMAL, .1, 0., False)


def check_case(probe, fixture, stabilization, tangent, results):
    neq = probe.model.neq
    trace = fixture.trace
    element = trace*np.linspace(.93, 1.07, neq)
    tau = np.diag(np.linspace(.12, .22, neq))
    if stabilization == 1:
        tau = np.diag(.11*np.arange(1, neq+1))
    ju, jh, jq, rhs, weights, target, diagnostic = probe.evaluate(fixture, element, tau, tangent)
    exact = reference_residual(probe, fixture, element, trace, probe.gradient, tau, tangent)
    expected_weights, expected_target = reference_target(probe.model, fixture, element, trace)
    reconstructed = ju @ element + jh @ trace + jq @ probe.gradient.ravel(order="F") - rhs
    results.residual = max(results.residual, max(abs(reconstructed-exact)))
    results.target = max(results.target, abs(target-expected_target))
    np.testing.assert_allclose(reconstructed, exact, rtol=2.e-11, atol=2.e-12)
    np.testing.assert_allclose(target, expected_target, rtol=2.e-11, atol=2.e-12)
    np.testing.assert_allclose(weights, expected_weights, atol=3.e-13, rtol=1.e-12)
    assert weights.min() >= 0 and weights.max() <= 1.+2.e-16
    assert abs(sum(weights)-1.) < 5.e-16
    if not tangent:
        np.testing.assert_allclose(sum(diagnostic), -exact[1], rtol=1.e-11, atol=2.e-12)
    results.states += 1
    derivatives = fixture.derivatives
    if probe.model.isothermal:
        derivatives = derivatives and abs(abs(trace[1]/trace[0])-1.) > .1
    if derivatives:
        check_jacobian(probe, fixture, element, tau, tangent, np.column_stack((ju, jh, jq)), results)


def check_jacobian(probe, fixture, element, tau, tangent, jacobian, results):
    neq = probe.model.neq
    x = np.r_[element, fixture.trace, probe.gradient.ravel(order="F")]

    def value(state):
        return reference_residual(probe, fixture, state[:neq], state[neq:2*neq],
                                  state[2*neq:].reshape((2, neq), order="F"), tau, tangent)

    for column in range(len(x)):
        # A sub-floor energy change is unresolvable in this whole-row difference.
        if column == neq+2 and fixture.trace[2] < 1.e-20:
            continue
        step = 1.e-7*max(abs(x[column]), .1)
        plus, minus = x.copy(), x.copy()
        plus[column] += step
        minus[column] -= step
        numeric = (value(plus)-value(minus))/(2.*step)
        results.derivative = max(results.derivative, np.max(abs(numeric-jacobian[:, column])))
        np.testing.assert_allclose(numeric, jacobian[:, column], rtol=2.e-6, atol=2.e-6,
                                   err_msg=f"column={column}, tangent={tangent}")
        results.columns += 1


def check_mach_endpoints(probe):
    for mach in [1., 1.+MACH_WIDTH]:
        trace = mach_state(probe.model.neq, 1., mach, density=1., cs=1.)
        fixture = Fixture('Mach endpoint tangent', trace, .2*NORMAL, .1, 0.)
        element = .97*trace
        tau = np.diag(np.linspace(.12, .22, probe.model.neq))
        analytic = probe.evaluate(fixture, element, tau, 0)[1][:, 1]
        errors = []
        # Crossing a clipped cubic's C1 endpoint gives first-order FD convergence.
        for step in [1.e-4, 1.e-5, 1.e-6]:
            plus, minus = trace.copy(), trace.copy()
            plus[1] += step
            minus[1] -= step
            fp = reference_residual(probe, fixture, element, plus, probe.gradient, tau, 0)
            fm = reference_residual(probe, fixture, element, minus, probe.gradient, tau, 0)
            errors.append(np.max(abs((fp-fm)/(2.*step)-analytic)))
        assert errors[1] < .2*errors[0] and errors[2] < .2*errors[1], errors


def check_prescribed_stabilization(probe, results):
    trace = mach_state(probe.model.neq, 1., 2.+MACH_WIDTH*.5)
    fixture = Fixture('prescribed per-equation tau', trace, .025*NORMAL, .1,
                      trace[2], derivatives=False)
    check_case(probe, fixture, 1, 0, results)


def check_model(model):
    results = Results()
    fixtures = [*transition_fixtures(model), *special_fixtures(model)]
    with compiled_probe(model) as probe:
        for fixture, tangent in product(fixtures, [0, 1]):
            try:
                check_case(probe, fixture, 5, tangent, results)
            except AssertionError as error:
                raise AssertionError(f"{fixture.label}; stab=5, tangent={tangent}") from error
        check_prescribed_stabilization(probe, results)
        if not model.isothermal:
            check_mach_endpoints(probe)
    configuration = "optimized" if model.optimized else "bounds-checked"
    print(f"Bohm PASS ({configuration}): {results.states} states, {results.columns} derivative columns; "
          f"max FD={results.derivative:.3g}, residual={results.residual:.3g}, "
          f"target={results.target:.3g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--isothermal", action="store_true")
    modes.add_argument("--neutral-gamma", action="store_true")
    args = parser.parse_args()
    check_model(Model(args.isothermal, args.neutral_gamma))
    if not args.isothermal:
        check_model(Model(args.isothermal, args.neutral_gamma, optimized=True))

if __name__ == "__main__":
    main()
