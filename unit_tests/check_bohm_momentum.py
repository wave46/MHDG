"""Compile extracted production Bohm stabilization/selector without a solver build.

Run: python3 unit_tests/check_bohm_momentum.py [--source path/to/hdg_BC.f90]
The optional source supports checking a pre-fix negative control.
"""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def routine(source, name):
    return re.search(
        rf"^\s*SUBROUTINE {name}\(.*?^\s*END\s*SUBROUTINE {name}\b[^\n]*",
        source, re.I | re.M | re.S,
    )[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "src/HDG/hdg_BC.f90")
    args = parser.parse_args()
    flags = ["-cpp", "-DNGAMMA", "-DTEMPERATURE", "-fdefault-real-8",
             "-fdefault-double-8", "-ffree-line-length-none", "-O0", "-fcheck=all",
             "-ffpe-trap=invalid,zero,overflow"]
    boundary = subprocess.check_output(
        ["gfortran", *flags, "-E", "-P", str(args.source.resolve())], text=True)
    assembly = routine(boundary, "assembly_bohm_bc")
    stabilization = assembly.split("! Case diagonal matrix stabilization", 1)[1]
    stabilization = stabilization.split("! Compute V(U^(k-1))", 1)[0]
    branch = routine(boundary, "set_Bohm_bc").split("! NEW BOHM BC", 1)[1]
    branch = branch.split("! Shape functions", 1)[0]
    physics = (ROOT / "src/Models/NGammaTiTe/physics.f90").read_text()
    helpers = "\n".join(routine(physics, name) for name in
                        ["softplus", "softplus_deriv", "compute_cs", "compute_dcs_du"])
    wrapper = """
module extracted_bohm
use,intrinsic::ieee_arithmetic
implicit none
type physics_type
 real*8::bohmth=.1d0,bohm_energy_thresh=0.d0
end type
type(physics_type)::phys
contains
""" + helpers + """
subroutine assemble(u,uhat,ds,s,delta,tau_input,full,ju,juhat)
real*8,intent(in)::u(5),uhat(5),ds(5),s,tau_input(5,5)
integer,intent(in)::delta
logical,intent(in)::full
real*8,intent(out)::ju(5,5),juhat(5,5)
type matrices
 real*8::All(5,5,1),Alu(5,5,1)
end type
type numer_type
 integer::stab
end type
type(matrices)::elMat
type(numer_type)::numer
real*8::ufg(5),uefg(5),dcs_du(5),setval,tau(5,5),NiNi(1,1)
integer::Neqstab,iel,i,j,ind(1),indi(1),indj(1),ind_asf(1),ind_ff(5),ind_fe(5)
elMat%All=0;elMat%Alu=0;numer%stab=1
if(full) numer%stab=6
ufg=uhat;uefg=u;dcs_du=ds;setval=s;tau=tau_input;NiNi=1
Neqstab=5;iel=1;ind_asf=0;ind_ff=[1,2,3,4,5];ind_fe=ind_ff
""" + stabilization + """
ju=elMat%Alu(:,:,1);juhat=elMat%All(:,:,1)
end subroutine
subroutine select_branch(uhat,bn,delta,setval)
real*8,intent(in)::uhat(5),bn
integer,intent(out)::delta
real*8,intent(out)::setval
real*8::ufg(1,5),upg(1,9),dcs_du(5),SoundSpeed,inc,sn,bohm_suppresion,energy_zero_threshold,bohm_velocity
integer::g
g=1;ufg(1,:)=uhat;upg=0;upg(1,2)=uhat(2)/uhat(1);inc=bn;energy_zero_threshold=1.d-20
call compute_cs(uhat,SoundSpeed)
call compute_dcs_du(uhat,dcs_du)
""" + branch + "\nend subroutine\nend module\n"
    with tempfile.TemporaryDirectory(prefix="bohm_momentum_") as directory:
        tmp = Path(directory)
        source = tmp / "extracted.f90"
        source.write_text(wrapper)
        subprocess.run(
            ["gfortran", *flags, str(source),
             str(ROOT / "unit_tests/fortran/test_bohm_momentum.f90"),
             "-o", str(tmp / "check")], cwd=tmp, check=True)
        subprocess.run([str(tmp / "check")], cwd=tmp, check=True)


if __name__ == "__main__":
    main()
