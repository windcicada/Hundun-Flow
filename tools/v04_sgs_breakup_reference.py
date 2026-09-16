#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Rebuild the fixed FP64 breakup fixture with complete original COAST routines.

Only the surrounding state modules and deterministic random draw are supplied
by this driver. Reference source bytes are hash-checked and remain read-only.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

EXPECTED = {'break_up.F': '6383b48430fc21f2449fc14b567526ffd7b1c2fb3072385048db58e0fbfa1782', 'cdf_spray.F': '0fd1549799064247883b54ae61c8d01c00b1c1e0eb01c18204e638105b500fde'}

MODULES = r"""module exchange
 logical :: master=.true.
 integer :: mydom=0
end module
module global
 real :: pi=3.1415926535897932384626433832795, big=1.e30
 integer :: scrn=6,mout=6
end module
module spray_module
 real :: part(22,1),xpart(1),ypart(1),zpart(1)
 integer :: index_cell(1),brokenp=0,nstations=20,uni_debug=6
 real :: dau_p_x(2),dau_p_y(2),dau_p_z(2),dau_p_u(2),dau_p_v(2),dau_p_w(2)
 real :: dau_p_d(2),dau_p_n(2),dau_p_t(2),dau_p_cell(2),dau_p_props(22,2)
 real :: cbr,csb,station(20),icdf(20),dstar(20),fstar(20)
 logical :: debug=.false.,stochabreak=.true.
end module
module nrtype
 integer,parameter :: sp=kind(1.0)
end module
module nr
 use nrtype
 real(sp) :: draw
 contains
 subroutine ran1(value)
 real(sp),intent(out)::value
 value=draw
 end subroutine
 function gasdev() result(value)
 real(sp)::value
 value=0
 end function
end module
"""

DRIVER = r"""program reference
 use spray_module
 use nr,only: draw
 implicit none
 integer :: code,flag,p
 real :: h1,h2,h3,h4,d,u,rg,rl,sigma,mu,epsilon,dt
 do
  read(*,*,iostat=code) h1,h2,h3,h4,p,d,u,rg,rl,sigma,mu,epsilon,dt,cbr,csb,flag,draw
  if(code/=0)exit
  part=0;part(1,1)=u;part(4,1)=d*1.e6;part(5,1)=1.;part(6,1)=400.
  part(9,1)=h4;part(11,1)=h3;part(12,1)=p;part(13,1)=h1;part(14,1)=h2
  xpart=0;ypart=0;zpart=0;index_cell=1;dau_p_d=0;dau_p_n=0;dau_p_props=0
  brokenp=0;stochabreak=flag==1
  call break_up(1,0.,0.,0.,rg,epsilon,d,rl,sigma,dt,mu)
  write(*,'(13(es26.17,1x))') part(13,1),part(14,1),part(11,1),part(9,1), &
      part(15,1),part(16,1),part(17,1),part(10,1)*1.e-6,real(brokenp), &
      dau_p_d(1)*1.e-6/d,dau_p_d(2)*1.e-6/d,dau_p_n(1),dau_p_n(2)
 end do
end program
"""


def requests():
    rows = []
    combinations = itertools.product((2e-5, 1e-4, 1e-3), (1., 1e4, 1e8),
                                     (0., 20.), (0, 1), (0, 1, 7), (.05, .5, .95))
    for i, (d, epsilon, speed, stochastic, poisson, quantile) in enumerate(combinations):
        dt = (1e-7, 1e-4, .01)[i % 3]
        age = (0., 1e-5, .001)[(i//3) % 3]
        rows.append([epsilon*.5, age, 1e3, age, poisson, d, speed,
                     (.4, 1., 4.)[(i//9) % 3], 750., .025, 2e-5, epsilon, dt,
                     1/math.sqrt(3), 2., stochastic, quantile])
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    for name, digest in EXPECTED.items():
        if hashlib.sha256((source/name).read_bytes()).hexdigest() != digest:
            raise ValueError("reference source identity: "+name)
    output.mkdir(parents=True, exist_ok=True)
    (output/"modules.f90").write_text(MODULES)
    (output/"driver.f90").write_text(DRIVER)
    rows = requests()
    inputs = "".join(" ".join(format(x, ".17g") for x in row)+"\n" for row in rows)
    (output/"input.dat").write_text(inputs)
    flags = ["-O2", "-ffp-contract=off", "-fdefault-real-8", "-fdefault-double-8",
             "-ffixed-line-length-none", "-fcheck=bounds"]
    with (output/"build.log").open("w") as log:
        subprocess.run([args.compiler]+flags+["-J"+str(output), str(output/"modules.f90"),
                       str(source/"break_up.F"), str(source/"cdf_spray.F"),
                       str(output/"driver.f90"), "-o", str(output/"reference")],
                       cwd=str(output), stdout=log, stderr=log, check=True)
    result = subprocess.run([str(output/"reference")], input=inputs.encode(),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    reference = result.stdout.decode()
    values = [[float(x) for x in row.split()] for row in reference.splitlines()]
    if len(values) != 324 or any(len(row) != 13 or not all(math.isfinite(x) for x in row) for row in values):
        raise ValueError("reference result layout")
    (output/"reference.dat").write_text(reference)
    fixture = str(len(rows))+"\n"+"".join(a+" "+b.strip()+"\n"
                for a, b in zip(inputs.splitlines(), reference.splitlines()))
    (output/"fixture.dat").write_text(fixture)
    report = dict(scope="complete break_up and cdf_spray; supplied state and random draw",
                  precision="FP64", requests=len(values), events=sum(row[8] == 2 for row in values),
                  source_sha256=EXPECTED, compiler_flags=flags,
                  compiler=subprocess.check_output([args.compiler, "--version"]).decode().splitlines()[0],
                  fixture_sha256=hashlib.sha256(fixture.encode()).hexdigest())
    (output/"reference.json").write_text(json.dumps(report, indent=2)+"\n")
    for name, digest in EXPECTED.items():
        if hashlib.sha256((source/name).read_bytes()).hexdigest() != digest:
            raise ValueError("reference source changed during execution: "+name)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
