#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
"""Compare native mixture faces with complete, read-only COAST VLS routines."""
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess


def run(command, directory, data=None):
    return subprocess.check_output(command, cwd=str(directory), input=data,
                                   universal_newlines=True)


def routine(source, name):
    pattern = r"(?im)^\s*subroutine " + name + r"\b[\s\S]*?^\s*end subroutine " + name + r"\b"
    found = re.search(pattern, source)
    if found is None:
        raise ValueError("COAST routine unavailable: " + name)
    return found.group(0) + "\n"


MODULES = """
module global
  implicit none
  integer, parameter :: lower=0,upper=15,mout=6
  real :: uround=epsilon(1.0),small=1.0e-30
end module
module arrays
  use global
  implicit none
  real :: ajc(lower:upper),x(lower:upper),y(lower:upper),z(lower:upper)
  integer :: imbflag(lower:upper)
end module
module imb_config
  implicit none
  logical :: ibm_config_present=.true.,ibm_immersed=.true.
end module
module coast_ibm_operators
  implicit none
  integer :: diffusion_limiter_adjustments=0
contains
"""

DRIVER = """
program face_probe
  use global
  use arrays
  implicit none
  integer,parameter :: imax=12
  integer :: io(0:imax),jo(0:imax),ko(0:imax),ibs(imax,imax),ibn(imax,imax)
  integer :: i,s,ns,thermal,allow,ios
  real :: b11(lower:upper),gi(lower:upper),w1(lower:upper),field(lower:upper)
  real :: samples(4,65),dependent(4),gam(imax),total,phi,physical
  io=[(i-1,i=0,imax)]; jo=0; ko=0; ibs=-100; ibn=-100
  ajc=1.0; x=[(real(i)/8.0,i=lower,upper)]; y=0.0; z=0.0
  b11=1.0; w1=0.5
  do
    read(*,*,iostat=ios) phi,physical,allow,ns,thermal
    if (ios/=0) exit
    do s=1,ns
      read(*,*) samples(:,s)
    enddo
    read(*,*) dependent
    if (thermal==1) read(*,*) samples(:,ns+2)
    samples(:,ns+1)=dependent
    gi=phi; imbflag=1
    if (allow==0) imbflag(5)=0
    total=physical
    do s=1,ns+1+thermal
      do i=lower,upper
        field(i)=samples(max(1,min(4,i-1)),s)
      enddo
      gam=physical
      call vls(field,0,b11,2,2,gi,ibs,ibn,8,w1,imax,io,imax,jo,imax,ko,gam)
      total=max(total,gam(3))
    enddo
    write(*,'(ES25.17E3)') total
  enddo
end program
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coast", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("check"))
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--dyn711", action="store_true", help="original gas VLS with upwind flat policy")
    parser.add_argument("--runner", default="", help="native executable launcher")
    args = parser.parse_args()
    coast, native, output = args.coast.resolve(), args.native.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    vls = coast / ("vls.F90" if args.dyn711 else "SRC.Coast/vls.F90")
    assets = [vls, native, Path(__file__).resolve()]
    source = output / "vls.F90"
    if args.dyn711:
        # The gas reference has no IBM selector argument. Keep its full routine
        # unchanged and adapt only the standalone caller and module stubs.
        modules = MODULES[:MODULES.index("module imb_config")]
        driver = DRIVER.replace("call vls(field,0,b11", "call vls(field,b11")
        source.write_text(modules + vls.read_text() + driver)
    else:
        ibm = coast / "SRC.Coast/physics/coast_ibm_operators.F90"
        original = ibm.read_text()
        wrapper = routine(original, "ibm_diffusion_limiter_allows_upwind_operator")
        limiter = routine(original, "imb_diffusion_limiter_allows_upwind")
        source.write_text(MODULES + wrapper + limiter + "end module\n" +
                          vls.read_text() + DRIVER)
        assets.append(ibm)
    assets.append(source)
    profiles = [
        ("closure", [[0.0, .2, .5, .5], [.6, .7, .5, .5]], []),
        ("linear", [[.1, .2, .3, .4], [.15, .2, .25, .3]], [[280, 290, 300, 310]]),
        ("extremum", [[.1, .3, .2, .15], [.6, .31, .4, .51]], []),
        ("constant_dependent", [[.1, .3, .2, .15], [.6, .4, .5, .55]], []),
        ("near_flat", [[.1, .3, .2, .15], [.6, .35, .45, .6]], []),
        ("thermal", [[.1, .2, .3, .4], [.15, .2, .25, .3]], [[300, 330, 310, 305]]),
        ("constant", [[.1, .2, .3, .4], [.2, .2, .2, .2]], [[300, 300, 300, 300]]),
        ("trace_roundoff", [[.1, .2, .3, .4], [1e-15, 2e-15, 1e-15, 2e-15]], []),
        ("trace_resolved", [[.1, .2, .3, .4], [1e-8, 2e-8, 1e-8, 2e-8]], []),
    ]
    if args.dyn711:
        profiles += [("signed", [[-2., -.5, .25, 1.5]], []),
                     ("signed_extremum", [[-1., .5, -.25, -2.]], [])]
    cases = []
    for name, fractions, thermal in profiles:
        for phi in (1.0, -1.0):
            for allow in ((1,) if args.dyn711 else (0, 1)):
                for physical in (0.0, .8):
                    cases.append((name, phi, physical, allow, fractions, thermal))
    records, reference_records = [], []
    for name, phi, physical, allow, fractions, thermal in cases:
        records.append("{} {} {} {} {}".format(phi, physical, allow,
                                               len(fractions), len(thermal)))
        records.extend(" ".join(map(str, field)) for field in fractions + thermal)
        # Feed identical canonical dependent fractions to the full reference
        # routine. FP32 quantization is then solely the selected REAL format.
        dependent = [float(Decimal(1)-sum(Decimal.from_float(field[i])
                      for field in fractions)) for i in range(4)]
        reference_records.append(records[-len(fractions)-len(thermal)-1])
        reference_records.extend(" ".join(map(str, field))
                                 for field in fractions + [dependent] + thermal)
    data = "\n".join(records) + "\n"
    reference_data = "\n".join(reference_records) + "\n"
    native_values = list(map(float, run(shlex.split(args.runner)+[str(native), "--flat-face" if args.dyn711 else "--face"], output, data).split()))
    if len(native_values) != len(cases):
        raise ValueError("Native face count mismatch")
    result = {"scope": ("dyn711 uniform Cartesian VLS; upwind flat policy; signed passive coordinates"
                         if args.dyn711 else "uniform Cartesian face VLS plus complete IBM gradient selector"),
              "source_directory": str(coast),
              "coast_head": (None if args.dyn711 else run(["git", "rev-parse", "HEAD"], coast).strip()),
              "compiler": run([args.compiler, "--version"], output).splitlines()[0],
              "inputs": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in assets}, "comparisons": []}
    passed = True
    for bits, flags, tolerance in ((32, [], 2e-6), (64, ["-fdefault-real-8"], 1e-11)):
        binary = output / ("vls" + str(bits))
        run([args.compiler, "-O2", "-ffree-line-length-none"] + flags +
            [str(source), "-o", str(binary)], output)
        values = list(map(float, run([str(binary)], output, reference_data).split()))
        if len(values) != len(cases):
            raise ValueError("COAST face count mismatch")
        maximum = 0.0
        exceptions, rounding = [], []
        for case, observed, expected in zip(cases, native_values, values):
            name, phi, physical, allow, _, _ = case
            if not args.dyn711 and (name.startswith("constant") or name=="trace_roundoff") and physical == 0.0:
                # Constant and sub-resolution coordinates leave the resolved
                # linear component at its central reconstruction.
                native_expected = .5 if name=="constant_dependent" and phi>0 and allow else 0.0
                valid = abs(observed-native_expected) <= 1e-13 and abs(expected-.5) <= tolerance
                passed = passed and valid
                exceptions.append({"case": name, "phi": phi, "allow": allow,
                                   "hundun": observed, "coast": expected, "passed": valid})
                continue
            error = abs(observed-expected)/max(1.0, abs(expected))
            if bits==32 and name=="near_flat":
                rounding.append({"case": name, "phi": phi, "allow": allow,
                                 "physical": physical, "hundun": observed,
                                 "coast": expected, "difference": error})
                passed = passed and 0.0<=observed<=max(physical,.5) and 0.0<=expected<=max(physical,.5)+tolerance
                continue
            maximum = max(maximum, error)
            passed = passed and error <= tolerance
        result["comparisons"].append({"real_bits": bits, "cases": len(cases),
                                      "equivalent_max_difference": maximum,
                                      "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
                                      "tolerance": tolerance, "constant_coordinate": exceptions,
                                      "fp32_flat_gradient_rounding": rounding})
    for name in ("global", "arrays", "imb_config", "coast_ibm_operators"):
        module = output / (name + ".mod")
        if module.exists():
            module.unlink()
    result["passed"] = passed
    destination = output / "vls.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
