#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare complete GTMC COAST step/mixer routines with the FP64 IEM row."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import subprocess

DRIVER = """
module arrays
  real :: ajc(64),f(64),gam_sgs(64),coef(2,64),gam(64),visc(64),pr(64),eps(64),dyn_LM(64)
  real :: fold(64),drhodt(64)
  integer :: jo(4),ko(4),nfo(8)
end module
module global
  integer :: l=2,m=2,n=2,nf=0,pc=1,bpc=2,nv=1,nvdp=2,lower=1,upper=64
  real :: pdf_mix_cd_value,dtim
end module
module chemistry
  integer :: isp=1,nsc=1,jfuel=1
  real :: kappa(64),temp_i(64),temp_k(64),tim_flow(64)
end module
module sgs_pdf
end module
program reference
  use arrays
  use global
  implicit none
  real :: volume,mu,mut,cd,rho_value,dt,mean,qold,rho(64),first_a,first_b,solution
  integer :: ios
  jo=0;ko=0;nfo=0
  do
    read(*,*,iostat=ios) volume,mu,mut,cd,rho_value,dt,mean,qold
    if(ios<0)exit
    if(ios/=0)stop 1
    ajc=volume;visc=mu;gam_sgs=mut;pdf_mix_cd_value=cd;rho=rho_value;dtim=dt
    f=mean;fold=qold;coef=0
    call step(rho)
    call mixer
    first_a=coef(1,2);first_b=coef(2,2)
    solution=first_b/first_a
    ! A repeated jstep starts from the same accepted history and new assembly.
    coef=0
    call step(rho)
    call mixer
    if(coef(1,2)/=first_a.or.coef(2,2)/=first_b)stop 2
    write(*,'(*(ES25.16E3,1X))') first_a,first_b,solution,first_a*solution-first_b
  end do
end program
"""


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="read-only GTMC SRC.Coast directory containing step.F90 and mixer.F90")
    parser.add_argument("--output", type=Path, default=Path("check/im"))
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--native-compiler",
                        help="also build the actual Hundun ESF workspace diagnostic using this C++ compiler")
    parser.add_argument("--equation-probe", type=Path,
                        help="compare registered IEM sources through the native species equation test executable")
    parser.add_argument("--c-z", type=float, default=1.)
    parser.add_argument("--schmidt", type=float, default=.7)
    args = parser.parse_args()
    if any(not math.isfinite(x) or x <= 0 for x in (args.c_z, args.schmidt)):
        raise ValueError("positive finite Cz and Schmidt numbers are required")
    source, output = args.source.resolve(), args.output.resolve()
    if (output / "evidence.json").exists():
        raise ValueError("choose a fresh output directory; this directory has frozen evidence")
    output.mkdir(parents=True, exist_ok=True)
    original = source / "mixer.F90"
    found = re.search(r"(?im)^\s*subroutine mixer\b[\s\S]*?^\s*end subroutine mixer\b", original.read_text())
    if found is None:
        raise ValueError("complete mixer routine is required")
    (output / "mix.f90").write_text(found.group(0) + "\n")
    (output / "driver.f90").write_text(DRIVER)
    command = [args.compiler, "-O2", "-fdefault-real-8", "-fdefault-double-8",
               "-J" + str(output), str(output / "driver.f90"),
               str(output / "mix.f90"), str(source / "step.F90"),
               "-o", str(output / "ref")]
    subprocess.check_call(command)
    rows = list(itertools.product([1e-12, 1e-9, 1e-6], [1e-5, 4e-5],
                [0., 1e-5, 1e-3], [0., 2.], [.1, 1., 5.],
                [1e-7, 1e-4], [.2, .8], [0., 1.]))
    data = "".join(" ".join(format(v, ".17g") for v in row) + "\n" for row in rows)
    (output / "input.dat").write_text(data)
    result = subprocess.check_output([str(output / "ref")], input=data,
                                     universal_newlines=True)
    (output / "coast.dat").write_text(result)
    lines = result.splitlines()
    if len(lines) != len(rows):
        raise ValueError("reference row count differs from input")
    error = [0., 0., 0., 0.]
    for row, line in zip(rows, lines):
        volume, mu, mut, cd, density, dt, mean, old = row
        beta = .5 * cd * (mu + mut) / volume**(2. / 3.)
        diagonal = density / dt + beta
        rhs = density / dt * old + beta * mean
        expected = [diagonal, rhs, rhs / diagonal]
        actual = list(map(float, line.split()))
        if len(actual) != 4 or any(not (-float("inf") < v < float("inf")) for v in actual):
            raise ValueError("nonfinite or incomplete reference row")
        for j in range(3):
            error[j] = max(error[j], abs(actual[j] - expected[j]) / max(abs(expected[j]), 1.))
        error[3] = max(error[3], abs(actual[3]) / max(abs(rhs), 1.))
    report = dict(states=len(rows), max_normalized_error=error,
                  scope="actual COAST step and mixer; independent analytic oracle",
                  precision="FP64 default REAL and literals", command=command,
                  compiler=subprocess.check_output([args.compiler, "--version"],
                                                   universal_newlines=True).splitlines()[0],
                  source={str(original): sha(original), str(source / "step.F90"): sha(source / "step.F90")},
                  reference_program_sha256=sha(output / "ref"))
    (output / "formula.json").write_text(json.dumps(report, indent=2) + "\n")
    if max(error) > 1e-11:
        raise ValueError("IEM row comparison exceeds 1e-11")
    print("IEM reference rows=%d normalized_error=%.17g" % (len(rows), max(error)))
    if args.native_compiler:
        compare_native(args, output, data, rows, lines)
    if args.equation_probe:
        compare_equations(args.equation_probe.resolve(), output, data, rows, lines)


def compare_equations(probe, output, data, rows, reference):
    command = [str(probe), "--iem"]
    result = subprocess.check_output(command, input=data, universal_newlines=True)
    (output / "rows.dat").write_text(result)
    native = result.splitlines()
    if len(native) != len(rows):
        raise ValueError("equation probe row count differs from input")
    error = [0.] * 4
    for line, ref in zip(native, reference):
        actual, expected = [list(map(float, s.split())) for s in (line, ref)]
        if len(actual) != 4 or any(not math.isfinite(v) for v in actual):
            raise ValueError("incomplete or nonfinite equation row")
        for j in range(3):
            error[j] = max(error[j], abs(actual[j] - expected[j]) / max(1., abs(expected[j])))
        error[3] = max(error[3], abs(actual[3]) / max(1., abs(expected[1])))
    report = dict(states=len(rows), max_normalized_error=error,
                  scope="native registered IEM contribution plus complete species assembly; uniform-field row reduction; repeated accepted history",
                  command=command, probe_sha256=sha(probe),
                  reference_sha256=sha(output / "formula.json"))
    (output / "rows.json").write_text(json.dumps(report, indent=2) + "\n")
    if max(error) > 1e-11:
        raise ValueError("native IEM equation comparison exceeds 1e-11")
    print("Native IEM equation rows=%d normalized_error=%.17g" % (len(rows), max(error)))


def compare_native(args, output, data, rows, reference):
    repo = Path(__file__).resolve().parents[1]
    source = repo / "versions/v0.4/src/models_esf.cpp"
    driver = repo / "tools/iem.cpp"
    executable = output / "native"
    command = [args.native_compiler, "-std=c++20", "-O2", "-ffunction-sections",
               "-fdata-sections", "-I" + str(repo / "versions/v0.4/include"),
               "-I" + str(source.parent), str(driver), str(source),
               "-Wl,--gc-sections", "-o", str(executable)]
    subprocess.check_call(command)
    run = [str(executable), str(args.c_z), str(args.schmidt)]
    result = subprocess.check_output(run, input=data, universal_newlines=True)
    (output / "native.dat").write_text(result)
    native = result.splitlines()
    if len(native) != len(rows):
        raise ValueError("native row count differs from input")
    worst = [{"absolute_update_gap": -1.} for _ in range(2)]
    mean_error = 0.
    for row, ref, line in zip(rows, reference, native):
        values = list(map(float, line.split()))
        if len(values) != 5 or any(not math.isfinite(v) for v in values):
            raise ValueError("incomplete or nonfinite native row")
        target = float(ref.split()[2])
        for j in range(2):
            gap = abs(values[2*j] - target)
            if gap > worst[j]["absolute_update_gap"]:
                worst[j] = dict(absolute_update_gap=gap, input=row,
                                coast_update=target, native_update=values[2*j],
                                native_relaxation=values[2*j+1])
        mean_error = max(mean_error, values[4])
    if mean_error > 1e-12:
        raise ValueError("native ensemble mean preservation check failed")
    report = dict(status="native_method_gap_measured", states=len(rows),
                  scope="spatially uniform two-field enthalpy mixing; actual Hundun Workspace.advance and complete COAST step/mixer",
                  material="unit Lewis, constant Pr=Sc; gamma=(mu+mut)/Sc",
                  c_z=args.c_z, schmidt=args.schmidt, chemistry="inactive",
                  current_scale=worst[0], matched_coast_rate=worst[1],
                  max_ensemble_mean_residual=mean_error,
                  command=command, run=run,
                  compiler=subprocess.check_output([args.native_compiler, "--version"],
                                                   universal_newlines=True).splitlines()[0],
                  sources={str(p): sha(p) for p in (source, driver, Path(__file__).resolve())},
                  program_sha256=sha(executable),
                  reference_sha256=sha(output / "formula.json"),
                  pending="implicit spatial ESF matrix, coupled mean authority, GTMC material properties and end-to-end acceptance")
    (output / "native.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Native IEM current-scale gap=%.9g matched-rate gap=%.9g" %
          (worst[0]["absolute_update_gap"], worst[1]["absolute_update_gap"]))


if __name__ == "__main__":
    main()
