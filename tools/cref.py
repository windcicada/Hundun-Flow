#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Build the external 40-layer periodic-cylinder COAST geometry probe.

The original objects are relinked read-only. The wrapper translates z coordinate
images after the original halo exchange and audits Cartesian volumes/metrics.
This is a reference preparation tool; full physics/performance acceptance is
recorded separately in docs/alg.md. Requires the original gfortran module ABI.
"""
import argparse
import hashlib
import json
import pathlib
import re
import subprocess

WRAPPER = r'''
! External geometry probe for the 40-layer, self-periodic z reference.
! Original COAST object files and all flow/scalar exchanges stay intact.
subroutine periodic_coordinate_probe(field, depth) bind(C, name='__wrap_mpi_halo_exchange_')
  use iso_c_binding, only: c_float, c_int
  use arrays, only: z, zv, jo, ko, ibl, ibr
  use global, only: lower, upper, l, m, n
  implicit none
  real(c_float), intent(inout) :: field(*)
  integer(c_int), intent(in) :: depth
  integer :: i, j, layer, ilower, iupper
  real(c_float) :: period
  interface
    subroutine original_exchange(values, halo) bind(C, name='__real_mpi_halo_exchange_')
      import c_float, c_int
      real(c_float) :: values(*)
      integer(c_int) :: halo
    end subroutine
  end interface
  call original_exchange(field, depth)
  if (.not. allocated(z)) return
  if (loc(field(1)) /= loc(z(lower))) return
  if (n /= 41 .or. depth < 1 .or. depth > 2) stop 9201
  if (any(ibl(2:l,2:m) /= -100) .or. any(ibr(2:l,2:m) /= -100)) stop 9202
  period = zv(2+jo(2)+ko(n))-zv(2+jo(2)+ko(1))
  if (abs(period-0.0314159265_c_float) > 1.0e-7_c_float) stop 9203
  ! Decomp/config exchanges exactly i=2:l, j=2:m for these z interfaces.
  ! Its received material/scalar values are periodic; coordinate images have
  ! the additional translation. The operation occurs only on the z array.
  do layer=1,depth
    do j=2,m
      do i=2,l
        ilower=i+jo(j)+ko(2-layer)-lower+1
        iupper=i+jo(j)+ko(n+layer)-lower+1
        field(ilower)=field(ilower)-period
        field(iupper)=field(iupper)+period
      enddo
    enddo
  enddo
end subroutine

subroutine periodic_geometry_audit() bind(C, name='__wrap_geom_')
  use iso_c_binding, only: c_float
  use arrays, only: x,y,z,xv,yv,zv,ajc,b11,b22,b33,jo,ko
  use global, only: l,m,n
  use exchange, only: master
  use mpi
  implicit none
  integer :: i,j,k,c,mpi_status
  real(c_float) :: dx,dy,dz,v,expected,errors(4),global_errors(4)
  interface
    subroutine original_geom() bind(C, name='__real_geom_')
    end subroutine
  end interface
  call original_geom()
  errors=0
  do k=2,n
    do j=2,m
      do i=2,l
        c=i+jo(j)+ko(k)
        dx=xv(c)-xv(c-1)
        dy=yv(c)-yv(c-jo(2))
        dz=zv(c)-zv(c-ko(2))
        v=dx*dy*dz
        errors(1)=max(errors(1),abs(ajc(c)/v-1))
        expected=dy*dz/(x(c+1)-x(c))
        errors(2)=max(errors(2),abs(b11(c)/expected-1))
        expected=dx*dz/(y(c+jo(2))-y(c))
        errors(3)=max(errors(3),abs(b22(c)/expected-1))
        expected=dx*dy/(z(c+ko(2))-z(c))
        errors(4)=max(errors(4),abs(b33(c)/expected-1))
      enddo
    enddo
  enddo
  call mpi_allreduce(errors,global_errors,4,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,mpi_status)
  if(master) write(*,'(a,4es16.7)') 'COAST_REFERENCE_GEOMETRY relative_error V,b11,b22,b33=',global_errors
  if(mpi_status /= MPI_SUCCESS .or. any(global_errors > 2.0e-5_c_float)) stop 9204
end subroutine

subroutine reference_boundary_audit() bind(C, name='__wrap_inlet_profile_apply_')
  use arrays, only: ibn,ibs,ibe,ibw,ibl,ibr
  use global, only: l,m,n
  use exchange, only: master
  use mpi
  implicit none
  integer :: local_counts(4),global_counts(4),status,phase
  interface
    subroutine original_profile() bind(C, name='__real_inlet_profile_apply_')
    end subroutine
  end interface
  do phase=1,2
    if(phase == 2) call original_profile()
    local_counts(1)=count(ibs(2:m,2:n) == -1)
    local_counts(2)=count(ibn(2:m,2:n) == -6)
    local_counts(3)=count(ibn(2:m,2:n) == -60)
    local_counts(4)=count(ibn(2:m,2:n) == -2)
    call mpi_allreduce(local_counts,global_counts,4,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    if(master) write(*,'(a,i1,a,4i10)') 'COAST_REFERENCE_BOUNDARY phase=',phase, &
      ' south(-1),north(-6,-60,-2)=',global_counts
    if(status /= MPI_SUCCESS) stop 9205
  enddo
end subroutine
'''


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coast", required=True, type=pathlib.Path)
    parser.add_argument("--output", default="check", type=pathlib.Path)
    parser.add_argument("--build", type=pathlib.Path,
                        help="external fresh build root containing x86_64 objects and modules")
    args = parser.parse_args()
    root, out = args.coast.resolve(), args.output.resolve()
    if root == out or root in out.parents:
        parser.error("choose an output directory outside the COAST checkout")
    src = root / "SRC.Coast"
    built = args.build.resolve() if args.build else src
    make = (src / "Makefile").read_text().replace("\\\n", " ")
    objects = []
    for name in ["MSRC", "FSRC"]:
        entries = re.search(r"^" + name + r"\s*=([^\n]+)", make, re.M).group(1)
        objects += [built / "x86_64" / re.sub(r"\.(F90|F|f)$", ".o", f)
                    for f in entries.split()]
    entries = re.search(r"^CPPOBJS\s*=([^\n]+)", make, re.M).group(1)
    objects += [built / "x86_64" / f for f in entries.split()]
    objects.append(built / "x86_64/imb_mesh.o")
    for path in objects:
        if not path.is_file():
            parser.error("original object required: " + str(path))
    out.mkdir(parents=True, exist_ok=True)
    driver, obj, binary = out / "pwrap.f90", out / "pwrap.o", out / "cref"
    driver.write_text(WRAPPER)
    compile_cmd = ["mpif90", "-O3", "-I" + str(built), "-c", str(driver), "-o", str(obj)]
    subprocess.check_call(compile_cmd, cwd=str(out))
    lib = root / "third_party/oce_sysroot/usr/lib/x86_64-linux-gnu"
    libraries = ("TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase TKXSBase TKShHealing "
                 "TKMesh TKTopAlgo TKBRep TKGeomAlgo TKGeomBase TKG3d TKG2d TKMath TKernel")
    link = (["mpif90", "-fno-lto"] + [str(path) for path in objects] +
            [str(obj), "-Wl,--wrap=mpi_halo_exchange_", "-Wl,--wrap=geom_", "-Wl,--wrap=inlet_profile_apply_",
             "-L" + str(lib), "-Wl,--disable-new-dtags,-rpath," + str(lib)] +
            ["-l" + name for name in libraries.split()] + ["-lstdc++", "-o", str(binary)])
    subprocess.check_call(link, cwd=str(out))
    report = {
        "purpose": "periodic geometry probe; physics/performance acceptance is separate",
        "scope": "40 cell self-periodic z, length 0.0314159265 m",
        "coast_source": str(src), "object_build": str(built),
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(root),
                                        universal_newlines=True).strip(),
        "compiler": subprocess.check_output(["mpif90", "--version"],
                                            universal_newlines=True).splitlines()[0],
        "compile": compile_cmd, "link": link,
        "objects": {str(path): sha(path) for path in objects},
        "wrapper_sha256": sha(driver), "executable_sha256": sha(binary),
        "makefile_sha256": sha(src / "Makefile"),
    }
    (out / "cref.json").write_text(json.dumps(report, indent=2) + "\n")
    print("linked {} original objects: {}".format(len(objects), binary))


if __name__ == "__main__":
    main()
