#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare Hundun's linked Smagorinsky kernel with the full COAST routine.

The external one-rank fixture supplies affine fields and geometry. Wall branches
are guarded; MPI halos are identity operations for this interior stencil test.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import shutil
import subprocess

FORTRAN = r'''
module global
 use mpi, only: mpi_real=>MPI_DOUBLE_PRECISION, mpi_comm_world, mpi_in_place, mpi_min, mpi_max, &
  mpi_allreduce, mpi_init, mpi_finalize, mpi_comm_size
 implicit none
 integer,parameter :: imax=5,jmax=5,kmax=5,l=4,m=4,n=4,lower=1,upper=125,nvu=1,nvv=2,nvw=3
 integer,parameter :: fsth=1,fnth=2,fwst=3,fest=4,flft=5,frht=6
 integer :: io(5),jo(5),ko(5),nfo(3),info
 real :: big=huge(1.),cs0,gamax,gamin
end module
module arrays
 implicit none
 real,target :: fbar(375)
 real :: x(125),y(125),z(125),ajc(125),rho(125),visc(125),gam_sgs(125),q(125),eps(125)
 integer :: ibs(5,5)=0,ibn(5,5)=0,ibw(5,5)=0,ibe(5,5)=0,ibl(5,5)=0,ibr(5,5)=0
end module
module extras
 use iso_c_binding
 implicit none
 contains
 function assign_pointer(first,low,high) result(p)
 real,target :: first
 integer :: low,high
 real,pointer :: p(:)
 call c_f_pointer(c_loc(first),p,[high-low+1])
 end function
end module
module exchange
 implicit none
 contains
 subroutine mpi_halo_exchange(values,nfield)
 use global,only:mpi_comm_size,mpi_comm_world
 real :: values(*)
 integer :: nfield,size,error
 call mpi_comm_size(mpi_comm_world,size,error)
 if(size/=1.or.nfield/=1) stop 2
 end subroutine
end module
real function slip_v(a,b,c)
 integer :: a,b,c
 stop 3
end function
real function wall_v(a,b,c,d,e,f)
 integer :: a,b,c,d,e,f
 stop 4
end function
program smag_reference
 use global
 use arrays
 implicit none
 integer :: sample,i,j,k,c,v,d
 real :: g(3,3),width(3),position(3),rr,cc
 call mpi_init(info)
 io=[0,1,2,3,4];jo=5*io;ko=25*io;nfo=[0,125,250]
 do sample=1,32
  read(*,*) width,rr,cc,((g(v,d),d=1,3),v=1,3)
  cs0=cc;rho=rr;visc=1.8e-5;ajc=product(width)
  do k=1,5
   do j=1,5
    do i=1,5
     c=i+jo(j)+ko(k);position=[real(i-3),real(j-3),real(k-3)]*width
     x(c)=position(1);y(c)=position(2);z(c)=position(3)
     do v=1,3
      fbar(c+nfo(v))=sum(g(v,:)*position)
     enddo
    enddo
   enddo
  enddo
  call gamma_smagorinsky
  write(*,'(I3,4(1X,ES25.17))') sample,gam_sgs(63),q(63),eps(63),eps(63)/rho(63)
 enddo
 call mpi_finalize(info)
end program
'''

CPP = r'''
#include "hundun/v04_physics.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
int main() {
  using namespace hundun::v04;
  for(int sample=1;sample<=32;++sample) {
    double dx,dy,dz,rho,cs,nu;VelocityGradient g;
    std::cin>>dx>>dy>>dz>>rho>>cs;
    for(auto& v:g.value)std::cin>>v;
    if(!std::cin || !smagorinsky_kinematic_viscosity(g,std::cbrt(dx*dy*dz),cs,nu))return 1;
    SgsState state;
    if(!sgs_state_from_viscosity(g,std::cbrt(dx*dy*dz),rho,1.8e-5,nu,state))return 2;
    std::cout<<sample<<' '<<std::setprecision(17)<<rho*nu<<' '
        <<state.kinetic_energy_m2_s2<<' '<<state.dissipation_w_m3<<' '
        <<state.specific_dissipation_m2_s3<<'\n';
  }
}
'''

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coast', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wrapper', type=Path,
                        help='shell wrapper providing the Hundun compiler/runtime environment')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parent.parent
    coast = args.coast.resolve()
    source = coast / 'SRC.Coast/gamma_smagorinsky.F90'
    out = args.output.resolve()
    if out == coast or coast in out.parents or out.exists():
        parser.error('choose a fresh output directory outside COAST')
    library = args.library.resolve()
    source_hash, library_hash = sha(source), sha(library)
    out.mkdir(parents=True)
    (out / 'driver.f90').write_text(FORTRAN)
    (out / 'probe.cpp').write_text(CPP)
    shutil.copyfile(str(source), str(out / 'gamma.f90'))
    rng = random.Random(713)
    rows = []
    for i in range(32):
        g = [rng.randint(-16,16)/8 for _ in range(9)]
        if i == 0: g = [0.] * 9
        if i == 1: g = [0.,2.,0.,-2.,0.,0.,0.,0.,0.]
        if i == 2: g = [1.,0.,0.,0.,1.,0.,0.,0.,1.]
        rows.append([.25,.5,.125,1.+(i%5)/4.,.125 if i<16 else .25]+g)
    inputs = ''.join(' '.join(map(str,row))+'\n' for row in rows)
    (out / 'input.dat').write_text(inputs)
    wrapper = ['bash',str(args.wrapper.resolve())] if args.wrapper else []
    fc = ['mpif90','-O2','-fdefault-real-8','-fdefault-double-8',
          '-ffree-line-length-none','driver.f90','gamma.f90','-o','coast']
    cxx = wrapper + ['env','OMPI_CXX=clang++','mpicxx','-std=c++17','-stdlib=libstdc++',
          '-D_GLIBCXX_USE_CXX11_ABI=1','-O3','-flto=thin','-fuse-ld=lld','-ffp-contract=off',
          '-I',str(repo/'versions/v0.4/include'),str(out/'probe.cpp'),str(library),'-o',str(out/'hundun')]
    with (out/'build.log').open('w') as log:
        subprocess.check_call(fc,cwd=str(out),stdout=log,stderr=subprocess.STDOUT)
        subprocess.check_call(cxx,cwd=str(repo),stdout=log,stderr=subprocess.STDOUT)
    reference = subprocess.check_output([str(out/'coast')],input=inputs.encode(),cwd=str(out)).decode()
    native = subprocess.check_output(wrapper+[str(out/'hundun')],input=inputs.encode(),cwd=str(repo)).decode()
    (out/'coast.dat').write_text(reference)
    (out/'hundun.dat').write_text(native)
    a = [list(map(float,line.split())) for line in reference.splitlines()]
    b = [list(map(float,line.split())) for line in native.splitlines()]
    assert len(a)==len(b)==32
    assert all(x[0]==y[0]==i+1 for i,(x,y) in enumerate(zip(a,b)))
    assert all(math.isfinite(v) for row in a+b for v in row)
    quantities = ['mu_sgs_pa_s','kinetic_energy_m2_s2','dissipation_w_m3','specific_dissipation_m2_s3']
    assert all(len(row)==5 for row in a+b)
    errors = {name:max(abs(x[j]-y[j])/max(1e-30,abs(x[j]),abs(y[j])) for x,y in zip(a,b))
              for j,name in enumerate(quantities,1)}
    error = max(errors.values())
    assert error<=1e-11,error
    assert sha(source)==source_hash and sha(library)==library_hash
    report = {'schema':'hundun_sgs_state_parity_v2', 'samples':32, 'passed':True,
        'scope':'one-rank interior affine gradients; full original routine; FP64',
        'relative_error_max':error,'quantity_relative_errors':errors,'threshold':1e-11,'coast_source_sha256':source_hash,
        'hundun_library_sha256':library_hash,'coast_source':str(source),'library':str(library),
        'compile':[fc,cxx], 'artifacts':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/'parity.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Smagorinsky FP64 parity: 32 samples, max relative error',error)


if __name__ == '__main__':
    main()
