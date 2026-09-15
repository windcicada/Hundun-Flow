#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Observe the original ICCG pressure row and density/flux update correspondence."""
import argparse
import json
from pathlib import Path
import subprocess
from ccost import sha

SOURCE = r'''
module parity_state
  implicit none
  double precision, allocatable, save :: row(:),previous(:),half_storage(:),compress_row(:)
  double precision, allocatable :: compress_delta(:,:)
  double precision :: metrics(9)
  logical :: pending=.false.
end module
subroutine observed_compress() bind(C,name='__wrap_compress_')
  use parity_state
  use arrays, only: coef
  use global, only: lower,upper,pc,sc,nc,wc,ec,lc,rc
  implicit none
  integer :: c
  interface
    subroutine original_compress() bind(C,name='__real_compress_')
    end subroutine
  end interface
  if(.not.allocated(compress_delta)) allocate(compress_delta(7,lower:upper))
  do c=lower,upper
    compress_delta(:,c)=dble([coef(pc,c),coef(sc,c),coef(nc,c), &
      coef(wc,c),coef(ec,c),coef(lc,c),coef(rc,c)])
  enddo
  call original_compress()
  do c=lower,upper
    compress_delta(:,c)=dble([coef(pc,c),coef(sc,c),coef(nc,c), &
      coef(wc,c),coef(ec,c),coef(lc,c),coef(rc,c)])-compress_delta(:,c)
  enddo
end subroutine

subroutine observed_update(u,v,w,dp) bind(C,name='__wrap_update_')
  use iso_c_binding, only: c_float
  use parity_state
  use arrays, only: jo,ko,coef,ajc,rho,rhold,gi,gj,gk,drhodp,imbflag_real
  use global, only: lower,upper,l,m,n,dtim,istep,pc,sc,nc,wc,ec,lc,rc,bpc
  use coast_screen_summary, only: coast_pressure_solver_backend
  use exchange, only: master
  use mpi
  implicit none
  real(c_float) :: u(lower:upper),v(lower:upper),w(lower:upper),dp(lower:upper)
  double precision :: volume,scale,before,divergence
  integer :: i,j,k,c,status
  interface
    subroutine original_update(u,v,w,dp) bind(C,name='__real_update_')
      import c_float
      real(c_float) :: u(*),v(*),w(*),dp(*)
    end subroutine
  end interface
  if(trim(coast_pressure_solver_backend)/='legacy_iccg') stop 9312
  if(.not.allocated(compress_delta)) stop 9313
  if (.not.allocated(row)) allocate(row(lower:upper),previous(lower:upper),half_storage(lower:upper),compress_row(lower:upper))
  metrics=0
  do k=2,n
    do j=2,m
      do i=2,l
        c=i+jo(j)+ko(k)
        if(imbflag_real(c)/=1.0) cycle
        volume=abs(dble(ajc(c)))
        scale=dble(rho(c))/dble(dtim)
        divergence=(dble(gi(c+1))-dble(gi(c))+dble(gj(c+jo(2)))-dble(gj(c)) &
          +dble(gk(c+ko(2)))-dble(gk(c)))/volume
        before=(dble(rho(c))-dble(rhold(c)))/dble(dtim)+divergence
        row(c)=(dble(coef(pc,c))*dble(dp(c)) &
          -dble(coef(sc,c))*dble(dp(c-1))-dble(coef(nc,c))*dble(dp(c+1)) &
          -dble(coef(wc,c))*dble(dp(c-jo(2)))-dble(coef(ec,c))*dble(dp(c+jo(2))) &
          -dble(coef(lc,c))*dble(dp(c-ko(2)))-dble(coef(rc,c))*dble(dp(c+ko(2))))/volume &
          -dble(coef(bpc,c))
        compress_row(c)=compress_delta(1,c)*dble(dp(c)) &
          -compress_delta(2,c)*dble(dp(c-1))-compress_delta(3,c)*dble(dp(c+1)) &
          -compress_delta(4,c)*dble(dp(c-jo(2)))-compress_delta(5,c)*dble(dp(c+jo(2))) &
          -compress_delta(6,c)*dble(dp(c-ko(2)))-compress_delta(7,c)*dble(dp(c+ko(2)))
        previous(c)=dble(rho(c))
        half_storage(c)=.5d0*dble(drhodp(c))*dble(dp(c))/dble(dtim)
        metrics(1)=max(metrics(1),abs(before+dble(coef(bpc,c)))/scale)
        metrics(4)=max(metrics(4),abs(before)/scale)
        metrics(6)=max(metrics(6),abs(row(c))/scale)
      enddo
    enddo
  enddo
  call original_update(u,v,w,dp)
  pending=.true.
end subroutine

subroutine observed_boundary(u,v,w) bind(C,name='__wrap_bndry3_')
  use iso_c_binding, only: c_float
  use parity_state
  use arrays, only: jo,ko,ajc,rho,rhold,gi,gj,gk,imbflag_real
  use global, only: lower,upper,l,m,n,dtim,istep
  use exchange, only: master
  use mpi
  implicit none
  real(c_float) :: u(lower:upper),v(lower:upper),w(lower:upper)
  double precision :: volume,scale,after,global_metrics(9)
  integer :: i,j,k,c,status
  interface
    subroutine original_boundary(u,v,w) bind(C,name='__real_bndry3_')
      import c_float
      real(c_float) :: u(*),v(*),w(*)
    end subroutine
  end interface
  if(pending) then
  do k=2,n
    do j=2,m
      do i=2,l
        c=i+jo(j)+ko(k)
        if(imbflag_real(c)/=1.0) cycle
        volume=abs(dble(ajc(c)))
        scale=previous(c)/dble(dtim)
        after=(dble(rho(c))-dble(rhold(c)))/dble(dtim) &
          +(dble(gi(c+1))-dble(gi(c))+dble(gj(c+jo(2)))-dble(gj(c)) &
          +dble(gk(c+ko(2)))-dble(gk(c)))/volume
        metrics(2)=max(metrics(2),abs(after-row(c))/scale)
        metrics(3)=max(metrics(3),abs(half_storage(c))/scale)
        metrics(5)=max(metrics(5),abs(after)/scale)
        if(i>2 .and. i<l .and. j>2 .and. j<m .and. k>2 .and. k<n) then
          if(imbflag_real(c-1)==1 .and. imbflag_real(c+1)==1 .and. &
             imbflag_real(c-jo(2))==1 .and. imbflag_real(c+jo(2))==1 .and. &
             imbflag_real(c-ko(2))==1 .and. imbflag_real(c+ko(2))==1) then
            metrics(8)=max(metrics(8),abs(after-row(c))/scale)
            metrics(9)=max(metrics(9),abs(after-row(c)+compress_row(c)-half_storage(c))/scale)
          endif
        endif
        metrics(7)=max(metrics(7),abs((dble(rho(c))-previous(c))/dble(dtim)-half_storage(c))/scale)
      enddo
    enddo
  enddo
  call MPI_Allreduce(metrics,global_metrics,9,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
  if(status/=MPI_SUCCESS) stop 9311
  if(master) write(*,'(A,I8,9(1X,ES22.14E3))') 'COAST_PRESSURE_PARITY ',istep,global_metrics
  pending=.false.
  endif
  call original_boundary(u,v,w)
end subroutine
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reference = json.loads(args.reference.read_text())
    out = args.output.resolve()
    coast = Path(reference['coast_source']).resolve().parent
    if out == coast or coast in out.parents:
        parser.error('output must be external to COAST')
    for path, digest in reference['objects'].items():
        if sha(path) != digest:
            parser.error('reference object changed: '+path)
    out.mkdir(parents=True, exist_ok=True)
    source, obj, binary = out/'par.f90', out/'par.o', out/'cpar'
    source.write_text(SOURCE)
    command = ['mpif90','-O2','-I'+reference['object_build'],'-c',str(source),'-o',str(obj)]
    subprocess.check_call(command,cwd=str(out))
    link=reference['link'][:]
    link[-1]=str(binary)
    link[2:2]=[str(obj),'-Wl,--wrap=update_','-Wl,--wrap=bndry3_','-Wl,--wrap=compress_']
    subprocess.check_call(link,cwd=str(out))
    report={'scope':__doc__, 'backend_required':'legacy_iccg',
        'scaling':'all maxima normalized by the pre-update rho/dt',
        'columns':['step','rhs_continuity_gap','row_update_gap','half_density_increment',
                   'continuity_before','continuity_after','matrix_residual','density_roundoff','interior_row_update_gap','interior_compress_remainder'],
        'reference_sha256':sha(args.reference),'source_sha256':sha(source),
        'program_sha256':sha(binary),'compile':command,'link':link}
    (out/'cpar.json').write_text(json.dumps(report,indent=2)+'\n')
    print(str(binary))

if __name__ == '__main__':
    main()
