#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Link read-only COAST kernels with per-rank wall timers and endpoint mass audit.

The observational executable is separate from the uninstrumented reference.
Kernel wall times include nested calls; the whole interval excludes output().
Use one final output step after the measured window and discard the first step.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


KERNELS = [('condif', '1+merge(0,1,nv<=nf)'),
           ('cgstab', '3+merge(0,1,nv<=nf)'),
           ('press', '5'), ('solve_pressure_cg', '6'),
           ('fieldpdf', '7'), ('viscos', '8'), ('gamma_vreman', '9'),
           ('source', '10')]
LABELS = ['momentum_transport', 'scalar_transport', 'momentum_linear',
          'scalar_linear', 'pressure_matrix', 'pressure_linear',
          'tpdf', 'viscosity', 'sgs_vreman', 'momentum_source']

STATE = r'''
module reference_cost
  implicit none
  double precision :: elapsed(10)=0, previous_output=0
  integer :: calls(10)=0
end module
'''

OUTPUT = r'''
subroutine measured_output() bind(C,name='__wrap_output_')
  use reference_cost
  use arrays, only: jo,ko,rho,rhold,ajc,gi,gj,gk,imbflag_real
  use global, only: l,m,n,dtim,istep,floin,flout
  use exchange, only: master
  use mpi
  implicit none
  double precision :: metrics(11),maximum(11),ledger(2),total(2),worst,global_worst
  double precision :: storage,divergence,volume,relative,begin_clock
  integer :: i,j,k,c,status,counts(10),minimum_counts(10)
  interface
    subroutine original_output() bind(C,name='__real_output_')
    end subroutine
  end interface
  begin_clock=MPI_Wtime()
  metrics(1)=0
  if(previous_output>0) metrics(1)=begin_clock-previous_output
  metrics(2:11)=elapsed
  ledger=0; worst=0
  do k=2,n
    do j=2,m
      do i=2,l
        c=i+jo(j)+ko(k)
        if(allocated(imbflag_real)) then
          if(imbflag_real(c)/=1.0) cycle
        endif
        volume=dble(abs(ajc(c)))
        storage=(dble(rho(c))-dble(rhold(c)))/dble(dtim)*volume
        divergence=dble(gi(c+1))-dble(gi(c)) &
          +dble(gj(c+jo(2)))-dble(gj(c)) &
          +dble(gk(c+ko(2)))-dble(gk(c))
        ledger(1)=ledger(1)+storage
        ledger(2)=ledger(2)+divergence
        worst=max(worst,abs(storage+divergence)/(volume*dble(rho(c))/dble(dtim)))
      enddo
    enddo
  enddo
  call MPI_Allreduce(metrics,maximum,11,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
  if(status/=MPI_SUCCESS) stop 9301
  call MPI_Allreduce(calls,counts,10,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,status)
  if(status/=MPI_SUCCESS) stop 9302
  call MPI_Allreduce(calls,minimum_counts,10,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,status)
  if(status/=MPI_SUCCESS .or. any(counts/=minimum_counts)) stop 9303
  call MPI_Allreduce(ledger,total,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
  if(status/=MPI_SUCCESS) stop 9304
  call MPI_Allreduce(worst,global_worst,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
  if(status/=MPI_SUCCESS) stop 9305
  relative=abs(sum(total))/max(abs(total(1)),abs(dble(floin))+abs(dble(flout)),1d-30)
  if(master) write(*,'(A,I8,11(1X,ES22.14E3),10I6,4(1X,ES22.14E3))') &
    'COAST_COST ',istep,maximum,counts,total,relative,global_worst
  call original_output()
  elapsed=0; calls=0
  previous_output=MPI_Wtime()
end subroutine
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wrapper(name, slot):
    imports = observation = ''
    if name == 'solve_pressure_cg':
        imports = '''  use coast_screen_summary, only: coast_pressure_solver_iterations, &
    coast_pressure_solver_initial_residual,coast_pressure_solver_final_residual, &
    coast_pressure_solver_backend,coast_pressure_solver_fallback
  use exchange, only: master
'''
        observation = '''  if(master) write(*,'(A,3I8,2(1X,ES22.14E3),1X,A,1X,L1)') &
    'COAST_PRESSURE_CALL ',istep,calls(slot),coast_pressure_solver_iterations, &
    coast_pressure_solver_initial_residual,coast_pressure_solver_final_residual, &
    trim(coast_pressure_solver_backend),coast_pressure_solver_fallback
'''
    elif name == 'cgstab':
        imports = '  use arrays, only: ncycl\n  use exchange, only: master\n'
        observation = '''  if(master) write(*,'(A,3I8)') 'COAST_TRANSPORT_CALL ',istep,nv,ncycl(nv)
'''
    return '''
subroutine measured_{name}() bind(C,name='__wrap_{name}_')
  use reference_cost
  use global, only: nv,nf,istep
{imports}\
  use mpi
  implicit none
  double precision :: begin_clock
  integer :: slot
  interface
    subroutine original_{name}() bind(C,name='__real_{name}_')
    end subroutine
  end interface
  slot={slot}
  begin_clock=MPI_Wtime()
  call original_{name}()
  elapsed(slot)=elapsed(slot)+MPI_Wtime()-begin_clock
  calls(slot)=calls(slot)+1
{observation}\
end subroutine
'''.format(name=name, slot=slot,imports=imports,observation=observation)


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
    source, obj, binary = out/'cost.f90', out/'cost.o', out/'ccost'
    source.write_text(STATE+''.join(wrapper(*k) for k in KERNELS)+OUTPUT)
    compile_cmd = ['mpif90', '-O2', '-I'+reference['object_build'], '-c',
                   str(source), '-o', str(obj)]
    subprocess.check_call(compile_cmd, cwd=str(out))
    link = reference['link'][:]
    link[-1] = str(binary)
    link[2:2] = [str(obj)]+['-Wl,--wrap='+name+'_' for name, _ in KERNELS]+[
        '-Wl,--wrap=output_']
    subprocess.check_call(link, cwd=str(out))
    report = {'scope': __doc__, 'columns': ['step', 'max_rank_interval_s']+
        [name+'_s' for name in LABELS]+[name+'_calls' for name in LABELS]+
        ['storage_kg_s', 'flux_divergence_kg_s', 'relative_mass_defect', 'max_cell_rho_dt_defect'],
        'reference_sha256': sha(args.reference), 'source_sha256': sha(source),
        'executable_sha256': sha(binary), 'compile': compile_cmd, 'link': link}
    (out/'ccost.json').write_text(json.dumps(report, indent=2)+'\n')
    print(str(binary))


if __name__ == '__main__':
    main()
