#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Audit complete original COAST condif/step rows against scalar identities.

Links an external driver with the fresh original objects, replacing only main.
Scope: all-fluid Cartesian interior, central faces, constant/variable diffusivity, BE.
The native row probe receives the reference face conductance explicitly; its
material interpolation policy is a separate operator comparison.
The production IBM, boundary and VLS operators have separate reference entries.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


DRIVER = r'''
program scalar_rows
  use arrays
  use global
  use imb_config, only: ibm_config_present,ibm_immersed
  implicit none
  integer :: i,j,k,c,t,sgn,constant,rr
  real :: residual
  lower=0; upper=511; l=3; m=3; n=3; imax=8; jmax=8; kmax=8
  nv=7; nf=6; dtim=0.125; big=huge(1.0); small=tiny(1.0)
  ibm_config_present=.false.; ibm_immersed=.false.
  allocate(io(0:8),jo(0:8),ko(0:8),nfo(7),tvd(7))
  do i=0,8
    io(i)=i-1; jo(i)=8*(i-1); ko(i)=64*(i-1)
  enddo
  nfo=0; tvd=.false.
  allocate(f(lower:upper),fold(lower:upper),rho(lower:upper), &
           ajc(lower:upper),gam(lower:upper),coef(8,lower:upper))
  allocate(gi(lower:upper),gj(lower:upper),gk(lower:upper), &
           w1(lower:upper),w2(lower:upper),w3(lower:upper))
  allocate(b11(lower:upper),b12(lower:upper),b13(lower:upper), &
           b21(lower:upper),b22(lower:upper),b23(lower:upper), &
           b31(lower:upper),b32(lower:upper),b33(lower:upper))
  allocate(gam_n(8),gam_e(8),gam_r(8))
  b12=0; b13=0; b21=0; b23=0; b31=0; b32=0; gam=0.125
  t=0
  do rr=1,3
    do constant=0,1
      do sgn=-1,1,2
        t=t+1
        do k=1,8
          do j=1,8
            do i=1,8
              c=i+jo(j)+ko(k)
              if (c > upper) cycle
              f(c)=0.25+(i+2*j+3*k)/64.0
              fold(c)=f(c)-0.03125
              if(constant==1) then
                f(c)=0.375; fold(c)=f(c)
              endif
              rho(c)=1.0+mod(i+j+k,3)/8.0
              ajc(c)=1.0+mod(i+j+k,2)/8.0
              gi(c)=sgn*(1.0+i/8.0)
              gj(c)=sgn*(-0.5+j/16.0)
              gk(c)=sgn*(0.25+k/32.0)
              w1(c)=0.25+mod(i,3)/16.0
              w2(c)=0.5+mod(j,3)/16.0
              w3(c)=0.375+mod(k,3)/16.0
              b11(c)=1.0+i/16.0
              b22(c)=1.0+j/32.0
              b33(c)=1.0+k/64.0
            enddo
          enddo
        enddo
        call condif
        call step(rho)
        do k=2,n
          do j=2,m
            do i=2,l
              c=i+jo(j)+ko(k)
              residual=coef(pc,c)*f(c)-coef(bpc,c) &
                -coef(nc,c)*f(c+1)-coef(sc,c)*f(c-1) &
                -coef(ec,c)*f(c+8)-coef(wc,c)*f(c-8) &
                -coef(rc,c)*f(c+64)-coef(lc,c)*f(c-64)
              write(*,'(4I5,9(1X,ES25.17E3))') t,i,j,k,coef(:,c),residual
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
end program
'''


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expected(t, i, j, k):
    constant = ((t-1)//2) % 2
    sign = -1 if t % 2 else 1
    ratio = ((t-1)//4 % 3+1)/2.0
    variable = t > 12
    gamma = lambda p: 0.125*(1+((p[0]+2*p[1]+3*p[2]) % 5)/8.0) if variable else 0.125
    position = [i, j, k]
    q = lambda p: 0.375 if constant else 0.25+(p[0]+2*p[1]+3*p[2])/64.0
    value = q(position)
    old = value if constant else value-0.03125
    density = 1+(i+j+k) % 3/8.0
    volume = 1+(i+j+k) % 2/8.0
    coefficients, neighbours, faces = [], [], []
    flux_sum = convection = diffusion = 0.0
    for axis in range(3):
        # COAST row order: plus, minus for each coordinate.
        for side in [1, -1]:
            p = position[:]
            p[axis] += side
            face = position[axis] if side == 1 else position[axis]-1
            flux_coordinate = position[axis]+1 if side == 1 else position[axis]
            flux = sign*([1, -0.5, 0.25][axis] +
                         flux_coordinate/[8.0, 16.0, 32.0][axis])
            weight = [0.25, 0.5, 0.375][axis]+face % 3/16.0
            face_gamma = (weight*gamma(position)+(1-weight)*gamma(p) if side==1 else
                          weight*gamma(p)+(1-weight)*gamma(position))
            transmissibility = face_gamma*(1+face/[16.0, 32.0, 64.0][axis])
            adjacent = q(p)
            face_value = weight*value+(1-weight)*adjacent if side == 1 else (
                weight*adjacent+(1-weight)*value)
            coefficient = (transmissibility-(1-weight)*flux if side == 1 else
                           transmissibility+weight*flux)/volume
            coefficients.append(coefficient)
            neighbours.append(adjacent)
            faces.append([transmissibility,weight,flux])
            flux_sum += side*flux/volume
            convection += side*flux*face_value/volume
            diffusion += transmissibility*(adjacent-value)/volume
    diagonal = sum(coefficients)+density/0.125
    rhs = density/0.125*old
    row_residual = diagonal*value-sum(a*b for a, b in zip(coefficients, neighbours))-rhs
    continuity = (ratio*density-density)/0.125+flux_sum
    conservative = (ratio*density*value-density*old)/0.125+convection-diffusion
    transformed = conservative-value*continuity
    # Native helper stores minus, plus for each axis.
    native_faces = [faces[i ^ 1] for i in range(6)]
    native = [volume,density,old,0.125]+[v for face in native_faces for v in face]
    return coefficients+[diagonal, rhs, row_residual], transformed, conservative, continuity, native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True,
                        help='fresh cref.json with complete object provenance')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--hundun', type=Path, required=True,
                        help='native species conservation test executable with --row')
    parser.add_argument('--variable-gamma', action='store_true',
                        help='append variable-material rows using original arithmetic face weights')
    args = parser.parse_args()
    reference = json.loads(args.reference.read_text())
    out = args.output.resolve()
    source = Path(reference['coast_source']).resolve()
    if out == source.parent or source.parent in out.parents:
        parser.error('output must be external to COAST')
    for name, digest in reference['objects'].items():
        if sha(Path(name)) != digest:
            parser.error('reference object changed: '+name)
    out.mkdir(parents=True, exist_ok=True)
    driver, obj, binary = out/'srow.f90', out/'srow.o', out/'srow'
    source_text = DRIVER
    if args.variable_gamma:
        source_text = source_text.replace('do rr=1,3', 'do rr=1,6').replace(
            '              rho(c)=',
            '              if (rr>3) gam(c)=0.125*(1.0+mod(i+2*j+3*k,5)/8.0)\n              rho(c)=')
    driver.write_text(source_text)
    compile_cmd = ['mpif90', '-O2', '-I'+reference['object_build'], '-c',
                   str(driver), '-o', str(obj)]
    subprocess.check_call(compile_cmd, cwd=str(out))
    link = []
    original = reference['link']
    for item in original:
        if item.endswith('/config_coast.o') or item.endswith('/pwrap.o') or '--wrap=' in item:
            continue
        link.append(item)
    link[-1] = str(binary)
    link.insert(2, str(obj))
    subprocess.check_call(link, cwd=str(out))
    trace = subprocess.check_output([str(binary)], cwd=str(out), universal_newlines=True)
    (out/'rows.dat').write_text(trace)
    rows = [line.split() for line in trace.splitlines() if line.strip()]
    expected_rows = 192 if args.variable_gamma else 96
    if len(rows) != expected_rows:
        raise ValueError('expected {} complete interior rows'.format(expected_rows))
    native_input = ''.join(' '.join(str(v) for v in expected(
        *[int(v) for v in row[:4]])[-1])+'\n' for row in rows)
    native_trace = subprocess.check_output([str(args.hundun.resolve()), '--row'],
        input=native_input, universal_newlines=True)
    (out/'native.dat').write_text(native_trace)
    native_rows = [[float(v) for v in line.split()] for line in native_trace.splitlines()]
    if len(native_rows) != len(rows) or any(len(row) != 8 for row in native_rows):
        raise ValueError('native probe must return eight coefficients per row')
    operator_error = identity_error = conservative_gap = constant_error = 0.0
    native_error = cross_error = 0.0
    for row, native_row in zip(rows,native_rows):
        key = [int(v) for v in row[:4]]
        actual = [float(v) for v in row[4:]]
        target, transformed, conservative, continuity, _ = expected(*key)
        if len(actual) != 9:
            raise ValueError('original probe must return eight coefficients and one residual')
        scale = max(1.0, max(abs(v) for v in target))
        native_error = max(native_error, max(abs(a-b)/scale for a,b in zip(native_row,target[:8])))
        cross_error = max(cross_error, max(abs(a-b)/scale for a,b in zip(native_row,actual[:8])))
        operator_error = max(operator_error, max(abs(a-b)/scale for a, b in zip(actual, target)))
        identity_error = max(identity_error, abs(target[-1]-transformed)/scale)
        conservative_gap = max(conservative_gap, abs(conservative-target[-1]))
        if ((key[0]-1)//2) % 2:
            constant_error = max(constant_error, abs(actual[-1])/scale)
    report = {'scope': __doc__, 'rows': len(rows), 'precision': 'original FP32',
              'variable_gamma': args.variable_gamma,
              'face_material_policy': 'original arithmetic weighted interpolation',
              'native_probe_scope': 'row assembly with explicitly supplied reference face conductance',
              'normalized_operator_error': operator_error,
              'normalized_row_identity_error': identity_error,
              'normalized_constant_error': constant_error,
              'normalized_native_FP64_error': native_error,
              'normalized_native_coast_FP32_difference': cross_error,
              'conservative_residual_gap': conservative_gap,
              'reference_sha256': sha(args.reference), 'driver_sha256': sha(driver),
              'executable_sha256': sha(binary), 'native_sha256': sha(args.hundun),
              'compile': compile_cmd, 'link': link}
    report['passed'] = (operator_error < 8e-7 and identity_error < 1e-14 and
                        constant_error < 8e-7 and conservative_gap > 0.1 and
                        native_error < 1e-11 and cross_error < 8e-7)
    (out/'srow.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ['passed', 'rows', 'normalized_operator_error',
          'normalized_row_identity_error', 'normalized_constant_error', 'conservative_residual_gap',
          'normalized_native_FP64_error','normalized_native_coast_FP32_difference']}))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
