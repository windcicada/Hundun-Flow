#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Decode frozen COAST input.d layouts into a traceable migration catalog."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

# These readers use nf=6. Hash admission binds record order, including the
# optional runtime regrid record and statistics_start_step fallback.
READERS = {
    '232f098855fa07031c7a28cfe33fc854f9f17838cc4b00ff2dfb495d330c5035': 'dyn711',
    'edd0671b800a670e436b0c00ae528cbe3828a00e9e712dd9f255ad8fd9b983f6': 'gtmc',
    'a755ae994ea7d740cafd617d152055b3283d1150ce41af1e63aa5c8193bba92a': '624cf',
}

class InputError(ValueError):
    pass

def tokens(record):
    """List-directed scalar subset, with quoted strings and / or ! comments.

    Null values and repetition syntax require a separately defined source
    default; reject them instead of inventing an initialized value.
    """
    result=[];i=0;after_comma=False
    while i<len(record):
        if record[i].isspace():i+=1;continue
        if record[i] in '/!':break
        if record[i]==',':
            if not result or after_comma:raise InputError('null list-directed value')
            after_comma=True;i+=1;continue
        if record[i] in "\"'":
            quote=record[i];i+=1;value='';closed=False
            while i<len(record):
                c=record[i];i+=1
                if c==quote:
                    if i<len(record) and record[i]==quote:value+=quote;i+=1
                    else:closed=True;break
                else:value+=c
            if not closed:raise InputError('unterminated quoted value')
            if i<len(record) and not record[i].isspace() and record[i] not in ',/!':
                raise InputError('missing separator after quoted value')
        else:
            start=i
            while i<len(record) and not record[i].isspace() and record[i] not in ',/!':i+=1
            value=record[start:i]
            if '*' in value:raise InputError('repeated list-directed value requires explicit expansion')
        result.append(value);after_comma=False
    return result

def logical(value):
    # The supported profiles use these conventional list-directed spellings.
    key=value.lower()
    if key in ('t','true','.t.','.true.'):return True
    if key in ('f','false','.f.','.false.'):return False
    raise InputError('expected logical value, got '+repr(value))

def integer(value):
    if not re.fullmatch(r'[+-]?\d+',value):raise InputError('expected integer, got '+repr(value))
    result=int(value)
    if not -2147483648<=result<=2147483647:raise InputError('source INTEGER range exceeded')
    return result

def real(value):
    try:result=float(value.replace('D','e').replace('d','e'))
    except ValueError as error:raise InputError('expected real, got '+repr(value)) from error
    if not math.isfinite(result):raise InputError('source real must be finite')
    return result

def comment(record):
    text=record.lstrip()
    return not text or text[0] in '!#' or bool(re.match(r'[Cc](?:\s|[-=.])',text))

class Records:
    def __init__(self,text):
        self.lines=text.splitlines();self.at=0;self.fields=[]
    def skip(self,n):
        if self.at+n>len(self.lines):raise InputError('truncated input at line '+str(self.at+1))
        self.at+=n
    def read(self,name,converters,units):
        if self.at>=len(self.lines):raise InputError(name+': missing line '+str(self.at+1))
        number=self.at+1;raw=self.lines[self.at];self.at+=1
        try:
            parts=tokens(raw)
            if len(parts)<len(converters):raise InputError('too few values')
            values=[convert(value) for convert,value in zip(converters,parts)]
        except InputError as error:raise InputError(f'{name}, line {number}: {error}') from error
        # Fortran discards trailing values on an ordinary list-directed read.
        self.fields.append(dict(name=name,line=number,raw=raw,values=values,units=units,
                                trailing_tokens=parts[len(converters):]))
        return values[0] if len(values)==1 else values
    def path(self,name):
        if self.at>=len(self.lines):raise InputError(name+': missing path')
        raw=self.lines[self.at];self.at+=1
        if not raw.strip() or len(raw)>200:raise InputError(name+': source path requires 1..200 characters')
        self.fields.append(dict(name=name,line=self.at,raw=raw,values=[raw.rstrip()],units=['path_suffix'],trailing_tokens=[]))
        return raw.rstrip()

def parse_text(text,profile):
    if profile not in READERS.values():raise InputError('unregistered source profile')
    r=Records(text);v={};runtime=profile!='dyn711'
    def read(name,converters,units):
        value=r.read(name,converters,units);v[name]=value;return value
    r.skip(7)
    for name in ('geometry_directory','restart_directory','statistics_directory'):v[name]=r.path(name)
    r.skip(3)
    read('formatted_io',[logical],['boolean'])
    read('restart',[logical,logical],['read','write'])
    if runtime and r.at<len(r.lines) and 'restart_regrid_replace_mesh' in r.lines[r.at]:
        read('restart_regrid_replace_mesh',[logical],['boolean'])
    else:
        v['restart_regrid_replace_mesh']=True if runtime else None
    read('steps',[integer],['source_lstep'])
    read('restart_interval',[integer],['steps'])
    read('output_interval',[integer],['steps'])
    r.skip(3)
    read('sgs_model',[str],['source_model_name'])
    read('sgs_restart',[logical],['boolean'])
    read('smagorinsky_coefficient',[real],['dimensionless'])
    # Runtime readers retry with two booleans if the optional integer fails.
    if runtime:
        save=r.at
        try:statistics=r.read('statistics',[logical,logical,integer],['collect','read','start_step'])
        except InputError:
            r.at=save;statistics=r.read('statistics',[logical,logical],['collect','read'])+[0]
        if statistics[2]<0:raise InputError('statistics start step must be nonnegative')
        v['statistics']=statistics
    else:read('statistics',[logical,logical],['collect','read'])
    read('legacy_inflow_flag',[logical],['boolean'])
    r.skip(3)
    read('outer_iterations',[integer],['iterations'])
    read('end_time',[real],['s'])
    read('screen_skip',[integer],['steps'])
    r.skip(3)
    read('cfl_window',[real,real,real],['minimum','maximum','target'])
    r.skip(3)
    read('wall_time_limit',[real],['hours'])
    r.skip(3)
    read('reference_properties',[real,real,real,real],['kg/m3','kg/m3','Pa s','Pa'])
    r.skip(3)
    for equation in ('u','v','w','pressure_correction','mixture_fraction','enthalpy'):
        read('equation_'+equation,[integer,real],['iterations_or_enabled','source_residual_tolerance'])
    r.skip(3)
    read('fuel',[str],['source_fuel_name'])
    read('mechanism',[str],['source_mechanism_name'])
    r.skip(1)
    read('fuel_dilution',[real,real],['CO2_volume_fraction','H2_volume_fraction'])
    r.skip(1)
    read('combustion',[logical,logical,str],['enabled','chemkin_format','integrator'])
    r.skip(3)
    read('pdf_enabled',[logical],['boolean'])
    read('pdf_fields',[integer],['fields'])
    read('pdf_restart',[logical,logical],['read','write'])
    read('noise_reduction',[logical],['boolean'])
    read('equation_pdf',[integer,real],['iterations','source_residual_tolerance'])
    for name,converters,units in (
        ('species_output',[str],['mass_or_mole_fraction']),
        ('read_heat_release',[logical],['boolean']),
        ('radiation',[logical],['boolean']),
        ('ignition',[logical],['boolean']),
        ('phase_averaging',[logical,logical],['collect','read'])):
        r.skip(3);read(name,converters,units)
    read('phase_frequency',[real,real],['Hz','start_time_s'])
    r.skip(3);read('compressible',[logical],['boolean'])
    r.skip(3);read('effective_boundary',[real,real],['sigma','length_m'])
    consumed=r.at
    extensions=[];group=None;group_lines=[];start=0
    for number,raw in enumerate(r.lines[consumed:],consumed+1):
        stripped=raw.strip()
        if group is not None:
            group_lines.append(raw)
            if stripped=='/' or stripped.lower() in ('&end','$end'):
                extensions.append(dict(kind='namelist',name=group,line=start,end_line=number,
                                       raw='\n'.join(group_lines),status='source_record_requires_model_binding'))
                group=None
            continue
        if comment(raw):continue
        match=re.fullmatch(r'[&$]([A-Za-z][A-Za-z0-9_]*)',stripped)
        if match:
            group=match[1].lower();start=number;group_lines=[raw];continue
        if stripped.startswith(('&','$')):
            raise InputError(f'namelist header requires its own record, line {number}')
        parts=tokens(raw)
        if parts:extensions.append(dict(kind='tokens',name=parts[0],line=number,raw=raw,
                                       values=parts[1:],status='source_record_requires_model_binding'))
    if group is not None:raise InputError(f'unterminated namelist {group}, line {start}')
    low,high,target=v['cfl_window']
    if not 0<low<=target<=high:raise InputError('invalid CFL minimum/target/maximum')
    if v['reference_properties'][3]<=0:raise InputError('thermodynamic reference pressure must be positive')
    if v['steps']<=0 or v['restart_interval']<0 or v['output_interval']<0:
        raise InputError('invalid run step counts')
    if v['pdf_enabled'] and (v['pdf_fields']<2 or v['pdf_fields']%2):
        raise InputError('paired PDF requires an even field count')
    lines={field['name']:field['line'] for field in r.fields}
    defaults=[]
    if runtime and 'restart_regrid_replace_mesh' not in lines:
        defaults.append(dict(name='restart_regrid_replace_mesh',value=True,
                             reason='registered runtime reader default for absent optional record'))
    statistics_field=next(field for field in r.fields if field['name']=='statistics')
    if runtime and len(statistics_field['values'])==2:
        defaults.append(dict(name='statistics_start_step',value=0,line=statistics_field['line'],
                             reason='registered runtime reader two-value fallback'))
    mappings=[
        dict(source='reference_properties[3],compressible',source_lines=[lines['reference_properties'],lines['compressible']],
             target='flow.thermodynamic_pressure_pa',value=None if v['compressible'] else v['reference_properties'][3],
             status='defined',semantics='null selects coupled EOS; positive p0 selects fixed thermodynamic pressure'),
        dict(source='cfl_window',source_lines=[lines['cfl_window']],target='time.convective_cfl_definition',
             value='directional_max',status='defined'),
        dict(source='cfl_window[2]',source_lines=[lines['cfl_window']],target='time.convective_cfl',value=target,status='defined'),
        dict(source='cfl_window',source_lines=[lines['cfl_window']],target='time.convective_cfl_margin',
             value=high-target if math.isclose(high-target,target-low,rel_tol=1e-12,abs_tol=1e-15) else None,
             status='defined' if math.isclose(high-target,target-low,rel_tol=1e-12,abs_tol=1e-15) else 'asymmetric_window_requires_policy'),
    ]
    for name,target in (('sgs_model','turbulence.model'),('pdf_fields','reaction.ensemble.fields'),
                        ('fuel','reaction assets'),('mechanism','reaction assets'),
                        ('noise_reduction','ESF dual-state closure'),('outer_iterations','solver reference_outer_iterations')):
        mappings.append(dict(source=name,source_lines=[lines[name]],target=target,value=v[name],
                             status='requires_native_model_policy'))
    return dict(schema='hundun_coast_input_catalog_v1',profile=profile,nf=6,
                numeric_representation='source decimal values; runtime REAL precision requires build manifest',
                core_end_line=consumed,source_line_count=len(r.lines),settings=v,fields=r.fields,
                extensions=extensions,mappings=mappings,reader_defaults=defaults,
                required_assets=['grid and decomposition','physical boundary definitions','thermophysics and mechanism',
                                 'synchronous restart and PDF','statistics and model histories'],
                conversion_status='source_catalog_ready_native_case_binding_pending')

def catalog(source,reader):
    raw=source.read_bytes();code=reader.read_bytes();sha=hashlib.sha256(code).hexdigest()
    if sha not in READERS:raise InputError('unregistered input reader SHA-256: '+sha)
    result=parse_text(raw.decode('utf-8'),READERS[sha])
    result['source']={'path':str(source.resolve()),'sha256':hashlib.sha256(raw).hexdigest()}
    result['reader']={'path':str(reader.resolve()),'sha256':sha}
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--reader',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    try:
        if args.output.resolve() in (args.source.resolve(),args.reader.resolve()):
            raise InputError('output must have its own path')
        result=catalog(args.source,args.reader)
        encoded=json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n'
        args.output.parent.mkdir(parents=True,exist_ok=True)
        temporary=args.output.with_name(args.output.name+'.tmp')
        with temporary.open('x',encoding='utf-8') as stream:stream.write(encoded)
        temporary.replace(args.output)
    except (InputError,OSError,UnicodeError) as error:parser.exit(2,str(error)+'\n')
    print(f"{result['profile']}: {len(result['fields'])} source fields, {len(result['extensions'])} extension records; {result['conversion_status']}")

if __name__=='__main__':main()
