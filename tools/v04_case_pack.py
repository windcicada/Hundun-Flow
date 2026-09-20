#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Bundle a native case/restart, source snapshot and matching Linux MPI runtime.

Run packaging in the program's build environment. This profile captures
x86-64 glibc, OpenMPI 4, its dlopen components, and the admitted Cantera SDK.
The relocated package starts through its private loader. Source rebuilding
uses the recorded Clang/libstdc++ ABI1 Release/ThinLTO toolchain contract.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import tempfile
import tarfile


def require(ok, message):
    if not ok:
        raise ValueError(message)


def command(args):
    r = subprocess.run([str(x) for x in args], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True, timeout=60)
    require(r.returncode == 0, '{}: {}'.format(args[0], r.stderr[-2000:]))
    return r.stdout


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def copy_tree(source, target):
    """Keep input hierarchy and only symlinks whose targets stay in the tree."""
    source = Path(source).resolve()
    require(source.is_dir(), 'input directory: '+str(source))
    for p in source.rglob('*'):
        if p.is_symlink():
            resolved = p.resolve()
            require(source in resolved.parents and resolved.exists(), 'external or broken input link: '+str(p))
        else:
            require(p.is_dir() or p.is_file(), 'special input file: '+str(p))
    shutil.copytree(str(source), str(target), symlinks=True)


def file_manifest(root):
    result = {}
    for p in sorted(root.rglob('*')):
        if p.is_symlink():
            result[str(p.relative_to(root))] = dict(link=os.readlink(str(p)))
        elif p.is_file():
            result[str(p.relative_to(root))] = dict(bytes=p.stat().st_size, mode=p.stat().st_mode & 0o777, sha256=sha(p))
    return result


def verify(root):
    root = Path(root).resolve()
    report = json.loads((root/'manifest.json').read_text())
    require(report['format'] == 'hundun_native_case_package_v1', 'package manifest format')
    for name, expected in report['files'].items():
        relative = Path(name)
        require(not relative.is_absolute() and '..' not in relative.parts, 'relative package path')
        p = root/relative
        require(root in p.resolve().parents, 'package path stays in root')
        if 'link' in expected:
            require(p.is_symlink() and os.readlink(str(p)) == expected['link'], 'package link: '+name)
        else:
            require(p.is_file() and not p.is_symlink(), 'package file: '+name)
            require(p.stat().st_size == expected['bytes'] and sha(p) == expected['sha256'] and
                    p.stat().st_mode & 0o777 == expected['mode'], 'package contents/mode: '+name)
    return len(report['files'])


COMMON = '''export OPAL_PREFIX="$root"
export OMPI_MCA_mca_base_component_path="$root/mpi"
export PMIX_MCA_mca_base_component_path="$root/pmix"
export HWLOC_PLUGINS_PATH="$root/hwloc"
export CANTERA_DATA="$root/ct/share/cantera/data"
export OMPI_MCA_pml="${OMPI_MCA_pml:-ob1}"
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,vader,tcp}"
export OMPI_MCA_osc="${OMPI_MCA_osc:-pt2pt}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
'''


def launcher(root, name, executable, parent=False, mpi=False):
    base = '$(dirname -- "$0")'+('/..' if parent else '')
    body = '#!/bin/sh\nset -eu\nroot=$(CDPATH= cd -- "'+base+'" && pwd)\n'+COMMON
    # OpenMPI prefix injection otherwise puts the bundled libc into an older
    # host /bin/sh before the wrapper can select its matching private loader.
    option = ' --noprefix' if mpi else ''
    body += 'exec "$root/lib/ld-linux-x86-64.so.2" --library-path "$root/lib" "$root/bin/'+executable+'"'+option+' "$@"\n'
    p = root/name
    p.write_text(body)
    p.chmod(0o755)


def pack(args):
    require(platform.machine() == 'x86_64', 'runtime profile requires Linux x86-64')
    output = Path(args.output).resolve()
    require(not output.exists(), 'output directory already exists')
    case, restart, source, sdk, build = [Path(getattr(args,n)).resolve()
        for n in ('case','restart','source','cantera','build')]
    program = Path(args.program).resolve()
    for p in (case, restart, sdk):
        require(p != output and p not in output.parents, 'output must be separate from copied input')
    require((case/'case.json').is_file() and (restart/'current').is_file(), 'native case/restart entry')
    current = (restart/'current').read_text().strip()
    require(re.fullmatch(r'generation-[0-9]+-[0-9]+',current) is not None and
            (restart/current/'manifest.bin').is_file(), 'native current generation')
    revision = command(['git','-C',source,'rev-parse','HEAD']).strip()
    mpi_version = command(['ompi_info','--version']).splitlines()[0]
    require('Open MPI v4.' in mpi_version, 'runtime profile requires OpenMPI 4')
    paths = {}
    for line in command(['ompi_info','--parsable','--path','all']).splitlines():
        if line.startswith('path:'):
            _, key, value = line.split(':',2)
            paths[key] = Path(value)
    snapshots = {name:file_manifest(p) for name,p in [('case',case),('restart',restart),('ct',sdk)]}
    output.parent.mkdir(parents=True,exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.pack-',dir=str(output.parent)))
    try:
        for name in ('bin','lib','mpi','pmix','hwloc','share','licenses','meta'):
            (stage/name).mkdir()
        queue, loaded, origins = [], set(), {}
        def copy_elf(src, dst):
            src, dst = Path(src), Path(dst)
            digest = sha(src)
            if dst.exists():
                require(sha(dst) == digest, 'runtime basename collision: '+dst.name)
                return
            shutil.copy2(str(src),str(dst))
            origins[str(dst.relative_to(stage))] = dict(path=str(src),sha256=digest)
            queue.append(dst)
        copy_elf(program,stage/'bin/hundun')
        copy_elf(paths['bindir']/'orterun',stage/'bin/mpiexec')
        copy_elf(paths['bindir']/'orted',stage/'bin/orted.bin')
        for src,name in [(paths['pkglibdir'],'mpi'),
                         (Path('/usr/lib/x86_64-linux-gnu/pmix2/lib/pmix'),'pmix'),
                         (Path('/usr/lib/x86_64-linux-gnu/hwloc'),'hwloc')]:
            require(src.is_dir(), 'runtime component directory: '+str(src))
            for p in sorted(src.glob('*.so')):
                copy_elf(p,stage/name/p.name)
        while queue:
            p = queue.pop()
            if p in loaded:
                continue
            loaded.add(p)
            deps = command(['ldd',p])
            require('not found' not in deps, 'unresolved ELF dependency: '+str(p))
            for dep in re.findall(r'(?:=>\s+|^\s*)(/[^\s]+)',deps,re.M):
                copy_elf(dep,stage/'lib'/Path(dep).name)
        require((stage/'lib/ld-linux-x86-64.so.2').is_file(), 'matching ELF loader')
        copy_tree(paths['pkgdatadir'],stage/'share/openmpi')
        copy_tree(case,stage/'case')
        copy_tree(restart,stage/'restart')
        copy_tree(sdk,stage/'ct')
        # License records follow every package-owned runtime component. The
        # Cantera SDK contains the admitted library and upstream notices.
        packages = set()
        for item in origins.values():
            path = item['path']
            r = subprocess.run(['dpkg-query','-S',path],stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,universal_newlines=True)
            if r.returncode == 0:
                for line in r.stdout.splitlines():
                    package = line.split(': ',1)[0].split(':',1)[0]
                    notice = Path('/usr/share/doc')/package/'copyright'
                    if notice.is_file():
                        shutil.copyfile(str(notice),str(stage/'licenses'/(package+'.txt')))
                        packages.add(package)
        for name in ('LICENSE','NOTICE'):
            p = source/name
            if p.is_file():
                shutil.copyfile(str(p),str(stage/name))
        command(['git','-C',source,'archive','--format=tar.gz','--output='+str(stage/'source.tar.gz'),revision])
        # Match the archived core inputs to the actual CMake build manifest.
        # Independent worktree documentation edits stay outside this snapshot.
        source_hashes = {}
        with tarfile.open(str(stage/'source.tar.gz')) as archive:
            for member in archive.getmembers():
                name = member.name
                if member.isfile() and name != 'versions/v0.4/src/app_main.cpp' and (name in ('VERSION','CMakeLists.txt','versions/v0.4/CMakeLists.txt') or
                    any(name.startswith(prefix) for prefix in
                        ('versions/v0.4/include/','versions/v0.4/src/','cmake/','third_party/yyjson/'))):
                    source_hashes[name] = hashlib.sha256(archive.extractfile(member).read()).hexdigest()
        canonical = 'HUNDUN_SOURCE_CONTENT_V1\n'+''.join(
            '{}={}\n'.format(name,source_hashes[name]) for name in sorted(source_hashes))
        source_digest = hashlib.sha256(canonical.encode()).hexdigest()
        core_manifest = (build/'versions/v0.4/generated/core.build-manifest.txt').read_text()
        expected = re.search(r'^core_source_content_sha256=(.+)$',core_manifest,re.M)
        app_manifest = (build/'versions/v0.4/generated/hundun.build-manifest.txt').read_text()
        entry = re.search(r'^entry_sha256=(.+)$',app_manifest,re.M)
        with tarfile.open(str(stage/'source.tar.gz')) as archive:
            app_sha = hashlib.sha256(archive.extractfile('versions/v0.4/src/app_main.cpp').read()).hexdigest()
        require(entry and entry.group(1) == app_sha, 'archived application entry differs from the build manifest')
        require(expected and expected.group(1) == source_digest,
                'archived core differs from the compiled source manifest')
        for name in ('core.build-manifest.txt','hundun.build-manifest.txt'):
            shutil.copyfile(str(build/'versions/v0.4/generated'/name),str(stage/'meta'/name))
        shutil.copyfile(str(build/'CMakeCache.txt'),str(stage/'meta/CMakeCache.txt'))
        launcher(stage,'run','hundun')
        launcher(stage,'mpirun','mpiexec',mpi=True)
        launcher(stage,'bin/orted','orted.bin',parent=True)
        (stage/'build').write_text('''#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
mkdir -p "$root/src"
tar -xzf "$root/source.tar.gz" -C "$root/src"
cmake -S "$root/src" -B "$root/b" -G Ninja \\
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \\
  -DCMAKE_C_FLAGS="-ffp-contract=off -flto=thin" \\
  -DCMAKE_CXX_FLAGS="-stdlib=libstdc++ -D_GLIBCXX_USE_CXX11_ABI=1 -ffp-contract=off -flto=thin" \\
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld -DCMAKE_SHARED_LINKER_FLAGS=-fuse-ld=lld \\
  -DHUNDUN_BUILD_TESTS=OFF -DHUNDUN_ENABLE_REACTING_CANTERA=ON \\
  -DHUNDUN_CANTERA_PACKAGE_ROOT="$root/ct"
cmake --build "$root/b" --parallel "${HUNDUN_BUILD_JOBS:-2}" --target hundun
''')
        (stage/'build').chmod(0o755)
        shutil.copyfile(__file__,str(stage/'pack.py'))
        (stage/'README.md').write_text('''# Hundun native case package

Linux x86-64 runtime with glibc, OpenMPI 4 components and Cantera.
The package selects its private ELF loader and libraries. Keep its folders together.
`case/` and `restart/` preserve the admitted inputs and original method history.

```sh
python3 pack.py --verify .
./run check case
./mpirun -n 2 --bind-to none ./run run case --restart restart --output out --steps 1 --restart-interval 1
./mpirun -n 2 --bind-to none ./run run case --restart out/Restart --output next --steps 1 --restart-interval 1
```

Launch from this directory. Each output directory receives a new run.
A fresh run uses the case boundary state, or an explicit `--initial-state`.
The manifest records source, program, input and dependency identities independently.

For source rebuilding, use CMake >=3.21, Ninja, Clang 15 with lld,
GCC 11 libstdc++ headers/libraries, and OpenMPI 4 development headers.
`./build` uses the supplied Cantera SDK, FP64, ABI1 and Release/ThinLTO.
The rebuilt executable is `b/versions/v0.4/hundun`.
Runtime validation and numerical acceptance are recorded by the producing task.
''')
        for name,p in [('case',case),('restart',restart),('ct',sdk)]:
            require(file_manifest(p) == snapshots[name], 'input changed while packaging: '+name)
            require(file_manifest(stage/name) == snapshots[name], 'copied input differs: '+name)
        for item in origins.values():
            require(sha(item['path']) == item['sha256'], 'runtime input changed while packaging')
        manifest = dict(format='hundun_native_case_package_v1', source_commit=revision,
            core_source_content_sha256=source_digest,
            program_sha256=sha(program), runtime_profile='linux-x86_64-glibc-openmpi4',
            mpi_version=mpi_version, elf_files=len(loaded), runtime_origins=origins,
            debian_license_packages=sorted(packages), files=file_manifest(stage))
        (stage/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
        require(not output.exists(), 'output directory appeared during packaging')
        os.rename(str(stage),str(output))
        return manifest
    finally:
        if stage.exists():
            shutil.rmtree(str(stage))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('program','case','restart','source','cantera','build','output'):
        p.add_argument('--'+name)
    p.add_argument('--verify',metavar='PACKAGE')
    args = p.parse_args()
    if args.verify:
        print('Verified package files:',verify(args.verify))
        return
    require(all(getattr(args,n) for n in ('program','case','restart','source','cantera','build','output')),
            'pack requires program, case, restart, source, cantera, build and output')
    r = pack(args)
    print(json.dumps({k:r[k] for k in ('source_commit','program_sha256','runtime_profile','elf_files')},indent=2))


if __name__ == '__main__':
    main()
