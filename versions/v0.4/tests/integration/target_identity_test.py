#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Configure and incremental reconfigure of the production identity functions."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmake', required=True)
    parser.add_argument('--module', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='hundun-identity-') as directory:
        root = Path(directory)
        (root / 'versions/v0.4/src').mkdir(parents=True)
        (root / 'versions/v0.4/src/model.cpp').write_text('int model;\n')
        (root / 'versions/v0.4/CMakeLists.txt').write_text('# product inputs\n')
        (root / 'core.c').write_text('int value;\n')
        (root / 'runner.cpp').write_text('int main() {}\n')
        (root / 'CMakeLists.txt').write_text('''cmake_minimum_required(VERSION 3.21)
project(identity_fixture NONE)
include("%s")
hundun_v04_identity_inputs("${CMAKE_CURRENT_SOURCE_DIR}" collected_inputs)
hundun_source_content_digest("${CMAKE_CURRENT_SOURCE_DIR}" collected ${collected_inputs})
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/collected" "${collected}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/core.c" "${CMAKE_CURRENT_SOURCE_DIR}/runner.cpp")
hundun_source_content_digest("${CMAKE_CURRENT_SOURCE_DIR}" content "${CMAKE_CURRENT_SOURCE_DIR}/core.c")
set(HUNDUN_BUILD_MANIFEST_CANONICAL "schema=V2\\ncore=${content}\\nc_flags=${C_FLAGS}\\ntests=${TESTS}\\nhypre=${HYPRE}\\nconfig=${CONFIG}\\n")
hundun_target_manifest(core none core)
file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/runner.cpp" entry)
hundun_target_manifest(runner "${entry}" runner)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ids" "${core}\\n${runner}\\n")
add_custom_target(identity_only ALL)
''' % args.module.resolve())
        build = root / 'build'

        def configure(*flags):
            subprocess.check_call([args.cmake, '-S', str(root), '-B', str(build)] + list(flags),
                                  stdout=subprocess.DEVNULL)
            return (build / 'ids').read_text().splitlines()

        baseline = configure('-DC_FLAGS=-O2', '-DTESTS=OFF', '-DHYPRE=OFF', '-DCONFIG=Release')
        assert len(baseline) == 2 and baseline[0] != baseline[1]
        content = (build / 'collected').read_text()
        scratch = root / 'another-build/CMakeFiles/CMakeScratch/transient'
        scratch.mkdir(parents=True)
        (scratch / 'CMakeLists.txt').write_text('throwaway compiler probe\n')
        assert configure() == baseline
        assert (build / 'collected').read_text() == content, 'build artifact entered source identity'
        (scratch / 'CMakeLists.txt').unlink()
        assert configure() == baseline
        assert (build / 'collected').read_text() == content
        (root / 'versions/v0.4/src/model.cpp').write_text('int model = 1;\n')
        configure()
        assert (build / 'collected').read_text() != content, 'real product input omitted'

        time.sleep(1.1)  # Makefile timestamp comparisons may have one-second resolution.
        (root / 'runner.cpp').write_text('int main() { return 0; }\n')
        subprocess.check_call([args.cmake, '--build', str(build)], stdout=subprocess.DEVNULL)
        runner = (build / 'ids').read_text().splitlines()
        assert runner[0] == baseline[0] and runner[1] != baseline[1], (baseline, runner)
        time.sleep(1.1)
        (root / 'core.c').write_text('int value = 1;\n')
        subprocess.check_call([args.cmake, '--build', str(build)], stdout=subprocess.DEVNULL)
        changed_c = (build / 'ids').read_text().splitlines()
        assert all(a != b for a, b in zip(runner, changed_c))
        previous = changed_c
        for flag in ('-DC_FLAGS=-O3', '-DTESTS=ON', '-DHYPRE=ON', '-DCONFIG=Debug'):
            current = configure(flag)
            assert all(a != b for a, b in zip(previous, current)), flag
            previous = current
        print('target identity: C-only, runner-only incremental, C flags, tests, HYPRE, config PASS')


if __name__ == '__main__':
    main()
