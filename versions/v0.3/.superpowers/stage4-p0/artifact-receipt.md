# Stage 4 P0 Cantera Artifact Receipt

- Recorded: `2026-08-09T23:20:29+08:00`
- Result: `PREFLIGHT_PASS`
- Result scope: `cantera_3_2_0_gcc11_release_artifact_candidate_only`
- Stage 4 product accepted: `false`
- Product changes: `none`
- Governance parent HEAD: `3ef96dde43d3ea1a86aa050fd5f470b9ec80b037`
- Input manifest SHA-256: `32622f5986e97903b9b07110d091c993a3c85e9b3a78ab160493050f38f1361e`
- Artifact manifest SHA-256: `efd2fbbc9f497b7b0f7212104497591464a5e0ae76d84527528e9586f8338296`
- External root: `/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0`
- Candidate artifact: `artifact/cantera-3.2.0-gcc11-release-v4`

This receipt accepts only the external Cantera Linux CPU artifact candidate produced by P0-2.
It does not accept HUNDUN CMake integration, a `ChemistryBackend`, reacting flow, a bundled
mechanism, relocation, or any Stage 4 scientific capability. P0-3 remains responsible for
workspace isolation, thread/MPI lifecycle and moved-prefix runtime evidence.

## Frozen builder and source identity

The candidate is built from Cantera `v3.2.0`, peeled commit
`4a8358eb80cfeb50474386b5f9ec0b3a83519889`, official release archive SHA-256
`a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b`.
The builder profile is Ubuntu 22.04, GCC 11/libstdc++, C++17, ABI=1 and generic x86-64:

```text
python_package=n
f90_interface=n
example_data=n
googletest=none
hdf_support=n
system_blas_lapack=n
system_eigen=n
system_fmt=n
system_sundials=n
system_yamlcpp=n
renamed_shared_libraries=yes
versioned_shared_library=yes
layout=standard
package_build=yes
use_rpath_linkage=n
CXX=g++-11
CC=gcc-11
boost_inc_dir=/p0/build/cantera-3.2.0-gcc11-release/boost-system-include-root
cxx_flags=-std=c++17 -D_GLIBCXX_USE_CXX11_ABI=1 -DEIGEN_MPL2_ONLY -march=x86-64 -mtune=generic
optimize=yes
prefix=/p0/install/cantera-3.2.0-gcc11-release
parallel_jobs=16
local_patches=0
```

| Evidence | SHA-256 / result |
|---|---|
| Outer builder `manifests/p0-2-run-cantera-build-retry8.sh` | `b8c68cb2d623caeb2169f3b116ce1f06d22e2d2948420da9284ec13128fe0ca3` |
| Inner builder `manifests/doxygen-builder-v1/build-retry8-inner.sh` | `67c83fa3abd7fd14d2ab4ac5d10ce81a336457f5a38a8157b63895d3953491bf` |
| Build log | `d9667c01580d9da019933f44893ee30272e19726f3460dc5dd26fb374cd49dce` |
| Build status | exit 0; `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| Build timing | wall 65.04 s; `e6300bb66f08536b01726e2caefa39c9734389fcf16a1794d4c51e8de44f9dc9` |
| Canonical retry8 non-data overlay manifest | `fa450a5ff9d489cf99ac166f7c67ea1734f5b679a38bbf2163cb97cede04f6f4` |

The timing wrapper observed a 5,240 KiB maximum RSS for the outer bubblewrap process. This is not
claimed as the recursive peak of all SCons/compiler children and is retained as a measurement
limitation rather than promoted to a package requirement.

The host network was not disconnected. A command/log policy scan found no executed
`curl`, `wget`, `git clone/fetch/pull`, `pip` or `apt download/install/update` in the canonical
retry8 build, projection or consumer path. Its script SHA-256 is
`17b16377d2e18f9ddef77284420e81bd21bec90ec704351bbb6c45cf2bb14fc4`;
the empty finding log SHA-256 is
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`;
the successful run log SHA-256 is
`1ccc15244dd5d3faf1f362c6c9ecc91ab850b269f9aabb0e2fc11bb4d094e696`.

## Artifact projection

The complete upstream install is retained as build evidence. The user-facing candidate is a
separate projection that removes the static library, pkg-config file, samples, manpages,
Python/Cython/test headers and all mechanisms. It adds the exact generated dependency headers and
separates the pure Cantera BSD license from the upstream combined dependency-license inventory.

| Field | Accepted value |
|---|---|
| Projection script SHA-256 | `c9d3c11c64ddb29c218eb7b5c13f87ca5b102482f96c2b7a2390f2378fefee26` |
| Projection log SHA-256 / exit | `119ad1836baeb7383875d8f70e233e29c8046d4f66df070ba88f8bf2c46a9090` / 0 |
| Regular files / symlinks | 653 / 2 |
| Projected bytes (`du -sb`) | 21,627,465 |
| Public Cantera headers | 194 |
| Generated dependency headers | 442 |
| License files | 15 |
| Data files / mechanisms | 1 / 0 |
| File-manifest SHA-256 | `5fd187220f23d11af133f43f4c1208f9744943bc793e89bc7543171363fa3004` |
| Symlink-manifest SHA-256 | `624d5d9488efc7bf8231d19ca6861e5118adada82b3cee4b1faaa512fa8b0269` |
| Size-manifest SHA-256 | `ddb5507d1dd6f29428aa4b5d1f4444ddf55006fff7414199ec685e5d9d7b3fa7` |

The two symlinks are internal and deterministic:

```text
lib/libcantera_shared.so   -> libcantera_shared.so.3.2.0
lib/libcantera_shared.so.3 -> libcantera_shared.so.3.2.0
```

The artifact contains only upstream `share/cantera/data/README.md`, SHA-256
`067c221d01413542325424c7b9c6708674e67c4b5a346c58dfdea891ea81754f`.
No mechanism inherits Cantera's code license by assumption. P0-3 must use a separately authored,
separately licensed and hash-bound synthetic mechanism.

## Binary and ABI result

| Field | Accepted value |
|---|---|
| Shared library | `lib/libcantera_shared.so.3.2.0` |
| SHA-256 / bytes | `093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760` / 10,164,736 |
| Unstripped SHA-256 / bytes | `54cdd1ffda2f1e0ba42c683c6c795d7dd1f96b40313fad5dbc94abf54589a38b` / 229,858,632 |
| Strip operation | Ubuntu binutils 2.38, `--strip-debug` |
| Build ID before and after strip | `ff123ec43dc02b7a2452eab03f049d8a4347a0ba` |
| SONAME | `libcantera_shared.so.3` |
| NEEDED | `libc.so.6`, `libgcc_s.so.1`, `libm.so.6`, `libstdc++.so.6` |
| RPATH / RUNPATH | absent / absent |
| Maximum GLIBC / GLIBCXX / CXXABI | `GLIBC_2.34` / `GLIBCXX_3.4.29` / `CXXABI_1.3.9` |
| ABI=1 evidence | 4,388 `std::__cxx11` symbols |
| Exceptions / RTTI evidence | 4 `__cxa_throw`; 1,218 typeinfo symbols |
| Jammy-rootfs `ldd -r` | exit 0; no missing library or unresolved symbol |

The accepted audit script SHA-256 is
`e677358c7e24c210764ef5c93d4a54d5519a3abce1aae2045a9a90f7bf5f95e8`;
its log SHA-256 is
`d964cdd5ada3b6222c7680b46e22c053b97c01ba23d9fbb1f145ca5753cc452b`.
The complete public Cantera-header projection and generated dependency-header projection are
byte-identical to their filtered retry8 sources. Their paired manifest SHA-256 values are
`f9b6242471396803f4dda58a5270db12fc0d0e3d6dbd9e389793e7582284bad3`
and `5454b93287db726805dc1e9a6f690de08cdaf693944323fedc1f941b42be7925`.

## Dependency and license closure

The main agent independently checked the final closure tables, not only the worker summary:

| Closure class | Count | Bad hash/status rows |
|---|---:|---:|
| C++ translation units | 207 | 0 nonzero preprocess return codes |
| Preprocessor dependency pairs | 150,402 | 0 |
| Boost headers consumed | 886 | 0 |
| Eigen headers consumed | 251 | 0 |
| Restricted Eigen files identified | 3 | 0 consumed |
| fmt/yaml-cpp/SUNDIALS compiled source files | 80 | 0 |

The aggregate closure manifest SHA-256 is
`afb65aa79948a2eed9c624d9a4cef8cff1fae0d6274a610d0da674e947e1a7ea`.
Every C++ command carries `EIGEN_MPL2_ONLY`. The generated header payload retains Eigen's
three-line `NonMPL2.h` guard, but the compiled dependency closure consumes no restricted Eigen
implementation. The full upstream Eigen license texts remain in the bundle as conservative
notices; their presence is not represented as GPL/LGPL implementation consumption.

The artifact includes the pure Cantera BSD-3-Clause license, SHA-256
`e92980b9712ce20e73898a97b0116889e84e07f548d6be8591e87dcad79c41bb`,
the upstream combined dependency-license inventory, SHA-256
`25bacb1e19380fc7beaff24973da5386c44fdb414188970b6d8fd0681135dca2`,
and individual fmt, yaml-cpp, SUNDIALS, Eigen and Boost license/package-copyright evidence.
No endorsement right is inferred.

All of these payload/runtime counts are zero:

```text
Boost headers/libraries and NEEDED
Doxygen
Clang/LLVM
Xapian
libxml2
typing_extensions
Python files and libpython
static libraries
pkg-config files
```

## Header/link smoke and retained failures

The minimal GCC 11/C++17/ABI=1 header/link consumer compiled and ran against the projected bundle:

| Evidence | Result |
|---|---|
| Source SHA-256 | `7f62d9fe76b3ce9a87a093072cc41b928981463b3d13aea187b4643145193b99` |
| Executable SHA-256 | `2317fa9f78470fe61696ed8ca942a2252e54671a1f75b6440bf47174eb157274` |
| Compile / run exit | 0 / 0 |
| Fixed output | `3.2.0` |

Three rejected projections and one rejected audit assertion remain preserved:

1. v1 used the wrong tar member prefix and retained Python/test headers;
2. v2 omitted the generated SUNDIALS configuration header;
3. v3 mislabeled the upstream combined license inventory as the pure Cantera license;
4. audit v1 required a complete compile-time data path that linker string merging does not preserve.

The rejected file-manifest SHA-256 values are
`849159de889189370c900c8e983e6a4a0730c253e43cf69bcacda39b763e5f63`,
`10264943faaabd5ebe4a687030d6d9bda92285af0388487b825cb4335df3ff0a`, and
`82b2b1f8cee141ada054908225fb69df53a9af8fa89707a414712e7dd30a2311`.
The failed audit evidence-manifest SHA-256 is
`b6136408e536c75ea6148ce680b62f52af948f6b576b069cc582db15765e45bd`.

## Manifest verification

The final artifact manifest was revalidated after `PREFLIGHT_PASS` was selected:

| Evidence | SHA-256 / result |
|---|---|
| Validator | `dcb40515d6affd292909761df48254d22576ffab5009ee7a299efcc2879a431d` |
| Final GREEN log | `48a601589d1d7f25675c80839acbda9fd4a6bec24fabd70f29085c2d2a701cf2`; exit 0 |
| Mutation runner v3 | `bc65534e0bdea34be3452ad74c849e71ab98ae1c75edb5e3de0e5e729821de56` |
| Mutation summary | `ad63a54673cbbd0dffb1ed787e819e021d175cef4837d0ff928b5f679c89874f` |
| Mutation result | 31 expected rejects; 0 unexpected passes |

The bounded manifest worker did not provide verifiable runtime `turn_context` before interruption,
so its model identity is not counted as independent evidence. The main agent independently reviewed
and amended every manifest field, re-hashed every listed external evidence file, and owns this
P0-2 acceptance decision.

## Deferred P0-3 proof

The stripped shared library contains one merged string fragment rooted at the retry8 install
prefix. This is not an ELF RPATH/RUNPATH or a runtime dependency. It is nevertheless not waived:

```text
relocation_status=pending_p0_3
required=p0_3_moved_prefix_runtime_nonconsumption_and_data_directory_override
```

P0-3 must prove that the standalone consumer uses only the moved bundle and the separately supplied
mechanism/data path, with no runtime reference to the original build/source/install roots. Until
that passes, this receipt does not call the complete package relocatable.
