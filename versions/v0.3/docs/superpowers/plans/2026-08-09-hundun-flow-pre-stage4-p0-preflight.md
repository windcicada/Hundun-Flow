# HUNDUN-FLOW Pre-Stage-4 P0 Preflight Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development`
> (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking. A worker receives exactly one task and file allowlist,
> does not commit and does not add DCO. The main agent owns all architecture, ABI, mathematical,
> provenance, complete-diff and acceptance decisions.

**Goal:** 在不修改 HUNDUN-FLOW 产品代码或活跃 Stage 3 工作树的前提下，形成可复现的
Cantera Linux CPU artifact 候选、独立 C++ 边界证据、双代理燃料来源候选、公开数学测试
向量和 Stage 3→4 intake 清单，从而缩短 Stage 4 正式关键路径。

**Architecture:** governance 独立 worktree 只保存计划、manifest、协议和 receipt；所有
第三方源码、rootfs、构建树、二进制和 standalone spike 放在 Git 外的内容寻址目录。
Ubuntu 22.04/GCC 11 builder 通过宿主已有的无特权 `bubblewrap` 和 Canonical 签名 rootfs
建立，不安装 Docker、不改宿主系统或网络。P0 证据始终标记为 preflight candidate，
Stage 3 接受后由正式 `4F-0` 和 `4P-1..4P-4` 逐哈希复核。

**Tech Stack:** Git、Markdown/JSON、CMake 3.21+ script mode、SHA-256、GPGV、curl、
Canonical Ubuntu 22.04 minimal rootfs、bubblewrap、GCC 11/libstdc++、C++17、Open MPI、
Cantera C++ 3.2.0、SCons（仅 maintainer builder）、readelf/ldd、Apache-2.0/DCO。

## Global Constraints

- Stage 4 产品实现仍要求 Stage 3 正式接受和用户明确启动；P0 不是 Stage 4 start。
- 不修改、暂存、清理或提交 product repo 和任一 Stage 3 worktree。
- 不访问 COAST、COAST-2、BOFFIN、私有研究源码、算例、数据或研究进程。
- 不关闭或重配宿主网络；network-independent 通过本地输入与 trace 验证。
- 正常 HUNDUN configure/build/install/test/runtime 不得要求 Python；Python/SCons 只允许在
  外部 maintainer builder rootfs 内。
- 目标 profile 固定为 Linux x86_64、Ubuntu 22.04/glibc 2.35+、GCC 11/libstdc++、
  C++17、`_GLIBCXX_USE_CXX11_ABI=1`、generic x86-64 ISA。
- Clang/libc++ 不与 GCC/libstdc++ Cantera artifact 混链；禁止 `_GLIBCXX_DEBUG`。
- Cantera、传递依赖、机理和液体数据分别登记来源、revision、SHA-256 和许可证。
- 不复制、翻译或机械改写第三方实现到 HUNDUN `src/`；P0 没有产品源码文件。
- 任意时刻最多一个 third-party build；默认 `-j16`，无资源冲突时最多 `-j32`。
- standalone spike 最多 2 MPI ranks；不运行网格、反应流、ESF/TCR、喷雾或长算例。
- 外部目录只追加；未知文件不删除。worker 不提交、不添加 DCO、不 push。

---

## 1. File and Ownership Map

### Tracked governance files

| Path | Owner | Responsibility |
|---|---|---|
| `AGENTS.md` | main agent | 登记 P0 是唯一 pre-Stage-4 非产品例外 |
| `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md` | main agent | 在状态机前增加 P0 候选输入边 |
| `.superpowers/stage4-p0/baseline-receipt.md` | main agent | exact planning baseline、宿主与隔离证据 |
| `.superpowers/stage4-p0/input-manifest.json` | provenance worker draft；main review | Cantera、submodule、rootfs、工具链和许可证 identity |
| `.superpowers/stage4-p0/validate-manifest.cmake` | bounded worker | JSON schema、必填值和磁盘 SHA 验证 |
| `.superpowers/stage4-p0/provenance-receipt.md` | main agent | 来源、许可、patch 和 binary-policy 结论 |
| `.superpowers/stage4-p0/artifact-manifest.json` | build worker draft；main review | 构建命令、环境、headers/data/shared libraries 和 hashes |
| `.superpowers/stage4-p0/artifact-receipt.md` | main agent | artifact 构建和可复现性结论 |
| `.superpowers/stage4-p0/linkage-receipt.md` | main agent | C++、thread、MPI、ABI、RPATH 和 relocation 证据 |
| `docs/references/2026-08-09-hundun-flow-stage4-p0-fuel-data-candidates.md` | provenance worker draft；main review | 双代理燃料来源、许可和可再分发状态 |
| `docs/numerics/2026-08-09-hundun-flow-stage4-6-p0-oracle-vectors.md` | main agent | 公开方程、单位、输入输出向量和 mutation |
| `.superpowers/stage4-p0/intake-dry-run.md` | main agent | accepted Stage 3→4 只读命令和 collision checklist |
| `.superpowers/stage4-p0/final-receipt.md` | main agent | P0 exact seal、复用/未证明边界和后台进程状态 |

### External generated tree

```text
/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/
  inputs/ubuntu/
  inputs/cantera/
  inputs/dependencies/
  source/rootfs-jammy/
  source/cantera-3.2.0/
  build/cantera-3.2.0-gcc11-release/
  install/cantera-3.2.0-gcc11-release/
  spikes/cantera-cxx-v1/
  logs/
  manifests/
```

This tree is never staged. Every command uses the explicit root above, never `$HOME`, `~`, a globbed
delete target or a shared Stage 3 build directory.

---

### Task P0-0: Activate the Isolated Preflight Lane

**Files:**
- Modify: `AGENTS.md`
- Modify: `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md`
- Create: `.superpowers/stage4-p0/baseline-receipt.md`

**Interfaces:**
- Consumes: approved design commit `aea268dcd8e04afb466f40b59f8668ae34af4e7b` and this plan.
- Produces: branch/worktree identity, P0 authority text, external-root allowlist and exact Stage 3
  non-interference baseline consumed by every later task.

- [ ] **Step 1: Create the execution worktree from the signed plan commit.** Main agent verifies
  `.worktrees` is the established location and the branch does not exist, then runs:

  ```bash
  p0_plan_head="$(git rev-parse HEAD)"
  git worktree add \
    /home/wyf/code_dev/.worktrees/hundun-flow-stage4-p0-preflight \
    -b coast/stage4-p0-preflight "$p0_plan_head"
  ```

  Immediately verify the new worktree HEAD equals `p0_plan_head`, then record `git rev-parse HEAD`,
  `git rev-parse HEAD^{tree}`, branch, absolute Git dir, common dir and
  `git status --porcelain=v1`.

- [ ] **Step 2: Write the policy RED before changing authority text.** Run these assertions and
  require both to fail because the approved P0 spec is not yet in `AGENTS.md` or integration state:

  ```bash
  rg -q 'pre-Stage-4 P0' AGENTS.md
  rg -q 'P0 preflight candidate' \
    docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md
  ```

- [ ] **Step 3: Add the minimal authority amendment.** In `AGENTS.md`, add the P0 design and plan
  to Required Reading, state that only this non-product lane is authorized before Stage 3 acceptance,
  and preserve the prohibition on `4F-0`/product implementation. In the integration plan, add:

  ```text
  approved P0 planning commit
      -> external provenance/artifact/oracle candidates
      -> hash revalidation after accepted Stage 3
      -> 4F-0 and formal 4P tasks
  ```

  It must explicitly say P0 is not an ancestor requirement for accepted Stage 3 and does not alter
  the Stage 4 version or capability state.

- [ ] **Step 4: Create the baseline receipt with fresh evidence.** Record planning HEAD/tree,
  product status, governance-main status, the three Stage 3 worktree HEAD/branch/porcelain counts,
  `uname -m`, glibc, GCC, CMake, MPI, bubblewrap path and this successful probe:

  ```bash
  bwrap --ro-bind / / --dev-bind /dev /dev --proc /proc \
    --unshare-user --uid 0 --gid 0 /bin/true
  ```

  Record only HUNDUN worktree state; do not enumerate or signal unrelated processes.

- [ ] **Step 5: Verify the governance-only boundary.** Run:

  ```bash
  git diff --check
  git diff --name-only | rg -v \
    '^(AGENTS\.md|docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration\.md|\.superpowers/stage4-p0/baseline-receipt\.md)$'
  ```

  Expected: the second command prints nothing. Re-run the two Step 2 assertions; both PASS. Confirm
  product `main` is clean and do not require active Stage 3 status to remain unchanged while its
  other agent is working.

- [ ] **Step 6: Main-agent review and commit.** Review the complete diff, then create a DCO commit:

  ```bash
  git add AGENTS.md \
    docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md \
    .superpowers/stage4-p0/baseline-receipt.md
  git commit -s -m 'docs: activate Stage 4 P0 preflight'
  ```

---

### Task P0-1: Lock Official Inputs, Submodules and Licenses

**Files:**
- Create: `.superpowers/stage4-p0/input-manifest.json`
- Create: `.superpowers/stage4-p0/validate-manifest.cmake`
- Create: `.superpowers/stage4-p0/provenance-receipt.md`

**Interfaces:**
- Consumes: P0-0 external root, official public HTTPS/Git sources and signature-verified
  distribution package metadata.
- Produces: schema `hundun.stage4_p0.inputs.v1`, verified local archive paths and dependency/license
  records consumed by P0-2, P0-3 and P0-4. The lock distinguishes four Cantera source-tree
  archives from builder-only header packages discovered by the frozen P0-2 configuration.

- [ ] **Step 1: Acquire signed Ubuntu metadata and exact Cantera identities.** Create only the
  explicit external directories with `mkdir -p`. Use the Canonical base URL
  `https://cloud-images.ubuntu.com/minimal/releases/jammy/release/` to download `SHA256SUMS`, its
  signature and `ubuntu-22.04-minimal-cloudimg-amd64-root.tar.xz`; verify:

  ```bash
  gpgv --keyring /usr/share/keyrings/ubuntu-cloudimage-keyring.gpg \
    /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/inputs/ubuntu/SHA256SUMS.gpg \
    /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/inputs/ubuntu/SHA256SUMS
  (cd /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/inputs/ubuntu && \
    rg 'ubuntu-22.04-minimal-cloudimg-amd64-root\.tar\.xz$' SHA256SUMS | sha256sum -c -)
  ```

  Verify Cantera annotated tag and peeled commit with `git ls-remote`; require peeled commit
  `4a8358eb80cfeb50474386b5f9ec0b3a83519889`. Download the official release asset from
  `https://github.com/Cantera/cantera/releases/download/v3.2.0/cantera-3.2.0.tar.gz` and require
  SHA-256 `a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b`。Do not substitute
  `https://github.com/Cantera/cantera/archive/refs/tags/v3.2.0.tar.gz`; that generated archive has
  a different identity (`f01e25e33f9d5e37db7ababe5af36b60caabff52dba04bb221d53e44735f60ec`).

- [ ] **Step 2: Freeze required submodule commits and archives.** Read `.gitmodules` and the Git tree
  at `v3.2.0`; record at least these exact commits before downloading their official archives:

  ```text
  ext/fmt       https://github.com/fmtlib/fmt.git
                a33701196adfad74917046096bf5a2aa0ab0bb50
  ext/yaml-cpp  https://github.com/jbeder/yaml-cpp.git
                0579ae3d976091d7d664aa9d2527e0d0cff25763
  ext/sundials  https://github.com/LLNL/sundials
                887af4374af2271db9310d31eaa9b5aeff49e829
  ext/eigen     https://gitlab.com/libeigen/eigen.git/
                3147391d946bb4b6c68edd901f2add6ac1f31f8c
  ```

  `googletest`, HighFive and example-data remain excluded unless P0-2 proves the selected C++-only
  configuration consumes them. If consumed, add their exact Cantera-tree commit, archive hash and
  license before rebuilding; never silently fetch a missing submodule.

  The first current-option P0-2 configuration proved that Cantera 3.2.0 also requires Boost
  headers. Freeze exact Jammy package `libboost1.74-dev=1.74.0-14ubuntu3` (source package
  `boost1.74`, `amd64`) as `builder_header_only`. The generic `libboost-dev` meta package is
  resolution evidence only and must not be installed or recorded as a consumed input. This package
  does not become a Cantera source-tree archive or an allowed artifact/runtime dependency.
  Canonical's Jammy apt source resolves the package through an HTTP retrieval URL, but apt verifies
  its size and SHA-256 against signed repository metadata; record both that resolved URL and the
  exact package digest instead of misrepresenting it as an unauthenticated source download.

  The isolated-Boost retry subsequently proved that Cantera 3.2's generated CLib source graph
  unconditionally invokes Doxygen even when `doxygen_docs=no`. Before another build, freeze this
  exact Jammy builder-tool closure as `builder_tool_dependencies`:

  ```text
  doxygen          1.9.1-2ubuntu2
  libclang-cpp14   1:14.0.0-1ubuntu1.1
  libclang1-14     1:14.0.0-1ubuntu1.1
  libllvm14        1:14.0.0-1ubuntu1.1
  libxapian30      1.4.18-4
  libxml2          2.9.13+dfsg-1ubuntu0.12
  ```

  Preserve the exact signed-apt URL, size, SHA-256 and complete Ubuntu package copyright inventory
  for every binary package. These are maintainer-builder tools, not Cantera source inputs.
  Doxygen/Clang/LLVM/Xapian are forbidden from the install payload and runtime closure. `libxml2`
  remains `audit_pending` until the finished artifact's `readelf`/`ldd` inventory proves whether it
  is only a Doxygen transitive dependency; do not infer that result before P0-2. Do not infer the
  generated or installed artifact's license from any builder package: retain the full package
  inventories and audit the actual Cantera/sourcegen outputs separately.

  The first Doxygen-enabled build may also discover Python modules used only by Cantera's external
  maintainer source generator. Freeze each such module before another build. For Python 3.10,
  Cantera 3.2.0 sourcegen imports `typing_extensions.Self`, while its sourcegen
  `pyproject.toml` declares only Jinja2 and Jammy's `python3-typing-extensions=3.10.0.2-1` does not
  provide the required API. Use the official PyPI `typing_extensions=4.15.0` pure-Python wheel as
  an exact `builder_pythonpath_only` input. Freeze the PyPI JSON metadata, wheel, source archive,
  byte sizes, SHA-256 values and matching PSF-2.0 license files. Do not install the wheel into the
  rootfs, run `pip`, or add it to the artifact: expose only the exact wheel through `PYTHONPATH`
  during the maintainer build and require zero installed/runtime Python payload afterward.

- [ ] **Step 3: Audit license and mechanism separation.** For Cantera and each consumed dependency,
  save the upstream license in the external `inputs/dependencies/` tree and record its SHA-256,
  copyright holder, SPDX identifier and redistribution obligations. For an OS builder package,
  retain both its exact package copyright inventory and any primary upstream license evidence.
  Do not collapse a multi-license package inventory into one SPDX conclusion: Boost remains
  `candidate_consumed_header_scope_audit_required` until P0-2 identifies the headers actually
  included or compiled. Record `h2o2.yaml` or any other spike mechanism as a separate asset with
  its own path/hash/source/licensing conclusion. A missing or ambiguous license sets the component
  status to `candidate_user_supplied` or rejects it; it is never inferred from Cantera's BSD
  license.

- [ ] **Step 4: Write the manifest and failing validator mutation.** Use `apply_patch` to create
  `input-manifest.json` with actual observed rootfs and dependency archive hashes. The immutable
  fields must include:

  ```json
  {
    "schema": "hundun.stage4_p0.inputs.v1",
    "result": "PREFLIGHT_PARTIAL",
    "stage4_product_accepted": false,
    "product_changes": "none",
    "target": {
      "os": "Ubuntu 22.04",
      "glibc_floor": "2.35",
      "arch": "x86_64",
      "compiler": "GCC 11",
      "stdlib": "libstdc++",
      "cxx_standard": 17,
      "glibcxx_cxx11_abi": 1,
      "isa": "x86-64"
    },
    "cantera": {
      "tag": "v3.2.0",
      "commit": "4a8358eb80cfeb50474386b5f9ec0b3a83519889",
      "source_sha256": "a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b"
    }
  }
  ```

  Create `validate-manifest.cmake` using `file(READ)`, `string(JSON ... GET)`, `file(REAL_PATH)`
  and `file(SHA256)`. Canonicalize the external root and every selected file before hashing；reject
  an absolute path, `..`, `file://` or a root-internal symlink resolving outside the root. A legal
  symlink resolving to another file inside the same root remains valid.
  First invoke it with
  `-DEXPECTED_CANTERA_SHA=0000000000000000000000000000000000000000000000000000000000000000`；
  expected failure is
  `Cantera archive SHA mismatch`.

- [ ] **Step 5: Validate real inputs and scan policy.** Run:

  ```bash
  cmake \
    -DMANIFEST=.superpowers/stage4-p0/input-manifest.json \
    -DEXTERNAL_ROOT=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0 \
    -DEXPECTED_CANTERA_SHA=a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b \
    -P .superpowers/stage4-p0/validate-manifest.cmake
  git diff --check
  ```

  Also scan the manifest for `file://`, private COAST/BOFFIN paths, unpinned `main/master/latest`,
  missing license status and any download command under product CMake. Expected: none.
  Add two path mutations: a root-internal symlink to the exact installed Ubuntu keyring must fail
  with `resolves outside EXTERNAL_ROOT`; a root-internal symlink to the signed checksum file inside
  the same P0 root must pass. Do not replace a real archive or write the symlink target.
  For every newly frozen builder input, add mutation-sensitive rejection of a wrong hash, version,
  role, package identity, bundle/runtime status, missing record and external-root symlink escape.
  Also reject any mutation that prematurely marks a pending consumed-file license audit complete.

- [ ] **Step 6: Main-agent provenance review and commit.** Verify every URL is official, every
  consumed binary has a source/license path, no endorsement claim exists and no mechanism inherits
  Cantera licensing by assumption. Update `result` to `PREFLIGHT_PASS` only for the input lock, write
  the receipt, then commit the three tracked files with subject
  `docs: lock Stage 4 P0 third-party inputs` and a valid DCO trailer.

---

### Task P0-2: Build the Ubuntu 22.04/GCC 11 Release Artifact Candidate

**Files:**
- Create: `.superpowers/stage4-p0/artifact-manifest.json`
- Create: `.superpowers/stage4-p0/artifact-receipt.md`
- Modify: `.superpowers/stage4-p0/input-manifest.json` only if P0-2 discovers a consumed dependency
  not already frozen by P0-1; re-run all P0-1 checks and commit the amended lock before proceeding.

**Interfaces:**
- Consumes: P0-1 verified rootfs and source archives.
- Produces: immutable external install root and schema `hundun.stage4_p0.artifact.v1` for P0-3.

- [ ] **Step 1: Materialize but do not trust the builder root.** Extract the verified rootfs with
  `tar --no-same-owner` into `source/rootfs-jammy/`. Launch it with `bubblewrap`, preserving host
  networking while isolating filesystem and PID state. The command must bind only the explicit P0
  external root and required kernel filesystems, for example:

  ```bash
  bwrap --unshare-user --uid 0 --gid 0 --unshare-pid --share-net \
    --bind /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/source/rootfs-jammy / \
    --dev-bind /dev /dev --proc /proc --ro-bind /sys /sys \
    --bind /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0 /p0 \
    --tmpfs /tmp --chdir /p0 /usr/bin/env -i \
    PATH=/usr/sbin:/usr/bin:/sbin:/bin LC_ALL=C.UTF-8 /bin/bash --noprofile --norc
  ```

  Verify inside the rootfs: x86_64, glibc 2.35, GCC/G++ major 11 after installing only recorded
  Jammy packages. Cache every acquired `.deb`, record version and SHA-256, and do not modify the host.
  Install the exact frozen `libboost1.74-dev` archive only after P0-1 validates it; do not install
  the `libboost-dev` meta package or any suggested Boost component package.
  Likewise install the six exact Doxygen closure archives only after the amended P0-1 validator
  verifies every package and copyright hash. Use `--no-install-recommends`; do not install suggested
  Doxygen documentation, GUI, LaTeX, Graphviz or Xapian tools. Record rootfs package state before and
  after installation and never install these packages on the host.

- [ ] **Step 2: Run the expected pre-build RED.** Compile a one-line C++ program including
  `<cantera/base/Solution.h>` against the empty install root. Expected: compilation fails because the
  Cantera headers do not exist. Record command, exit `!=0` and log SHA in `logs/`.

- [ ] **Step 3: Assemble exact source and build C++-only shared libraries.** Expand the verified
  top-level source and only the four accepted source dependency archives into their pinned `ext/`
  paths. Boost remains an independently frozen builder-root header package and is never copied into
  the Cantera source overlay. Preserve that full extract unchanged, then create a separately named
  canonical build overlay.
  The overlay keeps only upstream `data/README.md` under `data/`; all top-level YAML, including
  mechanisms and optional support databases, plus `data/example_data` remain in the immutable source
  extract but must not enter the installed
  artifact. This is a packaging filter, not a license inference or an upstream source patch. P0-3
  supplies a separately licensed HUNDUN-authored synthetic mechanism as an external input.

  Cantera 3.2.0 unconditionally creates a shared library in `src/SConscript`; it no longer exposes
  the historical `shared_libs` SCons option. `renamed_shared_libraries=yes` keeps the shared and
  static basenames distinct, while `versioned_shared_library=yes` controls the versioned SONAME and
  symlinks. Do not pass the absent `shared_libs` option.

  Do not set `boost_inc_dir=/usr/include`: SCons converts it to an early `-isystem /usr/include`,
  which reorders GCC's standard headers and breaks `cmath`/`#include_next <math.h>`. Instead create
  an otherwise empty build-only include root at
  `/p0/build/cantera-3.2.0-gcc11-release/boost-system-include-root` and add a bwrap read-only bind:

  ```text
  <verified-rootfs>/usr/include/boost
    -> /p0/build/cantera-3.2.0-gcc11-release/boost-system-include-root/boost
  ```

  Require a pre-build A/B diagnostic in the same bwrap environment: `<cmath>` and
  `<boost/version.hpp>` must both compile with the isolated root, the resolved Boost header must be
  byte-identical to the frozen package member, and `/usr/include` must remain the final standard
  include directory. Run SCons inside the builder with explicit settings:

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
  CXX=g++-11
  CC=gcc-11
  boost_inc_dir=/p0/build/cantera-3.2.0-gcc11-release/boost-system-include-root
  cxx_flags=-std=c++17 -D_GLIBCXX_USE_CXX11_ABI=1 -DEIGEN_MPL2_ONLY -march=x86-64 -mtune=generic
  optimize=yes
  prefix=/p0/install/cantera-3.2.0-gcc11-release
  ```

  Before building, run `scons help --options` and fail if any current named option is absent or has
  changed semantics. Also require `shared_libs` to be absent and independently verify the frozen
  3.2.0 source graph still constructs a `SharedLibrary`; if either expectation changes, stop for a
  new plan revision. Reject an absent Boost read-only bind, an explicit `/usr/include` Boost path or
  any Boost header whose content differs from the frozen `.deb`. Use `nice -n 10`, `-j16` by default
  and capture the complete command/environment/log.

  Retain the default generated CLib path. Do not use `clib_legacy=yes`, custom partial targets or
  copied generated sources to avoid Doxygen: the legacy switch changes the public API, while custom
  build/install paths violate the official-install/no-patch boundary. The accepted hypothesis is an
  exact builder-only Doxygen closure and exact builder-only sourcegen `PYTHONPATH`, followed by the
  upstream `build install` target. Do not use `pip` or modify the Cantera source to repair an
  undeclared sourcegen import.

- [ ] **Step 4: Inventory installed files and dynamic dependencies.** Hash every regular file under
  the install root in bytewise sorted relative-path order. Record symlink targets separately. Run
  `readelf -d`, `readelf --version-info`, `ldd`, `strings` policy scans and `g++-11 -v`. Reject:

  ```text
  NEEDED dependency outside the allowlist
  Python/Conda path or libpython
  build/rootfs absolute RPATH
  GLIBC requirement newer than 2.35
  GLIBCXX requirement unavailable in GCC 11 libstdc++
  -march=native or non-generic ISA
  missing EIGEN_MPL2_ONLY or a NonMPL2.h compile failure
  installed `include/boost/**`
  installed `lib*/cmake/Boost*/**` or `lib*/cmake/boost_headers*/**`
  installed `libboost*.so*` or `libboost*.a`
  any `libboost*.so` dynamic dependency
  static-only Cantera
  missing Cantera data directory
  any installed data file other than the exact upstream data/README.md
  ```

  Require `bundled_data_file_count=1`, an exact installed README hash and
  `bundled_mechanism_count=0`. The installed Cantera data directory contains only the upstream
  disclaimer README at P0; high-pressure/real-gas support-database capability is not claimed.
  Formal Stage 4 may add an asset only after its independent provenance/license gate. An empty
  mechanism set does not weaken the C++ linkage spike because P0-3 consumes its synthetic mechanism
  by exact absolute path.

  Independently inventory the complete transitive Boost header closure selected by compiler
  depfiles or the SCons dependency graph. Hash every consumed header and map it to the frozen Ubuntu
  package copyright inventory. If closure completeness cannot be proved, conservatively apply the
  complete package license inventory, keep the audit pending and leave the artifact result
  `PREFLIGHT_PARTIAL` or `PREFLIGHT_REJECT`. A missing or incompatible file-level license conclusion
  is a P0-2 blocker; it is not permission to relabel the package as uniformly BSL-1.0. Require
  `boost_payload_file_count=0` in the install tree, `boost_runtime_needed_count=0`, and no basename
  beginning with `libboost` in the `readelf`/`ldd` runtime closure.

  Independently require `doxygen_payload_file_count=0`, `clang_payload_file_count=0`,
  `llvm_payload_file_count=0` and `xapian_payload_file_count=0`. Their runtime-needed counts must also
  be zero. Record `libxml2_payload_file_count` and `libxml2_runtime_needed_count` from the actual
  artifact; if either is nonzero, stop for provenance/packaging analysis instead of silently
  reclassifying the frozen builder package as a product dependency.
  Require `typing_extensions_payload_file_count=0`, `python_payload_file_count=0`,
  `libpython_runtime_needed_count=0` and no Python/Conda path in installed files or dynamic metadata.

- [ ] **Step 5: Write and validate the artifact manifest.** Record source/input manifest SHA,
  builder rootfs SHA, package versions, exact SCons configuration, wall time/RSS, installed file
  hashes, SONAMEs, NEEDED libraries, data path, `bundled_data_file_count`,
  `bundled_mechanism_count`, packaging-filter hashes, license paths and logs. Record the exact Boost
  builder package identity, actual consumed-header inventory/license conclusion, and the negative
  Boost install/runtime scans. Keep the P0-1 input manifest's
  `consumed_header_audit_status=pending_p0_2` immutable; place the completed audit result only in the
  P0-2 artifact manifest and receipt. Set
  `stage4_product_accepted=false`, `product_changes=none`, and initially
  `result=PREFLIGHT_PARTIAL`. Independently re-hash all paths from the manifest and require equality.

- [ ] **Step 6: Main-agent build review and commit.** Inspect logs and source/options without copying
  upstream implementation into Git. If build and inventory pass, set artifact result to
  `PREFLIGHT_PASS`, write the receipt and commit only the tracked manifest/receipt with subject
  `build: record Stage 4 P0 Cantera artifact candidate`. If not, record `PREFLIGHT_REJECT` or
  `PREFLIGHT_PARTIAL`; do not switch to package strategy B/C.

---

### Task P0-3: Prove C++ Linkage, Workspace Isolation and Relocation

**Files:**
- Create: `.superpowers/stage4-p0/linkage-receipt.md`
- Modify: `.superpowers/stage4-p0/artifact-manifest.json` only to append immutable spike log/binary
  hashes; never change artifact input/build identity.

**Interfaces:**
- Consumes: P0-2 install root and one separately licensed/hash-verified tiny mechanism.
- Produces: standalone numeric summary schema `hundun.stage4_p0.cantera_spike.v1` and P0-3 verdict.

- [ ] **Step 1: Write the standalone spike in the external tree.** Use `apply_patch` to create
  `spikes/cantera-cxx-v1/main.cpp`. It must accept `<mechanism> <phase>`, construct one immutable
  configuration per rank and one complete mutable `Solution`/thermo/kinetics/transport workspace per
  OpenMP thread, then print in fixed key order:

  ```text
  schema
  cantera_version
  mechanism_sha256
  rank_count
  thread_count
  temperature_K
  pressure_Pa
  density_kg_m3
  cp_mass_J_kg_K
  viscosity_Pa_s
  net_production_rate_hash
  reactor_final_state_hash
  workspace_alias_count
  ```

  Use no HUNDUN headers or source. `workspace_alias_count` must be zero; pointer values themselves
  never cross ranks or enter the numeric equality hash.

- [ ] **Step 2: Compile and run the single-rank Release smoke.** Use `g++-11 -std=c++17
  -D_GLIBCXX_USE_CXX11_ABI=1`, the artifact include/lib paths and an explicit temporary RUNPATH.
  Run thermo, mixture-averaged transport and a fixed short 0D interval. Require finite positive
  temperature/density/cp/viscosity, finite rates, exact field presence and exit 0.

- [ ] **Step 3: Run thread and 1/2-rank isolation.** Compile with the builder's Open MPI wrapper and
  OpenMP. Run fixed one-thread, two-thread, `mpiexec -n 1` and `mpiexec -n 2` cases. Require:

  ```text
  all ranks exit 0
  workspace_alias_count=0
  fixed lane-order numeric hashes equal across repeated runs
  rank-local numeric summaries equal for identical inputs
  no Python process or module is launched
  no collective mismatch or timeout
  ```

  This is an interface spike, not a chemistry performance or decomposition test.

- [ ] **Step 4: Test ABI rejection and ordinary Debug consumption.** A mutation compiled with
  `_GLIBCXX_USE_CXX11_ABI=0` must be rejected by the preflight manifest policy before execution.
  Compile an ordinary `-O0 -g` GCC 11 consumer with ABI=1 and matching exceptions/RTTI; it must link
  the same Release artifact. Scan compile commands for `_GLIBCXX_DEBUG`, `clang`, `libc++` and
  unexpected host GCC; expected none.

- [ ] **Step 5: Move the complete install prefix and re-run.** Copy the bundle to a new explicit
  directory under `install/relocated/`, never edit individual `.so` files, then run using only
  relative `$ORIGIN`-style search paths or a wrapper-local loader path. `readelf -d`, `ldd`, runtime
  data discovery and numeric hashes must not reference the original build/source root. Record both
  tree hashes and package size.

- [ ] **Step 6: Main-agent review and commit.** Recompute binary/log hashes, verify no product or
  Stage 3 file changed, state exactly what the standalone spike does not prove, and commit the
  receipt plus allowed manifest append with subject `test: prove Stage 4 P0 Cantera boundary`.

---

### Task P0-4: Audit Dual-Surrogate Mechanism and Liquid-Property Candidates

**Files:**
- Create: `docs/references/2026-08-09-hundun-flow-stage4-p0-fuel-data-candidates.md`

**Interfaces:**
- Consumes: public official/project/paper sources only and P0-1 licensing vocabulary.
- Produces: separate n-dodecane and iso-octane candidate rows with statuses
  `bundled_candidate`, `user_supplied_only`, `rejected` or `unresolved`, consumed by formal
  Stage 4/6 provenance.

- [ ] **Step 1: Freeze evaluation fields before searching.** Start the document with one table row
  per asset and these required columns:

  ```text
  surrogate_family | asset_role | public_name | official_url | release_or_revision |
  sha256 | citation_or_doi | copyright_holder | license_or_terms_url |
  redistribution_status | Cantera_parse_status | intended_low_cost_check | excluded_claims
  ```

  The two families are `n_dodecane_kerosene_surrogate` and `iso_octane_gasoline_surrogate`.

- [ ] **Step 2: Research gas-mechanism candidates from primary sources.** Prefer upstream mechanism
  authors, institutional repositories or an explicitly licensed Cantera-distributed example. Do not
  use a search-result mirror as authority. Verify revision and SHA locally. If terms permit use but
  not redistribution, mark `user_supplied_only`; do not copy the mechanism into Git.

- [ ] **Step 3: Research pure-liquid property sources independently.** Record critical temperature,
  boiling point, molecular weight, density, viscosity, conductivity, heat capacity, latent heat and
  surface-tension correlation sources with units and validity ranges. Separate public numerical facts
  from copyrighted correlation text/code. Reject a candidate missing enough data for the planned
  Stage 6 service rather than filling values from COAST.

- [ ] **Step 4: Run low-cost external parse/property checks.** For each legally accessible mechanism,
  use the P0-3 C++ spike to load the named phase and report species/reaction counts, thermo state and
  mechanism SHA. For liquid data, independently evaluate two in-range points and one out-of-range
  failure in a standalone scratch calculation. No HUNDUN source or long flame/spray case is allowed.
  This step depends on the accepted P0-3 artifact/spike boundary. The public provenance/legal
  subcluster may be committed earlier with `PREFLIGHT_PARTIAL`, but that commit does not complete
  P0-4 and must list every parse/property check as deferred.

- [ ] **Step 5: Write capability limitations.** State that two surrogate families prove only
  interface generality and absence of fuel-name hardcoding. They do not validate real aviation
  kerosene/gasoline, ignition delay, flame speed, spray combustion or COAST similarity. State that
  later COAST `EXEC/Fuels` comparison still requires the user's exact-path confirmation.

- [ ] **Step 6: Main-agent source/legal review and commits.** Verify every substantive claim against
  a primary source, every `bundled_candidate` has explicit redistribution permission, and every
  uncertain asset is user-supplied, unresolved or rejected. The early public-provenance commit uses
  subject `docs: audit Stage 4 P0 fuel data candidates` and DCO while retaining
  `PREFLIGHT_PARTIAL`. After P0-3, append the hash-bound parse/property evidence in a separate
  validation commit before marking P0-4 complete.

---

### Task P0-5: Freeze Public Mathematical Oracle Vectors

**Files:**
- Create: `docs/numerics/2026-08-09-hundun-flow-stage4-6-p0-oracle-vectors.md`

**Interfaces:**
- Consumes: approved Stage 4--6 equations/reference catalog and public papers only.
- Produces: unit-tagged decimal/hex test vectors and mutation expectations; no executable product
  oracle and no COAST-derived values.

- [ ] **Step 1: Define a common vector schema and tolerances by arithmetic class.** Every vector has:

  ```text
  vector_id
  public_equation_and_citation
  SI_input
  operation_order
  expected_decimal
  expected_binary64_hex_when_exact
  comparison_class: bitwise | ulp_bounded | relative_absolute
  mutation_and_expected_failure
  capability_not_proved
  ```

  Exact identities use bitwise or hexadecimal expectations; transcendental/reference-backend values
  use a stated absolute/relative tolerance justified by the operation, never a product science gate.

- [ ] **Step 2: Add Stage 4 vectors.** Include at least: all-species sum and element matrix identity,
  ideal-gas density, total thermochemical enthalpy including formation terms, integrated two-species
  zero-mass chemistry delta, failed interval canonical zero delta, and Strang reversibility on a
  commuting analytic source/transport pair. Mutations must catch dependent-species storage,
  sensible-only enthalpy, endpoint rate substitution and partial failed-state publication.

- [ ] **Step 3: Add Stage 5 vectors.** Record published Philox4x32-10 vectors verbatim within source
  quotation limits, N=2 and N=4 antithetic Wiener sums, IEM exact exponential relaxation, laminar
  zero-stochastic coefficient and simplex/element/mean preservation. TCR rows contain only the two
  cited papers' independently derived algebra and are labelled `candidate_before_COAST_oracle`.

- [ ] **Step 4: Add Stage 6 vectors.** Include equal-and-opposite parcel/gas momentum, total
  thermochemical-energy source sign, D-squared evaporation over one fixed interval, drag relaxation,
  heat/mass-transfer limiting behavior, normalized deposition weights and rejected-trial zero publish.
  Provide both n-dodecane-like and iso-octane-like property identities without embedding production
  coefficients not yet licensed.

- [ ] **Step 5: Independently recompute each vector.** Use two methods where practical: closed-form
  hand algebra plus a tiny standalone C++17 calculation, or published vector plus C++ calculation.
  Record source/compiler hash and output hash outside Git. Deliberately apply every named mutation and
  verify it changes the expected result or status.

- [ ] **Step 6: Main-agent mathematical review and commit.** Check units, signs, ordering, exact-zero
  semantics, retry clock and capability limitations. Scan for copied source/control flow and private
  paths. Commit with subject `docs: freeze Stage 4-6 P0 oracle vectors` and DCO.

---

### Task P0-6: Prepare the Accepted Stage 3 Intake Dry-Run

**Files:**
- Create: `.superpowers/stage4-p0/intake-dry-run.md`

**Interfaces:**
- Consumes: tracked Stage 3 plans and current repository topology, never dirty file contents.
- Produces: exact read-only command list and expected inventory fields for formal `4F-0`.

- [ ] **Step 1: Write the intake identity block.** Require code/product/governance HEAD, tree, parent,
  version, accepted receipt, DCO, branch, remote, dirty/untracked status and linked `.git` pointer.
  Mark every value `execute_after_stage3_acceptance`; do not populate it from current active heads.

- [ ] **Step 2: Freeze the public/build inventory commands.** Include exact commands using `git`,
  `rg`, `cmake --build --target help`, `cmake --install`, `nm -D`, `readelf`, `ldd` and header list.
  Inventory current accepted naming domains and detect collisions with planned `chem_`, `comb_`,
  `spray_`, `rt_`, `hundun_chemistry`, `hundun_combustion` and `hundun_spray`.

- [ ] **Step 3: Freeze schema/persistence/diagnostics inventory.** Search tracked accepted files for
  schema v1--v3, Checkpoint v1--v3 IDs, diagnostics kinds, field identities, retry/rollback reports,
  final-flux authority, PISO corrector count and capability ledger. The dry-run must distinguish
  absent symbols from renamed symbols and forbid creation of a second authority.

- [ ] **Step 4: Add process and resource checks.** List only processes, systemd user units, commands
  and working directories identifiable as HUNDUN build/test/MPI jobs. The checklist says never stop a
  job during intake and never enumerate or infer unrelated research jobs.

- [ ] **Step 5: Mutation-review the checklist.** Remove one identity class at a time—Checkpoint IDs,
  diagnostics IDs, CMake targets, final flux, rollback or linked-worktree pointer—and verify the
  completeness table reports the omission. Run `git diff --check` and a private-path scan.

- [ ] **Step 6: Main-agent review and commit.** Confirm the document contains no current dirty
  Stage 3 value presented as accepted, no modification command and no inferred Stage 4 parent. Commit
  with subject `docs: prepare Stage 3 to Stage 4 intake` and DCO.

---

### Task P0-7: Seal the Preflight Candidate and Stop at the Stage Boundary

**Files:**
- Create: `.superpowers/stage4-p0/final-receipt.md`
- Modify: `.superpowers/stage4-p0/input-manifest.json` only to record final P0 result without changing
  immutable identities.
- Modify: `.superpowers/stage4-p0/artifact-manifest.json` only to record final P0 result without
  changing immutable identities.

**Interfaces:**
- Consumes: P0-0..P0-6 commits, external artifacts/logs and current read-only repository states.
- Produces: exact P0 seal with `PREFLIGHT_PASS`, `PREFLIGHT_PARTIAL` or `PREFLIGHT_REJECT`.

- [ ] **Step 1: Freeze the candidate.** Read `p0_plan_head` from P0-0's baseline receipt and record
  HEAD, parent, tree, P0 design/plan hashes, `p0_plan_head..HEAD` diff SHA-256, worktree status and
  DCO for every P0 commit. Do not modify tracked
  files consumed by P0-1..P0-6 after this identity is selected.

- [ ] **Step 2: Revalidate all hashes and low-cost evidence.** Re-run manifest validation, archive
  SHA, installed tree hash, binary/log SHA, `readelf`, `ldd`, single/thread/1/2-rank spike and moved
  prefix smoke. Do not rebuild merely because the receipt changes; rebuild only if immutable input,
  builder or artifact validation fails.

- [ ] **Step 3: Complete one consolidated review.** Review the whole P0 diff once for requirements,
  quality, source/license independence, host/product/Stage 3 boundaries, ABI, claims and task scope.
  Do not rescan the same diff under separate requirements/quality/complete-diff labels.

- [ ] **Step 4: Write reuse and non-proof tables.** At minimum state:

  ```text
  reusable: verified source/license hashes, builder profile, artifact candidate, standalone vectors,
            intake command templates
  not proved: HUNDUN CMake integration, ChemistryBackend, reacting transport, Restart v4,
              diagnostics v4, Stage 4 scientific acceptance, COAST equivalence, real-fuel validity
  ```

- [ ] **Step 5: Verify repository and process boundaries.** Require P0 worktree clean; product `main`
  clean; governance-main retains its unrelated Stage 7 untracked file; active Stage 3 worktrees are
  not cleaned or staged by P0. Confirm no P0 build, MPI or detached process remains. Do not inspect or
  stop unrelated research processes.

- [ ] **Step 6: Main-agent seal commit and stop.** Commit the final receipt/result-only update with
  DCO subject `docs: seal Stage 4 P0 preflight`. Report exact result and, if Stage 3 is still
  unaccepted, stop. Never start `4F-0`, modify product or cherry-pick P0 into the active Stage 3 lane.

---

## 2. Parallel Dispatch and Review Protocol

After P0-1 passes, the permitted overlap is:

```text
main agent:    P0-2 build supervision -> P0-3 ABI/linkage -> P0-5 math -> P0-7 seal
worker A:      P0-4 public source/license draft
worker B:      bounded P0-2 log/hash inventory or P0-6 mechanical command inventory
```

Only one worker may modify tracked files at a time. A build-log watcher owns no tracked file. Before
dispatch, main agent freezes the task allowlist and expected evidence. After return, main agent reads
the actual diff/logs, reruns the task gate and either accepts or rejects the worker result. Cross-stage
math, licenses, complete diff and final P0 result are never delegated.

No model/reasoning override is supplied to subagents. For a bounded worker task use the configured
`agent_type="luna_worker"`; the user has selected it over `ds_worker`. If runtime turn-context
identity is unavailable, record that fact and do not count the report as independently verified
model evidence. The main agent still owns provenance, mathematics, complete diff and acceptance.

## 3. Expected Schedule and Stop Conditions

| Work | Expected wall interval | Can overlap |
|---|---:|---|
| P0-0/P0-1 governance and provenance | 0.5--1 day | no, establishes inputs |
| P0-2 Jammy builder and Cantera artifact | 0.5--2 days | P0-4 research |
| P0-3 linkage/thread/MPI/relocation | 0.5--1 day | P0-6 dry-run drafting |
| P0-4 dual-surrogate provenance | 1--2 days | P0-2/P0-3 |
| P0-5 public vectors | 1--2 days | read-only review only |
| P0-6/P0-7 intake and seal | 0.5--1 day | no final overlap |

Stop immediately with a precise receipt if:

- official Cantera/rootfs/archive identity cannot be verified；
- a consumed component or mechanism lacks acceptable source/license evidence；
- the artifact requires a newer glibc, wrong libstdc++ ABI or non-generic ISA；
- normal consumption launches Python/Conda or requires an online fetch；
- rootless Jammy builder cannot be isolated without changing the host；
- P0 would need product/Stage 3/central-authority changes；
- packaging strategy A fails and proceeding would silently select B/C；
- a file or process cannot be confidently identified as belonging to P0.

These are P0 result conditions, not prompts for an unplanned workaround. Record
`PREFLIGHT_PARTIAL` or `PREFLIGHT_REJECT`, preserve evidence and continue only with safe independent
tasks. User authorization for this plan removes routine confirmation pauses but does not authorize
scope expansion, destructive cleanup, private-source access, push or publication.
