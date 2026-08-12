# Stage 4 P0 Official Input Provenance Receipt

- Recorded: 2026-08-09 (Asia/Shanghai)
- Schema: `hundun.stage4_p0.inputs.v1`
- Worker result: `PREFLIGHT_PARTIAL`
- Main-agent input-lock result: `PREFLIGHT_PASS`
- Result scope: `official_input_identity_lock_only`
- Stage 4 product accepted: `false`
- Product changes: `none`
- Base commit: `21b5a95bfe801428e339eed81ed3adaaeeff6267`
- This amendment parent HEAD: `3139044cb5310fddfba87f0093cb95e67dd917bd`
- External root: `/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0`
- Original P0-1 input file hash list: `manifests/p0-1-input-files.sha256`
- Original P0-1 input file hash-list SHA-256: `0dad0dd4d5e252dbd62cb799cdbf5bf7cff0d8dcc878bfb1eec8767f75c2f604`
- Current input manifest SHA-256: `32622f5986e97903b9b07110d091c993a3c85e9b3a78ab160493050f38f1361e`
- Current validator SHA-256: `0425f2e236d3200ee9f077dae078f24f4b6e85a55091c437cf2d9aca06c1523b`

This receipt locks official public inputs for main-agent review. It does not accept a
Stage 4 product, record any bundled binary, or authorize an ambiguous mechanism asset.
No local source patch was applied.

The main-agent result accepts only the immutable public-source, archive, Git-ref, checksum and
license-evidence lock required to start P0-2. It does not accept a binary, mechanism, package or
Stage 4 capability. Eigen remains subject to the P0-2 consumed-file license gate, and
`h2o2.yaml` remains blocked; neither status was relaxed to obtain `PREFLIGHT_PASS`.

## Immutable Target

| Field | Locked value |
|---|---|
| OS | Ubuntu 22.04 |
| glibc floor | 2.35 |
| architecture | x86_64 |
| compiler | GCC 11 |
| C++ standard library | libstdc++ |
| C++ standard | C++17 |
| libstdc++ C++11 ABI | 1 |
| ISA | x86-64 |

## Canonical Ubuntu Rootfs

Official base URL:
`https://cloud-images.ubuntu.com/minimal/releases/jammy/release/`

| Object | Bytes | SHA-256 |
|---|---:|---|
| `inputs/ubuntu/SHA256SUMS` | 1,327 | `58e962350c29fcec5f75dcd4b2f1da44c69a155526eecad4c12d5af2b75da8c1` |
| `inputs/ubuntu/SHA256SUMS.gpg` | 833 | `13d142a3675f07698752ce50bef280ebde455b9e296b1af9bd99f5a823532315` |
| `inputs/ubuntu/ubuntu-22.04-minimal-cloudimg-amd64-root.tar.xz` | 117,988,744 | `1bd7959194d3fa7b9660b169f010ff222109a885559108028c79e244d649234a` |

- Canonical metadata `Last-Modified`: `2026-08-07T02:44:02Z`.
- Verification keyring: `/usr/share/keyrings/ubuntu-cloudimage-keyring.gpg`.
- Keyring SHA-256: `2ddbc33fdd3acfa0715914e3970b6a033faade0de25985eb17995b3aa85f455e`.
- Signing key fingerprint: `D2EB44626FDDC30B513D5BB71A5D6C4C7DB87C81`.
- Signing identity reported by `gpgv`: `UEC Image Automatic Signing Key <cdimage@ubuntu.com>`.
- Signature creation time reported by `gpgv`: 2026-08-07 10:44:02 CST.
- `gpgv` exit: 0; log: `logs/ubuntu-gpgv.log`; log SHA-256:
  `cc525f5b8094dde86cc5d8e9b0f255ad3c61577beef6f4585eae5c08a9e95aea`.
- Signed rootfs checksum exit: 0; log: `logs/ubuntu-rootfs-sha256-check.log`;
  log SHA-256: `65aa15e2b89e0dc3b2686f1e2d282d2079cabcbebf8a024579a03fe9bbc7197f`.

The signed Canonical checksum record uses GNU binary-mode syntax
`<sha256> *ubuntu-22.04-minimal-cloudimg-amd64-root.tar.xz`. At P0-1 start, signed plan commit
`813670e` still used a literal-space sample filter; its expected formatting-only failure is
captured in `logs/ubuntu-rootfs-brief-regex.log` (SHA-256
`e4444568454c14d86a485e4eb17694ec4a9784bd2dac709185fb51f8476f7f7c`). Maintenance commit
`dd85bf1` changed the current command to a filename-suffix filter, which selects the same signed
`*filename` record and passes `sha256sum -c -`. The current command and cryptographic verification
both pass.

## Cantera 3.2.0 Identity

| Field | Locked value |
|---|---|
| Repository | `https://github.com/Cantera/cantera.git` |
| Tag | `v3.2.0` |
| Tag kind | annotated |
| Tag object | `3f552702cbcc1f74b473ff8543e4da2939342217` |
| Peeled commit | `4a8358eb80cfeb50474386b5f9ec0b3a83519889` |
| Official release asset | `https://github.com/Cantera/cantera/releases/download/v3.2.0/cantera-3.2.0.tar.gz` |
| Local path | `inputs/cantera/cantera-3.2.0.tar.gz` |
| Bytes | 3,073,054 |
| Required and observed SHA-256 | `a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b` |

The official Git ref and peeled commit were freshly observed with `git ls-remote`; log
`logs/cantera-v3.2.0-ls-remote.log` has SHA-256
`9956bf3ef81faa50029d49f115c2ac79285cd28e4d17682067b1b89aac84d541`. The tag object
contains the expected commit and release annotation but no cryptographic signature, so the
manifest records `not_present_tag_identity_verified_by_official_git_ref` and does not claim a
signed Cantera tag. Release-asset checksum log
`logs/cantera-release-asset-sha256-check.log` has SHA-256
`e6ceb0203718d3e5e03773bf3006264d487a7f29d04a57160002b17465900d66`.

The generated GitHub tag archive was not substituted:
`https://github.com/Cantera/cantera/archive/refs/tags/v3.2.0.tar.gz` has the distinct known
identity `f01e25e33f9d5e37db7ababe5af36b60caabff52dba04bb221d53e44735f60ec`.
The release asset's `.gitmodules`, `License.txt`, and `data/h2o2.yaml` match the corresponding
files read from the peeled Git tree; this is identity evidence for those selected files, not
a claim that the release asset and generated Git archive are byte-identical.

## Cantera Gitlinks and Official Archives

The exact `v3.2.0` tree was read before dependency download. The complete gitlink log is
`logs/cantera-v3.2.0-gitlinks.log`, SHA-256
`3ee8efac24c15a7ddf07caa29dbe6aedacf7eebc64ddededff99070065f67863`.

| Required path | Exact commit | Official tag/archive | Archive SHA-256 |
|---|---|---|---|
| `ext/fmt` | `a33701196adfad74917046096bf5a2aa0ab0bb50` | `https://github.com/fmtlib/fmt/archive/refs/tags/9.1.0.tar.gz` | `5dea48d1fcddc3ec571ce2058e13910a0d4a6bab4cc09a809d8b1dd1c88ae6f2` |
| `ext/yaml-cpp` | `0579ae3d976091d7d664aa9d2527e0d0cff25763` | `https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.7.0.tar.gz` | `43e6a9fcb146ad871515f0d0873947e5d497a1c9c60c58cb102a97b47208b7c3` |
| `ext/sundials` | `887af4374af2271db9310d31eaa9b5aeff49e829` | `https://github.com/LLNL/sundials/archive/refs/tags/archive-v5.3.0.tar.gz` | `fa9ed1c3751714fccd262f8d088261a54790ec89ae5a524399b6f06b950fe80a` |
| `ext/eigen` | `3147391d946bb4b6c68edd901f2add6ac1f31f8c` | `https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz` | `8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72` |

Official tag-ref logs and their SHA-256 values are:

| Dependency | Tag-ref result | Evidence log SHA-256 |
|---|---|---|
| fmt `9.1.0` | exact gitlink commit | `bb3d80528d3c01911ab01856d27ccfc5d95413ff7ce72c2b4bb77db28706079a` |
| yaml-cpp `yaml-cpp-0.7.0` | exact gitlink commit | `8b1c530d385a725a1db8782f190d13efbf2dd0ab3b86a2dbe3f850f0ff8cc966` |
| SUNDIALS `archive-v5.3.0` | exact gitlink commit | `d171b874d4809cbbca35e3ee531c30a2f47087cb02e37c9838420ef139603d81` |
| Eigen `3.4.0` | exact gitlink commit | `82b9586ea27165792e743fa42f002116ca8547e6ad828d10f89cfa1831243ad0` |

SUNDIALS requires the exact `archive-v5.3.0` tag: upstream `v5.3.0` points to the different
commit `cfaae739076128ad8a34f1d18cf480ecf256e24a`. It was not substituted.

The following Cantera gitlinks are recorded but excluded unless P0-2 proves consumption:

| Path | Commit | Status |
|---|---|---|
| `ext/googletest` | `e2239ee6043f73722e7aa812a459f54a28552929` | excluded unless consumed |
| `ext/HighFive` | `5513f28dcced33872a3e40a63e28d49272da20fc` | excluded unless consumed |
| `data/example_data` | `b9e0731611bc6e8d33b021332c7ade094bec1196` | excluded unless consumed |
| `ext/doxygen-awesome-css` | `df83fbf22cfff76b875c13d324baf584c74e96d0` | excluded documentation-only |

No missing submodule may be fetched silently. If P0-2 proves one of the first three is consumed,
its official archive hash and license must be added before rebuilding.

## P0-2-Discovered Builder-Only Boost Header Input

The first build using only current Cantera 3.2.0 options stopped during configuration with
`Boost could not be found`; it produced no Cantera library/product object or installed file, though
SCons did create compiler/configuration probe objects. This is a useful dependency discovery RED,
not a successful artifact build:

- Build log: `logs/p0-2-cantera-build-install.log`; exit 1; SHA-256
  `6489a62ad45f3479f13e4216f12c68d782f1ccc26062018ea26aa0386d370b4a`.
- Timing log: `logs/p0-2-cantera-build-install.time`; elapsed 1.57 s; maximum RSS 3,484 KiB;
  SHA-256 `072f541aef6673870e99a021708b63db16f46ddbf20b40e594aeeae85d573b32`.

P0-1 was therefore amended before another build. The added record is intentionally separate from
the four Cantera `ext/` source archives:

| Field | Locked value |
|---|---|
| Role | `builder_header_only` |
| Binary package | `libboost1.74-dev` |
| Source package | `boost1.74` |
| Version / architecture | `1.74.0-14ubuntu3` / `amd64` |
| Repository | Canonical Ubuntu Jammy main |
| Resolved apt URL | `http://archive.ubuntu.com/ubuntu/pool/main/b/boost1.74/libboost1.74-dev_1.74.0-14ubuntu3_amd64.deb` |
| Local path | `inputs/dependencies/boost/libboost1.74-dev_1.74.0-14ubuntu3_amd64.deb` |
| Bytes / SHA-256 | 9,608,510 / `4d9c90e43f0d25db6280d1ee326771cbb76462f73b9430f06bac1de8d05b7a78` |
| Builder dependency | already frozen `libstdc++-11-dev`; no additional required package |
| Consumption status | `required_builder_header_only_for_p0_2` |
| Consumed-header audit | `pending_p0_2` |
| Artifact bundle status | `forbidden_from_artifact_install` |
| Runtime status | `forbidden_runtime_dependency` |

The HTTP URL is the exact retrieval location recorded by Jammy apt. Trust comes from apt's signed
repository metadata and the locked package size/SHA-256, not from transport being mislabeled as
HTTPS. `p0-2-clean-apt-update.log` records the canonical apt update with SHA-256
`800c7670f55e3583ec06db5c26a02ee57e50dfc48a8b71249871a88b35ed3c89`.

The generic `libboost-dev=1.74.0.3ubuntu7` meta package was observed only to confirm resolution to
`libboost1.74-dev`; it is not a consumed manifest input and must not be installed. Its downloaded
resolution-evidence archive has SHA-256
`580eed276638d71a707d615e46747543807028545f433ea323c14b6ed30e0eca` and 3,490 bytes. Before the
amendment, both package queries confirmed neither `libboost1.74-dev` nor `libboost-dev` was
installed; logs have SHA-256 `978db882057624f227e6e4ec5485ea2bbe1560a22cafcdae13aa0d2a987fb73c`
and `2ac5aaca8ae26a510206c2d8ec65f8f35c51f1d82a0a85b545840efdd4313eb8`, respectively.

The package is not treated as uniformly BSL-1.0. Its Debian package copyright inventory contains
multiple file-scoped licenses and is frozen separately:

| Evidence | Path | SHA-256 | Status |
|---|---|---|---|
| Ubuntu package copyright member | `inputs/dependencies/boost/copyright` | `17369eeac3938acb31085c8a1d4f1a40dc88a518a9389958a12c0592bb2d5766` | `NOASSERTION`; `verified_package_license_inventory` |
| Official Boost Software License 1.0 | `inputs/dependencies/boost/LICENSE_1_0.txt` | `c9bff75738922193e67fa726fa225535870d2aa1059f91452c411736284ad566` | `BSL-1.0`; `verified_primary_license` |

P0-2 must derive and hash the complete transitive Boost-header closure from compiler depfiles or
the SCons dependency graph and map it to that package inventory. If closure completeness cannot be
proved, the artifact cannot claim a completed file-scope license audit. Independently,
`include/boost/**`, Boost CMake package files, `libboost*.so*`, `libboost*.a`, and every runtime
`NEEDED` basename starting with `libboost` remain forbidden in the install artifact.

Acquisition and license evidence logs are:

| Log | SHA-256 |
|---|---|
| `logs/p0-2-boost-builder-package-download.log` | `03f8d7afe66bb0430b8c3729217549dbc5070ed80c24dec6dc391e029124fc12` |
| `logs/p0-2-boost-builder-package-simulate.log` | `78f2f35bf538b4eb79e4999e6a2bbc012406983519ec8e7887999b47efac4d10` |
| `logs/p0-2-boost-builder-package-metadata-license-v2.log` | `866b83f7510bd5677037d58d0dae29fd2412b9a2d3e7190f0a87663e352be23d` |

### Rejected `/usr/include` injection and isolated include-root

The first post-install configuration retry followed the then-current plan literally and set
`boost_inc_dir=/usr/include`. SCons injected `-isystem /usr/include` ahead of GCC's C++ standard
headers. The `<cmath>` probe then found `/usr/include/c++/11/cmath`, but its
`#include_next <math.h>` could not reach the duplicate `/usr/include` entry that GCC had removed.
The retry stopped before a Cantera product object or install file was produced:

| Evidence | Result | SHA-256 |
|---|---|---|
| `logs/p0-2-cantera-build-install-retry2.log` | exit 1; `cmath` configuration failure | `a4833c5b4abd6827d24bb320ecef2cc655a3cbfea94bd0b55f907bcc3223ccaa` |
| `logs/p0-2-cantera-build-install-retry2.time` | 0.27 s; maximum RSS 3,364 KiB | `71c7c0c9d2744d6419b79f9273bc0cc7d2959f8be21c908dc8e57de936651933` |
| `logs/p0-2-boost-inc-dir-usr-include-diagnostic.log` | same TU: default search exit 0, injected search exit 1 | `89ba7a290cd42050279e8455ae6f00593526e51f91047b6b34c38ec2898abeae` |

The accepted configuration hypothesis uses a build-only include root whose `boost/` child is a
bwrap read-only bind of the verified rootfs `/usr/include/boost`. It exposes only the frozen Boost
header subtree and leaves GCC's final `/usr/include` search position unchanged. Before another
SCons build, the same GCC 11 flags compiled both `<cmath>` and `<boost/version.hpp>` successfully.
The resolved `boost/version.hpp` at the package, rootfs and bwrap views is 1,117 bytes with SHA-256
`90e046b8e3138a61c692abdd9bc2e45c1a95996cc5a8031cce1f110de5e64a70`:

| Evidence | Result | SHA-256 |
|---|---|---|
| `logs/p0-2-boost-isolated-include-root-diagnostic-v2.log` | `cmath_exit=0`, `boost_version_exit=0`, resolved include trace PASS | `0192b5591fd27b4035be91881dffb07c442f0b7c01db1288663c840af1e408d4` |
| `logs/p0-2-boost-version-header-hash-mapping-v2.log` | `.deb` member, rootfs and bind view byte-identical | `397e1861ae6474b5a13b3f8fb69886726e172549559e638ad577f9f7d802ce35` |
| `logs/p0-2-boost-isolated-include-root-diagnostic-v2.time` | diagnostic timing/status | `2985da18279f8728ed6dfbcbb6600169c70716d3a707526e6a18191b3c48063c` |

An earlier diagnostic invocation with an incorrect bind destination is preserved as failed
evidence and is not used for the GREEN conclusion. The isolated include root is an ephemeral build
view; it remains forbidden from the install tree and runtime closure.

## P0-2-Discovered Doxygen Builder-Tool Closure

The isolated-Boost retry reached Cantera 3.2.0's generated CLib source graph and stopped because
`doxygen` was absent. The upstream graph creates the CLib tag file and source-generation dependency
even when `doxygen_docs=no`; switching to `clib_legacy=yes` would select a different API, and a
hand-written partial install target would no longer be the official `build install` product. Both
alternatives are rejected. The retained source-graph audit is
`logs/p0-2-doxygen-builder-dependency-source-audit.log`, SHA-256
`db611d6250cb2aa444522b8e51fc9dbeadf903edce47e3c7ec07dd02fd16ffc5`.

The signed Jammy apt metadata resolves one exact six-package closure with no recommends and no
removals. The dry-run log is `logs/p0-2-doxygen-builder-package-simulate.log`, SHA-256
`17aad01f96e14dd9b293b9449605a2a62f9d28cd4f8d626841e0ad028df4760b`; its timing receipt has
SHA-256 `b70423d9e87827ac91d6e9ec9138624fd31b82d390b6ff2361259b46ed5bdc63`.

| Binary package | Source package | Version | Bytes | Package SHA-256 | Copyright-inventory SHA-256 | Runtime policy |
|---|---|---:|---:|---|---|---|
| `doxygen` | `doxygen` | `1.9.1-2ubuntu2` | 4,620,198 | `d3fe7f77f505262db0872ba810b691feb4a607887ceb368e17a1da66ee0387f7` | `286b6aadd5b010bb3df2e226dc67fa8cf5819af8870b3ac093442370214c50b9` | forbidden |
| `libclang-cpp14` | `llvm-toolchain-14` | `1:14.0.0-1ubuntu1.1` | 12,053,266 | `462546f6149fbef99b47f3e1e01e0743d256ae279337c7f89a80a976dfd02175` | `dc2a171c75ca72818f7f4312602f27cb0e4bf6d010a86e7073a3797cb228ddff` | forbidden |
| `libclang1-14` | `llvm-toolchain-14` | `1:14.0.0-1ubuntu1.1` | 6,792,182 | `b75b743f5d5effaab97790c1379fb1855d1a20bd5432a2387bf4dc82d86d45e3` | `dc2a171c75ca72818f7f4312602f27cb0e4bf6d010a86e7073a3797cb228ddff` | forbidden |
| `libllvm14` | `llvm-toolchain-14` | `1:14.0.0-1ubuntu1.1` | 23,967,046 | `9044b614a6c7fb6262e7cbeb13dc731fc0c92bed96281c1a3920dd706442ee8e` | `dc2a171c75ca72818f7f4312602f27cb0e4bf6d010a86e7073a3797cb228ddff` | forbidden |
| `libxapian30` | `xapian-core` | `1.4.18-4` | 700,928 | `5cfe52f4ca570e85efa828efda6d6831ceb0f667a32faf5438887cfaf528b7c2` | `d530c1aa427b7e55f170f48b43f3657f847255438496847796eb8a0496bdd5d6` | forbidden |
| `libxml2` | `libxml2` | `2.9.13+dfsg-1ubuntu0.12` | 764,660 | `b3678e6e4b166bc0e4226fb118d489ab51802c772914421397c9dcb2dd0e0d2b` | `ee746b96cfa5be73c3ea3e4cfb1285e9b315d4c9267f99b2ee9c5d911d9fe3f4` | `audit_pending_p0_2` |

All six records use `role=builder_tool_only`, are forbidden from the artifact install, and retain
their complete Ubuntu package copyright inventory rather than a single inferred SPDX label.
No artifact-license conclusion is inferred from these builder packages. The actual generated and
installed Cantera/sourcegen files remain subject to their own P0-2 inventory. Doxygen, LLVM/Clang
and Xapian remain build-time-only; `libxml2` is kept pending until the finished artifact's
`readelf`/`ldd` closure proves it is not consumed at runtime.

The no-install capture left the Cantera install target at zero regular files and zero symlinks.
The external evidence set is sealed as follows:

| Evidence | SHA-256 |
|---|---|
| `manifests/doxygen-builder-v1/packages.tsv` | `fc28203716b2682da5f47be409853469a1b68e23c732bf51957ee7bec8a74849` |
| `manifests/doxygen-builder-v1/report.md` | `44f1f6b1d0dddb871fa22ac29b32014fa3eeb96ae93efad8781fd1cd1c6c6224` |
| `manifests/doxygen-builder-v1/checksums.tsv` | `1e6f9ea7551773ce55b7d44ba12fdeefa961428794243f98474eeaa4d62b728a` |
| `logs/doxygen-builder-v1/checksums-verify-main-v2.log` | `e171eb57b2b713d772600386b07cd86023a902f7d6e5ac052ac9e07ad97d0b5e` |

`checksums.tsv` seals the files that existed before it was written and intentionally excludes
itself. `checksums-verify-main-v2.log` was then generated by a strict main-agent verification of
that checksum list, so it cannot appear in the earlier list without creating a second sealing
generation. Its separately recorded SHA-256 above is the post-seal verification receipt, not a
claim that it is a member of `checksums.tsv`.

## P0-2-Discovered Cantera Sourcegen Python Input

With the Doxygen closure installed, retry5 reached Cantera's official generated-CLib sourcegen and
then stopped before linking because Python 3.10 could not import `typing_extensions.Self`. The
build log `logs/p0-2-cantera-build-install-retry5.log` has SHA-256
`afae6799d7af07d864f37bb469373a07ec656616ca207eda22f3b3e5a393a533`; its timing receipt has
SHA-256 `50ec7108ee92389008eaced8f2c75862c0257e3028d6709288acc1a09a9f0dfb`. SCons had copied
252 regular files and zero symlinks into the install root before the sourcegen failure. That tree is
a failed partial-install witness, not an artifact candidate; it is isolated under
`install/failed-retry5-missing-typing-extensions`.

The frozen Cantera 3.2.0 source conditionally imports `typing_extensions.Self` on Python versions
below 3.11, while `interfaces/sourcegen/pyproject.toml` declares Python `>=3.10` and lists only
Jinja2 as a dependency. Jammy offers `python3-typing-extensions=3.10.0.2-1`; the frozen upstream
changelog records that runtime support for `Self` was added only in 4.0.0. That changelog has
SHA-256 `f61a6f3540b43f1545d42866718b3dc19b674afa4745d873d7e2827152409c8b`. The source and
metadata audit is
`logs/typing-extensions-builder-v1/source-and-metadata-audit.log`, SHA-256
`d368a764cfad4d8d48045d9c7433052eafa7008cf4263c85e35b6f6886b97e68`.

The selected input is the official PyPI pure-Python release `typing_extensions=4.15.0`:

| Field | Locked value |
|---|---|
| Role | `builder_pythonpath_only` |
| PyPI metadata | `inputs/dependencies/typing-extensions-4.15.0/pypi-metadata.json`; `e97e0b1087254aa1c7e8b2074c3796124dfd7d26e0f54ffcdc3a975b53047938` |
| Wheel | 44,614 bytes; `f0fa19c6845758ab08074a0cfa8b7aecb71c999ca73d62883bc25cc018c4e548` |
| Source archive | 109,391 bytes; `0cea48d173cc12fa28ecabc3b837ea3cf6f38c6d1136f85cbaaf598984861466` |
| Python requirement | `>=3.9` |
| `Self` API evidence | upstream sdist changelog; `f61a6f3540b43f1545d42866718b3dc19b674afa4745d873d7e2827152409c8b` |
| License | `PSF-2.0`; wheel and sdist license SHA-256 both `3b2f81fe21d181c499c59a256c8e1968455d6689d269aa85373bfb6af41da3bf` |
| Injection | exact wheel through `PYTHONPATH` |
| Rootfs install | forbidden |
| Artifact bundle/runtime | forbidden / forbidden |
| Artifact audit | `pending_p0_2` |

No `pip` command is used, the wheel is not installed into the rootfs, and no Cantera source is
patched. The plain Python 3.10 import fails as required; its log has SHA-256
`9d9aa164ccc6f4b3f3c97f64a61b7a847181d2b5acdba651c601274cb77770be`. The same interpreter
with the exact wheel on `PYTHONPATH` imports `Self`; its GREEN log has SHA-256
`8df0d02488ab88cf56b94bc590bf5a204539fdc8f3db3513cd85a02c9043e62f`.

## P0-2-Discovered Cantera CLib Sourcegen YAML Inputs

Retry6 injected the accepted `typing_extensions` wheel and passed the generated-CLib import that
blocked retry5. It then stopped in the same official sourcegen command because
`interfaces/sourcegen/src/sourcegen/_helpers.py` could not import `ruamel.yaml`. The failed build
log `logs/p0-2-cantera-build-install-retry6.log` has SHA-256
`67de7bf26a5c2153c92741829d0e045da9c905283a008cd44470708d29bd9d10`; the timing and exit-status
receipts have SHA-256 `9cdbf856ead400fa24b815161bc68b993b6820f2a29eb7fff1488c61262c9509`
and `53c234e5e8472b6ac51c1ae1cab3fe06fad053beb8ebfd8977b010655bfdd3c3`, respectively. Before the
sourcegen failure, SCons again copied 252 regular files and zero symlinks. The complete file-hash
manifest has SHA-256 `571dcb3883e96b09f371ecd6a46d71394b8ef9b728c3127cf4123df0ced7acff`;
its zero-row symlink manifest has SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`. The whole tree was moved to
`install/failed-retry6-missing-ruamel-yaml`, and the canonical install root was recreated empty.
It is failed evidence, not an artifact candidate.

The dependency is a build-system gap in the frozen upstream source, not a HUNDUN runtime need:

- `SConstruct` specifies `ruamel.yaml>=0.17.21,<1` and dispatches the CLib sourcegen even with
  `python_package=n`;
- sourcegen imports `ruamel.yaml` in exactly two files and uses only `YAML(typ="safe")` plus
  `load()` to read the frozen generator inputs;
- `interfaces/sourcegen/pyproject.toml` requires Python `>=3.10` but lists only Jinja2, while the
  top-level install-time check is bypassed after the explicit `python_package=n` selection;
- no Cantera source patch, legacy CLib fallback or custom generated header is permitted.

The frozen baseline uses the minimum versions expressly admitted by Cantera's own specifier. This
keeps the builder behavior surface narrow and is independently exercised under the exact Ubuntu
22.04 CPython 3.10 profile. Newer `0.18.x` releases are not needed by the two consumed APIs and are
therefore not introduced merely for recency. The selected main wheel explicitly requires the CLib
wheel on CPython below 3.11, so both are locked rather than silently omitting declared closure:

| Field | `ruamel.yaml` | `ruamel.yaml.clib` |
|---|---|---|
| Version / role | `0.17.21`; `builder_pythonpath_only` | `0.2.6`; `builder_pythonpath_only` |
| PyPI metadata SHA-256 | `e140839f8da85e5ba08ae2ce1d875bca7b4797a780dbf85c4319633ad33540c3` | `677b9e08e0143864232cccc1122275d0bf33f79af9aaf464ba19fc6789953c88` |
| Wheel bytes / SHA-256 | 109,478 / `742b35d3d665023981bd6d16b3d24248ce5df75fdb4e2924e93a05c1f8b61ca7` | 519,289 / `221eca6f35076c6ae472a531afa1c223b9c29377e62936f61bc8e6e8bdc5f9e7` |
| Source bytes / SHA-256 | 128,123 / `8b7ce697a2f212752a35c1ac414471dc16c424c9573be4926b56ff3f5d23b7af` | 180,695 / `4ff604ce439abb20794f05613c374759ce10e3595d1867764dd1ae675b85acbd` |
| Wheel tag | `py3-none-any` | `cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.manylinux_2_24_x86_64` |
| Python requirement | `>=3` | `>=3.5` |
| Declared dependency | `ruamel.yaml.clib>=0.2.6` for CPython `<3.11` | none |
| License | `MIT`; matching wheel/sdist SHA-256 `ab837b032c5aae84503fc0c733a116a26fd272e90dc4402fa68d3c9e51aed3b0` | `MIT`; matching wheel/sdist SHA-256 `16174d2cf8c2ee4b900bf8573106bc61c1ce0092da4fa781cc0ee81047a46539` |
| Injection | exact read-only wheel extraction through `PYTHONPATH` | exact read-only wheel extraction through `PYTHONPATH` |
| Rootfs / artifact / runtime | forbidden / forbidden / forbidden | forbidden / forbidden / forbidden |

The CLib builder wheel needs only `libc.so.6` and `libpthread.so.0`, with maximum referenced glibc
symbol `GLIBC_2.14`, below the frozen `2.35` floor. The two extracted wheel trees contain exactly
36 and 7 regular files, no symlinks, and are sealed by full file manifests with SHA-256
`ce5fef34f3fb24dd65ee21902b0d5c2b33123333c3ed95b6bd6241644702460f` and
`e1d31c44b7d6e12418793fcf1340f39b1a020efd7c095a381b9c77202415ff8a`.

The plain isolated import RED has SHA-256
`67986322bb7cf8b0aacae9b476e85b099308ab0ad8a1b01c3303c978c17b9032`. A first attempted GREEN
using the main namespace wheel directly as a zip produced the same expected import failure and is
retained under `logs/ruamel-yaml-builder-v1/import-green.log`; namespace-wheel injection is therefore
rejected. The exact extracted-wheel injection loads the compiled closure, parses a nested safe-YAML
fixture and reports `__with_libyaml__=True`; its GREEN log has SHA-256
`11002e5dd21d5ff5b1838c4b644bae6da1f4581702c849ff004fc63e1bae0d37`. The final source, metadata,
tree, ABI and no-install audit log has SHA-256
`65499efb7c76647361d87fe78f62a7ec1002b87fc855c0d4b917e1c73fb665c2`. No `pip` command or rootfs
installation is used.

## License Lock

| Component | SPDX / status | License evidence SHA-256 | Holder and redistribution obligation |
|---|---|---|---|
| Cantera | `BSD-3-Clause`; verified redistributable | `e92980b9712ce20e73898a97b0116889e84e07f548d6be8591e87dcad79c41bb` | Caltech, Sandia Corporation, Cantera Developers; retain/reproduce notice, conditions and disclaimer; no endorsement |
| fmt | `MIT`; verified redistributable | `825c9324e70f8c839c8ba910543dd4a7daee243b86ef960594c11381a19980b8` | Victor Zverovich; include copyright and permission notice in copies/substantial portions; upstream includes an optional compiled-object exception |
| yaml-cpp | `MIT`; verified redistributable | `aa6fcc27be034e41e21dd832f9175bfe694a48491d9e14ff0fa278e19ad14f1b` | Jesse Beder; include copyright and permission notice in copies/substantial portions |
| SUNDIALS | `BSD-3-Clause`; verified redistributable | license `fe9b98949f9ac131af1130eb3301b8433d23145367ccebaf23204f2322c4d882`; notice `5cb9c6b0ca7c09b26fe2d1f9ef45231086c9ea2cdb32eb5f9dce7012e2364bc2` | Lawrence Livermore National Security and Southern Methodist University; retain/reproduce license and disclaimer, preserve notice, no endorsement |
| Eigen | primary `MPL-2.0`; `candidate_consumed_file_scope_audit_required` | guide `c83230b770f17ef1386ea1fd3681271dd98aa93646bdbfb5bff3a1b7050fff9d`; MPL-2.0 `fab3dd6bdab226f1c08630b1dd917e11fcb4ec5e1e020e2c16f83a0a13863e85` | Eigen contributors and file-specific third-party holders; preserve file notices and make corresponding MPL-covered source available when applicable; P0-2 must audit files actually included/compiled |
| Boost builder headers | package inventory `NOASSERTION`; primary `BSL-1.0`; `candidate_consumed_header_scope_audit_required` | package inventory `17369eeac3938acb31085c8a1d4f1a40dc88a518a9389958a12c0592bb2d5766`; primary license `c9bff75738922193e67fa726fa225535870d2aa1059f91452c411736284ad566` | Boost contributors and file-specific holders; P0-2 must prove and audit the complete transitive consumed-header closure; no Boost payload/runtime is allowed in the artifact |
| Doxygen builder-tool closure | package inventories retained as `NOASSERTION`; builder-only | Doxygen `286b6aadd5b010bb3df2e226dc67fa8cf5819af8870b3ac093442370214c50b9`; LLVM packages `dc2a171c75ca72818f7f4312602f27cb0e4bf6d010a86e7073a3797cb228ddff`; Xapian `d530c1aa427b7e55f170f48b43f3657f847255438496847796eb8a0496bdd5d6`; libxml2 `ee746b96cfa5be73c3ea3e4cfb1285e9b315d4c9267f99b2ee9c5d911d9fe3f4` | Preserve each package's file-scoped inventory; never bundle Doxygen/LLVM/Clang/Xapian; decide libxml2 runtime status only from the P0-2 artifact closure |
| typing_extensions builder wheel | `PSF-2.0`; verified matching wheel/sdist license | `3b2f81fe21d181c499c59a256c8e1968455d6689d269aa85373bfb6af41da3bf` | Python Software Foundation and historical licensors named in `LICENSE`; retain exact license with builder-source provenance; wheel is forbidden from artifact/runtime payload |
| ruamel.yaml builder wheel | `MIT`; verified matching wheel/sdist license | `ab837b032c5aae84503fc0c733a116a26fd272e90dc4402fa68d3c9e51aed3b0` | Copyright 2014–2022 Anthon van der Neut, Ruamel bvba; retain copyright and permission notice; builder-only and forbidden from artifact/runtime payload |
| ruamel.yaml.clib builder wheel | `MIT`; verified matching wheel/sdist license | `16174d2cf8c2ee4b900bf8573106bc61c1ce0092da4fa781cc0ee81047a46539` | Copyright 2014–2022 Anthon van der Neut, Ruamel bvba; retain copyright and permission notice; builder-only and forbidden from artifact/runtime payload |

The complete Eigen license pack is retained because upstream documents file-scoped BSD,
MINPACK, Apache-2.0, GPL and LGPL material in addition to its primary MPL-2.0 coverage:

| File | SHA-256 |
|---|---|
| `COPYING.BSD` | `51928dce36213c5333ba3172e847d735d4c6e9b7ff2722a326c49067155b82eb` |
| `COPYING.MINPACK` | `c87b7f8ee88f6195e91743820c00354833583aef091b72e2d4a49c8e28e798a0` |
| `COPYING.APACHE` | `03379001a7b12a2ec997a25554247d985270b353c10d5bafee9ac8d6519820b7` |
| `COPYING.GPL` | `8ceb4b9ee5adedde47b31e975c1d90c73ad27b6b165a1dcd80c7c545eb65b903` |
| `COPYING.LGPL` | `dc626520dcd53a22f727af3ee42c770e56c97a64fe3adb063799d8ab032fe551` |

Presence of a license text is not evidence that its file-scoped code will be consumed. Binary
acceptance remains pending P0-2's compiled-file and runtime inventory.

P0-2 must compile with `-DEIGEN_MPL2_ONLY`. Eigen's own
`Eigen/src/Core/util/NonMPL2.h` turns an attempted non-MPL2 include into a compile error under that
definition. Such a failure is a design blocker, not permission to remove the definition. The
post-build consumed-file audit remains required even when this compile-time guard passes.

## Mechanism Asset Separation

`h2o2.yaml` is a separate candidate asset, not Cantera code-license material:

| Field | Value |
|---|---|
| Official source container | Cantera release asset at `data/h2o2.yaml` |
| Local path | `inputs/dependencies/mechanisms/h2o2/h2o2.yaml` |
| SHA-256 | `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4` |
| Upstream description | Hydrogen-Oxygen submechanism extracted from GRI-Mech 3.0 and modified to include N2 |
| SPDX | `NOASSERTION` |
| License status | `candidate_user_supplied` |
| Consumption status | `blocked_until_asset_specific_redistribution_permission` |

The accompanying Cantera data disclaimer has SHA-256
`067c221d01413542325424c7b9c6708674e67c4b5a346c58dfdea891ea81754f` and says the
distributed input files are illustrative. It does not supply an asset-specific redistribution
grant. Cantera's BSD license is therefore not inherited by this mechanism.

## Validator Mutation and Policy Scan

RED command:

```bash
cmake \
  -DMANIFEST=.superpowers/stage4-p0/input-manifest.json \
  -DEXTERNAL_ROOT=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0 \
  -DEXPECTED_CANTERA_SHA=0000000000000000000000000000000000000000000000000000000000000000 \
  -P .superpowers/stage4-p0/validate-manifest.cmake
```

- Observed exit: 1.
- Required diagnostic: `Cantera archive SHA mismatch`.
- Log: `logs/validate-manifest-red.log`.
- Log SHA-256: `cec64733e82eb00f606d69bbf30f06e746533105778a0930914dd522cbe5e93f`.

GREEN command:

```bash
cmake \
  -DMANIFEST=.superpowers/stage4-p0/input-manifest.json \
  -DEXTERNAL_ROOT=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0 \
  -DEXPECTED_CANTERA_SHA=a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b \
  -P .superpowers/stage4-p0/validate-manifest.cmake
```

- Observed exit: 0.
- Diagnostic: `Stage 4 P0 input manifest validated`.
- Log: `logs/validate-manifest-green.log`.
- Log SHA-256: `56cb61c6be6d2559a8cf7a4ee92d089404b2de6ecdad59f53422eb8831b91c94`.

### Canonical-path boundary regression

`verify_external_file` canonicalizes both `EXTERNAL_ROOT` and every selected external file with
`file(REAL_PATH)`. It compares the resolved file against a canonical root prefix ending in `/`, so
`/path/root-other` cannot pass as a child of `/path/root`. Hashing occurs only after this
component-safe containment check.

The fixture generator command was:

```bash
cmake \
  -DSOURCE_MANIFEST="$PWD/.superpowers/stage4-p0/input-manifest.json" \
  -DEXTERNAL_ROOT=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0 \
  -DOUTPUT_MANIFEST=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests/validate-symlink-escape-manifest.json \
  -DINSIDE_OUTPUT_MANIFEST=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests/validate-symlink-inside-manifest.json \
  -P /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests/p0-1-make-symlink-escape-fixture.cmake
```

- Generator exit: 0.
- Generator SHA-256: `f2a9da703d95ed62f5689dca4516f602ee2f67409eac88a4c5cea87ff9925e5a`.
- Generator log: `logs/validate-manifest-symlink-fixture.log`; log SHA-256:
  `63e9a6042f617cc9cbacca98facd3393a2921e6a35364509b028d50e9051e8d7`.
- Escape manifest SHA-256: `7f8334901a76efb9a72fb5864c45222f4386c72b4574a0dd5956e15f50f49470`.
- Inside-root manifest SHA-256: `56e3da62958b397123a372230604b3c7f8ac730c6858aca73368ad303233619a`.

The escape manifest changes only the signed-sums path/hash to a root-internal symlink:

```text
inputs/validator-fixture/outside-root-keyring.gpg
  -> /usr/share/keyrings/ubuntu-cloudimage-keyring.gpg
```

The target is the already-installed, read-only system keyring; no object outside the P0 external
root was written, replaced, or moved. Before the fix, the validator followed this symlink, matched
the keyring SHA and returned exit 0. The TDD assertion therefore failed as intended:

- Pre-fix raw validator exit: 0.
- Diagnostic: `TDD_RED: expected outside-root symlink rejection, but old validator accepted the manifest`.
- Log: `logs/validate-manifest-symlink-escape-before-fix.log`; log SHA-256:
  `8fc339ca853fa8d135cafe2471700f108d123355836c4ed1ca2d79a89d25119f`.

After the fix, the same validator command and manifest produce:

- Raw validator exit: 1.
- Diagnostic: `Ubuntu signed_sums resolves outside EXTERNAL_ROOT`, including both the manifest
  path and resolved `/usr/share/keyrings/ubuntu-cloudimage-keyring.gpg` target.
- Log: `logs/validate-manifest-symlink-escape-red.log`; log SHA-256:
  `701bb81649f9bc019cedd3399fad3d25eccb902d0da9b48700f8e6cdb1e67651`.

A second manifest selects the root-internal symlink
`inputs/validator-fixture/inside-root-SHA256SUMS`, which resolves to the real signed checksum file
under `inputs/ubuntu/`. It remains legal and validates with exit 0:

- Log: `logs/validate-manifest-symlink-inside-green.log`; log SHA-256:
  `56cb61c6be6d2559a8cf7a4ee92d089404b2de6ecdad59f53422eb8831b91c94`.

The manifest's system keyring remains intentionally outside the external-root file helper. The
validator requires the manifest path to equal
`/usr/share/keyrings/ubuntu-cloudimage-keyring.gpg` and the manifest SHA to equal
`2ddbc33fdd3acfa0715914e3970b6a033faade0de25985eb17995b3aa85f455e`, then applies
`file(REAL_PATH)` and hashes that resolved system file. The real-manifest GREEN run proves this
separate check was not incorrectly subjected to the external-root rule.

### Exact identity and schema lock mutations

The review mutation generator ran with exit 0:

```bash
cmake -DSOURCE_MANIFEST="$PWD/.superpowers/stage4-p0/input-manifest.json" \
  -DOUTPUT_DIRECTORY=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests \
  -P /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests/p0-1-make-lock-mutation-fixtures.cmake
```

Generator log `logs/validate-lock-mutation-fixtures.log` has SHA-256
`c6f4b2843298974f2cd8b3fc5e1bb0d6944bfbcb3efcce844bca6218cbca93bb`. Before the exact locks,
all four manifests incorrectly validated with raw exit 0, providing the TDD RED. After the fix,
the same command template

```bash
cmake -DMANIFEST=<mutation.json> \
  -DEXTERNAL_ROOT=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0 \
  -DEXPECTED_CANTERA_SHA=a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b \
  -P .superpowers/stage4-p0/validate-manifest.cmake
```

produces these required failures:

| Mutation manifest SHA-256 | Exit / diagnostic | Log SHA-256 |
|---|---|---|
| keyring path+SHA to `/etc/os-release`: `a0919dbf752e0c987afe1cb37151f8c816a83909e1641d42d69038974f744c9e` | 1 / `Ubuntu keyring path mismatch` | `31fb962beba1bb76f1aae55d32f5aeee75d4b2c97c96c1038dcc1af1d9e2bb75` |
| Eigen license status to `verified_redistributable`: `66309784c082cd5122f5340be341e4a7ba705f8696b001f743620b730ba113be` | 1 / `eigen license status mismatch` | `5b00e16e7fa77a62a5a607909acb5a206dcab1f5430c524d0478e306b7bb9399` |
| duplicate second mechanism: `d34db35bb8b81469973277aa3f527a1eb5a39d8dfb2e62a3ce85bf409ccd17c9` | 1 / `Mechanism count mismatch` | `5e8ef7f986500aad4616ce9a256d202fa198dbcd8ac199f5f6af23ffb455b993` |
| h2o2 consumption status to `unblocked`: `7fe14375b615a48c86d0e138113eb13f6fe7db67735d55b90d985ace3ac871a5` | 1 / `Mechanism consumption status mismatch` | `2ec7499a4715c683348c437020e43333d7632d0668caca7e086c719bcc2fdf6d` |

The validator also exactly locks Cantera/fmt/yaml-cpp/SUNDIALS license status, Eigen's candidate
license/consumption/bundle statuses and required `EIGEN_MPL2_ONLY` obligation, and the sole h2o2
mechanism's name, phase, SPDX, blocked status, empty license fields, content hash and terms hash.

### Boost builder-input schema RED and mutations

After the Boost record was first added, the previous validator still returned exit 0 because it did
not inspect the new array. This is the mutation-sensitive schema RED; log
`logs/validate-boost-schema-before-fix.log` has SHA-256
`48f2f12bf8289348d7d52b474f0faf7a6e9d0a5446001e11dfe8498efb3480ec`. After adding exact
builder-input checks, the same canonical manifest returns exit 0; log
`logs/validate-boost-input-green-v2.log` has the same hash because the validator's one-line success
output is intentionally unchanged.

The external fixture generator has SHA-256
`b3227206cd3aa5561fc69275070997b34a2a68cd4897603207a991041513fa6e`; its execution log has
SHA-256 `49d0b6edd91c07a1c4d2cf819d2f8bee93848c8a8ab4691f1a18e2f84dd21396`. All ten mutations are
rejected:

| Mutation manifest SHA-256 | Exit / diagnostic | Log SHA-256 |
|---|---|---|
| wrong Boost archive SHA: `1eccd3f39e5f97aff3c1b7b57445ea34ed32ef07e00bc6d95f3ed938cd44051f` | 1 / `Boost archive manifest SHA mismatch` | `ceb8e757c86d027012353595d5bd106cdd7cae1855a5b0a27378514321a5ff04` |
| wrong version: `7bda543e05f9e1e30000fe02a6d833107ba189051bc9ffca7106ea8eda67c1a6` | 1 / `Boost version mismatch` | `7fdb749a1a04f2c4856deb5ebed15a3557674c0107aa5645a97fbc6db751c383` |
| wrong role: `8c93a39303f1ac1a16d9c197fe4335adb564dcc1c3fe3756a62f1d2553daf155` | 1 / `Boost builder dependency role mismatch` | `295dcb6979f88e0804017c31cb903d64f326f75303083d93d5f7c13387c4ebb8` |
| bundle allowed: `e16d8604085820c4ac60bf718a770a07891c673824457dd8a9b574e504e98a92` | 1 / `Boost bundle status mismatch` | `4902f111c8ef6a05e23cdc6e13ba3976f09d5435ec51f73b1a7272501cda4af0` |
| runtime allowed: `ae68a8aed9d8d2bb754b398f4d8d7b53a80876588a2992314443532cc72fd4fd` | 1 / `Boost runtime status mismatch` | `5c570b1bb287f4d55478dcef9658584c6ebfbd117b2573082d27e0bd03283b9d` |
| missing builder array: `53497de4df5a01bdf6b1f56c1c71e4c9b00cf3d3a07ea7f9deb88068abe58549` | 1 / `Invalid builder-header dependency array` | `40d2968d8dce8ed3b849010e11e4ecca46935aec92170d1b02a57207717010d2` |
| meta package substituted: `705f9689e41a089349ee702698489d241552740e4f161cc1c2c677af398ebe56` | 1 / `Boost builder package mismatch` | `2b5f7cf2bdc327f32106db663b379831eb7cbf8fd6a992d58984122d756378a5` |
| Boost path symlink escapes external root: `f6fb16bbf5d4eeacac66ac9da7df5c3b3be031b52f9cfa70618d46c25d3ce3ae` | 1 / `Boost builder package resolves outside EXTERNAL_ROOT` | `b342eab25977fe764c93ffeb6daff32406f84181f3d1233442cb68d4ed88f06f` |
| license scope prematurely called verified: `2692163cbfc38e9f0d9deaceabe3a50e66358c52bd9e464ad1f3f211d522f833` | 1 / `Boost license status mismatch` | `e4578988364fc90e1db8ec8df984bfb9774b22c6a74e3b8986684e8ca18504db` |
| consumed-header audit prematurely called verified: `0a76fa19afc5011fc88b9c2008a72abcdea61ee04fe01ed7ceb85093efca3dad` | 1 / `Boost consumed-header audit status mismatch` | `84c5051284dbe7088e78bb2b732768638f4ee32d3cc1d3eebe3a13b8ce62ddc5` |

### Doxygen builder-tool schema RED and mutations

After the six-package array was added, the pre-amendment validator returned exit 0 without reading
it. This is the schema RED; `logs/validate-doxygen-builder-schema-before-fix.log` has SHA-256
`48f2f12bf8289348d7d52b474f0faf7a6e9d0a5446001e11dfe8498efb3480ec`. The amended validator
locks array size and order, package/source/version/repository/URL/path/size/hash, copyright
inventory, role, consumption, audit, bundle/runtime and empty patch state. The final canonical
GREEN log `logs/validate-doxygen-builder-green-v3.log` has the same one-line-output SHA-256
`48f2f12bf8289348d7d52b474f0faf7a6e9d0a5446001e11dfe8498efb3480ec`.

The fixture generator has SHA-256
`36627d79b6edcc72ff544c4836d6f2600f365bf8b03866310304f049878c123a`; the runner has SHA-256
`3e8a98cf4968736429e1b9114a0a0a8cb4a2947be51a221879c0fb89a0ac813c`. The final mutation
summary `logs/doxygen-builder-v1/validator-mutations-v3/summary.tsv` has SHA-256
`f54691035560ebc550ed2632b2c4fb8a2154a2cfeafc166a2fe68a5d1937d604` and records
`failure_count=0` for all 14 rejection cases:

| Mutation | Exit / required diagnostic |
|---|---|
| wrong Doxygen archive SHA | 1 / `doxygen archive manifest SHA mismatch` |
| wrong libclang-cpp14 version | 1 / `libclang-cpp14 version mismatch` |
| wrong libclang1-14 role | 1 / `libclang1-14 builder dependency role mismatch` |
| wrong libllvm14 source package | 1 / `libllvm14 source package mismatch` |
| wrong libxapian30 package | 1 / `libxapian30 builder package mismatch` |
| Doxygen bundle allowed | 1 / `doxygen bundle status mismatch` |
| libclang-cpp14 runtime allowed | 1 / `libclang-cpp14 runtime status mismatch` |
| libxml2 runtime decided prematurely | 1 / `libxml2 runtime status mismatch` |
| Doxygen license status collapsed | 1 / `doxygen license status mismatch` |
| builder-tool array removed | 1 / `Invalid builder-tool dependency array` |
| Doxygen copyright inventory removed | 1 / missing `copyright_inventory.spdx` |
| Doxygen copyright SHA changed | 1 / `doxygen copyright inventory manifest SHA mismatch` |
| archive path replaced by outside-root fixture | 1 / `doxygen archive path mismatch` |
| artifact audit marked complete before P0-2 | 1 / `doxygen artifact audit status mismatch` |

The archive-path mutation is rejected by the stronger exact-path identity lock before filesystem
resolution. The earlier canonical-path regression independently proves that
`verify_external_file` rejects a root-internal symlink resolving outside `EXTERNAL_ROOT` while
accepting an internal symlink whose target remains inside the root.

### typing_extensions builder-input schema RED and mutations

After the builder-Python record was added, the previous validator still returned exit 0 without
reading it. The schema RED log `logs/validate-typing-extensions-schema-before-fix.log` has SHA-256
`48f2f12bf8289348d7d52b474f0faf7a6e9d0a5446001e11dfe8498efb3480ec`. The amended validator
locks the one-row array, PyPI identity, metadata/wheel/sdist paths, sizes and hashes, Python
requirement, matching license evidence, injection mode, no-rootfs/no-bundle/no-runtime policy,
pending artifact audit and empty patch set. Its canonical GREEN log
`logs/validate-typing-extensions-green-v3.log` has the same one-line-output SHA-256
`48f2f12bf8289348d7d52b474f0faf7a6e9d0a5446001e11dfe8498efb3480ec`.

The fixture generator has SHA-256
`790c2cf733498d1c0ed06cdd3df76d57efcccc693a9e13153e27d601db78a04f`; the runner has SHA-256
`e410ae6a217c10bfb5d37295fd50947f7e1dc3c032015371409ecd3f1b0d73c2`. The final 16-case summary
`logs/typing-extensions-builder-v1/validator-mutations-v3/summary.tsv` has SHA-256
`926fb59c8463fe8ffd21a7b461b0db595615127150183f3b5d04ecda34ce4632` and records
`failure_count=0`. It rejects a missing array; wrong name, version or role; wrong metadata, wheel,
sdist or license hash; wrong Python requirement or SPDX; `pip` injection; rootfs install; bundle or
runtime allowance; premature audit completion; and a wheel path resolving outside the external
root.

Because this amendment expands the builder-Python array from one row to three, the same 16 fixtures
were regenerated from the current manifest and rerun against the current validator. The three-row
schema summary `logs/typing-extensions-builder-v1/validator-mutations-v4-three-row-schema/summary.tsv`
has SHA-256 `9e6ffea5f17c0e396d9b980c9b704f1282379e022d59279dbf13a9222247f6a8` and again records
`failure_count=0`; the v3 result remains historical evidence rather than being overwritten.

### ruamel.yaml builder-input schema RED and mutations

After the two YAML builder records were added, the previous validator rejected the three-row array
because it knew only the earlier single-row schema. This is the schema RED;
`logs/validate-ruamel-schema-before-fix-v2.log` has SHA-256
`8fc08abcfe3d8398c55b1f456cb2e5ec77165bf532ebd5b076751712a13fd3a7`. A prior insertion-location
mistake was also preserved in `logs/validate-ruamel-schema-before-fix.log`; it is not cited as the
semantic RED.

The amended validator locks all three builder-Python rows in order, all ruamel PyPI and artifact
identities, exact relative paths, wheel/sdist sizes and hashes, wheel METADATA, matching licenses,
the upstream dependency expression, CPython/manylinux ABI, no-rootfs/no-bundle/no-runtime policy,
empty patch sets and the two complete extracted-wheel trees. The main-agent-reviewed canonical
GREEN log `logs/validate-ruamel-schema-after-main-review-v4.log` has SHA-256
`48f2f12bf8289348d7d52b474f0faf7a6e9d0a5446001e11dfe8498efb3480ec`.

The mutation generator has SHA-256
`1656b53e3aae2d2970eb1492850869b383a981206ae21b4ed9b46cdbc60c226f`; the runner has SHA-256
`620920d1d6c25fe4bec8264cc5ced4bbd7b7b62f57b9bdde8c72364dcdb23ac1`. The 25-case summary
`logs/ruamel-yaml-builder-v1/validator-mutations-v1/summary.tsv` has SHA-256
`e837ee60237ec4ef703bd3896c5308d1b61131518b2f6e0bdb80df1befb72bb9` and records
`failure_count=0`. It rejects missing schema, altered names/versions/roles, wrong paths or hashes,
lost declared closure, changed license, non-approved injection, rootfs/bundle/runtime allowance,
premature audit completion, tree-manifest changes, local patches and CLib ABI drift.

The policy scan found no `file://` URL, private COAST/BOFFIN token, unpinned
`main`/`master`/`latest` reference, missing required license status, or executable product
download command. It used the repository's syntax-aware
`hundun_assert_public_dependency_policy` function; quoted detector regular expressions in the
guard itself were not treated as executable commands. The policy log is
`logs/p0-1-policy-scan.log`, SHA-256
`3f5b604c2e819285848acee8ceae182c351dd4390db82a10a9c7d92f7b7cd3b8`. The same scan was rerun
after the Doxygen amendment; `logs/p0-1-policy-scan-doxygen-amendment.log` has SHA-256
`e0bc7c9bae8cf310bc9419fc5bd4e887a5a5bd6ff81fd0121dcca2c4aeb5a240`. The same one-line result
was obtained after the typing_extensions amendment; its log
`logs/p0-1-policy-scan-typing-extensions-amendment.log` has the same SHA-256.
The ruamel amendment produced the same product-policy result in
`logs/p0-1-policy-scan-ruamel-amendment.log`, again SHA-256
`e0bc7c9bae8cf310bc9419fc5bd4e887a5a5bd6ff81fd0121dcca2c4aeb5a240`. Its separate input-manifest
boundary scan has SHA-256 `8008af35367ede4b7e51bb012ac09934b969653008358420ccb8b9ae7fc5ca1e` and
confirms no file URL, private source/data path, unpinned branch token, `pip` injection or
rootfs/bundle/runtime allowance.

## Boundary and Review Conclusion

- Original P0-1 writes were confined to `inputs/`, `source/`, `logs/`, and `manifests/` under the
  Stage 4 P0 external root. The initial P0-2 Boost and Doxygen dependency-discovery attempts entered
  their isolated `build/` paths and stopped before a successful install; those probes remain failure
  evidence. Retries 5 and 6 later reached the upstream install actions before sourcegen failed,
  each producing a separately recorded 252-file partial install. Both partial trees are rejected and
  isolated under explicit failed-attempt paths; the canonical install root is empty. No `spikes/`
  path was entered by this amendment.
- No product source or product CMake file was changed.
- No accepted Cantera artifact binary exists yet. The Boost and Doxygen-closure `.deb` files,
  typing_extensions wheel and ruamel wheel pair are immutable builder inputs, not accepted bundle or
  runtime components, so this receipt makes no artifact/runtime closure claim. The retry5 and retry6
  partial installs are explicitly rejected as candidates.
- P0-0 supplied only the external-root and policy boundary. Rootfs, Git refs, archives, hashes,
  licenses and mechanism conclusions above were independently observed for P0-1.
- The worker handed off `PREFLIGHT_PARTIAL`; the main agent completed the original provenance
  review, committed its three tracked artifacts at
  `7fb02ed28ce47d3e762b02950495e4047fc57e5b`, and set only the official input identity lock to
  `PREFLIGHT_PASS`. This receipt now extends that input lock with the P0-2-discovered Boost header,
  Doxygen builder-tool closure and exact sourcegen Python inputs before any further retry. Retry6
  proved the YAML dependency gap and its rejected partial install has now also been isolated. Binary,
  mechanism, package and Stage 4 capability acceptance remain false.
- Downstream gates remain explicit: P0-2 must audit actually consumed Eigen files and any newly
  consumed excluded gitlink; it must also prove the complete transitive Boost-header closure,
  retain the file-scope audit in the artifact evidence, and require both
  `boost_payload_file_count=0` and `boost_runtime_needed_count=0`. It must also prove zero
  Doxygen/LLVM/Clang/Xapian payload/runtime entries and resolve `libxml2` from the finished binary
  closure. It must prove zero Python/typing_extensions payload and no `libpython` dynamic
  dependency. The accepted P0-2 artifact must additionally prove zero ruamel.yaml/CLib payload and
  no dependency on the builder-only `_ruamel_yaml` extension. `h2o2.yaml` remains blocked absent
  asset-specific permission or a user-supplied replacement.
