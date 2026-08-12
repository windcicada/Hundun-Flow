# Source and contribution policy

HUNDUN-FLOW is independently implemented.
BOFFIN is used only as a private legal comparison baseline.
COAST and COAST-2 are private scientific references, not product source
ancestors. Stages 1--3 reuse zero files from those programs. Stage 5/6 may use
the narrowly controlled, separate-process oracle described below; no private
source or legacy adapter is accepted into Git, the product or an installation.

Project source is Apache-2.0 and contributions require a valid DCO
`Signed-off-by` trailer. The vendored yyjson 0.12.0 files retain their MIT
license and recorded upstream commit in `third_party/yyjson/UPSTREAM.json` and
`LICENSES/yyjson-MIT.txt`.

Cantera and its transitive dependencies are third-party components rather than
HUNDUN original source. Each bundled source archive or binary must retain its
upstream identity, copyright, license and disclaimer, and record the official
URL, release/tag, commit, SHA-256, local patches, ABI profile and transitive
licenses. A Cantera mechanism file is a separate copyright object; it may be
redistributed only under its own recorded permission.

Do not copy, translate, mechanically rewrite, imitate, or adapt source,
control flow, ABI, array layout, comments, messages, legacy input, Decomp,
Restart, or compatibility surfaces from an earlier program. Algorithms and
scientific decisions are derived from published sources and verified by
project-owned tests.

The Stage 5 ESF/TCR and Stage 6 fuel-provenance tasks permit one limited COAST
oracle workflow after the user confirms the exact current source realpath and
version. The governance manifest must identify an allowlist of pure
mathematical modules. A runner may copy those modules into an untracked,
generated directory outside the repository and compile them with a
HUNDUN-owned standalone driver for synthetic differential tests. The oracle
runs in a separate process and exposes only numeric input/output records. Its
source, build tree, cases, data, executable, messages, ABI and control flow do
not enter Git, the product, installed artifacts or public tests. A standalone
test driver may use Fortran where required to call the private oracle, but it
must not make Fortran a product build or runtime dependency.

Private comparison evidence remains outside the public repository. The normal
HUNDUN configure, source build, install, test, formal acceptance and runtime
paths have no Python dependency and do not fetch source dependencies online.
A maintainer-only Cantera artifact producer may use the frozen upstream
Python/SCons toolchain inside its isolated release-builder environment when
upstream requires it. Normal CMake never invokes that producer, and neither its
toolchain nor Python enters the HUNDUN package or acceptance runtime.

“Offline” and “network-independent” mean that the build uses already pinned
local inputs and performs no configure-time fetch. Tests must not disconnect
the host network, because doing so can interrupt the development service.
