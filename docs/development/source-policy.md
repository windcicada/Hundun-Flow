# Source and contribution policy

HUNDUN-FLOW is independently implemented.
BOFFIN is used only as a private legal comparison baseline.
COAST and COAST-2 are black-box scientific references, not source ancestors.
Stage 1 reuses zero files from those programs.
No Fortran source or legacy adapter is accepted.

Project source is Apache-2.0 and contributions require a valid DCO
`Signed-off-by` trailer. The vendored yyjson 0.12.0 files retain their MIT
license and recorded upstream commit in `third_party/yyjson/UPSTREAM.json` and
`LICENSES/yyjson-MIT.txt`.

Do not copy, translate, mechanically rewrite, imitate, or adapt source,
control flow, ABI, array layout, comments, messages, legacy input, Decomp,
Restart, or compatibility surfaces from an earlier program. Algorithms and
scientific decisions are derived from published sources and verified by
project-owned tests.

Private comparison evidence remains outside the public repository. The public
configure, build, test, and runtime paths have no Python dependency and do not
fetch source dependencies online.
