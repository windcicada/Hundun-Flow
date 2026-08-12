# Cantera local patch ledger

HUNDUN's pinned Cantera 3.2.0 artifact is built with **zero local source
patches**. The source identity is the official `v3.2.0` release asset recorded
in `UPSTREAM.json`.

Build flags, projection rules, debug stripping, and relocation metadata are
packaging operations. They do not alter upstream source and are not recorded
as source patches. Any future source patch requires a new ledger entry, hash,
license review, artifact identity, and Stage acceptance.
