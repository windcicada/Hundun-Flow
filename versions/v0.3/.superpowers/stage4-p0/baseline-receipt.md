# Stage 4 P0 Baseline Receipt

Recorded at: `2026-08-09T17:18:52+08:00`

## Authority

```text
plan_head=813670efc2a0ce6adb1a033fb98b7582b05a0fce
plan_tree=a9f56697f5313bb94fab8f7d5eca820dee1275c2
branch=coast/stage4-p0-preflight
worktree=/home/wyf/code_dev/.worktrees/hundun-flow-stage4-p0-preflight
git_dir=/home/wyf/code_dev/hundun-flow-governance/.git/worktrees/hundun-flow-stage4-p0-preflight
git_common=/home/wyf/code_dev/hundun-flow-governance/.git
external_root=/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0
external_root_allowlist=inputs,source,build,install,spikes,logs,manifests
other_external_write_roots=forbidden
stage4_product_accepted=false
product_changes=none
```

The execution worktree was clean at intake. The two pre-amendment authority
queries both returned exit `1` as the intended RED:

```text
rg -q 'pre-Stage-4 P0' AGENTS.md
rg -q 'P0 preflight candidate' docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md
```

## Host and isolation probe

```text
architecture=x86_64
host_glibc=2.31
host_gcc=gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0
cmake=3.31.12
mpi=OpenRTE 2.1.1
bubblewrap=/usr/bin/bwrap
bubblewrap_unprivileged_user_probe_exit=0
```

The host is not the frozen Ubuntu 22.04/glibc 2.35/GCC 11 release profile.
No artifact built directly against the host can satisfy P0-2. The approved
builder route is a Canonical signed Jammy rootfs launched with unprivileged
bubblewrap; it must not alter the host or disconnect its network.

## Repository state at intake

```text
product_path=/home/wyf/code_dev/hundun-flow
product_branch=main
product_status=clean

governance_main_path=/home/wyf/code_dev/hundun-flow-governance
governance_main_branch=main
governance_main_status=user_untracked_stage7_plan_only
governance_main_untracked=docs/plans/2026-07-21-hundun-flow-stage7-neural-warm-start.md
```

The unrelated Stage 7 file is user-owned and must remain untouched.

## Active Stage 3 observations

These values are a non-interference snapshot, not accepted Stage 3 identities
and not Stage 4 interface authority. Another agent may legitimately change them.

| Worktree | Branch | HEAD | Modified | Untracked | Total |
|---|---|---|---:|---:|---:|
| `/home/wyf/code_dev/.worktrees/hundun-flow-stage3` | `stage3-static-lfp-gcibm-wale` | `64e3a70f49da4ccb3585d0eb6afa9966fb82282d` | 52 | 258 | 310 |
| `/home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework` | `coast/stage3-framework-completion` | `cd7f00c0d6749d40c200048f09144237de790244` | 15 | 2 | 17 |
| `/home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure` | `coast/stage3-infrastructure-lane` | `3735742b65d8372d54d6ac517eac932281f37fe2` | 6 | 1 | 7 |

P0 must not edit, stage, clean, build in or infer accepted interfaces from
these worktrees. It also does not inspect or signal unrelated research
processes.
