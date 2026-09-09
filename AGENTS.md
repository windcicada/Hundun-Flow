## Agent skills

### HUNDUN-FLOW workflow

- As requested on 2026-09-05, suspend use of `ponytail` for HUNDUN-FLOW until the user explicitly re-enables it.
- Before exploring unfamiliar code, check the active checkout with `codegraphf status --json .`; use `codegraphf sync .` when the index is stale. Start with `context`, `query`, `callers`, `callees`, or `impact`, then use `rg` and source reads to verify results. State when CodeGraphF is unavailable and use the fallback.
- Index the actual active worktree, not another checkout of the same repository. Use the CLI if the current session does not expose the configured MCP tools.
- Use short file and directory names and a shallow layout for cases, builds and debugging. Reuse existing flat working directories such as `check` and `cases/g`.
- Remove obsolete generated source and debug artifacts when retiring an experiment, without making backups. Preserve files belonging to active work.

### Issue tracker

Issues and specs live as GitHub issues, operated with the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Default triage vocabulary is in use: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context layout: one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.

### Commit identity

- Use `WANG YUDONG <wangyudong@buaa.edu.cn>` for this project's automated author, committer and DCO signature. Add a co-author only when the user explicitly identifies one.
- Before integrating or pushing commits, check author, committer and trailers in every incoming commit. Correct the previously misconfigured automation identity `oooo <vvvvmarisa@163.com>` (GitHub `NotOmee`) to the identity above; preserve actual third-party authorship.
