## Agent skills

### HUNDUN-FLOW workflow

- As requested on 2026-09-05, suspend use of `ponytail` for HUNDUN-FLOW until the user explicitly re-enables it.
- Before exploring unfamiliar code, check the active checkout with `codegraphf status --json .`; use `codegraphf sync .` when the index is stale. Start with `context`, `query`, `callers`, `callees`, or `impact`, then use `rg` and source reads to verify results. State when CodeGraphF is unavailable and use the fallback.
- Index the actual active worktree, not another checkout of the same repository. Use the CLI if the current session does not expose the configured MCP tools.

### Issue tracker

Issues and specs live as GitHub issues, operated with the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Default triage vocabulary is in use: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context layout: one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.
