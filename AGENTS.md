## Agent skills

### HUNDUN-FLOW workflow

- As requested on 2026-09-05, suspend use of `ponytail` for HUNDUN-FLOW until the user explicitly re-enables it.
- Before exploring unfamiliar code, check the active checkout with `codegraphf status --json .`; use `codegraphf sync .` when the index is stale. Start with `context`, `query`, `callers`, `callees`, or `impact`, then use `rg` and source reads to verify results. State when CodeGraphF is unavailable and use the fallback.
- Index the actual active worktree, not another checkout of the same repository. Use the CLI if the current session does not expose the configured MCP tools.
- Use short file and directory names and a shallow layout for cases, builds and debugging. Reuse existing flat working directories such as `check` and `cases/g`.
- Remove obsolete generated source and debug artifacts when retiring an experiment, without making backups. Preserve files belonging to active work.

### COAST 数值算法基线

- Hundun-Flow 优先采用 COAST 的成熟方法，以 Hundun 风格的 C++ 实现，继承已验证有效的修正和优化；COAST 源码作为只读参考。
- 每项迁移先固定 COAST 的版本、编译配置及算例实际运行后端，再对齐方程系数、IBM、物性、边界条件、时间推进、耦合顺序、线性求解和收敛判据。逐模块记录对应关系及验收状态。
- 替代算法成为默认配置须提供同网格、同场、同时间步、同硬件的完整方法对照，统一原方程残差与物理精度标准，并证明确切收益；收益相当时采用 COAST 方法。
- 通用替代工作按 `docs/alg.md` 的阶段和验收门槛推进。单个算例或单个算法组件的测试，按实际覆盖范围报告。
- 用户于 2026-09-13 将湍流模型范围确定为现有 Vreman 的通用验收。动态 Lilly／Piomelli、额外 SGS 模型及 COAST 模型目录复制已退出本期开发范围；保留已有有效实现，优先公共方程、选定物理组合和算例能力。
- 燃烧与蒸发接入先读 `docs/rc.md`：用户指定 rsfz-143 的 GTMC 和 admin-253 的 624CF 为续算调试参考。保持参考长测和检查点，遇到影响物理定义或验收标准的冲突时暂停依赖工作并请用户决策。

### Issue tracker

Issues and specs live as GitHub issues, operated with the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Default triage vocabulary is in use: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context layout: one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.

### Commit identity

- Use `WANG YUDONG <wangyudong@buaa.edu.cn>` for this project's automated author, committer and DCO signature. Add a co-author only when the user explicitly identifies one.
- Before integrating or pushing commits, check author, committer and trailers in every incoming commit. Correct the previously misconfigured automation identity `oooo <vvvvmarisa@163.com>` (GitHub `NotOmee`) to the identity above; preserve actual third-party authorship.
