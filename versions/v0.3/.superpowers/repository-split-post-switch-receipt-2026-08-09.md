# HUNDUN-FLOW 仓库拆分与目录切换完成记录

状态：`COMPLETE`

## 最终路径

- governance：`/home/wyf/code_dev/hundun-flow-governance`
- product：`/home/wyf/code_dev/hundun-flow`
- Stage 3 branch：`coast/stage3-framework-completion`
- Stage 3 worktree：`/home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework`

## 身份

- governance split seal HEAD：`ee4d2b18d0c68b3080edcc2e132045175961cfb8`
- governance split seal tree：`fc3c3f170d17200938a317f5a4360fd1c0174c3c`
- product HEAD：`ae3d08bbb220d1d3b28ec070d1cba9c33fb85877`
- product tree：`833b633939765365a07fc8d49e34802399954399`
- product version：`0.1.0`
- product commit count：1
- product remote count：0

## 切换结果

旧 mixed repository 通过同文件系统目录改名成为 governance；product candidate 通过目录改名占用原产品路径。两个动作均可通过反向改名恢复，没有删除仓库或历史。

8 个 linked worktree 的 `.git` 绝对指针已从旧 product 路径改到 governance：

- `hundun-flow-flat-layout`
- `hundun-flow-stage3`
- `hundun-flow-task11-exact-candidate`
- `hundun-task11-0884-ab`
- `hundun-task11-directional-widest-ab`
- `hundun-task11-force-decomp`
- `hundun-task11-rhie-face-ab`
- `hundun-task11-weighted-ab`

逐个 `rev-parse --git-common-dir` 均返回 `/home/wyf/code_dev/hundun-flow-governance/.git`。新的 Stage 3 工作树也使用该 common Git directory。

## 保留内容

- 旧 Stage 3 工作树 HEAD：`64e3a70f49da4ccb3585d0eb6afa9966fb82282d`
- 切换后旧 Stage 3 dirty/untracked 条目：310
- governance 根原有未跟踪文件：`docs/plans/2026-07-21-hundun-flow-stage7-neural-warm-start.md`

上述内容未清理、未覆盖、未提交。旧 Stage 3 工作树不再作为开发基线。

## 进程和安全边界

切换前后均未发现 HUNDUN-FLOW build、CTest、MPI 或数值进程。没有访问私有参考源码、研究数据或研究进程；没有 push 或发布。

## 后续权威

Stage 3 后续由以下文件控制：

- `docs/superpowers/specs/2026-08-08-hundun-flow-post-task11-semantic-port-architecture-design.md`
- `docs/superpowers/specs/2026-08-09-hundun-flow-stage3-compact-scientific-design.md`
- `docs/superpowers/plans/2026-08-09-hundun-flow-stage3-framework-completion.md`

下一产品同步只允许在 Task 21 exact-HEAD 接受后进行，并把产品版本提升为 `0.2.0`。
