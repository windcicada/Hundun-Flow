# HUNDUN-FLOW Stage 3→4 只读接收演练

日期：2026-08-09
schema：`hundun.stage4_p0.intake_dry_run.v1`
status：`template_only`
execution_state：`execute_after_stage3_acceptance`

本文是正式 Stage 4 `4F-0` 的检查模板，不执行 `4F-0`，也不把当前活跃 Stage 3 工作树
当成已接受基线。P0 期间禁止用当前 HEAD、dirty 文件、二进制或进程输出填充下列表格。

## 1. 启动条件和身份块

只有同时满足下列条件，主 agent 才能复制本模板并执行命令：

1. Stage 3 governance receipt 明确写出正式接受结果、tested code HEAD 和 product HEAD；
2. receipt 对应提交已进入治理历史，且产品/测试子树在 tested code HEAD 后没有变化；
3. 用户已明确启动 Stage 4；
4. 主 agent 已取得 accepted governance、code worktree、product repo 和 accepted build/install
   的绝对路径；
5. 未接受 worktree、计划编写 HEAD 和 P0 branch 均未被当作 Stage 4 parent。

正式 `4F-0` receipt 必须填满下表；任一 `execute_after_stage3_acceptance` 尚未替换即阻断：

| identity class | required value |
|---|---|
| accepted Stage 3 governance repo realpath | `execute_after_stage3_acceptance` |
| accepted Stage 3 governance HEAD / tree / parent | `execute_after_stage3_acceptance` |
| accepted Stage 3 tested code worktree realpath | `execute_after_stage3_acceptance` |
| accepted Stage 3 tested code HEAD / tree / parent | `execute_after_stage3_acceptance` |
| accepted Stage 3 product repo realpath | `execute_after_stage3_acceptance` |
| accepted Stage 3 product HEAD / tree / parent | `execute_after_stage3_acceptance` |
| accepted receipt path / SHA-256 / receipt commit | `execute_after_stage3_acceptance` |
| code/product version and banner | `execute_after_stage3_acceptance` |
| code/product branch and remote | `execute_after_stage3_acceptance` |
| governance/code/product dirty and untracked status | `execute_after_stage3_acceptance` |
| code/product DCO result | `execute_after_stage3_acceptance` |
| linked worktree `.git` pointer / Git dir / common dir | `execute_after_stage3_acceptance` |
| Stage 3 deferred-science/capability limitations | `execute_after_stage3_acceptance` |

### 1.1 只读身份命令

执行者先把 receipt 中的绝对路径逐字写入 task-specific variables。不得使用 `~`、`$HOME`、
glob 或从 dirty worktree 猜路径：

```bash
stage4_intake_governance=/absolute/accepted/governance/repo
stage4_intake_code=/absolute/accepted/tested/code/worktree
stage4_intake_product=/absolute/accepted/product/repo
stage4_intake_receipt=/absolute/accepted/stage3/receipt
stage4_intake_build=/absolute/accepted/stage3/build
stage4_intake_install=/absolute/accepted/stage3/install
```

对 governance、code、product 三个根分别执行：

```bash
realpath -e "$stage4_intake_governance"
git -C "$stage4_intake_governance" rev-parse HEAD HEAD^{tree} HEAD^
git -C "$stage4_intake_governance" branch --show-current
git -C "$stage4_intake_governance" remote -v
git -C "$stage4_intake_governance" status --porcelain=v1 --untracked-files=all
git -C "$stage4_intake_governance" log -1 --format=fuller

realpath -e "$stage4_intake_code"
git -C "$stage4_intake_code" rev-parse HEAD HEAD^{tree} HEAD^
git -C "$stage4_intake_code" branch --show-current
git -C "$stage4_intake_code" remote -v
git -C "$stage4_intake_code" status --porcelain=v1 --untracked-files=all
git -C "$stage4_intake_code" log -1 --format=fuller

realpath -e "$stage4_intake_product"
git -C "$stage4_intake_product" rev-parse HEAD HEAD^{tree} HEAD^
git -C "$stage4_intake_product" branch --show-current
git -C "$stage4_intake_product" remote -v
git -C "$stage4_intake_product" status --porcelain=v1 --untracked-files=all
git -C "$stage4_intake_product" log -1 --format=fuller

realpath -e "$stage4_intake_receipt"
sha256sum "$stage4_intake_receipt"
```

DCO 只验证现有 trailer，不添加或伪造 sign-off：

```bash
git -C "$stage4_intake_code" show -s --format='%H%n%P%n%an <%ae>%n%cn <%ce>%n%s%n%(trailers:key=Signed-off-by,valueonly)' HEAD
git -C "$stage4_intake_product" show -s --format='%H%n%P%n%an <%ae>%n%cn <%ce>%n%s%n%(trailers:key=Signed-off-by,valueonly)' HEAD
```

linked worktree 指针必须是普通文本 `.git` 文件，且其目标与 Git 自报结果相同：

```bash
test -f "$stage4_intake_code/.git"
sed -n '1p' "$stage4_intake_code/.git"
git -C "$stage4_intake_code" rev-parse --absolute-git-dir
git -C "$stage4_intake_code" rev-parse --git-common-dir
git -C "$stage4_intake_code" rev-parse --show-toplevel
```

若 code 是主 worktree 而不是 linked worktree，`test -f` 可按 receipt 中已登记的拓扑改为
`test -d`；这种差异必须写入 receipt，不能静默跳过 pointer class。

## 2. public/build/install inventory

### 2.1 公开头和目标图

```bash
rg --files "$stage4_intake_code/include" | LC_ALL=C sort
rg -n 'add_library|add_executable|add_subdirectory|install\(|target_link_libraries|target_include_directories|EXPORT|ALIAS' \
  "$stage4_intake_code/CMakeLists.txt" \
  "$stage4_intake_code/src" \
  "$stage4_intake_code/tests"
cmake --build "$stage4_intake_build" --target help
```

保存并核对：所有 public headers、安装组件、library/executable/test target、export set、
standalone-header test 和 tests-off 路径。Stage 4 planned names 的碰撞扫描为：

```bash
rg -n -i '(^|[^a-z0-9])(chem_|comb_|spray_|rt_|hundun_chemistry|hundun_combustion|hundun_spray)([^a-z0-9]|$)' \
  "$stage4_intake_code/include" \
  "$stage4_intake_code/src" \
  "$stage4_intake_code/tests" \
  "$stage4_intake_code/CMakeLists.txt"
```

每个匹配项标成 `same_authority_reuse`、`name_collision` 或 `unrelated_text`。发现同职责现有
authority 时更新 Stage 4 文件表，禁止新增第二个 registry/driver/flux authority。

### 2.2 已接受安装和 ABI

先只读检查 receipt 指定的 accepted install。不得从当前 P0 artifact 或系统路径替代它：

```bash
realpath -e "$stage4_intake_install"
find "$stage4_intake_install" -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum
find "$stage4_intake_install" -type l -printf '%P -> %l\n' | LC_ALL=C sort
find "$stage4_intake_install" -type f -perm -0100 -print
find "$stage4_intake_install" -type f \( -name '*.so' -o -name '*.so.*' \) -print
```

对上一步确认的每个 HUNDUN executable/shared library 逐个执行以下模板；不对未知系统或研究
二进制运行：

```bash
nm -D --defined-only /absolute/accepted/install/path/to/hundun-library
readelf -h -d --version-info /absolute/accepted/install/path/to/hundun-library
ldd /absolute/accepted/install/path/to/hundun-library
```

正式 `4F-0` 还需验证一次 isolated scratch install。它只写主 agent 为该 task 新建的空目录，
不写 accepted repo/build/install，P0 不执行此命令：

```bash
stage4_intake_scratch_install=/absolute/new/stage4-4F-0-scratch-install
cmake --install "$stage4_intake_build" --prefix "$stage4_intake_scratch_install"
find "$stage4_intake_scratch_install" -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum
```

scratch 路径必须在 `4F-0` receipt 中逐字登记且开始时不存在。模板不包含删除命令；是否保留或
移入治理证据区由正式 task 决定。

## 3. schema、persistence、diagnostics 和 flow authority

### 3.1 精确搜索

以下搜索只针对 accepted tracked files：

```bash
git -C "$stage4_intake_code" grep -n -E 'schema[_ -]?version|schema v[123]|schema_version'
git -C "$stage4_intake_code" grep -n -E 'Checkpoint|checkpoint|Restart|restart|CRC|publish[-_ ]last|presence'
git -C "$stage4_intake_code" grep -n -E 'diagnostic|diagnostics|provider|counter|capability|ledger'
git -C "$stage4_intake_code" grep -n -E 'rollback|retry|collective|attempt|accepted[_ -]?step'
git -C "$stage4_intake_code" grep -n -E 'final[_ -]?(mass[_ -]?)?flux|FaceMassFlux|PISO|corrector'
git -C "$stage4_intake_code" grep -n -E 'authority|registry|descriptor|dispatch|run_.*case'
git -C "$stage4_intake_code" grep -n -E 'field[_ -]?identity|fingerprint|serialization|section[_ -]?id|kind[_ -]?id'
```

每类结果必须给出 `symbol`、`declaring_file`、`defining_file`、`callers`、`tests`、
`serialized_id`、`restart_impact`、`Stage4_action`。零匹配只能标为 `absent_after_exact_search`，
不能解释为已重命名；候选重命名必须用 declaration、definition 和 caller 三方证据确认。

### 3.2 必填 authority 表

| authority class | accepted Stage 3 value | Stage 4 decision |
|---|---|---|
| root executable/case dispatch | `execute_after_stage3_acceptance` | reuse / extend / collision blocker |
| schema v1/v2/v3 loader and IDs | `execute_after_stage3_acceptance` | v4 additive owner |
| Checkpoint v1/v2/v3 dispatcher and IDs | `execute_after_stage3_acceptance` | v4 section owner |
| diagnostics registry and stable IDs | `execute_after_stage3_acceptance` | provider extension point |
| capability ledger root | `execute_after_stage3_acceptance` | Stage 4 rows only |
| retry/rollback transaction owner | `execute_after_stage3_acceptance` | chemistry state participant |
| collective failure owner | `execute_after_stage3_acceptance` | backend failure mapping |
| accepted final mass-flux authority | `execute_after_stage3_acceptance` | reacting scalars consume same flux |
| PISO corrector count authority | `execute_after_stage3_acceptance` | remain exactly two |
| IBM wall/force authority | `execute_after_stage3_acceptance` | reacting boundary extension only |
| WALE effective-viscosity authority | `execute_after_stage3_acceptance` | transport service consumption |
| field identity/fingerprint owner | `execute_after_stage3_acceptance` | composition/mechanism identity extension |

缺少任一行是 intake blocker，不允许 worker凭计划中的旧文件名补值。

## 4. HUNDUN 进程和资源检查

只查看命令行明确包含 accepted HUNDUN repo/worktree/build/install 路径或已登记 HUNDUN unit
前缀的作业。不得列举、推断、检查、进入或停止其他研究进程：

```bash
ps -eo pid=,ppid=,etimes=,stat=,args= | \
  rg -F \
    -e "$stage4_intake_governance" \
    -e "$stage4_intake_code" \
    -e "$stage4_intake_product" \
    -e "$stage4_intake_build" \
    -e "$stage4_intake_install"
systemctl --user list-units --type=service --all --no-pager | rg -i 'hundun|stage[0-9].*hundun'
systemctl --user list-unit-files --type=service --no-pager | rg -i 'hundun|stage[0-9].*hundun'
```

不得把上述 exact-path filter 改成模糊的 `mpirun|python|solver` 全机扫描。结果记录 PID、unit、
exact command、cwd（只在 cwd 已由 HUNDUN receipt/command 识别时）、elapsed 和 owner。
intake 只报告，不发送 signal、不 stop unit、不等待任务完成。

## 5. completeness 与 mutation review

### 5.1 完整性断言

正式 receipt 必须逐项出现以下稳定键：

```text
accepted_governance_head
accepted_governance_tree
accepted_code_head
accepted_code_tree
accepted_product_head
accepted_product_tree
accepted_receipt_sha256
accepted_version
dco_result
worktree_git_pointer
public_header_inventory_sha256
exported_symbol_inventory_sha256
cmake_target_inventory_sha256
schema_v1_v2_v3_inventory
checkpoint_v1_v2_v3_inventory
diagnostics_inventory
capability_inventory
retry_rollback_authority
collective_failure_authority
final_flux_authority
piso_corrector_authority
background_hundun_jobs
deferred_scientific_validation
stage4_product_changes=none
```

### 5.2 mutation 表

主 agent 在 Git 外复制 receipt candidate，每次只删一类，再运行 completeness checker：

| removed class | expected result |
|---|---|
| code/product/governance identity | `REJECT missing accepted identity` |
| linked `.git` pointer | `REJECT missing worktree topology` |
| CMake targets/public headers/exported symbols | `REJECT missing build/API inventory` |
| schema v1--v3 | `REJECT missing schema authority` |
| Checkpoint v1--v3 | `REJECT missing persistence authority` |
| diagnostics/capability IDs | `REJECT missing observable capability authority` |
| retry/rollback/collective failure | `REJECT missing transaction authority` |
| final flux/PISO | `REJECT missing numerical authority` |
| HUNDUN background jobs | `REJECT missing resource audit` |
| deferred capability limitations | `REJECT capability overclaim risk` |

checker 还必须拒绝：任一 `execute_after_stage3_acceptance` 占位、dirty tree 被写成 accepted、
P0 HEAD 被写成 Stage 4 parent、无法验证的 DCO、第二个同职责 authority 或
`stage4_product_changes` 非 `none`。

## 6. P0 自检与停止边界

P0 只运行以下模板检查，不执行上述 accepted-state 命令：

```bash
test "$(rg -o 'execute_after_stage3_acceptance' .superpowers/stage4-p0/intake-dry-run.md | wc -l)" -ge 10
rg -n 'accepted_governance_head|accepted_code_head|accepted_product_head|checkpoint_v1_v2_v3_inventory|final_flux_authority|piso_corrector_authority' \
  .superpowers/stage4-p0/intake-dry-run.md
git diff --check -- .superpowers/stage4-p0/intake-dry-run.md
```

P0 结束时本文仍为 `template_only`。正式 `4F-0` 必须从 accepted Stage 3 重新执行全部身份和
inventory 命令，并生成新的 receipt；不能把本模板提交 SHA 当成接受证据。
