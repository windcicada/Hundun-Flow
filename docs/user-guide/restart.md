# Restart（默认 versions/v0.4）

```sh
hundun run /path/to/case --restart /path/to/old-run/Restart \
  --output /path/to/new-run --steps 10 --output-interval 0 --restart-interval 10
```

输出必须使用独立目录。保留 Restart 根目录的 `current`、所指 generation 的
`manifest.bin` 和全部 rank 文件，不修改或覆盖来源文件。

默认严格校验输入/计划、字段目录、几何身份、文件大小与内容哈希。
当前格式2保存两层已接受场、两层非对流速率和面通量、压力参考及时间控制历史。
旧格式1缺少完整历史，读取后需要一阶恢复步；这与下面的MG存储兼容是两件事。
格式2支持读取到不同MPI分区，但局部网格仍须满足产品编译约束；不是任意rank数均可用。

## 旧 MG 外层 ghost 布局的显式兼容

`0b0eff1` 消除了 MG 原始线性 bundle 外层重复的 ghost 分配。
历史字段目录及依赖它的边界/产品指纹因而不同，但 checkpoint 不包含 MG 工作区。
对于仅有这项存储布局差异、物理与数值配置保持相同的格式2文件，可显式请求：

```sh
hundun run /path/to/case --restart /path/to/old-run/Restart \
  --restart-storage-compatibility mg-bundle-ghost-v1 \
  --output /path/to/new-run --steps 1 --output-interval 0 --restart-interval 1
```

自定义 `v04_thin_domain_runner` 使用同名兼容选项和原有 `--restart-root` 参数。
它仍验证旧 checkpoint 完成标记、statistics 和 accumulator 的来源身份；
下一次输出使用新产品身份，不改写历史统计文件。

兼容路径根据当前物理配置重建明确的历史目录、边界和产品指纹，要求三者相符，
并保留全部几何、字段、完整性及数值检查。不能用这个选项加载物性、边界、网格、
数值策略或其他历史版本不匹配的文件，也不能通过修改manifest绕过拒绝。
镜像保留来源 plan/schema 与 manifest SHA；运行报告显式记录
`restart_storage_migrated`。物理场和BDF历史恢复后重新建立工作区；
新写出的checkpoint使用当前身份，后续无需兼容选项。

已验证：公共接口1/2/4 ranks，原`69d8eee`二进制生成的小算例旧文件续算，
以及转换后的checkpoint在默认严格路径下2→4 ranks续算。
这不代表所有历史Re3900文件已经迁移或长期稳定。

旧文档中的 `restart.read`、Checkpoint v3/presence 等属于其他版本接口，
不是这里默认v0.4命令行的配置项。
