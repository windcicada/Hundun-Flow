# R4历史工具附件

这里的`.snapshot`是不可执行的历史证据副本，不是当前工具源码。
文件名为内容的SHA-256，字节来自提交
`59985dbf2d10dee167eddfdd38ae5c2752051605`，与R4凭据逐字节一致。
不要格式化、补水印或修正这些快照；改动会破坏历史身份。

| 快照SHA-256 | 原始文件 |
| --- | --- |
| `7336f5dc1c697ea514b9b8a3011a14e3ad2a2465c9957bb10e7b96dfd714ead6` | `tools/v04_literature_extract.py` |
| `c73fb4eb663c450f2b288edbe23d7c4a4233c94f13570ebb790834a986c3dc25` | `tools/v04_parnaudeau_vector_extract.py` |

发布提交`aa2b6219e265ff3978dda29cb4d0fc974f56191a`给工具增加了水印，
因此当前工具与历史工具不再具有相同哈希。校验器的历史模式检查快照字节，
同时单独报告当前校验器哈希；不会执行历史脚本，也不会用当前哈希替换R4记录。

从仓库根目录执行：

```sh
python3 tools/v04_literature_extract.py receipt-validate \
  docs/verification/v0.4-literature-data-receipt-r4-partial.json \
  --historical --relocate-root "/home/wyf/code_dev/hundun-flow=$PWD" \
  --artifact-store docs/verification/literature-r4-artifacts --json
```

PDF及原始数字化资料仍在凭据指定的外部目录，没有复制进仓库。
更换机器时，需要持有这些原始附件，并用额外的`--relocate-root OLD=NEW`
指定其新根目录。缺失或哈希不符会失败，不会按同名文件猜测位置。
本地CTest的R4校验也依赖这些附件；便携的CLI回归另外使用临时合成附件，
二者不能相互替代。

路径迁移只解决“到哪里找原来的字节”，不补充实验数据。
R4及本次新生成的当前凭据均保持`complete=false`，
`--require-complete`仍拒绝未完成的文献依据。
