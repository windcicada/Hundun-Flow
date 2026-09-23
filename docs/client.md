# Hundun Client 0.1.0

## 安装与首次配置

安装 Python 3.11+，运行 `python -m pip install hundun_client-0.1.0-py3-none-any.whl`，随后执行 `hundun-client serve`。浏览器访问 `http://127.0.0.1:8765`。Windows 与 Linux 使用相同入口，网页已随 wheel 打包。源码开发采用 `python tools.py build` 构建网页和安装包。

状态默认位于平台用户状态目录；`hundun-client serve --state <目录> --port 8766` 可选择独立目录。SQLite 保存主机、输入版本、作业请求、AI 任务及事件。目录中的 `hosts`、`ai`、`cache` 保存运行记录、会话和图件。

“连接与模型”设置主机、算例目录、外部求解器入口、MPI 入口、进程及单次步数上限。配置工作目录即可读取已有结果；连接 Linux 求解器后启用原生检查、新算、场输出、保存停止和续算。每次新算及续算使用独立运行目录。

## SSH 计算端

SSH 主机需要 Python 3.11+、可运行的外部求解器和 MPI。高级配置支持 `key_file`、`known_hosts`、`password_ref` 和 `python`。SSH 校验主机公钥，使用用户已有 known_hosts 或指定文件。密码引用采用 `env:变量名` 或 `keyring:条目名`。计算端工作目录应为用户可写的绝对路径。

客户端以 SFTP 部署按内容哈希命名的辅助程序至 `<工作目录>/.client`，通过标准输入传递 JSON 固定操作。MPI 作业在独立会话运行，浏览器、客户端或 SSH 连接关闭后继续计算。重新连接时核对节点、启动标识、PID、参数、输出目录和原生控制回执。待核实状态单独显示。

首次远端绘图在独立 `.client/plot` 环境安装绘图依赖；计算节点需能获取这些 Python 包，也可由管理员预装该环境。大场文件保留在计算端，客户端获取图像、摘要及来源。超过 256 MiB 的归档提供远端路径，使用 SFTP 下载。

## 输入、结果与归档

“输入版本”从模板或已有 `case.json` 目录创建独立副本；完整 JSON 编辑形成新版本。版本记录保存来源、修改理由、配置差异和资产 SHA256。原生检查结果记录在对应版本。主机高级配置的 `tools` 可登记几何、网格工具，例如 `{"mesh":["/opt/tools/mesh","{case}"]}`。

各页面采用物理时间、局部 CFL、单步耗时、方程指标、收支与实际模型字段。截面图支持已发布的 VTK/Visit 场、坐标轴、位置和变量，保留物理单位、时间及来源。原生短文件名、多进程拼接及重复点插值均由公共绘图组件处理。

“报告与归档”导出原始运行记录，或打包输入、资产、运行日志、已接受检查点、SHA256 清单和恢复说明。完整检查点归档在保存停止之后进行。需要携带程序时选择外部发行包目录；包内保留程序本身的依赖和许可。

## 内置 AI

模型设置支持兼容端点和模型名称。凭据使用启动进程的环境变量或系统凭据库，例如 `env:HUNDUN_MODEL_KEY`。也可通过 `POST /api/credentials` 写入系统凭据库；设置保存引用。模型待配置时，任务显示“配置模型后继续”。

固定 SDK 版本为 nanobot-ai 0.3.5，适配集中在 `app/agent_worker.py`。每个任务独立保存会话；MCP 工具调用与网页共用 `app/service.py`、`app/ops.py`。工具覆盖能力检查、输入版本、原生检查、计算、控制、绘图、报告和归档。工具定义与 SDK 来源见 [官方 Python SDK 文档](https://github.com/HKUDS/nanobot/blob/v0.3.5/docs/python-sdk.md)。

任务在指定主机的工作目录和资源上限内执行。缺少物理定义时，助手通过任务提问；补充消息后继续。程序每 5 秒核对关联作业，在完成、失败或连续 300 秒进度保持相同时触发分析；每种状态事件触发一次。任务重启后先读取已持久化的请求与作业状态。单任务自动分析上限为 30 轮，可通过补充目标建立下一任务。

“取消 AI 任务”结束后续助手操作，已派发计算保持自身状态。“保存停止计算”向求解器写入原生停止请求，由求解器提交检查点并返回回执。两项操作具有独立状态和按钮。

## 接口与迁移

所有原 `/api/runs` 入口保留，`host_id` 查询参数选择主机。新增 `/api/hosts`、`/api/model`、`/api/revisions`、`/api/operations`、`/api/tasks`。任务事件使用 `/api/tasks/<id>/events?stream=true`，支持 SSE 事件游标和重新连接。

计算及输入变更携带 `request_id`，同一身份和参数返回已持久化结果；同一身份对应不同参数返回冲突。派发前登记请求，连接中断后的未决请求进入核对状态。

旧 `jobs.json` 的一次性迁移入口：`POST /api/migrate`，JSON 为 `{"host_id":"local","state":"旧状态目录"}`。迁移在原目录建立 SQLite，保留原运行路径、检查点和旧文件，并登记已有输入目录。后续重启读取 SQLite。

## 开发与后续

`app/` 提供 Python 后端、适配器与 AI；`ui/` 提供 React 页面；`tests/` 提供接口、SSH、SDK和绘图检查。`python -m unittest discover -s tests -t .` 执行常规检查。实际求解器的可选测试入口为 `python -m tests.smoke --case <小算例> --program <程序> --mpi <MPI> --output <测试目录>`。

本期采用单用户本机服务。多用户部署、队列调度、桌面外壳和其他 Agent 运行时按后续分支计划推进。
