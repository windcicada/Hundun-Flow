# Hundun Client

独立的仿真工作台与内置 AI。Windows 和 Linux 客户端通过本机命令或 SSH，连接 Linux 上的 Hundun-Flow 求解器。

## 安装与启动

安装发行 wheel 后运行：

```sh
python -m pip install hundun_client-0.1.0-py3-none-any.whl
hundun-client serve
```

浏览器打开 http://127.0.0.1:8765 。首次在“连接与模型”中设置工作目录、计算程序和模型服务。Python 3.11+ 承载客户端；前端资源随安装包提供。

客户端可独立查看原生记录、绘制截面、管理输入版本、导出归档。配置求解器后启用新算和恢复；配置模型后启用 AI 任务。Windows 客户端通过 SSH 连接 Linux 计算节点。

## 开发

```sh
python -m venv .venv
python -m pip install -e . build
npm ci --prefix ui
python tools.py build
hundun-client serve --port 8766
```

使用方式与接口见 [客户端说明](docs/client.md)，验证范围见 [验证记录](docs/check.md)。

## 许可

Apache License 2.0。第三方依赖保留各自许可与归属。
