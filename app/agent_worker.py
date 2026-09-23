"""Pinned framework adapter, isolated from the HTTP lifecycle and credentials."""

import asyncio
import json
import sys
import signal
from pathlib import Path
from .store import Store


def emit(value):
    print(json.dumps(value, ensure_ascii=False, default=str), flush=True)


async def run(state, key, message):
    if sys.platform != "win32":
        asyncio.get_running_loop().add_signal_handler(
            signal.SIGTERM, asyncio.current_task().cancel
        )
    from nanobot import Nanobot
    from loguru import logger

    logger.remove()
    store = Store(state)
    task = store.get("task", key)
    model = store.get("model", "default")
    work = Path(state) / "ai" / key / "work"
    work.mkdir(parents=True, exist_ok=True)
    config = {
        "agents": {
            "defaults": {
                "workspace": str(work),
                "model": model["model"],
                "provider": "custom",
                "maxToolIterations": 30,
            }
        },
        "providers": {
            "custom": {"apiBase": model["base_url"], "apiKey": "${HUNDUN_MODEL_KEY}"}
        },
        "tools": {
            "exec": {"enable": False},
            "file": {"enable": False},
            "web": {"enable": False},
            "restrictToWorkspace": True,
            "mcpServers": {
                "simulation": {
                    "command": sys.executable,
                    "args": ["-m", "app.mcp", str(state), key],
                    "env": {"PYTHONPATH": str(Path(__file__).resolve().parent.parent)},
                    "toolTimeout": 360,
                }
            },
        },
    }
    path = work.parent / "model.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    system = (
        "你是 Hundun Client 模拟助手。仅通过 simulation 工具操作模拟；依据真实结果回答。"
        "对已有输入创建新版本，通过原生检查后运行。任务授权仅适用于所选主机与资源额度。"
        "缺少物理定义时提出明确问题。普通监测由客户端执行，作业运行中可结束当前回复等待事件。"
        "完成用户全部目标后调用 finish_task；需要用户信息时调用 ask_user。"
        f"当前授权：主机={task['host_id']}，MPI上限={task['max_ranks']}，单次推进步数上限={task['max_steps']}。"
        f"原始目标：{task['message']}\n本次消息：{message}"
    )
    async with Nanobot.from_config(config_path=path, workspace=work) as bot:
        # This version-specific seam keeps generic framework tools outside the product tool surface.
        await bot._mcp_provider.connect()
        for name in list(bot._loop.tools.tool_names):
            if not name.startswith("mcp_simulation_"):
                bot._loop.tools.unregister(name)
        if not bot._loop.tools.tool_names:
            raise RuntimeError("模拟工具连接尚未就绪")
        async for e in bot.stream(system, session_key="task:" + key):
            kind = str(e.type)
            if kind.endswith("text.delta"):
                emit({"type": "text", "text": e.delta})
            elif kind.endswith("tool.started"):
                emit({"type": "tool", "name": e.name, "arguments": e.arguments})
            elif kind.endswith("tool.completed"):
                emit({"type": "tool_done", "name": e.name})
            elif kind.endswith("run.completed"):
                if e.result.stop_reason in ("error", "tool_error"):
                    emit({"type": "failed", "message": e.result.content})
                else:
                    emit(
                        {
                            "type": "answer",
                            "text": e.result.content,
                            "usage": e.result.usage,
                            "stop_reason": e.result.stop_reason,
                        }
                    )
            elif kind.endswith("run.failed"):
                emit({"type": "failed", "message": e.error})


if __name__ == "__main__":
    try:
        asyncio.run(run(sys.argv[1], sys.argv[2], json.load(sys.stdin)["message"]))
    except Exception as exc:
        emit({"type": "failed", "message": str(exc)[:1200]})
        sys.exit(1)
