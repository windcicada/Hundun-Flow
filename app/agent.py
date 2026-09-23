"""Persistent AI tasks; cancelling assistance is separate from stopping CFD."""

import asyncio
import json
import os
import sys
import time
import uuid
from pathlib import Path
from .jobs import ApiError
from .service import secret


class Tasks:
    def __init__(self, service):
        self.service = service
        self.store = service.store
        self.workers = {}
        self.running = {}

    def create(self, p):
        request = p.get("request_id") or uuid.uuid4().hex
        with self.store.exclusive():
            previous = self.store.get("task_request", request)
            if previous:
                if previous["input"] != p:
                    raise ApiError("request_conflict", "请求编号已用于其他目标", 409)
                return self.store.get("task", previous["task"])
            result = self._create(p)
            self.store.put("task_request", request, {"input": p, "task": result["id"]})
            return result

    def _create(self, p):
        host = self.store.get("host", p.get("host_id", "local"))
        if not host:
            raise ApiError("host_missing", "请选择计算主机")
        text = str(p.get("message", "")).strip()
        if not text:
            raise ApiError("message", "请输入模拟目标")
        key = uuid.uuid4().hex[:16]
        value = {
            "id": key,
            "host_id": host["id"],
            "message": text,
            "created": time.time(),
            "updated": time.time(),
            "status": "queued",
            "jobs": [],
            "max_ranks": min(
                int(p.get("max_ranks", host.get("max_ranks", 1))),
                host.get("max_ranks", 1),
            ),
            "max_steps": min(
                int(p.get("max_steps", 100)), host.get("max_steps", 10000000)
            ),
            "turns": 0,
            "pending": [],
        }
        if value["max_ranks"] < 1 or value["max_steps"] < 1:
            raise ApiError("task_budget", "任务资源额度需为正整数")
        self.store.put("task", key, value)
        self.store.event(key, {"type": "created", "message": text})
        self.schedule(key, text)
        return value

    def schedule(self, key, message):
        if key in self.running and not self.running[key].done():
            return
        self.running[key] = asyncio.create_task(self.run(key, message))

    def update(self, key, **fields):
        task = self.store.get("task", key)
        if not task:
            raise ApiError("task_missing", "任务不存在", 404)
        if task["status"] == "cancelled" and fields.get("status") != "queued":
            return task
        task.update(fields, updated=time.time())
        return self.store.put("task", key, task)

    async def cancel(self, key):
        task = self.update(key, status="cancelled")
        self.store.event(
            key,
            {"type": "cancelled", "message": "AI 后续操作已停止，计算作业保持自身状态"},
        )
        process = self.workers.get(key)
        if process and process.returncode is None:
            process.terminate()
        return task

    def followup(self, key, message):
        task = self.store.get("task", key)
        if not task:
            raise ApiError("task_missing", "任务不存在", 404)
        if key in self.running and not self.running[key].done():
            task.setdefault("pending", []).append(message)
            self.store.put("task", key, task)
        else:
            self.update(key, status="queued")
            self.schedule(key, message)
        self.store.event(key, {"type": "user", "message": message})
        return self.store.get("task", key)

    async def run(self, key, message):
        task = self.store.get("task", key)
        if task["status"] == "cancelled":
            return
        model = self.store.get("model", "default")
        if not model or not model.get("model") or not model.get("base_url"):
            self.update(key, status="needs_model")
            self.store.event(
                key,
                {
                    "type": "needs_model",
                    "message": "请在连接与模型页面配置模型端点后恢复任务",
                },
            )
            return
        if task.get("turns", 0) >= 30:
            self.update(key, status="needs_input")
            self.store.event(
                key,
                {
                    "type": "needs_input",
                    "message": "本任务已完成30轮分析，请补充下一阶段目标",
                },
            )
            return
        self.update(key, status="running", turns=task.get("turns", 0) + 1)
        env = os.environ.copy()
        env["PYTHONUTF8"] = "1"
        env["PYTHONPATH"] = str(Path(__file__).resolve().parent.parent)
        try:
            credential = secret(model.get("credential_ref"))
            if model.get("credential_ref") and not credential:
                raise ApiError("credential_missing", "模型凭据引用尚未解析")
            env["HUNDUN_MODEL_KEY"] = credential or "local"
            worker = await asyncio.create_subprocess_exec(
                sys.executable,
                "-m",
                "app.agent_worker",
                str(self.service.state),
                key,
                stdin=asyncio.subprocess.PIPE,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.DEVNULL,
                env=env,
            )
            self.workers[key] = worker
            worker.stdin.write(json.dumps({"message": message}).encode())
            await worker.stdin.drain()
            worker.stdin.close()
            failed = False
            async for line in worker.stdout:
                try:
                    event = json.loads(line)
                except ValueError:
                    continue
                if credential:
                    event = json.loads(
                        json.dumps(event).replace(
                            json.dumps(credential)[1:-1], "[credential]"
                        )
                    )
                self.store.event(key, event)
                if event.get("type") == "failed":
                    failed = True
            code = await worker.wait()
            if self.store.get("task", key)["status"] in (
                "cancelled",
                "completed",
                "needs_input",
            ):
                return
            current = self.store.get("task", key)
            self.update(
                key,
                status="failed"
                if code or failed
                else ("waiting" if current.get("jobs") else "needs_input"),
            )
        except Exception as exc:
            self.update(key, status="failed")
            self.store.event(key, {"type": "failed", "message": str(exc)[:1000]})
        finally:
            self.workers.pop(key, None)
            current = self.store.get("task", key)
            if current and current.get("pending") and current["status"] != "cancelled":
                message = current["pending"].pop(0)
                self.store.put("task", key, current)
                asyncio.get_running_loop().call_later(0.1, self.schedule, key, message)

    async def monitor(self):
        for t in self.store.all("task"):
            if t["status"] in ("running", "queued"):
                self.update(t["id"], status="waiting")
        while True:
            for t in self.store.all("task"):
                if t["status"] != "waiting":
                    continue
                try:
                    reports = []
                    for job in t.get("jobs", []):
                        d = await self.service.execute(
                            t["host_id"], "detail", {"run_id": job}
                        )
                        r = d["run"]
                        reports.append(
                            {
                                "id": job,
                                "step": r["step"],
                                "status": r["status"],
                                "alive": r["alive"],
                            }
                        )
                    signature = json.dumps(reports, sort_keys=True)
                    now = time.time()
                    if signature != t.get("progress"):
                        self.update(t["id"], progress=signature, progress_time=now)
                    stalled = (
                        bool(reports)
                        and signature == t.get("progress")
                        and now - t.get("progress_time", now) >= 300
                    )
                    terminal = (
                        not reports
                        or all(r["alive"] is False for r in reports)
                        or any(r["status"] == "failed" for r in reports)
                        or stalled
                    )
                    if terminal and signature != t.get("notified"):
                        self.update(t["id"], notified=signature, status="queued")
                        self.schedule(
                            t["id"],
                            "任务恢复核对结果："
                            + signature
                            + "。先核对既有操作与原目标，分析结果并完成后处理。完成全部目标后调用 finish_task。",
                        )
                except Exception as exc:
                    if t.get("monitor_error") != str(exc):
                        self.update(t["id"], monitor_error=str(exc))
                        self.store.event(
                            t["id"], {"type": "connection", "message": str(exc)[:500]}
                        )
            await asyncio.sleep(5)
