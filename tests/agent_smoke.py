"""Opt-in deterministic-model check of the complete SDK/MCP/native workflow."""

import argparse
import asyncio
import json
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from app.agent import Tasks
from app.service import Service
from app.jobs import Registry
from app.store import Store


async def main(args):
    root = Path(args.output).resolve()
    service = Service(root)
    profile = {
        "id": "local",
        "kind": "local",
        "roots": [str(Path(args.case).resolve())],
        "program": str(Path(args.program).resolve()),
        "mpi": str(Path(args.mpi).resolve()),
        "max_ranks": 2,
        "max_steps": 1,
    }
    service.save_host(profile)
    calls = []

    def next_call():
        i = len(calls)
        store = Store(root / "hosts/local")
        revision = (store.all("revision") or [{}])[-1]
        reg = Registry(
            profile["roots"], root / "hosts/local", profile["program"], profile["mpi"]
        )
        jobs = reg.database()["jobs"]
        latest = max(jobs, key=lambda k: jobs[k]["created_at"]) if jobs else ""
        if i in (3, 4, 5, 6) and latest:
            until = time.monotonic() + 40
            while time.monotonic() < until:
                r = reg.summary(latest)
                if r["status"] == "completed" and r["alive"] is False:
                    break
                time.sleep(0.2)
        actions = [
            (
                "case_prepare",
                {"source": profile["roots"][0], "request_id": "ai-prepare-001"},
            ),
            ("case_check", {"case_id": revision.get("id", ""), "ranks": 2}),
            (
                "job_start",
                {
                    "case_id": revision.get("id", ""),
                    "steps": 1,
                    "ranks": 2,
                    "request_id": "ai-start-001",
                },
            ),
            (
                "job_start",
                {
                    "run_id": latest,
                    "steps": 1,
                    "ranks": 2,
                    "request_id": "ai-resume-001",
                },
            ),
            (
                "result_render",
                {"run_id": latest, "options": {"variable": "temperature"}},
            ),
            ("result_report", {"run_id": latest}),
            ("case_pack", {"run_id": latest, "request_id": "ai-pack-001"}),
            (
                "finish_task",
                {"summary": "独立输入、检查、新算、续算、绘图、报告与归档完成"},
            ),
        ]
        return actions[i] if i < len(actions) else None

    class Model(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_POST(self):
            json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            action = next_call()
            calls.append(action)
            if action:
                name, params = action
                message = {
                    "role": "assistant",
                    "tool_calls": [
                        {
                            "id": "call-" + str(len(calls)),
                            "type": "function",
                            "function": {
                                "name": "mcp_simulation_" + name,
                                "arguments": json.dumps(params),
                            },
                        }
                    ],
                }
            else:
                message = {"role": "assistant", "content": "流程完成。"}
            chunk = {
                "id": "workflow",
                "object": "chat.completion.chunk",
                "created": 1,
                "model": "deterministic",
                "choices": [{"index": 0, "delta": message, "finish_reason": None}],
            }
            if action:
                message["tool_calls"][0]["index"] = 0
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.end_headers()
            self.wfile.write(("data: " + json.dumps(chunk) + "\n\n").encode())
            chunk["choices"] = [
                {
                    "index": 0,
                    "delta": {},
                    "finish_reason": "tool_calls" if action else "stop",
                }
            ]
            self.wfile.write(
                ("data: " + json.dumps(chunk) + "\n\ndata: [DONE]\n\n").encode()
            )
            self.wfile.flush()

    server = ThreadingHTTPServer(("127.0.0.1", 0), Model)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        service.store.put(
            "model",
            "default",
            {
                "base_url": f"http://127.0.0.1:{server.server_port}/v1",
                "model": "deterministic",
            },
        )
        tasks = Tasks(service)
        task = tasks.create(
            {
                "message": "执行已登记小算例的准备、运行、续算、绘图和归档",
                "max_ranks": 2,
                "max_steps": 1,
            }
        )
        await asyncio.wait_for(tasks.running[task["id"]], 240)
        events = service.store.events(task["id"])
        operations = [e["action"] for e in events if e["type"] == "operation"]
        assert operations == [
            "prepare",
            "case_check",
            "start",
            "resume",
            "plot",
            "report",
            "pack",
        ], events
        assert service.store.get("task", task["id"])["status"] == "completed", events
        (root / "result.json").write_text(
            json.dumps(
                {"task": task["id"], "operations": operations, "events": events},
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        print("PASS SDK -> MCP -> native prepare/check/start/resume/plot/report/pack")
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ("case", "program", "mpi", "output"):
        p.add_argument("--" + name, required=True)
    asyncio.run(main(p.parse_args()))
