"""Exercise the pinned SDK and real MCP subprocess with a deterministic model endpoint."""

import asyncio
import json
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from app.service import Service
from app.agent import Tasks


class AgentTests(unittest.IsolatedAsyncioTestCase):
    async def test_sdk_mcp_finish_and_recover(self):
        calls = []

        class Model(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_POST(self):
                body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
                tools = [t["function"]["name"] for t in body.get("tools", [])]
                calls.append(tools)
                if len(calls) == 1:
                    message = {
                        "role": "assistant",
                        "content": None,
                        "tool_calls": [
                            {
                                "id": "c1",
                                "type": "function",
                                "function": {
                                    "name": "mcp_simulation_workspace_list",
                                    "arguments": '{"kind":"probe"}',
                                },
                            }
                        ],
                    }
                elif len(calls) == 2:
                    message = {
                        "role": "assistant",
                        "content": None,
                        "tool_calls": [
                            {
                                "id": "c2",
                                "type": "function",
                                "function": {
                                    "name": "mcp_simulation_finish_task",
                                    "arguments": '{"summary":"连接检查完成"}',
                                },
                            }
                        ],
                    }
                else:
                    message = {"role": "assistant", "content": "已检查客户端能力。"}
                out = {
                    "id": "test",
                    "object": "chat.completion",
                    "created": 1,
                    "model": "test",
                    "choices": [
                        {
                            "index": 0,
                            "message": message,
                            "finish_reason": "tool_calls"
                            if message.get("tool_calls")
                            else "stop",
                        }
                    ],
                    "usage": {
                        "prompt_tokens": 1,
                        "completion_tokens": 1,
                        "total_tokens": 2,
                    },
                }
                self.send_response(200)
                if body.get("stream"):
                    self.send_header("Content-Type", "text/event-stream")
                    self.end_headers()
                    delta = dict(message)
                    if delta.get("tool_calls"):
                        delta["tool_calls"] = [
                            {**t, "index": i} for i, t in enumerate(delta["tool_calls"])
                        ]
                    chunk = {
                        **out,
                        "object": "chat.completion.chunk",
                        "choices": [
                            {"index": 0, "delta": delta, "finish_reason": None}
                        ],
                    }
                    self.wfile.write(("data: " + json.dumps(chunk) + "\n\n").encode())
                    chunk["choices"] = [
                        {
                            "index": 0,
                            "delta": {},
                            "finish_reason": out["choices"][0]["finish_reason"],
                        }
                    ]
                    self.wfile.write(
                        ("data: " + json.dumps(chunk) + "\n\ndata: [DONE]\n\n").encode()
                    )
                    self.wfile.flush()
                else:
                    self.send_header("Content-Type", "application/json")
                    self.end_headers()
                    self.wfile.write(json.dumps(out).encode())

        server = ThreadingHTTPServer(("127.0.0.1", 0), Model)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            with tempfile.TemporaryDirectory() as root:
                service = Service(root)
                service.store.put(
                    "model",
                    "default",
                    {
                        "base_url": f"http://127.0.0.1:{server.server_port}/v1",
                        "model": "test",
                    },
                )
                manager = Tasks(service)
                task = manager.create({"message": "检查连接能力"})
                await asyncio.wait_for(manager.running[task["id"]], 60)
                events = service.store.events(task["id"])
                self.assertEqual(
                    service.store.get("task", task["id"])["status"], "completed", events
                )
                self.assertTrue(any(e["type"] == "operation" for e in events), events)
                self.assertTrue(calls)
                self.assertTrue(
                    all(
                        all(t.startswith("mcp_simulation_") for t in names)
                        for names in calls
                    )
                )
                self.assertEqual(
                    Service(root).store.get("task", task["id"])["summary"],
                    "连接检查完成",
                )
        finally:
            server.shutdown()
            server.server_close()
            thread.join()
