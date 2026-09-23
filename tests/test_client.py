import json
import tempfile
import unittest
from pathlib import Path
from fastapi.testclient import TestClient
from app.backend import create_app
from app.service import Service
from app.agent import Tasks
from app.jobs import ApiError, Registry


class ClientTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.service = Service(self.root / "state")
        self.api = TestClient(
            create_app(service=self.service), base_url="http://localhost"
        )

    def tearDown(self):
        self.temp.cleanup()

    def test_solver_absent(self):
        self.assertFalse(self.api.get("/api/health").json()["program_available"])
        self.assertEqual(self.api.get("/api/runs").json(), {"runs": []})
        self.assertEqual(
            self.api.post(
                "/api/model",
                json={
                    "base_url": "http://localhost:1234/v1",
                    "model": "local",
                    "credential_ref": "env:TEST_KEY",
                },
            ).status_code,
            200,
        )

    def test_revision_idempotence_and_restart(self):
        source = self.root / "input"
        source.mkdir()
        (source / "case.json").write_text('{"mesh":{}}')
        self.service.save_host(
            {"id": "local", "kind": "local", "roots": [str(source)], "max_ranks": 2}
        )
        p = {
            "host_id": "local",
            "action": "prepare",
            "params": {
                "source": str(source),
                "config": {"mesh": {"n": 4}},
                "request_id": "revision-123",
            },
        }
        a = self.api.post("/api/operations", json=p)
        self.assertEqual(a.status_code, 200, a.text)
        b = self.api.post("/api/operations", json=p)
        self.assertEqual(a.json(), b.json())
        self.assertEqual(json.loads((source / "case.json").read_text()), {"mesh": {}})
        again = TestClient(
            create_app(service=Service(self.root / "state")),
            base_url="http://localhost",
        )
        self.assertEqual(len(again.get("/api/revisions").json()["revisions"]), 1)
        p["params"]["config"] = {"other": 1}
        self.assertEqual(self.api.post("/api/operations", json=p).status_code, 409)

    def test_legacy_migration_once(self):
        state = self.root / "old"
        state.mkdir()
        value = {"jobs": {}, "requests": {"previous": {"result": {"ok": True}}}}
        (state / "jobs.json").write_text(json.dumps(value))
        self.assertEqual(Registry(state=state).database(), value)
        (state / "jobs.json").write_text("{}")
        self.assertEqual(Registry(state=state).database(), value)

    def test_workspace_escape(self):
        response = self.api.post(
            "/api/operations",
            json={
                "action": "prepare",
                "params": {"source": str(self.root), "request_id": "escape-123"},
            },
        )
        self.assertEqual(response.status_code, 400)

    def test_node_bundle_stdlib_probe(self):
        import subprocess, sys

        bundle, _ = self.service.bundle()
        response = subprocess.run(
            [sys.executable, "-S", str(bundle)],
            input=json.dumps(
                {"profile": {}, "state": str(self.root / "node"), "action": "probe"}
            ),
            text=True,
            capture_output=True,
        )
        value = json.loads(response.stdout)
        self.assertTrue(value["result"]["capabilities"]["read"])
        self.assertFalse(value["result"]["capabilities"]["compute"])


class TaskTests(unittest.IsolatedAsyncioTestCase):
    async def test_missing_model_cancel_and_scope(self):
        with tempfile.TemporaryDirectory() as root:
            service = Service(root)
            manager = Tasks(service)
            t = manager.create({"message": "check inputs"})
            await manager.running[t["id"]]
            self.assertEqual(
                service.store.get("task", t["id"])["status"], "needs_model"
            )
            await manager.cancel(t["id"])
            with self.assertRaises(ApiError):
                await service.execute("local", "runs", task_id=t["id"])
            self.assertEqual(service.store.all("asset"), [])


if __name__ == "__main__":
    unittest.main()
