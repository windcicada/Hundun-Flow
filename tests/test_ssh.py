import tempfile
import unittest
from pathlib import Path
from app.service import Service
from .ssh import start


class SSHTests(unittest.IsolatedAsyncioTestCase):
    async def test_deploy_disconnect_reconnect_and_identity(self):
        with tempfile.TemporaryDirectory() as root:
            server, p = await start(root)
            try:
                service = Service(Path(root) / "state")
                p.update(
                    id="node", name="SSH test", roots=[root], max_ranks=2, max_steps=3
                )
                service.save_host(p)
                first = await service.execute("node", "probe")
                self.assertFalse(first["program_available"])
                source = Path(root) / "input"
                source.mkdir()
                (source / "case.json").write_text('{"mesh":{}}')
                params = {"source": str(source), "request_id": "ssh-copy-123"}
                one = await service.execute("node", "prepare", params)
                service = Service(Path(root) / "state")
                two = await service.execute("node", "prepare", params)
                self.assertEqual(one, two)
                self.assertEqual(
                    len((await service.execute("node", "revisions"))["revisions"]), 1
                )
                p["known_hosts"] = str(Path(root) / "empty")
                Path(p["known_hosts"]).write_text("")
                service.save_host(p)
                with self.assertRaises(Exception):
                    await service.execute("node", "probe")
            finally:
                server.close()
                await server.wait_closed()
