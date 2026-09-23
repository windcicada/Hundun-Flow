"""Isolated real SSH/SFTP transport for integration checks."""

import asyncio
import os
import shlex
from pathlib import Path
import asyncssh


class Server(asyncssh.SSHServer):
    def begin_auth(self, username):
        return False


async def execute(process):
    child = await asyncio.create_subprocess_exec(
        *shlex.split(process.command),
        stdin=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    output, error = await child.communicate(
        (await process.stdin.read()).encode()
        if process.command.endswith(".pyz")
        else b""
    )
    process.stdout.write(output.decode())
    process.stderr.write(error.decode())
    process.exit(child.returncode)


async def start(root):
    key = asyncssh.generate_private_key("ssh-ed25519")
    server = await asyncssh.create_server(
        Server,
        "127.0.0.1",
        0,
        server_host_keys=[key],
        process_factory=execute,
        sftp_factory=asyncssh.SFTPServer,
    )
    port = server.get_port()
    known = Path(root) / "known_hosts"
    known.write_text(
        "[127.0.0.1]:" + str(port) + " " + key.export_public_key().decode()
    )
    return server, {
        "kind": "ssh",
        "hostname": "127.0.0.1",
        "username": "test",
        "port": port,
        "known_hosts": str(known),
        "work_dir": str(Path(root) / "remote"),
        "python": os.sys.executable,
    }
