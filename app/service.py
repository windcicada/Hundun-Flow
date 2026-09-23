"""Host-independent task interface shared by the browser and AI tools."""

import asyncio
import hashlib
import json
import os
import re
import shlex
import uuid
import zipfile
from pathlib import Path
from .jobs import ApiError
from .ops import Operations
from .store import Store

ACTIONS = {
    "health",
    "probe",
    "runs",
    "detail",
    "import_run",
    "cases",
    "start",
    "resume",
    "control",
    "case_check",
    "prepare",
    "revisions",
    "fields",
    "plot",
    "report",
    "pack",
    "external",
}
MUTATIONS = {"start", "resume", "control", "prepare", "pack", "external"}


def secret(ref):
    if not ref:
        return None
    if ref.startswith("env:"):
        return os.environ.get(ref[4:])
    if ref.startswith("keyring:"):
        import keyring

        return keyring.get_password("hundun-client", ref[8:])
    raise ApiError("credential_ref", "凭据使用环境变量或系统凭据库引用")


class Service:
    def __init__(self, state):
        self.store = Store(state)
        self.state = self.store.root
        if not self.store.get("host", "local"):
            self.store.put(
                "host",
                "local",
                {
                    "id": "local",
                    "name": "本机",
                    "kind": "local",
                    "roots": [],
                    "program": None,
                    "mpi": None,
                    "max_ranks": os.cpu_count() or 1,
                    "max_steps": 10000000,
                },
            )

    def hosts(self):
        return self.store.all("host")

    def save_host(self, p):
        kind = p.get("kind", "local")
        if kind not in ("local", "ssh"):
            raise ApiError("host_kind", "主机类型为 local 或 ssh")
        key = p.get("id") or uuid.uuid4().hex[:12]
        if not re.fullmatch("[a-zA-Z0-9_-]{1,40}", key):
            raise ApiError("host_id", "主机编号格式错误")
        roots = p.get("roots", [])
        if not isinstance(roots, list) or any(not isinstance(x, str) for x in roots):
            raise ApiError("roots", "工作目录需为路径列表")
        if kind == "ssh" and (
            not p.get("hostname")
            or not p.get("username")
            or not str(p.get("work_dir", "")).startswith("/")
        ):
            raise ApiError("ssh_config", "SSH 需要主机、用户和绝对工作目录")
        for k in ("max_ranks", "max_steps"):
            if type(p.get(k, 1)) != int or p.get(k, 1) < 1:
                raise ApiError("resource_limit", "资源额度需为正整数")
        keys = {
            "name",
            "kind",
            "roots",
            "program",
            "mpi",
            "max_ranks",
            "max_steps",
            "hostname",
            "username",
            "port",
            "key_file",
            "password_ref",
            "known_hosts",
            "work_dir",
            "python",
            "tools",
            "legacy_state",
        }
        out = {k: v for k, v in p.items() if k in keys}
        out.update(id=key, kind=kind)
        if kind == "local":
            out["roots"] = [str(Path(v).expanduser().resolve()) for v in roots]
        return self.store.put("host", key, out)

    async def execute(self, host, action, params=None, task_id=None):
        if action not in ACTIONS:
            raise ApiError("unknown_action", "未登记的任务操作")
        profile = self.store.get("host", host)
        if not profile:
            raise ApiError("host_missing", "请先登记计算主机", 404)
        params = params or {}
        task = self.store.get("task", task_id) if task_id else None
        if task_id:
            if not task or task["host_id"] != host or task["status"] == "cancelled":
                raise ApiError("task_scope", "任务授权已结束或主机不同", 409)
            if action in MUTATIONS and task.get("status") not in (
                "running",
                "queued",
                "waiting",
            ):
                raise ApiError("task_scope", "该任务当前处于暂停状态", 409)
            if action in ("start", "resume"):
                if (
                    params.get("ranks", 1) > task["max_ranks"]
                    or params.get("steps", 0) > task["max_steps"]
                ):
                    raise ApiError("task_budget", "操作超过本任务资源额度")
        if profile["kind"] == "local":
            result = await asyncio.to_thread(
                Operations(
                    profile, profile.get("legacy_state") or self.state / "hosts" / host
                ).dispatch,
                action,
                params,
            )
        else:
            result = await self.remote(profile, action, params)
        if action in ("plot", "report", "pack"):
            path = Path(result["path"])
            if profile["kind"] == "ssh":
                try:
                    path = await self.fetch(profile, result["path"])
                except ApiError as exc:
                    if exc.code != "asset_size":
                        raise
                    path = None
                    result = {
                        **result,
                        "host_id": host,
                        "storage": "remote",
                        "download": "SFTP",
                    }
            if path is not None:
                token = (
                    hashlib.sha256(str(path).encode()).hexdigest()[:24] + path.suffix
                )
                self.store.put("asset", token, {"path": str(path), "host_id": host})
                result = {**result, "url": "/api/assets/" + token}
        if task_id:
            if action in ("start", "resume"):
                current = self.store.get("task", task_id)
                current["jobs"] = list(
                    dict.fromkeys(current.get("jobs", []) + [result["run"]["id"]])
                )
                self.store.put("task", task_id, current)
            self.store.event(
                task_id,
                {
                    "type": "operation",
                    "action": action,
                    "host_id": host,
                    "result": result,
                },
            )
        return result

    def bundle(self):
        out = self.state / "node.pyz"
        source = Path(__file__).parent
        files = [
            "__init__.py",
            "store.py",
            "data.py",
            "jobs.py",
            "ops.py",
            "node.py",
            "plot.py",
        ]
        digest = hashlib.sha256(
            b"".join((source / f).read_bytes() for f in files)
        ).hexdigest()
        if not out.exists() or self.store.get("meta", "node") != digest:
            tmp = out.with_suffix(".tmp")
            with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
                z.writestr("__main__.py", "from app.node import main\nmain()\n")
                for f in files:
                    z.write(source / f, "app/" + f)
            os.replace(tmp, out)
            self.store.put("meta", "node", digest)
        return out, digest

    async def connection(self, p):
        import asyncssh

        kw = {
            "host": p["hostname"],
            "port": int(p.get("port", 22)),
            "username": p["username"],
            "connect_timeout": 15,
        }
        if p.get("key_file"):
            kw["client_keys"] = [str(Path(p["key_file"]).expanduser())]
        if p.get("password_ref"):
            kw["password"] = secret(p["password_ref"])
        if p.get("known_hosts"):
            kw["known_hosts"] = str(Path(p["known_hosts"]).expanduser())
        # Library defaults verify the user's known_hosts; no trust bypass.
        return await asyncssh.connect(**kw)

    async def remote(self, p, action, params):
        import posixpath

        bundle, digest = self.bundle()
        work = p["work_dir"].rstrip("/") + "/.client"
        async with await self.connection(p) as conn:
            async with conn.start_sftp_client() as sftp:
                await sftp.makedirs(work, exist_ok=True)
                target = posixpath.join(work, "node-" + digest[:12] + ".pyz")
                if not await sftp.exists(target):
                    temp = target + "." + uuid.uuid4().hex + ".tmp"
                    await sftp.put(str(bundle), temp)
                    await sftp.rename(temp, target)
            python = p.get("python", "python3")
            if action in ("fields", "plot"):
                venv = work + "/plot"
                check = await conn.run(
                    shlex.quote(venv + "/bin/python")
                    + " -c "
                    + shlex.quote("import pyvista, matplotlib"),
                    check=False,
                )
                if check.exit_status:
                    cmd = (
                        shlex.quote(python)
                        + " -m venv "
                        + shlex.quote(venv)
                        + " && "
                        + shlex.quote(venv + "/bin/python")
                        + " -m pip install "
                        + shlex.quote("pyvista==0.49.0")
                        + " "
                        + shlex.quote("matplotlib==3.11.2")
                    )
                    setup = await conn.run(cmd, check=False, timeout=300)
                    if setup.exit_status:
                        raise ApiError(
                            "plot_setup",
                            "远端绘图环境准备失败：" + setup.stderr[-2000:],
                            409,
                        )
                python = venv + "/bin/python"
            profile = {
                k: p.get(k)
                for k in ("roots", "program", "mpi", "max_ranks", "max_steps", "tools")
                if p.get(k) is not None
            }
            request = {
                "profile": profile,
                "state": work + "/state",
                "action": action,
                "params": params,
            }
            result = await conn.run(
                shlex.quote(python) + " " + shlex.quote(target),
                input=json.dumps(request),
                check=False,
                timeout=360,
            )
            try:
                value = json.loads(result.stdout)
            except ValueError:
                raise ApiError(
                    "remote_protocol",
                    "远端辅助程序返回异常：" + result.stderr[-1500:],
                    409,
                )
            if "error" in value:
                raise ApiError(**value["error"])
            return value["result"]

    async def fetch(self, p, remote):
        root = p["work_dir"].rstrip("/") + "/.client/state/"
        if not remote.startswith(root) or "/.." in remote:
            raise ApiError("asset_path", "远端产物路径异常")
        folder = self.state / "cache"
        folder.mkdir(exist_ok=True)
        target = folder / (
            hashlib.sha256((p["id"] + remote).encode()).hexdigest()[:24]
            + Path(remote).suffix
        )
        async with await self.connection(p) as conn:
            async with conn.start_sftp_client() as sftp:
                info = await sftp.stat(remote)
                if info.size > 256 * 1024**2:
                    raise ApiError(
                        "asset_size", "大型归档保存在计算端，请通过文件传输下载"
                    )
                temp = target.with_suffix(".tmp")
                await sftp.get(remote, str(temp))
                os.replace(temp, target)
        return target
