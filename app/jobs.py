"""Local process identity, durable idempotency and native control requests."""

import json
import os
import re
import subprocess
import socket
import time
import uuid
from pathlib import Path
from .store import Store

try:
    from platformdirs import user_state_dir
except ImportError:

    def user_state_dir(name):
        return str(Path.home() / ".local" / "state" / name)


try:
    from .data import ident, is_run, has_fields, read_json, rows, walk, history
except ImportError:
    from data import ident, is_run, has_fields, read_json, rows, walk, history


class ApiError(Exception):
    def __init__(self, code, message, status=400):
        self.code, self.message, self.status = code, message, status
        super().__init__(message)


def node_identity():
    boot = Path("/proc/sys/kernel/random/boot_id")
    return {
        "hostname": socket.gethostname(),
        "boot": boot.read_text(encoding="utf-8").strip() if boot.is_file() else None,
    }


def process(pid):
    if os.name != "posix" or not Path("/proc").is_dir():
        return None
    try:
        root = Path("/proc") / str(pid)
        if root.stat().st_uid != os.getuid():
            return None
        fields = (root / "stat").read_text(encoding="utf-8").rsplit(")", 1)[1].split()
        if fields[0] == "Z":
            return None
        args = (root / "cmdline").read_bytes().split(b"\0")
        environ = (root / "environ").read_bytes().split(b"\0")
        size = next(
            (
                v.split(b"=", 1)[1].decode()
                for v in environ
                if v.startswith(b"OMPI_COMM_WORLD_SIZE=") or v.startswith(b"PMI_SIZE=")
            ),
            "",
        )
        return {
            "pid": int(pid),
            "starttime": fields[19],
            "argv": [a.decode(errors="replace") for a in args if a],
            "cwd": (root / "cwd").resolve(),
            "mpi": int(size) if size.isdigit() else None,
        }
    except (OSError, ValueError, IndexError):
        return None


class Registry:
    def __init__(self, roots=None, state=None, program=None, mpi=None):
        self.roots = [
            Path(p).expanduser().resolve()
            for p in (
                roots
                if roots is not None
                else os.environ.get("HUNDUN_CASE_ROOTS", "").split(",")
            )
            if str(p).strip()
        ]
        self.state = (
            Path(
                state
                or os.environ.get(
                    "HUNDUN_CLIENT_STATE", user_state_dir("hundun-client")
                )
            )
            .expanduser()
            .resolve()
        )
        self.state.mkdir(parents=True, exist_ok=True, mode=0o700)
        self.program = Path(program).expanduser().resolve() if program else None
        self.mpi = Path(mpi).expanduser().resolve() if mpi else None
        self.store = Store(self.state)
        self._migrate()
        self.runs, self.cases = {}, {}
        self._processes = []
        self._process_time = 0
        self._discovered_at = 0
        self.refresh()

    def _migrate(self):
        if self.store.get("registry", "main") is not None:
            return
        old = self.state / "jobs.json"
        value = read_json(old) if old.exists() else {"jobs": {}, "requests": {}}
        if (
            not isinstance(value, dict)
            or not isinstance(value.get("jobs"), dict)
            or not isinstance(value.get("requests"), dict)
        ):
            raise ApiError(
                "state_invalid", "旧任务登记格式异常，保留原文件以便核查", 409
            )
        self.store.put("registry", "main", value)

    def database(self):
        return self.store.get("registry", "main", {"jobs": {}, "requests": {}})

    def save(self, db):
        self.store.put("registry", "main", db)

    def refresh(self, force=False):
        if not force and time.monotonic() - self._discovered_at < 10:
            return
        self.runs, self.cases = {}, {}
        for path in walk(self.roots):
            if is_run(path):
                self.runs[ident(path)] = path
            if (path / "case.json").is_file():
                self.cases[ident(path)] = path
        for raw in self.database().get("imports", []):
            path = Path(raw).resolve()
            if any(
                root == path or root in path.parents for root in self.roots
            ) and is_run(path):
                self.runs[ident(path)] = path
        for raw in self.database().get("case_imports", []):
            candidate = Path(raw).resolve()
            if self.state in candidate.parents and (candidate / "case.json").is_file():
                self.cases[ident(candidate)] = candidate
        for key, job in self.database()["jobs"].items():
            path = Path(job["path"])
            if self.state in path.resolve().parents:
                self.runs[key] = path
                case = Path(job["case_path"])
                if any(root == case or root in case.parents for root in self.roots):
                    self.cases[ident(case)] = case
        self._discovered_at = time.monotonic()

    def resolve_run(self, key):
        if key not in self.runs:
            self.refresh(force=True)
        if key not in self.runs:
            raise ApiError("run_not_found", "未找到已登记运行目录", 404)
        path = self.runs[key].resolve()
        if ident(path) != key or not any(
            root == path or root in path.parents for root in self.roots + [self.state]
        ):
            raise ApiError("path_changed", "运行目录已移出允许范围", 409)
        return path

    def native_process(self, info, path):
        if not info:
            return False
        args = info["argv"]
        if self.program is None:
            return False
        trusted = {str(self.program), str(self.program.parent / "bin/hundun")}
        if not any(arg in trusted for arg in args):
            return False
        return self.output_process(info, path)

    def output_process(self, info, path):
        if not info:
            return False
        args = info["argv"]
        try:
            output = Path(args[args.index("--output") + 1])
            if not output.is_absolute():
                output = info["cwd"] / output
            return "run" in args and output.resolve() == path
        except (ValueError, IndexError):
            return False

    def observed(self, key):
        """An output-matching process is evidence of life, never permission."""
        path = self.resolve_run(key)
        if time.monotonic() - self._process_time > 1:
            self._processes = [
                info
                for entry in (Path("/proc").iterdir() if Path("/proc").is_dir() else [])
                if entry.name.isdigit()
                for info in [process(entry.name)]
                if info
            ]
            self._process_time = time.monotonic()
        for info in self._processes:
            if self.output_process(info, path):
                current = process(info["pid"])
                if (
                    current
                    and current["starttime"] == info["starttime"]
                    and self.output_process(current, path)
                ):
                    return current
        return None

    def life(self, key, phase, job, trusted=None):
        if job.get("node") and job["node"] != node_identity():
            return None, "unknown", None
        observed = trusted or self.observed(key)
        if observed:
            return True, "trusted" if trusted else "external", observed
        if job.get("pid") and job.get("starttime"):
            launcher = process(job["pid"])
            if launcher and launcher["starttime"] == job["starttime"]:
                return True, "owned_launcher", launcher
            return False, "exited", None
        if phase in ("completed", "stopped", "failed"):
            return False, "exited", None
        return None, "unknown", None

    def live(self, key):
        path = self.resolve_run(key)
        job = self.database()["jobs"].get(key)
        if job and job.get("node") and job["node"] != node_identity():
            return None
        # Inspect the native process as well as its registered launcher. PID reuse
        # alone never authorizes a command against an unrelated process.
        if job and job.get("pid"):
            launcher = process(job.get("pid", -1))
            if not launcher or launcher["starttime"] != job.get("starttime"):
                return None
        if time.monotonic() - self._process_time > 1:
            self._processes = [
                info
                for entry in (Path("/proc").iterdir() if Path("/proc").is_dir() else [])
                if entry.name.isdigit()
                for info in [process(entry.name)]
                if info
            ]
            self._process_time = time.monotonic()
        for info in self._processes:
            if self.native_process(info, path):
                # Re-read the selected process to reject exit/PID reuse between scans.
                current = process(info["pid"])
                if (
                    current
                    and current["starttime"] == info["starttime"]
                    and self.native_process(current, path)
                ):
                    return current
        return None

    def case_for(self, key):
        path = self.resolve_run(key)
        job = self.database()["jobs"].get(key)
        if job:
            candidate = Path(job["case_path"])
            return candidate if ident(candidate) in self.cases else None
        live = self.live(key)
        if live:
            args = live["argv"]
            try:
                candidate = Path(args[args.index("run") + 1])
                if not candidate.is_absolute():
                    candidate = live["cwd"] / candidate
                candidate = candidate.resolve()
                if ident(candidate) in self.cases:
                    return candidate
            except (ValueError, IndexError):
                pass
        candidates = {
            p.resolve()
            for p in (path, path / "case", path.parent, path.parent / "case")
            if ident(p) in self.cases
        }
        return next(iter(candidates)) if len(candidates) == 1 else None

    def summary(self, key):
        path = self.resolve_run(key)
        state = read_json(path / "status.json", {})
        if not isinstance(state, dict):
            raise ApiError("invalid_status", "运行状态文件格式错误", 409)
        _, meta, values, evidence = history(path, 65536)
        last = values[-1] if values else {}
        ev = evidence[-1] if evidence else {}
        job = self.database()["jobs"].get(key, {})
        live = self.live(key)
        case_path = self.case_for(key)
        case = read_json(case_path / "case.json", {}) if case_path else {}
        if not isinstance(case, dict):
            raise ApiError("invalid_case", "算例配置格式错误", 409)
        mesh = case.get("mesh", {}).get("exact_cells") if case else None
        if isinstance(mesh, dict):
            mesh = [mesh.get(a) for a in ("x", "y", "z")]
        phase = state.get("phase", "unknown")
        alive, process_status, observed = self.life(key, phase, job, live)
        if alive is False:
            # Process inspection can span the final publication of status.json.
            state = read_json(path / "status.json", state)
            phase = state.get("phase", phase)
        if job and phase == "unknown" and alive is not None:
            phase = "starting" if alive else "failed"
        if alive is None:
            phase = "unverified"
        elif alive is False and phase in ("ready", "writing", "solving", "retrying"):
            phase = "inactive"
        mpi = job.get("ranks") or (observed or {}).get("mpi")
        if mpi is None:
            mpi = next(
                (
                    ev[k]
                    for k in ("mpi_size", "mpi_ranks", "ranks")
                    if type(ev.get(k)) is int and ev[k] > 0
                ),
                None,
            )
        seconds = [
            r["seconds"]
            for r in values
            if isinstance(r.get("seconds"), (float, int)) and r["seconds"] > 0
        ]
        has_restart = (path / "Restart/current").exists()
        available = self.program_available()
        return {
            "id": key,
            "name": path.name,
            "path": str(path),
            "case_path": str(case_path) if case_path else None,
            "status": phase,
            "step": state.get("step", last.get("step", ev.get("step"))),
            "case_binding": "registered"
            if job
            else "live_process"
            if live and case_path
            else "directory_candidate"
            if case_path
            else "unknown",
            "time": state.get("time", last.get("time", ev.get("time"))),
            "dt": state.get("dt", last.get("dt")),
            "mpi": mpi,
            "mesh": mesh,
            "scheme": case.get("time", {}).get("scheme") if case else None,
            "coupling": ev.get(
                "coupling", case.get("solver", {}).get("coupling") if case else None
            ),
            "models": [
                str(case[k].get("model", case[k].get("mode", k)))
                for k in ("reaction", "spray", "turbulence")
                if isinstance(case.get(k), dict)
            ],
            "seconds_per_step": sum(seconds) / len(seconds) if seconds else None,
            "timing_window": meta,
            "target_steps": job.get("target_step", job.get("steps")),
            "requested_steps": job.get("steps"),
            "alive": alive,
            "process_status": process_status,
            "updated_at": state.get("updated_unix_ms"),
            "has_restart": has_restart,
            "has_fields": has_fields(path),
            "capabilities": {
                "output": bool(live),
                "pause": bool(live),
                "resume": bool(
                    alive is False and has_restart and case_path and available
                ),
                "start": False,
            },
            "pending_controls": [
                a
                for a in ("stop", "output")
                if (path / a).exists() or (path / (a + ".pending")).exists()
            ],
        }

    def import_run(self, raw):
        path = Path(raw).expanduser().resolve()
        if not any(
            path == root or root in path.parents for root in self.roots
        ) or not is_run(path):
            raise ApiError("invalid_run", "请选择已配置根目录内的原生运行目录")
        with self.store.exclusive():
            db = self.database()
            db["imports"] = sorted(set(db.get("imports", []) + [str(path)]))
            self.save(db)
        self.runs[ident(path)] = path
        return self.summary(ident(path))

    def once(self, request_id, signature, callback):
        if not re.fullmatch(r"[A-Za-z0-9_-]{8,80}", request_id):
            raise ApiError("request_id", "请求编号应为8至80位字母、数字或连字符")
        with self.store.exclusive():
            db = self.database()
            old = db["requests"].get(request_id)
            if old:
                if old["signature"] != signature:
                    raise ApiError("request_conflict", "请求编号已用于其他操作", 409)
                if old.get("result"):
                    return old["result"]
                if old.get("error"):
                    raise ApiError(**old["error"])
                raise ApiError(
                    "request_pending",
                    "该请求已登记，请查看运行状态；重复启动已阻止",
                    409,
                )
            db["requests"][request_id] = {
                "signature": signature,
                "created": time.time(),
            }
            self.save(db)
            try:
                result = callback(db)
            except ApiError as exc:
                db["requests"][request_id]["error"] = {
                    "code": exc.code,
                    "message": exc.message,
                    "status": exc.status,
                }
                self.save(db)
                raise
            db["requests"][request_id]["result"] = result
            self.save(db)
            return result

    def control(self, key, action, request_id):
        if action not in ("pause", "output"):
            raise ApiError("invalid_action", "支持保存并停止或输出场数据")

        def perform(db):
            path = self.resolve_run(key)
            if not self.live(key):
                raise ApiError(
                    "not_running", "未找到与程序及输出目录匹配的活跃进程", 409
                )
            name = "stop" if action == "pause" else "output"
            if (path / (name + ".pending")).is_symlink() or (
                (path / (name + ".pending")).exists()
                and not (path / (name + ".pending")).is_file()
            ):
                raise ApiError("invalid_control", "控制文件路径异常", 409)
            if not (path / (name + ".pending")).exists():
                try:
                    fd = os.open(
                        str(path / name),
                        os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW,
                        0o600,
                    )
                    try:
                        os.fsync(fd)
                    finally:
                        os.close(fd)
                    directory_fd = os.open(str(path), os.O_RDONLY)
                    try:
                        os.fsync(directory_fd)
                    finally:
                        os.close(directory_fd)
                except FileExistsError:
                    if (path / name).is_symlink() or not (path / name).is_file():
                        raise ApiError("invalid_control", "控制文件路径异常", 409)
            return {
                "ok": True,
                "request_id": request_id,
                "pending": True,
                "control": {
                    "action": action,
                    "native_request": name,
                    "completion": "control.jsonl",
                    "published": "Restart/current" if name == "stop" else "Visit",
                },
            }

        return self.once(request_id, [key, action], perform)

    def command(self, argv, timeout=120):
        try:
            result = subprocess.run(
                argv,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                universal_newlines=True,
                cwd=str(self.state),
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise ApiError("preflight_failed", str(exc), 409)
        if result.returncode:
            raise ApiError("preflight_failed", result.stdout[-6000:], 409)
        return result.stdout

    def program_available(self):
        return os.name == "posix" and all(
            p is not None and p.is_file() and os.access(str(p), os.X_OK)
            for p in (self.program, self.mpi)
        )

    def launch(self, key, steps, ranks, request_id, resume=False):
        if (
            type(steps) is not int
            or not 1 <= steps <= 10000000
            or type(ranks) is not int
            or not 1 <= ranks <= min(1024, os.cpu_count() or 1)
        ):
            raise ApiError("invalid_limits", "步数或MPI进程数超出允许范围")

        def perform(db):
            self.refresh(force=True)
            source = self.resolve_run(key) if resume else None
            if resume and self.summary(key)["alive"] is not False:
                raise ApiError(
                    "process_unverified",
                    "源运行仍存活或进程状态待核实，请在确认结束后恢复",
                    409,
                )
            case = self.case_for(key) if resume else self.cases.get(key)
            if not case:
                raise ApiError("case_unknown", "运行目录尚未绑定唯一的已登记算例", 409)
            for existing in db["jobs"].values():
                if existing.get("case_path") != str(case):
                    continue
                prior = process(existing.get("pid", -1))
                if prior and prior["starttime"] == existing.get("starttime"):
                    raise ApiError(
                        "case_running",
                        "该算例已有登记任务运行中，请等待结束或保存并停止",
                        409,
                    )
                if not existing.get("pid"):
                    raise ApiError(
                        "launch_unresolved",
                        "该算例有尚未确认派发结果的任务，请先核查任务登记与日志",
                        409,
                    )
            if not self.program_available():
                raise ApiError(
                    "program_missing", "请配置可执行的HUNDUN_PROGRAM和HUNDUN_MPI", 409
                )
            preflight = {
                "check": self.command(
                    [
                        str(self.mpi),
                        "-n",
                        str(ranks),
                        str(self.program),
                        "check",
                        str(case),
                    ]
                )[-6000:]
            }
            start_step = 0
            if resume:
                if not (source / "Restart/current").exists():
                    raise ApiError(
                        "restart_missing", "源运行尚无完整Restart/current", 409
                    )
                source_evidence, _ = rows(source / "evidence.jsonl", 65536)
                source_case = (
                    source_evidence[-1].get("case") if source_evidence else None
                )
                checked_case = re.search(r"\bVALID case=(\d+)", preflight["check"])
                if (
                    source_case is not None
                    and checked_case
                    and str(source_case) != checked_case.group(1)
                ):
                    raise ApiError(
                        "case_identity_mismatch",
                        "所选算例配置与原运行记录的算例身份不同；恢复已停止",
                        409,
                    )
                preflight["restart"] = self.command(
                    [str(self.program), "restart-info", str(source / "Restart")]
                )[-6000:]
                preflight["restore_compatibility"] = "checked_on_native_load"
                try:
                    info = json.loads(preflight["restart"].strip().splitlines()[-1])
                    start_step = int(info["step"])
                except (ValueError, KeyError, TypeError, IndexError):
                    raise ApiError(
                        "restart_info_invalid", "原生检查点信息缺少有效步号", 409
                    )
            output = self.state / "runs" / ("r" + uuid.uuid4().hex[:10])
            output.mkdir(parents=True, exist_ok=False)
            argv = [
                str(self.mpi),
                "-n",
                str(ranks),
                str(self.program),
                "run",
                str(case),
                "--output",
                str(output),
                "--steps",
                str(steps),
                "--monitor-interval",
                "1",
                "--diagnostics-interval",
                "1",
                "--output-interval",
                str(min(steps, 20)),
                "--restart-interval",
                str(min(steps, 20)),
            ]
            if resume:
                argv += ["--restart", str(source / "Restart")]
            job = {
                "node": node_identity(),
                "path": str(output),
                "case_path": str(case),
                "steps": steps,
                "start_step": start_step,
                "target_step": start_step + steps,
                "ranks": ranks,
                "argv": argv,
                "created_at": time.time(),
                "request_id": request_id,
                "source_run": key if resume else None,
                "preflight": preflight,
            }
            new_key = ident(output)
            db["jobs"][new_key] = job
            self.save(db)  # A crash after dispatch must not cause a duplicate launch.
            with (
                (output / "stdout.log").open("ab") as stdout,
                (output / "stderr.log").open("ab") as stderr,
            ):
                child = subprocess.Popen(
                    argv,
                    stdout=stdout,
                    stderr=stderr,
                    stdin=subprocess.DEVNULL,
                    start_new_session=True,
                    cwd=str(output),
                )
            info = process(child.pid)
            job.update(
                {"pid": child.pid, "starttime": info["starttime"] if info else None}
            )
            self.save(db)
            self.runs[new_key] = output
            return {
                "ok": True,
                "request_id": request_id,
                "run": {"id": new_key, **job},
                "pending": True,
            }

        return self.once(
            request_id, [key, "resume" if resume else "start", steps, ranks], perform
        )
