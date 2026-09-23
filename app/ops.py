"""Shared native operations used by HTTP, MCP and the SSH node helper."""

import difflib
import hashlib
import json
import os
import shutil
import sys
import tarfile
import time
import uuid
from pathlib import Path
from . import data
from .jobs import Registry, ApiError


def sha(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


class Operations:
    def __init__(self, profile, state):
        self.profile = profile
        self.registry = Registry(
            roots=profile.get("roots", []),
            state=state,
            program=profile.get("program"),
            mpi=profile.get("mpi"),
        )
        self.store = self.registry.store

    def allowed(self, raw, exists=True):
        p = Path(raw).expanduser().resolve()
        roots = self.registry.roots + [self.registry.state]
        if not any(p == r or r in p.parents for r in roots) or (
            exists and not p.exists()
        ):
            raise ApiError("outside_workspace", "路径应位于已配置的工作目录")
        return p

    def copy(self, source, target, links=False):
        # Copy only ordinary assets; keep each snapshot self-contained.
        source = self.allowed(source)

        def ignore(folder, names):
            skipped = []
            for name in names:
                item = Path(folder) / name
                if item.is_dir() and (
                    data.is_run(item) or item.resolve() == target.resolve()
                ):
                    skipped.append(name)
                    continue
                if (
                    item.is_symlink()
                    and links
                    and item.is_file()
                    and source in item.resolve().parents
                ):
                    continue
                if item.is_symlink() or not (item.is_file() or item.is_dir()):
                    raise ApiError(
                        "asset_type", "资产包含链接或特殊文件，请准备普通文件目录"
                    )
            return skipped

        shutil.copytree(source, target, ignore=ignore)

    def once(self, p, action, fn):
        def work(db):
            result = fn()
            db.update(self.registry.database())
            return result

        return self.registry.once(
            p["request_id"],
            [action, {k: v for k, v in p.items() if k != "request_id"}],
            work,
        )

    def dispatch(self, action, p=None):
        p = p or {}
        r = self.registry
        if action == "health":
            return {
                "ok": True,
                "program": str(r.program) if r.program else None,
                "program_available": r.program_available(),
                "case_roots": list(map(str, r.roots)),
                "state": str(r.state),
                "api_version": 1,
                "platform": sys.platform,
            }
        if action == "probe":
            result = self.dispatch("health")
            result["cpu_count"] = os.cpu_count()
            result["protocol"] = "native_v1"
            if r.program_available():
                result["version"] = r.command([str(r.program), "--version"]).strip()
            result["capabilities"] = {
                "read": True,
                "plot": True,
                "compute": r.program_available(),
                "prepare": True,
            }
            return result
        if action == "runs":
            r.refresh(force=True)
            items = []
            for key, path in sorted(r.runs.items()):
                try:
                    items.append(data.safe(r.summary(key)))
                except Exception as exc:
                    items.append(
                        {
                            "id": key,
                            "name": path.name,
                            "path": str(path),
                            "status": "unavailable",
                            "alive": None,
                            "models": [],
                            "capabilities": {},
                            "error": {
                                "code": getattr(exc, "code", "run_read_failed"),
                                "message": str(getattr(exc, "message", exc)),
                            },
                        }
                    )
            return {"runs": items}
        if action == "detail":
            path = r.resolve_run(p["run_id"])
            case = r.case_for(p["run_id"])
            return data.detail(
                path,
                r.summary(p["run_id"]),
                data.read_json(case / "case.json") if case else None,
            )
        if action == "import_run":
            return {"run": r.import_run(p["path"])}
        if action == "cases":
            r.refresh(force=True)
            return {
                "cases": [
                    {
                        "id": k,
                        "name": v.name,
                        "path": str(v),
                        "config": data.read_json(v / "case.json"),
                        "capabilities": {"start": r.program_available()},
                    }
                    for k, v in sorted(r.cases.items())
                ]
            }
        if action == "control":
            return r.control(p["run_id"], p["action"], p["request_id"])
        if action in ("start", "resume"):
            ranks = p.get("ranks", 1)
            if type(ranks) != int or ranks > self.profile.get(
                "max_ranks", os.cpu_count() or 1
            ):
                raise ApiError("resource_limit", "MPI 进程数超过该主机授权额度")
            if p.get("steps", 0) > self.profile.get("max_steps", 10000000):
                raise ApiError("resource_limit", "推进步数超过该主机授权额度")
            return r.launch(
                p.get("run_id") if action == "resume" else p["case_id"],
                p["steps"],
                ranks,
                p["request_id"],
                resume=action == "resume",
            )
        if action == "case_check":
            if type(p.get("ranks", 1)) != int or not 1 <= p.get(
                "ranks", 1
            ) <= self.profile.get("max_ranks", os.cpu_count() or 1):
                raise ApiError("resource_limit", "检查进程数超过授权额度")
            r.refresh(force=True)
            case = r.cases.get(p["case_id"])
            if not case:
                raise ApiError("case_unknown", "请选择已登记算例")
            if not r.program_available():
                raise ApiError("program_missing", "请先配置求解器和 MPI")
            result = {
                "ok": True,
                "output": r.command(
                    [
                        str(r.mpi),
                        "-n",
                        str(p.get("ranks", 1)),
                        str(r.program),
                        "check",
                        str(case),
                        "--dry-plan",
                    ]
                ),
            }
            revision = self.store.get("revision", p["case_id"])
            if revision:
                revision["check"] = {**result, "time": time.time()}
                self.store.put("revision", p["case_id"], revision)
            return result
        if action == "prepare":
            return self.once(p, action, lambda: self.prepare(p))
        if action == "revisions":
            return {"revisions": self.store.all("revision")}
        if action == "fields":
            from .plot import discover_fields

            return discover_fields(r.resolve_run(p["run_id"]))
        if action == "plot":
            from .plot import render_slice

            out = r.state / "plots"
            out.mkdir(exist_ok=True)
            return render_slice(r.resolve_run(p["run_id"]), p.get("options", {}), out)
        if action == "report":
            d = self.dispatch("detail", p)
            out = r.state / "exports"
            out.mkdir(exist_ok=True)
            path = out / (p["run_id"] + ".json")
            path.write_text(
                json.dumps(data.safe(d), ensure_ascii=False, indent=2), encoding="utf-8"
            )
            return {"path": str(path), "run_id": p["run_id"]}
        if action == "pack":
            return self.once(p, action, lambda: self.pack(p))
        if action == "external":
            return self.once(p, action, lambda: self.external(p))
        raise ApiError("unknown_action", "该操作未在客户端接口中登记")

    def prepare(self, p):
        r = self.registry
        target = r.state / "cases" / ("c" + uuid.uuid4().hex[:10])
        target.parent.mkdir(exist_ok=True)
        before = {}
        if p.get("source"):
            source = self.allowed(p["source"])
            before = data.read_json(source / "case.json", {})
            self.copy(source, target)
        elif p.get("case_id"):
            r.refresh(force=True)
            source = r.cases.get(p["case_id"])
            if not source:
                raise ApiError("case_unknown", "算例编号不存在")
            before = data.read_json(source / "case.json", {})
            self.copy(source, target)
        else:
            if not r.program_available():
                raise ApiError(
                    "program_missing", "新模板需要连接求解器，已有输入可直接导入"
                )
            r.command([str(r.program), "init-case", "--output", str(target)])
        if "config" in p:
            if not isinstance(p["config"], dict):
                raise ApiError("invalid_config", "配置必须为 JSON 对象")
            (target / "case.json").write_text(
                json.dumps(p["config"], ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
        config = data.read_json(target / "case.json")
        if not isinstance(config, dict):
            raise ApiError("invalid_config", "输入需要有效的 case.json")
        assets = {
            str(f.relative_to(target)): sha(f) for f in target.rglob("*") if f.is_file()
        }
        key = data.ident(target)
        diff = "\n".join(
            difflib.unified_diff(
                json.dumps(before, indent=2, sort_keys=True).splitlines(),
                json.dumps(config, indent=2, sort_keys=True).splitlines(),
                fromfile="source",
                tofile="revision",
            )
        )
        revision = {
            "id": key,
            "path": str(target),
            "source": str(p.get("source") or p.get("case_id") or "native_template"),
            "created": time.time(),
            "reason": p.get("reason", "用户输入"),
            "assets": assets,
            "diff": diff,
            "config": config,
        }
        self.store.put("revision", key, revision)
        db = r.database()
        db.setdefault("case_imports", []).append(str(target))
        r.save(db)
        r.refresh(force=True)
        return {"case": {"id": key, "path": str(target)}, "revision": revision}

    def pack(self, p):
        r = self.registry
        source = r.resolve_run(p["run_id"])
        case = r.case_for(p["run_id"])
        if not case:
            raise ApiError("case_unknown", "归档需要唯一算例身份")
        # Immutable accepted checkpoints only; native checksum inspection is retained.
        if p.get("checkpoint", True) and (source / "Restart/current").exists():
            if r.summary(p["run_id"])["alive"] is not False:
                raise ApiError("run_active", "请在保存停止后归档完整检查点")
            if r.program_available():
                r.command([str(r.program), "restart-info", str(source / "Restart")])
        out = r.state / "exports"
        out.mkdir(exist_ok=True)
        stage = out / ("b" + uuid.uuid4().hex[:10])
        stage.mkdir()
        self.copy(case, stage / "case")
        (stage / "run").mkdir()
        for name in (
            "status.json",
            "monitor.jsonl",
            "diagnostics.jsonl",
            "evidence.jsonl",
            "control.jsonl",
        ):
            f = source / name
            if f.is_file() and not f.is_symlink():
                shutil.copy2(f, stage / "run" / name)
        if p.get("checkpoint", True) and (source / "Restart").is_dir():
            # Native current may be a relative symlink; resolve while preserving file data.
            restart = (source / "Restart").resolve()
            for f in restart.rglob("*"):
                if f.is_symlink() and restart not in f.resolve().parents:
                    raise ApiError("asset_type", "检查点链接超出目录")
            shutil.copytree(restart, stage / "run/Restart", symlinks=False)
        for name in p.get("results", []):
            result = (source / name).resolve()
            if (
                source not in result.parents
                or not result.exists()
                or Path(name).is_absolute()
            ):
                raise ApiError("result_path", "结果应位于所选运行目录")
            dest = stage / "run" / result.relative_to(source)
            if dest.exists():
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            if result.is_dir():
                self.copy(result, dest)
            else:
                shutil.copy2(result, dest)
        if p.get("runtime"):
            self.copy(self.allowed(p["runtime"]), stage / "runtime", links=True)
        files = {
            str(f.relative_to(stage)): {"bytes": f.stat().st_size, "sha256": sha(f)}
            for f in stage.rglob("*")
            if f.is_file()
        }
        (stage / "manifest.json").write_text(
            json.dumps(
                {
                    "format": "hundun-client-case-v1",
                    "program": str(r.program),
                    "files": files,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        (stage / "RUN.md").write_text(
            "恢复：在 Linux 计算端配置原生程序及 MPI，执行\n\nmpirun -n <ranks> hundun run case --restart run/Restart --output resumed --steps <steps>\n",
            encoding="utf-8",
        )
        archive = out / (stage.name + ".tar.gz")
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(stage, arcname=stage.name)
        shutil.rmtree(stage)
        return {"path": str(archive), "sha256": sha(archive)}

    def external(self, p):
        tool = self.profile.get("tools", {}).get(p.get("tool"))
        if not isinstance(tool, list) or not tool:
            raise ApiError("tool_missing", "请在主机配置中登记几何或网格工具命令")
        case = self.allowed(p["path"])
        argv = [str(x).replace("{case}", str(case)) for x in tool]
        return {"output": self.registry.command(argv), "path": str(case)}
