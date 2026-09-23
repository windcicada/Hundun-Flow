"""Loopback HTTP interface for the shared client service."""

import asyncio
import contextlib
import json
from pathlib import Path
from urllib.parse import urlsplit
from contextlib import asynccontextmanager
from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.exceptions import RequestValidationError
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, StrictInt, StrictStr
from platformdirs import user_state_dir
from .jobs import ApiError
from .ops import Operations
from .service import Service
from .agent import Tasks


class Control(BaseModel):
    action: StrictStr
    request_id: StrictStr


class Launch(BaseModel):
    steps: StrictInt
    ranks: StrictInt
    request_id: StrictStr


class ImportRun(BaseModel):
    path: StrictStr


def create_app(registry=None, service=None, state=None):
    service = service or Service(state or user_state_dir("hundun-client"))
    tasks = Tasks(service)

    @asynccontextmanager
    async def lifespan(api):
        monitor = asyncio.create_task(tasks.monitor()) if registry is None else None
        yield
        if monitor:
            monitor.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await monitor
        for worker in list(tasks.workers.values()):
            if worker.returncode is None:
                worker.terminate()
                try:
                    await asyncio.wait_for(worker.wait(), 5)
                except asyncio.TimeoutError:
                    worker.kill()
                    await worker.wait()

    api = FastAPI(
        title="Hundun Client", docs_url=None, redoc_url=None, lifespan=lifespan
    )
    api.state.registry = registry
    api.state.service = service
    api.state.tasks = tasks

    async def call(action, p=None, host="local"):
        if registry is not None:
            ops = Operations.__new__(Operations)
            ops.registry = registry
            ops.store = registry.store
            ops.profile = {}
            return await asyncio.to_thread(ops.dispatch, action, p)
        return await service.execute(host, action, p)

    @api.exception_handler(ApiError)
    async def api_error(request, exc):
        return JSONResponse(
            {"error": {"code": exc.code, "message": exc.message}},
            status_code=exc.status,
        )

    from .plot import PlotError

    @api.exception_handler(PlotError)
    async def plot_error(request, exc):
        return JSONResponse(
            {"error": {"code": exc.code, "message": exc.message}}, status_code=400
        )

    @api.exception_handler(RequestValidationError)
    async def invalid(request, exc):
        return JSONResponse(
            {
                "error": {
                    "code": "invalid_input",
                    "message": "请求字段或数值类型不符合接口要求",
                }
            },
            status_code=422,
        )

    @api.exception_handler(KeyError)
    @api.exception_handler(ValueError)
    @api.exception_handler(TypeError)
    async def malformed(request, exc):
        return JSONResponse(
            {
                "error": {
                    "code": "invalid_input",
                    "message": "请检查操作所需字段、数值和配置格式",
                }
            },
            status_code=422,
        )

    @api.exception_handler(OSError)
    async def io_error(request, exc):
        return JSONResponse(
            {"error": {"code": "io_error", "message": str(exc)}}, status_code=409
        )

    import asyncssh

    @api.exception_handler(asyncssh.Error)
    async def ssh_error(request, exc):
        return JSONResponse(
            {"error": {"code": "ssh_connection", "message": str(exc)}}, status_code=409
        )

    @api.middleware("http")
    async def local_only(request, call_next):
        host = request.headers.get("host", "")
        try:
            hostname = urlsplit("http://" + host).hostname
        except ValueError:
            hostname = None
        if hostname not in ("localhost", "127.0.0.1", "::1"):
            return JSONResponse(
                {"error": {"code": "invalid_host", "message": "仅允许本机访问"}},
                status_code=403,
            )
        origin = request.headers.get("origin")
        if request.method not in ("GET", "HEAD", "OPTIONS") and (
            (origin and origin.rstrip("/") not in ("http://" + host, "https://" + host))
            or request.headers.get("sec-fetch-site") == "cross-site"
        ):
            return JSONResponse(
                {
                    "error": {
                        "code": "invalid_origin",
                        "message": "写操作须由同源界面发起",
                    }
                },
                status_code=403,
            )
        return await call_next(request)

    @api.get("/api/health")
    async def health(host_id: str = "local"):
        return await call("health", host=host_id)

    @api.get("/api/runs")
    async def runs(host_id: str = "local"):
        return await call("runs", host=host_id)

    @api.get("/api/cases")
    async def cases(host_id: str = "local"):
        return await call("cases", host=host_id)

    @api.get("/api/revisions")
    async def revisions(host_id: str = "local"):
        return await call("revisions", host=host_id)

    @api.post("/api/runs/import")
    async def import_run(body: ImportRun, host_id: str = "local"):
        return await call("import_run", body.model_dump(), host_id)

    @api.get("/api/runs/{key}")
    async def detail(key: str, host_id: str = "local"):
        return await call("detail", {"run_id": key}, host_id)

    @api.get("/api/runs/{key}/fields")
    async def fields(key: str, host_id: str = "local"):
        return await call("fields", {"run_id": key}, host_id)

    @api.post("/api/runs/{key}/control")
    async def control(key: str, body: Control, host_id: str = "local"):
        return await call("control", {"run_id": key, **body.model_dump()}, host_id)

    @api.post("/api/runs/{key}/resume")
    async def resume(key: str, body: Launch, host_id: str = "local"):
        return await call("resume", {"run_id": key, **body.model_dump()}, host_id)

    @api.post("/api/cases/{key}/start")
    async def start(key: str, body: Launch, host_id: str = "local"):
        return await call("start", {"case_id": key, **body.model_dump()}, host_id)

    @api.post("/api/runs/{key}/plot")
    async def plot(key: str, body: dict, host_id: str = "local"):
        return await call("plot", {"run_id": key, "options": body}, host_id)

    @api.get("/api/hosts")
    async def hosts():
        return {"hosts": service.hosts()}

    @api.post("/api/hosts")
    async def save_host(body: dict):
        return service.save_host(body)

    @api.post("/api/hosts/{key}/probe")
    async def probe(key: str):
        return await service.execute(key, "probe")

    @api.post("/api/operations")
    async def operation(body: dict):
        return await call(
            body.get("action", ""), body.get("params", {}), body.get("host_id", "local")
        )

    @api.post("/api/migrate")
    async def migrate(body: dict):
        profile = service.store.get("host", body.get("host_id", "local"))
        if not profile or profile["kind"] != "local":
            raise ApiError("migration_host", "请选择旧登记所在的本机连接")
        root = Path(body["state"]).expanduser().resolve()
        from .jobs import Registry

        old = Registry(state=root)
        roots = set(profile.get("roots", []))
        for job in old.database()["jobs"].values():
            roots.add(str(Path(job["case_path"]).resolve()))
        profile.update(legacy_state=str(root), roots=sorted(roots))
        service.save_host(profile)
        return {
            "imported": len(old.database()["jobs"]),
            "state": str(root),
            "format": "SQLite",
        }

    @api.get("/api/model")
    async def model():
        return service.store.get("model", "default") or {}

    @api.post("/api/model")
    async def save_model(body: dict):
        value = {
            k: body[k] for k in ("base_url", "model", "credential_ref") if k in body
        }
        if not isinstance(value.get("base_url", ""), str) or (
            value.get("base_url")
            and urlsplit(value["base_url"]).scheme not in ("http", "https")
        ):
            raise ApiError("model_url", "模型端点使用 HTTP 或 HTTPS 地址")
        return service.store.put("model", "default", value)

    @api.post("/api/credentials")
    async def credential(body: dict):
        import keyring

        try:
            keyring.set_password("hundun-client", body["name"], body["value"])
        except Exception:
            raise ApiError(
                "credential_store", "系统凭据库暂不可用，可使用 env:变量名 配置凭据引用"
            )
        return {"credential_ref": "keyring:" + body["name"]}

    @api.get("/api/tasks")
    async def task_list():
        return {"tasks": service.store.all("task")}

    @api.post("/api/tasks")
    async def task_create(body: dict):
        return tasks.create(body)

    @api.post("/api/tasks/{key}/cancel")
    async def task_cancel(key: str):
        return await tasks.cancel(key)

    @api.post("/api/tasks/{key}/followup")
    async def task_followup(key: str, body: dict):
        return tasks.followup(key, body["message"])

    @api.get("/api/tasks/{key}/events")
    async def events(key: str, request: Request, after: int = 0, stream: bool = False):
        if not stream:
            return {"events": service.store.events(key, after)}

        async def generate():
            cursor = max(after, int(request.headers.get("last-event-id", "0")))
            while not await request.is_disconnected():
                items = service.store.events(key, cursor)
                for item in items:
                    cursor = item["seq"]
                    yield (
                        "id: "
                        + str(cursor)
                        + "\ndata: "
                        + json.dumps(item, ensure_ascii=False)
                        + "\n\n"
                    )
                if not items:
                    yield ": keepalive\n\n"
                await asyncio.sleep(1)

        return StreamingResponse(
            generate(),
            media_type="text/event-stream",
            headers={"Cache-Control": "no-cache"},
        )

    @api.get("/api/assets/{name}")
    async def asset(name: str):
        record = service.store.get("asset", name)
        if not record:
            raise ApiError("asset_missing", "产物不存在", 404)
        path = Path(record["path"])
        if not path.is_file():
            raise ApiError("asset_missing", "产物文件待恢复", 404)
        return FileResponse(str(path))

    web = Path(__file__).parent / "web"
    if not web.is_dir():
        web = Path(__file__).parent.parent / "ui/dist"
    if web.is_dir():
        api.mount("/", StaticFiles(directory=str(web), html=True), name="ui")
    return api


app = create_app()
