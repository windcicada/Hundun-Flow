"""Loopback-only native run API. Start with uvicorn ui.backend:app."""
import os
from pathlib import Path
from urllib.parse import urlsplit

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, StrictInt, StrictStr

try:
    from . import data
    from .jobs import ApiError, Registry
except ImportError:
    import data
    from jobs import ApiError, Registry


class Control(BaseModel):
    action: StrictStr
    request_id: StrictStr


class Launch(BaseModel):
    steps: StrictInt
    ranks: StrictInt
    request_id: StrictStr


class ImportRun(BaseModel):
    path: StrictStr


def create_app(registry=None):
    api = FastAPI(title='Hundun 工作台', docs_url=None, redoc_url=None)
    registry = registry or Registry()
    api.state.registry = registry

    @api.exception_handler(ApiError)
    async def api_error(request, exc):
        return JSONResponse({'error': {'code': exc.code, 'message': exc.message}}, status_code=exc.status)

    @api.exception_handler(RequestValidationError)
    async def invalid_input(request, exc):
        return JSONResponse({'error': {'code': 'invalid_input', 'message': '请求字段或数值类型不符合接口要求'}}, status_code=422)

    @api.exception_handler(OSError)
    async def io_error(request, exc):
        return JSONResponse({'error': {'code': 'io_error', 'message': str(exc)}}, status_code=409)

    @api.middleware('http')
    async def local_only(request: Request, call_next):
        host = request.headers.get('host', '')
        try:
            parsed = urlsplit('http://'+host)
        except ValueError:
            return JSONResponse({'error': {'code': 'invalid_host', 'message': '主机名无效'}}, status_code=403)
        if parsed.hostname not in ('localhost', '127.0.0.1', '::1'):
            return JSONResponse({'error': {'code': 'invalid_host', 'message': '仅允许本机访问'}}, status_code=403)
        if request.method not in ('GET', 'HEAD', 'OPTIONS'):
            origin = request.headers.get('origin')
            # Vite's proxy preserves its local Origin and Host. External origins
            # and cross-port browser writes are rejected. CLI writes omit Origin.
            if origin and origin.rstrip('/') not in ('http://'+host, 'https://'+host):
                return JSONResponse({'error': {'code': 'invalid_origin', 'message': '写操作须由同源界面发起'}}, status_code=403)
            if request.headers.get('sec-fetch-site') == 'cross-site':
                return JSONResponse({'error': {'code': 'invalid_origin', 'message': '写操作须由同源界面发起'}}, status_code=403)
        return await call_next(request)

    @api.get('/api/health')
    def health():
        return {'ok': True, 'program': str(registry.program), 'program_available': registry.program_available(),
            'case_roots': [str(p) for p in registry.roots], 'state': str(registry.state), 'api_version': 1}

    @api.get('/api/runs')
    def runs():
        registry.refresh()
        output = []
        for key, path in sorted(list(registry.runs.items())):
            try:
                output.append(registry.summary(key))
            except Exception as exc:
                output.append({'id': key, 'name': path.name, 'path': str(path), 'status': 'unavailable',
                    'case_path': None, 'case_binding': 'unknown', 'step': None, 'time': None, 'dt': None,
                    'mpi': None, 'mesh': None, 'scheme': None, 'coupling': None, 'models': [],
                    'seconds_per_step': None, 'target_steps': None, 'updated_at': None, 'has_restart': False,
                    'pending_controls': [], 'alive': None, 'process_status': 'unknown', 'capabilities': {'start': False, 'resume': False, 'pause': False, 'output': False},
                    'error': {'code': getattr(exc, 'code', 'run_read_failed'), 'message': str(getattr(exc, 'message', exc))}})
        return {'runs': data.safe(output)}

    @api.post('/api/runs/import')
    def import_run(body: ImportRun):
        return {'run': data.safe(registry.import_run(body.path))}

    @api.get('/api/runs/{key}')
    def run(key: str):
        path = registry.resolve_run(key)
        case = registry.case_for(key)
        return data.detail(path, registry.summary(key), data.read_json(case/'case.json') if case else None)

    @api.post('/api/runs/{key}/control')
    def control(key: str, body: Control):
        return registry.control(key, body.action, body.request_id)

    @api.post('/api/runs/{key}/resume')
    def resume(key: str, body: Launch):
        return registry.launch(key, body.steps, body.ranks, body.request_id, resume=True)

    @api.get('/api/cases')
    def cases():
        registry.refresh()
        return {'cases': [{'id': key, 'name': path.name, 'path': str(path), 'config': data.safe(data.read_json(path/'case.json')),
            'capabilities': {'start': registry.program_available()}} for key, path in sorted(registry.cases.items())]}

    @api.post('/api/cases/{key}/start')
    def start(key: str, body: Launch):
        return registry.launch(key, body.steps, body.ranks, body.request_id)

    def plot_module():
        try:
            if __package__:
                from . import plot
            else:
                import plot
            return plot
        except ImportError:
            raise ApiError('plot_unavailable', '绘图模块或依赖尚未就绪', 503)

    @api.get('/api/runs/{key}/fields')
    def fields(key: str):
        plot = plot_module()
        try:
            return plot.discover_fields(registry.resolve_run(key))
        except plot.PlotError as exc:
            raise ApiError(exc.code, str(exc))

    @api.post('/api/runs/{key}/plot')
    def render(key: str, body: dict):
        plot = plot_module()
        output = registry.state/'plots'
        output.mkdir(exist_ok=True)
        try:
            result = plot.render_slice(registry.resolve_run(key), body, output)
        except plot.PlotError as exc:
            raise ApiError(exc.code, str(exc))
        path = Path(result['path']).resolve()
        if path.parent != output.resolve():
            raise ApiError('invalid_asset', '绘图结果位置异常', 500)
        return {**result, 'url': '/api/assets/'+path.name}

    @api.get('/api/assets/{name}')
    def asset(name: str):
        root = (registry.state/'plots').resolve()
        path = (root/name).resolve()
        if path.parent != root or path.suffix.lower() not in ('.png', '.jpg', '.svg', '.pdf') or not path.is_file():
            raise ApiError('asset_missing', '图像不存在', 404)
        return FileResponse(str(path))

    dist = Path(__file__).resolve().parent/'dist'
    if dist.is_dir():
        api.mount('/', StaticFiles(directory=str(dist), html=True), name='ui')
    return api


app = create_app()
