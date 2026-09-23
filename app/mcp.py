"""Typed simulation tools backed by the same operations as the browser."""

import sys
import uuid
from typing import Literal
from mcp.server.fastmcp import FastMCP
from .service import Service

service = Service(sys.argv[1])
task_id = sys.argv[2]
server = FastMCP("Hundun simulation")


async def call(action, params):
    task = service.store.get("task", task_id)
    try:
        return await service.execute(task["host_id"], action, params, task_id)
    except Exception as exc:
        return {"error": str(exc), "code": getattr(exc, "code", "operation_failed")}


def request(value):
    return value or uuid.uuid4().hex


@server.tool()
async def workspace_list(
    kind: Literal["cases", "runs", "revisions", "probe"] = "cases",
) -> dict:
    """List real registered cases, runs, input revisions or solver capabilities."""
    return await call(kind, {})


@server.tool()
async def case_prepare(
    source: str = "",
    case_id: str = "",
    config: dict | None = None,
    reason: str = "",
    request_id: str = "",
) -> dict:
    """Copy a native input directory or registered case to a new revision. Omit source and case_id to use the connected solver template. config replaces case.json as a complete object."""
    p = {
        "source": source,
        "case_id": case_id,
        "reason": reason,
        "request_id": request(request_id),
    }
    if config is not None:
        p["config"] = config
    return await call("prepare", p)


@server.tool()
async def case_check(case_id: str, ranks: int = 1) -> dict:
    """Ask the native solver to validate a registered input revision."""
    return await call("case_check", {"case_id": case_id, "ranks": ranks})


@server.tool()
async def job_start(
    case_id: str = "",
    run_id: str = "",
    steps: int = 1,
    ranks: int = 1,
    request_id: str = "",
) -> dict:
    """Start a new run, or resume run_id into a new directory. Use one of case_id and run_id."""
    return await call(
        "resume" if run_id else "start",
        {
            "case_id": case_id,
            "run_id": run_id,
            "steps": steps,
            "ranks": ranks,
            "request_id": request(request_id),
        },
    )


@server.tool()
async def job_status(run_id: str) -> dict:
    """Read actual accepted-step history, budgets, logs, process identity and checkpoint evidence."""
    return await call("detail", {"run_id": run_id})


@server.tool()
async def job_control(
    run_id: str, action: Literal["output", "pause"], request_id: str = ""
) -> dict:
    """Request field output or save checkpoint then stop. Completion is confirmed by native receipts."""
    return await call(
        "control",
        {"run_id": run_id, "action": action, "request_id": request(request_id)},
    )


@server.tool()
async def result_render(run_id: str, options: dict | None = None) -> dict:
    """Without options list real field frames/variables; otherwise render a slice with frame, variable, normal x/y/z and coordinate center or meters."""
    return await call(
        "fields" if options is None else "plot", {"run_id": run_id, "options": options}
    )


@server.tool()
async def result_report(run_id: str) -> dict:
    """Export the real run record and budgets; this is not a physical validation certificate."""
    return await call("report", {"run_id": run_id})


@server.tool()
async def case_pack(
    run_id: str,
    checkpoint: bool = True,
    runtime: str = "",
    request_id: str = "",
    results: list[str] | None = None,
) -> dict:
    """Archive inputs, evidence and accepted checkpoint, optionally with a selected external runtime directory."""
    return await call(
        "pack",
        {
            "run_id": run_id,
            "checkpoint": checkpoint,
            "runtime": runtime,
            "results": results or [],
            "request_id": request(request_id),
        },
    )


@server.tool()
async def geometry_tool(tool: str, path: str, request_id: str = "") -> dict:
    """Run a host-configured geometry/mesh command on an input directory; command definitions belong to host settings."""
    return await call(
        "external", {"tool": tool, "path": path, "request_id": request(request_id)}
    )


@server.tool()
def finish_task(summary: str) -> dict:
    """Mark the user's entire goal completed only after checking all required evidence."""
    t = service.store.get("task", task_id)
    if t["status"] == "cancelled":
        return {"status": "cancelled"}
    t.update(status="completed", summary=summary)
    service.store.put("task", task_id, t)
    service.store.event(task_id, {"type": "completed", "message": summary})
    return {"status": "completed"}


@server.tool()
def ask_user(question: str) -> dict:
    """Pause when required physical definitions or user inputs are missing."""
    t = service.store.get("task", task_id)
    if t["status"] == "cancelled":
        return {"status": "cancelled"}
    t.update(status="needs_input")
    service.store.put("task", task_id, t)
    service.store.event(task_id, {"type": "question", "message": question})
    return {"status": "needs_input"}


if __name__ == "__main__":
    server.run(transport="stdio")
