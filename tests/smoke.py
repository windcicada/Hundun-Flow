"""Opt-in native smoke: local and SSH, bounded to an independent small case."""

import argparse
import asyncio
import json
import sys
import uuid
from pathlib import Path
from app.service import Service
from tests.ssh import start


async def check(service, host, source):
    async def op(operation, **p):
        if operation in ("prepare", "start", "resume", "control", "pack"):
            p.setdefault("request_id", uuid.uuid4().hex)
        return await service.execute(host, operation, p)

    prepared = await op("prepare", source=str(source))
    case = prepared["case"]["id"]
    launched = await op("start", case_id=case, ranks=2, steps=4)
    run = launched["run"]["id"]
    receipt = False
    controls = {}
    for _ in range(90):
        d = await op("detail", run_id=run)
        if d["run"]["capabilities"].get("output"):
            controls["output"] = await op("control", run_id=run, action="output")
            controls["pause"] = await op("control", run_id=run, action="pause")
            receipt = True
            break
        await asyncio.sleep(0.2)
    if not receipt:
        raise RuntimeError("Native control window missed: " + json.dumps(d["run"]))
    for _ in range(180):
        d = await op("detail", run_id=run)
        if d["run"]["alive"] is False and d["run"]["status"] == "stopped":
            break
        await asyncio.sleep(0.25)
    assert d["run"]["status"] == "stopped", d["run"]
    assert d["run"]["has_restart"], d["run"]
    evidence = Path(d["run"]["path"]) / "control.jsonl"
    controls["receipts"] = evidence.read_text() if evidence.exists() else ""
    restored = await op("resume", run_id=run, steps=1, ranks=2)
    new = restored["run"]["id"]
    for _ in range(180):
        nextd = await op("detail", run_id=new)
        if nextd["run"]["alive"] is False and nextd["run"]["status"] == "completed":
            break
        await asyncio.sleep(0.25)
    assert nextd["run"]["status"] == "completed", nextd["run"]
    assert nextd["run"]["step"] > d["run"]["step"]
    plot = await op("plot", run_id=new, options={"variable": "temperature"})
    report = await op("report", run_id=new)
    archive = await op("pack", run_id=new)
    return {
        "host": host,
        "case": case,
        "stopped": d["run"],
        "resumed": nextd["run"],
        "controls": controls,
        "plot": plot,
        "report": report,
        "archive": archive,
    }


async def main(args):
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    service = Service(root / "state")
    p = {
        "id": "local",
        "kind": "local",
        "roots": [str(Path(args.case).resolve())],
        "program": str(Path(args.program).resolve()),
        "mpi": str(Path(args.mpi).resolve()),
        "max_ranks": 2,
        "max_steps": 4,
    }
    service.save_host(p)
    results = [await check(service, "local", Path(args.case))]
    server, remote = await start(root)
    try:
        remote.update(p, id="ssh", kind="ssh")
        service.save_host(remote)
        # Use this check environment's installed plotting dependencies.
        plot = Path(remote["work_dir"]) / ".client/plot"
        plot.parent.mkdir(parents=True, exist_ok=True)
        plot.symlink_to(
            Path(sys.executable).absolute().parent.parent, target_is_directory=True
        )
        results.append(await check(service, "ssh", Path(args.case)))
    finally:
        server.close()
        await server.wait_closed()
    (root / "result.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    print(
        json.dumps(
            [{k: r[k] for k in ("host", "case", "plot")} for r in results],
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", required=True)
    p.add_argument("--program", required=True)
    p.add_argument("--mpi", required=True)
    p.add_argument("--output", required=True)
    asyncio.run(main(p.parse_args()))
