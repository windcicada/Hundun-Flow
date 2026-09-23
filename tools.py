"""Build precompiled frontend into the independently installable wheel."""

import pathlib
import shutil
import subprocess
import sys

root = pathlib.Path(__file__).parent.resolve()
if sys.argv[1:] == ["build"]:
    subprocess.run(
        [
            "npm.cmd" if sys.platform == "win32" else "npm",
            "run",
            "build",
            "--prefix",
            str(root / "ui"),
        ],
        check=True,
    )
    shutil.rmtree(root / "app/web", ignore_errors=True)
    shutil.copytree(root / "ui/dist", root / "app/web")
    subprocess.run([sys.executable, "-m", "build"], cwd=root, check=True)
else:
    raise SystemExit("Usage: python tools.py build")
