"""Inspect independently built client archives and product naming."""

import hashlib
import tarfile
import zipfile
from pathlib import Path

root = Path(__file__).resolve().parent.parent
terms = ("co" + "ast", "co" + "dex")
files = sorted((root / "dist").glob("*"))
checks = []
for path in files:
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as z:
            items = [(n, z.read(n)) for n in z.namelist() if not n.endswith("/")]
    elif path.name.endswith(".tar.gz"):
        with tarfile.open(path) as z:
            items = [(n.name, z.extractfile(n).read()) for n in z if n.isfile()]
    else:
        continue
    for name, body in items:
        value = name.lower() + "\n" + body.decode("utf-8", errors="ignore").lower()
        assert not any(t in value for t in terms), (path.name, name)
        assert "/node_modules/" not in name and "/check/" not in name, name
    assert any("/web/index.html" in n for n, _ in items)
    assert any(n.endswith("licenses.txt") for n, _ in items)
    checks.append(hashlib.sha256(path.read_bytes()).hexdigest() + "  " + path.name)
(root / "dist/SHA256SUMS").write_text("\n".join(checks) + "\n")
print("PASS client wheel/source naming, frontend, licenses, archive scope")
