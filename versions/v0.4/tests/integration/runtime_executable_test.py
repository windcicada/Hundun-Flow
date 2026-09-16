# SPDX-License-Identifier: Apache-2.0
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
binary = Path(sys.argv[1]).resolve()
expected = hashlib.sha256(binary.read_bytes()).hexdigest()
with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    moved = root/'program with spaces'
    shutil.copy2(str(binary),str(moved))
    for prefix in ([], ['/lib64/ld-linux-x86-64.so.2']):
        for executable in (str(binary),'./program with spaces'):
            r = subprocess.run(prefix+[executable],cwd=str(root),stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,universal_newlines=True,timeout=30)
            assert r.returncode==0 and r.stdout.strip()==expected, (prefix,executable,r.returncode,r.stdout,r.stderr)
print('direct and explicit-loader program hashes match the actual executable after relocation')
