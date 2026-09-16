#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import json
from pathlib import Path
import subprocess
import sys
import tempfile
binary,mpi=map(lambda p:str(Path(p).resolve()),sys.argv[1:3])
def call(*args,ok=True):
    r=subprocess.run([mpi,'--oversubscribe','--bind-to','none','-n','2',binary,*map(str,args)],
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,timeout=45)
    assert (r.returncode==0)==ok,r.stdout[-5000:]
    return r.stdout
with tempfile.TemporaryDirectory(prefix='hf-info-') as tmp:
    root=Path(tmp);case=root/'c';run=root/'r'
    call('init-case','--output',case)
    call('run',case,'--output',run,'--steps',1,'--max-dt','1e-5','--output-interval',0)
    info=json.loads(call('restart-info',run/'Restart'))
    assert info['step']==1 and info['source_ranks']==2 and info['rank_file_bytes']>0
    assert info['integrity']=='manifest_and_rank_checksums'
    generation=run/'Restart'/info['generation']
    path=next(generation.glob('rank-*.bin'))
    with path.open('r+b') as f:
        f.seek(-10,2);b=f.read(1);f.seek(-1,1);f.write(bytes([b[0]^1]))
    call('restart-info',run/'Restart',ok=False)
    assert json.loads(call('restart-info',run/'Restart','--metadata-only'))['integrity']=='manifest_checksum'
    with (generation/'manifest.bin').open('r+b') as f:
        f.seek(20);b=f.read(1);f.seek(20);f.write(bytes([b[0]^1]))
    call('restart-info',run/'Restart','--metadata-only',ok=False)
print('Restart inspect: metadata, streamed checksums and corruption detection PASS')
