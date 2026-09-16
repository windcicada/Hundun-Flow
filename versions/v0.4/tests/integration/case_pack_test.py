# SPDX-License-Identifier: Apache-2.0
"""Filesystem contract tests; run directly with the pack tool path."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile

spec = importlib.util.spec_from_file_location('pack',sys.argv[1])
pack = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pack)
with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    source = root/'source'
    source.mkdir()
    (source/'data').write_bytes(b'original checkpoint')
    (source/'current').symlink_to('data')
    (source/'bin').mkdir()
    (source/'bin/run').write_text('#!/bin/sh\nexit 0\n')
    (source/'bin/run').chmod(0o755)
    pack.copy_tree(source,root/'copy')
    assert pack.file_manifest(source) == pack.file_manifest(root/'copy')
    manifest = dict(format='hundun_native_case_package_v1',files=pack.file_manifest(source))
    (source/'manifest.json').write_text(json.dumps(manifest))
    assert pack.verify(source) == 3
    failures = 0
    def reject(call):
        global failures
        try:
            call()
        except ValueError:
            failures += 1
        else:
            raise AssertionError('invalid package accepted')
    (source/'data').write_bytes(b'changed checkpoint')
    reject(lambda:pack.verify(source))
    (source/'data').write_bytes(b'original checkpoint')
    (source/'bin/run').chmod(0o644)
    reject(lambda:pack.verify(source))
    (source/'bin/run').chmod(0o755)
    (source/'escape').symlink_to(root/'copy/data')
    reject(lambda:pack.copy_tree(source,root/'bad'))
    assert not (root/'bad').exists()
    (source/'escape').unlink()
    (source/'escape').symlink_to('missing')
    reject(lambda:pack.copy_tree(source,root/'bad'))
    (source/'escape').unlink()
    os.mkfifo(str(source/'fifo'))
    reject(lambda:pack.copy_tree(source,root/'bad'))
    (source/'fifo').unlink()
    altered = dict(format=manifest['format'],files={'../copy/data':manifest['files']['data']})
    (source/'manifest.json').write_text(json.dumps(altered))
    reject(lambda:pack.verify(source))
    (source/'manifest.json').write_text(json.dumps(manifest))
    assert pack.verify(source) == 3
    assert pack.file_manifest(root/'copy') == manifest['files']
    print('internal links, exact copies, executable modes, six failures and source preservation passed')
