#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Exercise historical provenance through the supported receipt CLI only."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    repository = Path(__file__).resolve().parents[4]
    tool = repository / "tools/v04_literature_extract.py"
    with tempfile.TemporaryDirectory(prefix="hundun-receipt-cli-") as temporary:
        root = Path(temporary)
        old = root / "not-present" / "old-checkout"
        bundle = root / "relocated"
        bundle.mkdir()
        store = root / "historical-data"
        store.mkdir()
        reference = bundle / "reference.json"
        shutil.copyfile(str(repository / "docs/references/cylinder-re3900-norberg.json"),
                        str(reference))
        payload = bundle / "payload.dat"
        payload.write_bytes(b"original attachment\n")
        # Historical scripts are inert data: running this snapshot must fail.
        historical_bytes = b"raise RuntimeError('never execute historical data')\n"
        historical_sha = digest(historical_bytes)
        historical = store / (historical_sha + ".snapshot")
        historical.write_bytes(historical_bytes)
        document = {
            "schema": "HUNDUN_V04_LITERATURE_RECEIPT_V1",
            "complete": False,
            "incomplete_references": [str(old / reference.name)],
            "artifacts": [
                {"kind": "reference", "path": str(old / reference.name),
                 "sha256": digest(reference.read_bytes())},
                {"kind": "attachment", "path": str(old / payload.name),
                 "sha256": digest(payload.read_bytes())},
                {"kind": "extractor", "path": str(old / "extractor.py"),
                 "sha256": historical_sha},
            ],
        }
        receipt = root / "receipt.json"
        receipt.write_text(json.dumps(document), encoding="utf-8")
        original_receipt = receipt.read_bytes()
        base = [sys.executable, str(tool), "receipt-validate", str(receipt)]
        historical_options = ["--historical", "--relocate-root",
                              "{}={}".format(old, bundle),
                              "--artifact-store", str(store)]
        checks = 0

        def run(arguments, success, expected_error=None):
            result = subprocess.run(base + arguments, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, universal_newlines=True)
            assert (result.returncode == 0) == success, result.stderr
            if expected_error:
                assert expected_error in result.stderr, result.stderr
            assert receipt.read_bytes() == original_receipt, "validator rewrote receipt"
            return result

        result = run(historical_options + ["--json"], True)
        report = json.loads(result.stdout)
        assert report["complete"] is False
        assert report["mode"] == "historical"
        assert report["bound_extractor_sha256"] == historical_sha
        assert report["validator_sha256"] == digest(tool.read_bytes())
        checks += 1
        run([], False, "missing artifact")
        run(historical_options + ["--require-complete"], False, "authority is incomplete")
        run(historical_options[1:], False, "require --historical")
        # A string-prefix sibling is not the same directory root.
        run(["--historical", "--relocate-root", "{}={}".format(str(old)[:-1], bundle),
             "--artifact-store", str(store)], False, "missing artifact")
        checks += 4
        for artifact in (reference, payload, historical):
            original = artifact.read_bytes()
            artifact.write_bytes(original + b"corrupted\n")
            run(historical_options, False, "hash mismatch")
            artifact.unlink()
            run(historical_options, False, "missing artifact")
            artifact.write_bytes(original)
            checks += 2
        run(historical_options + ["--relocate-root", "{}={}".format(old, store)],
            False, "duplicate relocation root")
        checks += 1
        # Current receipts retain the default strict validator binding.
        current = root / "current.json"
        subprocess.run([sys.executable, str(tool), "receipt", "--output", str(current),
                        "--reference", str(reference), "--attachment", str(payload)],
                       check=True)
        subprocess.run([sys.executable, str(tool), "receipt-validate", str(current)],
                       check=True)
        assert json.loads(current.read_text())["complete"] is False
        checks += 1
        print("receipt CLI checks={} complete=false history_unchanged=true".format(checks))


if __name__ == "__main__":
    main()
