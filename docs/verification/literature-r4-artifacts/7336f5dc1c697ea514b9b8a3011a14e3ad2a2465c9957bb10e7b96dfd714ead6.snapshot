#!/usr/bin/env python3
"""Validate frozen literature data and write an immutable hash receipt."""

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path


SCHEMA = "HUNDUN_V04_LITERATURE_V1"
PARNAUDEAU_DOI = "10.1063/1.2957018"
PARNAUDEAU_STATIONS = (1.06, 1.54, 2.02)
PARNAUDEAU_PROFILE_FIGURES = {
    "mean_u_over_uc": 11,
    "mean_v_over_uc": 12,
    "uu_over_uc2": 13,
    "vv_over_uc2": 14,
    "uv_over_uc2": 15,
}


class LiteratureError(ValueError):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def finite_number(value):
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def valid_sha256(value):
    return (isinstance(value, str) and len(value) == 64 and
            all(character in "0123456789abcdef" for character in value))


def validate_parnaudeau_profiles(path, document):
    source = document.get("source", {})
    if str(source.get("doi", "")).lower() != PARNAUDEAU_DOI:
        return
    if document.get("pending_profiles") is not None:
        return
    authority = document.get("profile_authority")
    if not isinstance(authority, dict) or authority.get("status") not in (
            "controlled_digitization_complete", "author_arrays_complete"):
        raise LiteratureError(
            "{}: completed Parnaudeau reference lacks profile authority".format(path))
    source_artifact = authority.get("source_artifact")
    if (not isinstance(source_artifact, dict) or
            source_artifact.get("kind") not in (
                "publisher_vor", "author_uploaded_primary",
                "institution_primary", "author_arrays") or
            not source_artifact.get("public_url") or
            not source_artifact.get("acquired_utc") or
            not source_artifact.get("provenance") or
            not valid_sha256(source_artifact.get("sha256")) or
            not isinstance(source_artifact.get("bytes"), int) or
            source_artifact["bytes"] <= 0):
        raise LiteratureError(
            "{}: invalid Parnaudeau primary source artifact".format(path))
    artifact_hashes = authority.get("receipt_attachment_sha256")
    if (not isinstance(artifact_hashes, list) or not artifact_hashes or
            any(not valid_sha256(value) for value in artifact_hashes) or
            len(set(artifact_hashes)) != len(artifact_hashes) or
            source_artifact["sha256"] not in artifact_hashes):
        raise LiteratureError(
            "{}: Parnaudeau digitization artifacts are not enumerated".format(path))
    if (not isinstance(authority.get("rule_document"), str) or
            not authority["rule_document"].strip() or
            not valid_sha256(authority.get("rule_document_sha256")) or
            authority["rule_document_sha256"] not in artifact_hashes):
        raise LiteratureError(
            "{}: controlled-digitization rule is not bound".format(path))
    roles = authority.get("artifact_roles")
    if not isinstance(roles, dict):
        raise LiteratureError(
            "{}: Parnaudeau artifact roles are missing".format(path))
    if (authority["status"] == "controlled_digitization_complete" and
            (source_artifact["kind"] == "author_arrays" or
             len(artifact_hashes) < 5)):
        raise LiteratureError(
            "{}: controlled digitization must bind all artifact roles".format(path))
    if authority["status"] == "controlled_digitization_complete":
        required_roles = {
            "source_primary", "rule", "calibration", "raw_trace",
            "extraction_script"}
        if (set(roles) != required_roles or
                any(not valid_sha256(value) for value in roles.values()) or
                len(set(roles.values())) != len(roles) or
                roles["source_primary"] != source_artifact["sha256"] or
                roles["rule"] != authority["rule_document_sha256"] or
                not set(roles.values()).issubset(set(artifact_hashes))):
            raise LiteratureError(
                "{}: controlled-digitization artifact roles are incomplete".format(path))
    if (authority["status"] == "author_arrays_complete" and
            source_artifact["kind"] != "author_arrays"):
        raise LiteratureError(
            "{}: author-array authority has the wrong source kind".format(path))
    if authority["status"] == "author_arrays_complete":
        required_roles = {"source_author_arrays", "rule", "field_dictionary"}
        if (set(roles) != required_roles or
                any(not valid_sha256(value) for value in roles.values()) or
                len(set(roles.values())) != len(roles) or
                roles["source_author_arrays"] != source_artifact["sha256"] or
                roles["rule"] != authority["rule_document_sha256"] or
                not set(roles.values()).issubset(set(artifact_hashes))):
            raise LiteratureError(
                "{}: author-array artifact roles are incomplete".format(path))
    interval = authority.get("common_y_over_d_interval")
    if (not isinstance(interval, list) or len(interval) != 2 or
            any(not finite_number(value) for value in interval) or
            interval[0] >= interval[1]):
        raise LiteratureError(
            "{}: invalid common Parnaudeau transverse interval".format(path))
    method = authority.get("digitization_method")
    required_method = (
        "axis_calibration", "curve_selection", "sampling_rule",
        "error_composition", "independent_repeat_rule")
    if (not isinstance(method, dict) or
            any(not isinstance(method.get(key), str) or not method[key].strip()
                for key in required_method)):
        raise LiteratureError(
            "{}: incomplete controlled-digitization method".format(path))
    profiles = authority.get("profiles")
    if not isinstance(profiles, list):
        raise LiteratureError("{}: missing Parnaudeau profiles".format(path))
    expected = set((station, quantity)
                   for station in PARNAUDEAU_STATIONS
                   for quantity in PARNAUDEAU_PROFILE_FIGURES)
    seen = set()
    for profile in profiles:
        if not isinstance(profile, dict):
            raise LiteratureError("{}: invalid Parnaudeau profile".format(path))
        station = profile.get("station_x_over_d")
        quantity = profile.get("quantity")
        key = (station, quantity)
        if key not in expected or key in seen:
            raise LiteratureError(
                "{}: duplicate or unexpected Parnaudeau profile".format(path))
        seen.add(key)
        if (profile.get("figure") != PARNAUDEAU_PROFILE_FIGURES[quantity] or
                not isinstance(profile.get("legend_identity"), str) or
                not profile["legend_identity"].strip() or
                not isinstance(profile.get("source_locator"), str) or
                not profile["source_locator"].strip() or
                not isinstance(profile.get("extraction_method"), str) or
                not profile["extraction_method"].strip()):
            raise LiteratureError(
                "{}: incomplete Parnaudeau profile locator".format(path))
        uncertainty = profile.get("extraction_uncertainty")
        if (not isinstance(uncertainty, dict) or
                uncertainty.get("status") not in (
                    "controlled_digitization_bound", "author_array_precision") or
                uncertainty.get("composition") not in (
                    "conservative_l1", "not_applicable_author_array") or
                not finite_number(uncertainty.get("coordinate_absolute")) or
                uncertainty["coordinate_absolute"] < 0 or
                not finite_number(uncertainty.get("value_absolute")) or
                uncertainty["value_absolute"] < 0 or
                uncertainty.get("experimental_uncertainty_status") !=
                "not_reported_separately"):
            raise LiteratureError(
                "{}: invalid Parnaudeau extraction uncertainty".format(path))
        data = profile.get("data")
        if not isinstance(data, list) or len(data) < 8:
            raise LiteratureError(
                "{}: Parnaudeau profile requires at least eight points".format(path))
        prior = None
        for row in data:
            if (not isinstance(row, list) or len(row) != 2 or
                    not finite_number(row[0]) or not finite_number(row[1])):
                raise LiteratureError(
                    "{}: invalid Parnaudeau profile row".format(path))
            if prior is not None and row[0] <= prior:
                raise LiteratureError(
                    "{}: Parnaudeau y/D coordinate not increasing".format(path))
            prior = row[0]
        if data[0][0] > interval[0] or data[-1][0] < interval[1]:
            raise LiteratureError(
                "{}: Parnaudeau profile misses common transverse interval".format(path))
    if seen != expected:
        raise LiteratureError(
            "{}: Parnaudeau profile matrix is incomplete".format(path))


def required_attachment_hashes(document):
    source = document.get("source", {})
    if (str(source.get("doi", "")).lower() != PARNAUDEAU_DOI or
            document.get("pending_profiles") is not None):
        return set()
    return set(document["profile_authority"]["receipt_attachment_sha256"])


def validate_reference(path):
    document = load(path)
    if not isinstance(document, dict) or document.get("schema") != SCHEMA:
        raise LiteratureError("{}: wrong schema".format(path))
    source = document.get("source")
    if not isinstance(source, dict):
        raise LiteratureError("{}: missing source".format(path))
    for key in ("title", "authors", "year", "primary_url", "locator"):
        if not source.get(key):
            raise LiteratureError("{}: source missing {}".format(path, key))
    observables = document.get("observables")
    if not isinstance(observables, list) or not observables:
        raise LiteratureError("{}: no observables".format(path))
    names = set()
    for observable in observables:
        if not isinstance(observable, dict):
            raise LiteratureError("{}: invalid observable".format(path))
        for key in ("name", "units", "definition", "source_locator",
                    "extraction_method", "uncertainty"):
            if key not in observable:
                raise LiteratureError(
                    "{}: observable missing {}".format(path, key))
        name = observable["name"]
        if not isinstance(name, str) or not name or name in names:
            raise LiteratureError("{}: duplicate observable".format(path))
        names.add(name)
        uncertainty = observable["uncertainty"]
        if not isinstance(uncertainty, dict) or uncertainty.get("status") not in (
                "reported", "not_reported", "derived"):
            raise LiteratureError("{}: invalid uncertainty authority".format(path))
        if uncertainty["status"] == "reported":
            if (not finite_number(uncertainty.get("value")) or
                    uncertainty["value"] < 0 or not uncertainty.get("kind")):
                raise LiteratureError("{}: invalid reported uncertainty".format(path))
        elif uncertainty.get("value") is not None:
            raise LiteratureError("{}: unavailable uncertainty must be null".format(path))
        has_value = finite_number(observable.get("value"))
        has_data = isinstance(observable.get("data"), list) and observable["data"]
        if has_value == bool(has_data):
            raise LiteratureError(
                "{}: observable requires exactly one value or data array".format(path))
        if has_data:
            prior = None
            for row in observable["data"]:
                if (not isinstance(row, list) or len(row) != 2 or
                        not finite_number(row[0]) or not finite_number(row[1])):
                    raise LiteratureError("{}: invalid data row".format(path))
                if prior is not None and row[0] <= prior:
                    raise LiteratureError("{}: data coordinate not increasing".format(path))
                prior = row[0]
    required_observables = document.get("required_observables", [])
    if (not isinstance(required_observables, list) or
            any(not isinstance(name, str) or not name
                for name in required_observables) or
            len(set(required_observables)) != len(required_observables)):
        raise LiteratureError("{}: malformed required_observables".format(path))
    required_observables = set(required_observables)
    pending_observables = document.get("pending_observables")
    pending_names = set()
    if pending_observables is not None:
        if not isinstance(pending_observables, list) or not pending_observables:
            raise LiteratureError("{}: malformed pending_observables".format(path))
        for pending_observable in pending_observables:
            if (not isinstance(pending_observable, dict) or
                    not isinstance(pending_observable.get("name"), str) or
                    not pending_observable["name"] or
                    not isinstance(pending_observable.get("status"), str) or
                    not pending_observable["status"] or
                    not isinstance(pending_observable.get("reason"), str) or
                    not pending_observable["reason"] or
                    pending_observable["name"] in pending_names):
                raise LiteratureError(
                    "{}: malformed pending_observables".format(path))
            pending_names.add(pending_observable["name"])
    missing_required = required_observables - names
    if pending_names != missing_required:
        raise LiteratureError(
            "{}: required observable coverage disagrees with pending_observables"
            .format(path))
    if names & pending_names:
        raise LiteratureError(
            "{}: observable is both complete and pending".format(path))
    pending = document.get("pending_profiles")
    if pending is not None:
        if (not isinstance(pending, dict) or
                not isinstance(pending.get("stations_x_over_d"), list) or
                not pending["stations_x_over_d"] or
                not isinstance(pending.get("quantities"), list) or
                not pending["quantities"] or
                not isinstance(pending.get("status"), str) or
                not isinstance(pending.get("reason"), str)):
            raise LiteratureError("{}: malformed pending_profiles".format(path))
    validate_parnaudeau_profiles(path, document)
    return document


def reference_complete(document):
    return (document.get("pending_profiles") is None and
            document.get("pending_observables") is None)


def import_csv(path):
    data = []
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        for line, row in enumerate(reader, 1):
            if not row or row[0].lstrip().startswith("#"):
                continue
            if len(row) != 2:
                raise LiteratureError("{}:{} requires two columns".format(path, line))
            values = [float(row[0]), float(row[1])]
            if any(not math.isfinite(value) for value in values):
                raise LiteratureError("{}:{} nonfinite value".format(path, line))
            if data and values[0] <= data[-1][0]:
                raise LiteratureError("{}:{} coordinate not increasing".format(path, line))
            data.append(values)
    if not data:
        raise LiteratureError("{}: empty digitization".format(path))
    return data


def write_receipt(output, references, attachments):
    if output.exists():
        raise LiteratureError("receipt is immutable: {} exists".format(output))
    records = []
    incomplete = []
    required_attachments = set()
    for path in references:
        document = validate_reference(path)
        if not reference_complete(document):
            incomplete.append(str(path.resolve()))
        required_attachments.update(required_attachment_hashes(document))
        records.append({"kind": "reference", "path": str(path.resolve()),
                        "sha256": sha256(path)})
    for path in attachments:
        if not path.is_file():
            raise LiteratureError("missing attachment: {}".format(path))
        records.append({"kind": "attachment", "path": str(path.resolve()),
                        "sha256": sha256(path)})
    supplied_attachments = set(record["sha256"] for record in records
                               if record["kind"] == "attachment")
    if not required_attachments.issubset(supplied_attachments):
        raise LiteratureError("completed references require bound attachments")
    script = Path(__file__).resolve()
    records.append({"kind": "extractor", "path": str(script),
                    "sha256": sha256(script)})
    document = {"schema": "HUNDUN_V04_LITERATURE_RECEIPT_V1",
                "complete": not incomplete,
                "incomplete_references": incomplete,
                "artifacts": records}
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(str(output), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        payload = (json.dumps(document, indent=2, sort_keys=True) + "\n").encode("utf-8")
        cursor = 0
        while cursor < len(payload):
            cursor += os.write(descriptor, payload[cursor:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def validate_receipt(path, require_complete=False):
    document = load(path)
    if (not isinstance(document, dict) or
            document.get("schema") != "HUNDUN_V04_LITERATURE_RECEIPT_V1"):
        raise LiteratureError("{}: wrong receipt schema".format(path))
    artifacts = document.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise LiteratureError("{}: empty receipt".format(path))
    incomplete = []
    required_attachments = set()
    supplied_attachments = set()
    extractor_seen = False
    seen = set()
    for artifact in artifacts:
        if (not isinstance(artifact, dict) or
                artifact.get("kind") not in ("reference", "attachment", "extractor") or
                not isinstance(artifact.get("path"), str)):
            raise LiteratureError("{}: malformed receipt artifact".format(path))
        artifact_path = Path(artifact["path"])
        canonical = str(artifact_path.resolve())
        if canonical in seen:
            raise LiteratureError("{}: duplicate artifact {}".format(
                path, artifact_path))
        seen.add(canonical)
        if artifact["kind"] == "extractor":
            if extractor_seen:
                raise LiteratureError("{}: duplicate extractor".format(path))
            current_extractor = Path(__file__).resolve()
            if artifact.get("sha256") != sha256(current_extractor):
                raise LiteratureError(
                    "{}: receipt binds a different extractor".format(path))
            extractor_seen = True
            continue
        if not artifact_path.is_file():
            raise LiteratureError("{}: missing artifact {}".format(
                path, artifact_path))
        if artifact.get("sha256") != sha256(artifact_path):
            raise LiteratureError("{}: artifact hash mismatch {}".format(
                path, artifact_path))
        if artifact["kind"] == "attachment":
            supplied_attachments.add(artifact["sha256"])
        if artifact["kind"] == "reference":
            reference = validate_reference(artifact_path)
            if not reference_complete(reference):
                incomplete.append(canonical)
            required_attachments.update(required_attachment_hashes(reference))
    if not extractor_seen:
        raise LiteratureError("{}: extractor is not bound".format(path))
    if not required_attachments.issubset(supplied_attachments):
        raise LiteratureError(
            "{}: completed reference attachments are not bound".format(path))
    expected_complete = not incomplete
    if document.get("complete") is not expected_complete:
        raise LiteratureError("{}: complete flag disagrees with references".format(path))
    if sorted(document.get("incomplete_references", [])) != sorted(incomplete):
        raise LiteratureError("{}: incomplete reference list mismatch".format(path))
    if require_complete and not expected_complete:
        raise LiteratureError("{}: literature authority is incomplete".format(path))
    return document


def self_test():
    with tempfile.TemporaryDirectory(prefix="hundun-v04-literature-") as root_text:
        root = Path(root_text)
        base = {
            "schema": SCHEMA,
            "source": {
                "title": "synthetic primary source",
                "authors": "test",
                "year": 2026,
                "primary_url": "https://example.invalid/primary",
                "locator": "table 1",
            },
            "observables": [{
                "name": "value",
                "value": 1.0,
                "units": "1",
                "definition": "synthetic",
                "source_locator": "table 1",
                "extraction_method": "direct",
                "uncertainty": {
                    "status": "reported", "value": 0.1, "kind": "absolute"
                },
            }],
        }
        complete_reference = root / "complete.json"
        complete_reference.write_text(json.dumps(base), encoding="utf-8")
        complete_receipt = root / "complete-receipt.json"
        write_receipt(complete_receipt, [complete_reference], [])
        validate_receipt(complete_receipt, True)

        relocated_document = json.loads(complete_receipt.read_text(
            encoding="utf-8"))
        extractor = next(item for item in relocated_document["artifacts"]
                         if item["kind"] == "extractor")
        extractor["path"] = str(root / "original-checkout-no-longer-present.py")
        relocated_receipt = root / "relocated-receipt.json"
        relocated_receipt.write_text(json.dumps(relocated_document),
                                     encoding="utf-8")
        validate_receipt(relocated_receipt, True)

        pending = json.loads(json.dumps(base))
        pending["pending_profiles"] = {
            "stations_x_over_d": [1.0],
            "quantities": ["mean_u"],
            "status": "pending",
            "reason": "primary array unavailable",
        }
        pending_reference = root / "pending.json"
        pending_reference.write_text(json.dumps(pending), encoding="utf-8")
        pending_receipt = root / "pending-receipt.json"
        write_receipt(pending_receipt, [pending_reference], [])
        validate_receipt(pending_receipt, False)
        try:
            validate_receipt(pending_receipt, True)
            raise AssertionError("incomplete literature receipt accepted")
        except LiteratureError:
            pass
        before = pending_receipt.read_bytes()
        try:
            write_receipt(pending_receipt, [pending_reference], [])
            raise AssertionError("immutable receipt overwritten")
        except LiteratureError:
            if pending_receipt.read_bytes() != before:
                raise AssertionError("failed overwrite changed receipt")

        required = json.loads(json.dumps(base))
        required["required_observables"] = ["value", "force_value"]
        required["pending_observables"] = [{
            "name": "force_value",
            "status": "pending_primary_authority",
            "reason": "direct force authority unavailable",
        }]
        required_path = root / "required-pending.json"
        required_path.write_text(json.dumps(required), encoding="utf-8")
        validate_reference(required_path)
        required_receipt = root / "required-pending-receipt.json"
        write_receipt(required_receipt, [required_path], [])
        validate_receipt(required_receipt, False)
        try:
            validate_receipt(required_receipt, True)
            raise AssertionError("pending required observable accepted")
        except LiteratureError:
            pass
        required_bypass = json.loads(json.dumps(required))
        del required_bypass["pending_observables"]
        required_bypass_path = root / "required-bypass.json"
        required_bypass_path.write_text(
            json.dumps(required_bypass), encoding="utf-8")
        try:
            validate_reference(required_bypass_path)
            raise AssertionError("required observable pending deletion accepted")
        except LiteratureError:
            pass

        parnaudeau = json.loads(json.dumps(base))
        parnaudeau["source"]["doi"] = PARNAUDEAU_DOI
        parnaudeau["pending_profiles"] = {
            "stations_x_over_d": list(PARNAUDEAU_STATIONS),
            "quantities": list(PARNAUDEAU_PROFILE_FIGURES),
            "status": "pending",
            "reason": "primary artifact unavailable",
        }
        parnaudeau_pending = root / "parnaudeau-pending.json"
        parnaudeau_pending.write_text(json.dumps(parnaudeau), encoding="utf-8")
        validate_reference(parnaudeau_pending)
        bypass = json.loads(json.dumps(parnaudeau))
        del bypass["pending_profiles"]
        bypass_path = root / "parnaudeau-bypass.json"
        bypass_path.write_text(json.dumps(bypass), encoding="utf-8")
        try:
            validate_reference(bypass_path)
            raise AssertionError("Parnaudeau pending field deletion accepted")
        except LiteratureError:
            pass

        primary = root / "primary-artifact.pdf"
        primary.write_bytes(b"%PDF-1.4\nsynthetic self-test artifact\n%%EOF\n")
        primary_sha = sha256(primary)
        calibration = root / "calibration.json"
        calibration.write_text('{"synthetic":true}', encoding="utf-8")
        trace = root / "raw-trace.csv"
        trace.write_text("1,2\n", encoding="utf-8")
        rule = root / "digitization-rule.md"
        rule.write_text("synthetic frozen rule\n", encoding="utf-8")
        extraction_script = root / "extract.py"
        extraction_script.write_text("# synthetic\n", encoding="utf-8")
        calibration_sha = sha256(calibration)
        trace_sha = sha256(trace)
        rule_sha = sha256(rule)
        extraction_script_sha = sha256(extraction_script)
        completed = json.loads(json.dumps(bypass))
        completed["profile_authority"] = {
            "status": "controlled_digitization_complete",
            "source_artifact": {
                "kind": "publisher_vor",
                "public_url": "https://example.invalid/primary.pdf",
                "acquired_utc": "2026-08-22T00:00:00Z",
                "provenance": "synthetic self-test publisher artifact",
                "sha256": primary_sha,
                "bytes": primary.stat().st_size,
            },
            "receipt_attachment_sha256": [
                primary_sha, calibration_sha, trace_sha, rule_sha,
                extraction_script_sha],
            "rule_document": "docs/verification/synthetic-rule.md",
            "rule_document_sha256": rule_sha,
            "artifact_roles": {
                "source_primary": primary_sha,
                "rule": rule_sha,
                "calibration": calibration_sha,
                "raw_trace": trace_sha,
                "extraction_script": extraction_script_sha,
            },
            "common_y_over_d_interval": [-1.0, 1.0],
            "digitization_method": {
                "axis_calibration": "independent affine tick fit",
                "curve_selection": "target PIV legend only",
                "sampling_rule": "published marker centroids",
                "error_composition": "conservative component sum",
                "independent_repeat_rule": "two independent traces",
            },
            "profiles": [],
        }
        rows = [[-1.4 + 0.4 * index, 0.01 * index]
                for index in range(8)]
        for station in PARNAUDEAU_STATIONS:
            for quantity, figure in PARNAUDEAU_PROFILE_FIGURES.items():
                completed["profile_authority"]["profiles"].append({
                    "station_x_over_d": station,
                    "quantity": quantity,
                    "figure": figure,
                    "legend_identity": "present PIV",
                    "source_locator": "synthetic figure panel",
                    "extraction_method": "controlled marker digitization",
                    "extraction_uncertainty": {
                        "status": "controlled_digitization_bound",
                        "composition": "conservative_l1",
                        "coordinate_absolute": 0.01,
                        "value_absolute": 0.01,
                        "experimental_uncertainty_status":
                        "not_reported_separately",
                    },
                    "data": rows,
                })
        completed_path = root / "parnaudeau-complete.json"
        completed_path.write_text(json.dumps(completed), encoding="utf-8")
        validate_reference(completed_path)
        for mutation_name, mutation in (
                ("missing-profile", lambda value:
                 value["profile_authority"]["profiles"].pop()),
                ("short-profile", lambda value:
                 value["profile_authority"]["profiles"][0].update(
                     {"data": rows[:7]})),
                ("missing-common-interval", lambda value:
                 value["profile_authority"].pop("common_y_over_d_interval"))):
            mutated = json.loads(json.dumps(completed))
            mutation(mutated)
            mutated_path = root / (mutation_name + ".json")
            mutated_path.write_text(json.dumps(mutated), encoding="utf-8")
            try:
                validate_reference(mutated_path)
                raise AssertionError(
                    "Parnaudeau {} mutation accepted".format(mutation_name))
            except LiteratureError:
                pass
        missing_attachment_receipt = root / "missing-attachment-receipt.json"
        try:
            write_receipt(missing_attachment_receipt, [completed_path], [])
            raise AssertionError("completed profile receipt omitted attachment")
        except LiteratureError:
            if missing_attachment_receipt.exists():
                raise AssertionError("failed receipt left partial output")
        profile_receipt = root / "profile-receipt.json"
        write_receipt(profile_receipt, [completed_path],
                      [primary, calibration, trace, rule, extraction_script])
        validate_receipt(profile_receipt, True)


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command")
    validate = commands.add_parser("validate")
    validate.add_argument("references", nargs="+", type=Path)
    csv_parser = commands.add_parser("csv")
    csv_parser.add_argument("path", type=Path)
    receipt = commands.add_parser("receipt")
    receipt.add_argument("--output", required=True, type=Path)
    receipt.add_argument("--reference", action="append", default=[], type=Path)
    receipt.add_argument("--attachment", action="append", default=[], type=Path)
    receipt_validate = commands.add_parser("receipt-validate")
    receipt_validate.add_argument("path", type=Path)
    receipt_validate.add_argument("--require-complete", action="store_true")
    commands.add_parser("self-test")
    arguments = parser.parse_args()
    if arguments.command is None:
        parser.error("a command is required")
    try:
        if arguments.command == "validate":
            for path in arguments.references:
                validate_reference(path)
        elif arguments.command == "csv":
            print(json.dumps(import_csv(arguments.path), separators=(",", ":")))
        elif arguments.command == "receipt":
            if not arguments.reference:
                raise LiteratureError("at least one --reference is required")
            write_receipt(arguments.output, arguments.reference,
                          arguments.attachment)
        elif arguments.command == "receipt-validate":
            validate_receipt(arguments.path, arguments.require_complete)
        else:
            self_test()
    except (LiteratureError, OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
