#!/usr/bin/env python3
"""Check every cross-link in the self-contained Paper 1 supplement.

This standard-library-only verifier binds state records, vertex rows,
no-design certificates, rational H-tables, and the primary-source page map.
The object-level symbolic and polyhedral verifiers perform the mathematical
checks; this script prevents correct objects from being relabeled, omitted, or
combined with data from another dataset.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import tomllib
from fractions import Fraction
from pathlib import Path


SYSTEM_COUNTS = {
    "(3,6)": 4,
    "(3,7)": 10,
    "(3,8)": 38,
    "(4,8)": 22,
    "(3,9)": 58,
    "(4,9)": 103,
    "(3,10)": 113,
    "(4,10)": 159,
    "(5,10)": 292,
}
SOURCE_SYSTEMS = ("(3,9)", "(4,9)", "(3,10)", "(4,10)", "(5,10)")
SOURCE_COUNTS = {
    "(3,9)": 52,
    "(4,9)": 60,
    "(3,10)": 93,
    "(4,10)": 125,
    "(5,10)": 161,
}
VERTEX_PATHS = {
    f"data/vertices/vertices_{system.strip('()').replace(',', '_')}.json"
    for system in SYSTEM_COUNTS
}
SCRIPT_PATHS = {
    "scripts/verify_constraint_source_standalone.py",
    "scripts/verify_interference_certificates_standalone.py",
    "scripts/verify_paper1_bundle_standalone.py",
    "scripts/verify_states_standalone.py",
    "scripts/verify_vertex_exhaustion_standalone.py",
}
EXPECTED_PATHS = frozenset(
    {
        ".python-version",
        "LICENSE",
        "README.md",
        "SHA256SUMS",
        "SOURCE.json",
        "data/constraint_source_map.json",
        "data/constraints.json",
        "data/interference_certificates.json",
        "data/source_rows.json",
        "data/states.jsonl",
        "pyproject.toml",
        "uv.lock",
    }
    | VERTEX_PATHS
    | SCRIPT_PATHS
)
STATE_KEYS = {
    "system",
    "index",
    "classified",
    "integer_form",
    "denominator",
    "closed_form",
}
CLOSED_FORM_KEYS = {"den", "pretty", "support_dets", "weights"}
SOURCE_IDENTITY_KEYS = {
    "schema_version",
    "artifact",
    "artifact_version",
    "source_repository",
    "source_commit",
    "source_tree_dirty",
    "python_version",
    "sympy_version",
    "pdftotext_version",
}
SOURCE_METADATA_KEYS = {
    "appendix_pdf_pages",
    "author",
    "bitstream_url",
    "method",
    "note",
    "pdf_sha256",
    "pdftotext_version",
    "repository_url",
    "title",
    "year",
}
PRINTED_MAP_KEYS = {
    "coeffs",
    "pdf_page",
    "rhs",
    "row",
    "row_on_page",
    "source_kind",
    "thesis_page",
}
RESTRICTION_MAP_KEYS = PRINTED_MAP_KEYS | {"source_row", "source_system"}
EXPECTED_PDF_SHA256 = (
    "42b29e4b27cbf3d38d66d08661be2dd1fef301d97f19a628cda27d31eea0e31d"
)


def _fractions(values) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in values)


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _archive_files(root: Path) -> set[str]:
    symlinks = [path for path in root.rglob("*") if path.is_symlink()]
    if symlinks:
        raise AssertionError(
            f"supplement must not contain symbolic links: "
            f"{[str(path.relative_to(root)) for path in symlinks]}"
        )
    return {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file()
    }


def _verify_checksums(root: Path) -> None:
    manifest = root / "SHA256SUMS"
    if not manifest.is_file():
        raise AssertionError("mandatory SHA256SUMS manifest is missing")
    archive_files = _archive_files(root)
    listed: set[str] = set()
    listed_order: list[str] = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, separator, relative = line.partition("  ")
        if separator != "  " or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise AssertionError(f"malformed checksum entry: {line!r}")
        if relative == "SHA256SUMS" or relative not in archive_files:
            raise AssertionError(f"invalid checksum path: {relative}")
        path = root / relative
        if relative in listed:
            raise AssertionError(f"duplicate checksum entry: {relative}")
        listed.add(relative)
        listed_order.append(relative)
        if not path.is_file():
            raise AssertionError(f"checksummed file is missing: {relative}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != digest:
            raise AssertionError(f"checksum mismatch: {relative}")
    if listed_order != sorted(listed_order):
        raise AssertionError("checksum entries are not in canonical path order")
    actual_files = archive_files - {"SHA256SUMS"}
    if listed != actual_files:
        raise AssertionError(
            f"checksum coverage mismatch: unlisted={sorted(actual_files - listed)}, "
            f"missing={sorted(listed - actual_files)}"
        )


def _verify_expected_paths(root: Path) -> None:
    actual = _archive_files(root)
    if actual != EXPECTED_PATHS:
        raise AssertionError(
            f"archive path set mismatch: extra={sorted(actual - EXPECTED_PATHS)}, "
            f"missing={sorted(EXPECTED_PATHS - actual)}"
        )


def _verify_source_identity(root: Path, source_map: dict) -> None:
    identity = _read_json(root / "SOURCE.json")
    if set(identity) != SOURCE_IDENTITY_KEYS:
        raise AssertionError("SOURCE.json field set drifted")
    expected_scalars = {
        "schema_version": 1,
        "artifact": "gpc-census-paper1-supplement",
        "artifact_version": "1.0.0",
        "source_repository": "https://github.com/theorlandog/gpc-census",
    }
    for key, expected in expected_scalars.items():
        if identity[key] != expected:
            raise AssertionError(f"SOURCE.json {key} drifted")
    if not re.fullmatch(r"[0-9a-f]{40}", identity["source_commit"]):
        raise AssertionError("SOURCE.json source_commit is not a full Git object id")
    if type(identity["source_tree_dirty"]) is not bool:
        raise AssertionError("SOURCE.json source_tree_dirty must be Boolean")
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", identity["python_version"]):
        raise AssertionError("SOURCE.json python_version is not a full version")
    if (root / ".python-version").read_text(encoding="utf-8").strip() != "3.12.13":
        raise AssertionError("supplement Python version must be 3.12.13")
    with (root / "pyproject.toml").open("rb") as stream:
        project = tomllib.load(stream)["project"]
    if project.get("name") != identity["artifact"]:
        raise AssertionError("SOURCE.json artifact disagrees with pyproject.toml")
    if project.get("version") != identity["artifact_version"]:
        raise AssertionError("SOURCE.json artifact version disagrees with pyproject.toml")
    sympy_pins = [
        value.removeprefix("sympy==")
        for value in project.get("dependencies", [])
        if value.startswith("sympy==")
    ]
    if sympy_pins != [identity["sympy_version"]]:
        raise AssertionError("SOURCE.json SymPy version disagrees with pyproject.toml")
    if identity["pdftotext_version"] != source_map["source"]["pdftotext_version"]:
        raise AssertionError("SOURCE.json Poppler version disagrees with the source map")


def _verify_printed_page_metadata(system: str, rows: list[dict], first: int, last: int) -> None:
    previous_page = first
    counts: dict[int, int] = {}
    for index, row in enumerate(rows, start=1):
        page = row.get("pdf_page")
        thesis_page = row.get("thesis_page")
        row_on_page = row.get("row_on_page")
        if type(page) is not int or not first <= page <= last:
            raise AssertionError(f"source-map PDF page invalid for {system} row {index}")
        if page < previous_page:
            raise AssertionError(f"source-map pages are not ordered for {system}")
        if thesis_page != page - 9:
            raise AssertionError(f"source-map thesis page invalid for {system} row {index}")
        counts[page] = counts.get(page, 0) + 1
        if row_on_page != counts[page]:
            raise AssertionError(f"source-map row-on-page invalid for {system} row {index}")
        previous_page = page


def _verify_source_concordance(constraints: dict, source_map: dict, source_rows: dict) -> None:
    if set(source_map) != {"schema_version", "source", "systems"}:
        raise AssertionError("source-map top-level fields drifted")
    if source_map["schema_version"] != 1:
        raise AssertionError("source-map schema version drifted")
    if set(source_map["source"]) != SOURCE_METADATA_KEYS:
        raise AssertionError("source-map metadata fields drifted")
    metadata = source_map["source"]
    if metadata["pdf_sha256"] != EXPECTED_PDF_SHA256:
        raise AssertionError("primary-source PDF identity drifted")
    if metadata["appendix_pdf_pages"] != [93, 118]:
        raise AssertionError("source-map Appendix A page range drifted")
    if not isinstance(metadata["pdftotext_version"], str) or not metadata["pdftotext_version"]:
        raise AssertionError("source-map pdftotext version is missing")
    if set(source_map["systems"]) != set(SOURCE_SYSTEMS):
        raise AssertionError("source-map system set drifted")
    if set(source_rows) != set(SOURCE_SYSTEMS):
        raise AssertionError("source-row system set drifted")

    first_page, last_page = metadata["appendix_pdf_pages"]
    for system in SOURCE_SYSTEMS:
        mapped = source_map["systems"][system]
        if set(mapped) != {"printed_rows", "rows", "table_rows"}:
            raise AssertionError(f"source-map system fields drifted for {system}")
        raw_rows = source_rows[system]
        table_rows = constraints[system.strip("()")]["inequalities"]
        expected_count = SOURCE_COUNTS[system]
        if not (
            mapped["table_rows"] == expected_count
            and len(mapped["rows"]) == expected_count
            and len(raw_rows) == expected_count
            and len(table_rows) == expected_count
        ):
            raise AssertionError(f"source concordance count disagrees for {system}")
        expected_printed = expected_count - 1 if system == "(3,9)" else expected_count
        if mapped["printed_rows"] != expected_printed:
            raise AssertionError(f"printed-row count disagrees for {system}")

        printed_mappings = mapped["rows"][:expected_printed]
        _verify_printed_page_metadata(system, printed_mappings, first_page, last_page)
        for index, (mapping, raw, normalized) in enumerate(
            zip(mapped["rows"], raw_rows, table_rows), start=1
        ):
            is_restriction = system == "(3,9)" and index == expected_count
            expected_keys = RESTRICTION_MAP_KEYS if is_restriction else PRINTED_MAP_KEYS
            if set(mapping) != expected_keys:
                raise AssertionError(f"source-map row fields drifted for {system} row {index}")
            if set(raw) != {"coeffs", "ineq", "rhs"} or not isinstance(raw["ineq"], str):
                raise AssertionError(f"source-row fields drifted for {system} row {index}")
            if mapping["row"] != index:
                raise AssertionError(f"source-map order disagrees for {system} row {index}")
            expected_kind = "restriction" if is_restriction else "printed"
            if mapping["source_kind"] != expected_kind:
                raise AssertionError(f"source kind disagrees for {system} row {index}")
            row_pair = (raw["coeffs"], int(raw["rhs"]))
            if row_pair != (normalized["coeffs"], int(normalized["rhs"])):
                raise AssertionError(
                    f"source rows disagree with normalized table for {system} row {index}"
                )
            if row_pair != (mapping["coeffs"], int(mapping["rhs"])):
                raise AssertionError(
                    f"source rows disagree with page map for {system} row {index}"
                )

    recovered = source_map["systems"]["(3,9)"]["rows"][51]
    source_310 = source_map["systems"]["(3,10)"]["rows"][81]
    raw_310 = source_rows["(3,10)"][81]
    if not (
        recovered["source_system"] == "(3,10)"
        and recovered["source_row"] == 82
        and recovered["pdf_page"] == source_310["pdf_page"]
        and recovered["thesis_page"] == source_310["thesis_page"]
        and recovered["row_on_page"] == source_310["row_on_page"]
        and raw_310["coeffs"][-1] == 0
        and recovered["coeffs"] == raw_310["coeffs"][:-1]
        and int(recovered["rhs"]) == int(raw_310["rhs"])
    ):
        raise AssertionError("the documented (3,9) restriction row drifted")


def verify(root: Path) -> None:
    _verify_checksums(root)
    _verify_expected_paths(root)
    data = root / "data"
    vertices_dir = data / "vertices"
    state_rows = [
        json.loads(line)
        for line in (data / "states.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(state_rows) != 799:
        raise AssertionError(f"expected 799 state records, got {len(state_rows)}")

    vertices: dict[tuple[str, int], dict] = {}
    for system, expected_count in SYSTEM_COUNTS.items():
        n, d = system.strip("()").split(",")
        rows = _read_json(vertices_dir / f"vertices_{n}_{d}.json")
        if len(rows) != expected_count:
            raise AssertionError(f"{system} vertex count is {len(rows)}, not {expected_count}")
        for index, row in enumerate(rows):
            vertices[(system, index)] = row

    states: dict[tuple[str, int], dict] = {}
    class_counts = {"DESIGN-INT": 0, "DESIGN-REAL": 0, "INTERFERENCE": 0}
    for record in state_rows:
        if set(record) != STATE_KEYS:
            raise AssertionError("state projection contains nonexact or missing fields")
        closed_form_keys = set(record["closed_form"])
        if closed_form_keys not in (CLOSED_FORM_KEYS, CLOSED_FORM_KEYS | {"family"}):
            raise AssertionError("closed-form state field set drifted")
        key = (record["system"], int(record["index"]))
        if key in states:
            raise AssertionError(f"duplicate state key: {key}")
        if key not in vertices:
            raise AssertionError(f"state key is not a vertex: {key}")
        vertex = vertices[key]
        n, d = (int(value) for value in record["system"].strip("()").split(","))
        denominator = record["denominator"]
        integer_form = record["integer_form"]
        if (
            type(denominator) is not int
            or denominator <= 0
            or len(integer_form) != d
            or any(type(value) is not int for value in integer_form)
            or integer_form != sorted(integer_form, reverse=True)
            or any(value < 0 or value > denominator for value in integer_form)
            or sum(integer_form) != n * denominator
            or math.gcd(denominator, *integer_form) != 1
        ):
            raise AssertionError(f"state canonical integer spectrum invalid at {key}")
        if record["integer_form"] != vertex["integer_form"]:
            raise AssertionError(f"state integer form disagrees at {key}")
        if int(record["denominator"]) != int(vertex["denominator"]):
            raise AssertionError(f"state denominator disagrees at {key}")
        if _fractions(vertex["spectrum"]) != tuple(
            Fraction(value, int(record["denominator"]))
            for value in record["integer_form"]
        ):
            raise AssertionError(f"vertex spectrum disagrees at {key}")
        classification = record["classified"]
        if classification not in class_counts:
            raise AssertionError(f"unknown class at {key}: {classification}")
        class_counts[classification] += 1
        states[key] = record
    if set(states) != set(vertices):
        raise AssertionError("state and vertex key sets differ")
    if class_counts != {"DESIGN-INT": 631, "DESIGN-REAL": 12, "INTERFERENCE": 156}:
        raise AssertionError(f"classification counts drifted: {class_counts}")

    certificates = _read_json(data / "interference_certificates.json")
    if certificates.get("record_count") != 12 or len(certificates.get("records", [])) != 12:
        raise AssertionError("expected 12 no-design certificates")
    seen_certificates: set[tuple[str, int]] = set()
    for certificate in certificates["records"]:
        system = f"({certificate['system'][0]},{certificate['system'][1]})"
        key = (system, int(certificate["vertex_index"]))
        if key in seen_certificates:
            raise AssertionError(f"duplicate certificate key: {key}")
        if key not in vertices or key not in states:
            raise AssertionError(f"certificate key is not in the atlas: {key}")
        vertex = vertices[key]
        state = states[key]
        if state["classified"] != "INTERFERENCE":
            raise AssertionError(f"certificate target is not INTERFERENCE: {key}")
        if certificate["integer_form"] != vertex["integer_form"]:
            raise AssertionError(f"certificate integer form disagrees at {key}")
        if int(certificate["denominator"]) != int(vertex["denominator"]):
            raise AssertionError(f"certificate denominator disagrees at {key}")
        if _fractions(certificate["spectrum"]) != _fractions(vertex["spectrum"]):
            raise AssertionError(f"certificate spectrum disagrees at {key}")
        seen_certificates.add(key)
    expected_certificates = {
        key
        for key, state in states.items()
        if state["classified"] == "INTERFERENCE"
        and (key[0] == "(3,8)" or key == ("(4,9)", 65))
    }
    if seen_certificates != expected_certificates:
        raise AssertionError("no-design certificate target set drifted")

    constraints = _read_json(data / "constraints.json")
    source_map = _read_json(data / "constraint_source_map.json")
    source_rows = _read_json(data / "source_rows.json")
    _verify_source_concordance(constraints, source_map, source_rows)
    _verify_source_identity(root, source_map)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="supplement root containing the exact checksummed archive members",
    )
    args = parser.parse_args()
    verify(args.root)
    print(
        "PASS: 799 state-to-vertex links, 12 certificate targets, five "
        "source maps, and checksum coverage agree"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
