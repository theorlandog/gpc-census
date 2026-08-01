#!/usr/bin/env python3
"""Verify the rank-9 and rank-10 constraint transcription against the thesis.

The script uses only the Python standard library and Poppler's ``pdftotext``.
It reads Appendix A directly from Murat Altunbulak's 2008 dissertation,
parses every displayed inequality, and compares the rows in order with the
machine-readable transcription.  It can also regenerate the compact page map
provided in the Paper 1 supplementary archive.

The PDF is not redistributed.  Download it from the Bilkent repository and run

    python3 scripts/verify_constraint_source_standalone.py thesis.pdf

Use ``--write-map PATH`` to regenerate the page-level provenance artifact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROWS = ROOT / "results/data/constraints/ak_thesis_rank9_10_constraints.json"
DEFAULT_NORMALIZED = ROOT / "src/gpc_census/data/constraints.json"
EXPECTED_PDF_SHA256 = "42b29e4b27cbf3d38d66d08661be2dd1fef301d97f19a628cda27d31eea0e31d"
SOURCE_URL = "https://repository.bilkent.edu.tr/items/dc2cf80a-82c0-4902-93a7-e83189a0d431"
BITSTREAM_URL = "https://repository.bilkent.edu.tr/bitstreams/4c5aaf14-4dbc-438e-8ef9-47e051dece53/download"
SYSTEMS = ("(3,9)", "(4,9)", "(3,10)", "(4,10)", "(5,10)")


def _parse_inequality(text: str, d: int) -> tuple[list[int], int]:
    lhs, rhs_text = text.split("≤", 1)
    rhs_match = re.match(r"\s*(\d+)", rhs_text)
    if rhs_match is None:
        raise ValueError(f"missing right-hand side: {text!r}")
    normalized = lhs.replace("−", "-").replace("–", "-").replace(" ", "")
    terms = re.findall(r"([+-]?)(\d*)λ(\d+)", normalized)
    if len(terms) != normalized.count("λ"):
        raise ValueError(f"could not parse every lambda term: {text!r}")
    coefficients = [0] * d
    for sign, magnitude, index_text in terms:
        coefficient = int(magnitude or "1")
        if sign == "-":
            coefficient = -coefficient
        index = int(index_text) - 1
        if not 0 <= index < d:
            raise ValueError(f"orbital index outside 1..{d}: {text!r}")
        coefficients[index] += coefficient
    return coefficients, int(rhs_match.group(1))


def _extract_rows(pdf: Path) -> dict[str, list[dict]]:
    with tempfile.TemporaryDirectory(prefix="gpc-source-") as directory:
        text_path = Path(directory) / "appendix-a.txt"
        subprocess.run(
            [
                "pdftotext",
                "-f",
                "93",
                "-l",
                "118",
                "-layout",
                str(pdf),
                str(text_path),
            ],
            check=True,
        )
        pages = text_path.read_text(encoding="utf-8").split("\f")

    result: dict[str, list[dict]] = {system: [] for system in SYSTEMS}
    current: str | None = None
    row_on_page: dict[tuple[str, int], int] = {}
    for pdf_page, page in enumerate(pages, start=93):
        for raw_line in page.splitlines():
            line = raw_line.strip()
            compact = line.replace(" ", "")
            if compact.startswith("A.2ν-representability"):
                current = None
                break
            marker = re.search(r"Forthesystem∧(\d+)H(\d+)", compact)
            if marker:
                current = f"({marker.group(1)},{marker.group(2)})"
                continue
            if current not in result or "≤" not in line or "λ" not in line:
                continue
            d = int(current.split(",")[1][:-1])
            coefficients, rhs = _parse_inequality(line, d)
            key = (current, pdf_page)
            row_on_page[key] = row_on_page.get(key, 0) + 1
            result[current].append(
                {
                    "coeffs": coefficients,
                    "rhs": rhs,
                    "pdf_page": pdf_page,
                    "thesis_page": pdf_page - 9,
                    "row_on_page": row_on_page[key],
                }
            )
    return result


def _row_pair(row: dict) -> tuple[tuple[int, ...], int]:
    return tuple(int(value) for value in row["coeffs"]), int(row["rhs"])


def _verify(
    extracted: dict[str, list[dict]],
    source_rows: dict[str, list[dict]],
    normalized: dict[str, dict],
    pdf_sha256: str,
    pdftotext_version: str,
) -> dict:
    mapping: dict[str, dict] = {}
    for system in SYSTEMS:
        printed = extracted[system]
        expected = source_rows[system]
        normalized_rows = normalized[system.strip("()")]["inequalities"]
        if [_row_pair(row) for row in expected] != [
            _row_pair(row) for row in normalized_rows
        ]:
            raise AssertionError(f"raw and normalized tables disagree for {system}")

        mapped_rows: list[dict] = []
        if system == "(3,9)":
            if len(printed) != 51 or len(expected) != 52:
                raise AssertionError("expected the documented 51/52 (3,9) page-break case")
            if [_row_pair(row) for row in printed] != [
                _row_pair(row) for row in expected[:51]
            ]:
                raise AssertionError("the 51 printed (3,9) rows disagree with the table")
            restricted = list(source_rows["(3,10)"][81]["coeffs"])
            if restricted[-1] != 0:
                raise AssertionError("(3,10) row 82 does not restrict to nine orbitals")
            if _row_pair(expected[51]) != (
                tuple(restricted[:-1]),
                source_rows["(3,10)"][81]["rhs"],
            ):
                raise AssertionError("the recovered (3,9) row is not the row-82 restriction")
            source_310 = extracted["(3,10)"][81]
            for index, row in enumerate(printed, start=1):
                mapped_rows.append({"row": index, "source_kind": "printed", **row})
            mapped_rows.append(
                {
                    "row": 52,
                    "source_kind": "restriction",
                    "source_system": "(3,10)",
                    "source_row": 82,
                    "pdf_page": source_310["pdf_page"],
                    "thesis_page": source_310["thesis_page"],
                    "row_on_page": source_310["row_on_page"],
                    "coeffs": expected[51]["coeffs"],
                    "rhs": expected[51]["rhs"],
                }
            )
        else:
            if [_row_pair(row) for row in printed] != [
                _row_pair(row) for row in expected
            ]:
                raise AssertionError(f"printed rows disagree with the table for {system}")
            mapped_rows = [
                {"row": index, "source_kind": "printed", **row}
                for index, row in enumerate(printed, start=1)
            ]
        mapping[system] = {
            "table_rows": len(expected),
            "printed_rows": len(printed),
            "rows": mapped_rows,
        }

    return {
        "schema_version": 1,
        "source": {
            "author": "Murat Altunbulak",
            "title": "The Pauli Principle, Representation Theory, and Geometry of Flag Varieties",
            "year": 2008,
            "repository_url": SOURCE_URL,
            "bitstream_url": BITSTREAM_URL,
            "pdf_sha256": pdf_sha256,
            "pdftotext_version": pdftotext_version,
            "appendix_pdf_pages": [93, 118],
            "method": "Fresh pdftotext -layout extraction followed by ordered coefficient comparison",
            "note": "Appendix A states that (3,9) has 52 independent inequalities but prints 51. Supplementary table row 52 equals printed (3,10) row 82 restricted to lambda_10=0; the identity does not establish historical intent.",
        },
        "systems": mapping,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path, help="Altunbulak dissertation PDF")
    parser.add_argument("--source-rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--normalized", type=Path, default=DEFAULT_NORMALIZED)
    parser.add_argument("--write-map", type=Path)
    parser.add_argument(
        "--expected-map",
        type=Path,
        help="require the regenerated page map to equal this JSON object",
    )
    parser.add_argument(
        "--allow-different-pdf-hash",
        action="store_true",
        help="parse another byte-level copy of the dissertation",
    )
    args = parser.parse_args()

    digest = hashlib.sha256(args.pdf.read_bytes()).hexdigest()
    if digest != EXPECTED_PDF_SHA256 and not args.allow_different_pdf_hash:
        raise SystemExit(
            f"source PDF hash differs: expected {EXPECTED_PDF_SHA256}, got {digest}"
        )
    source_rows = json.loads(args.source_rows.read_text(encoding="utf-8"))
    normalized = json.loads(args.normalized.read_text(encoding="utf-8"))
    version_process = subprocess.run(
        ["pdftotext", "-v"], capture_output=True, text=True, check=True
    )
    version_line = (version_process.stderr or version_process.stdout).splitlines()[0]
    mapping = _verify(
        _extract_rows(args.pdf),
        source_rows,
        normalized,
        digest,
        version_line,
    )
    if args.expected_map:
        expected = json.loads(args.expected_map.read_text(encoding="utf-8"))
        if mapping != expected:
            raise AssertionError(
                f"regenerated source map differs from {args.expected_map}"
            )
    if args.write_map:
        args.write_map.parent.mkdir(parents=True, exist_ok=True)
        args.write_map.write_text(
            json.dumps(mapping, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    counts = ", ".join(
        f"{system}: {mapping['systems'][system]['table_rows']}"
        for system in SYSTEMS
    )
    print(f"PASS: Appendix A constraint transcription verified row by row ({counts})")
    print("PASS: all 60 (4,9) and all 379 rank-10 rows match the printed tables")
    print(
        "PASS: supplied (3,9) row 52 equals printed (3,10) row 82 at "
        "lambda_10=0; this identity does not establish historical intent"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
