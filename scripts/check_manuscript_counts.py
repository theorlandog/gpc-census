#!/usr/bin/env python3
"""Check every structured census claim retained in Paper 1.

The checker reads the manuscript itself, the exact state atlas, the nine vertex
files, and the exact global no-design artifact. It fails if the Markdown table,
headline totals, evidence scope, or atlas-to-vertex bijection drifts.
"""

from __future__ import annotations

import json
import pathlib
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction


ROOT = pathlib.Path(__file__).resolve().parents[1]
STATES = ROOT / "results/data/states.jsonl"
VERTICES = ROOT / "results/data/vertices"
CERTS = ROOT / "results/data/interference_certificates.json"
SOURCE_MAP = ROOT / "results/data/constraints/altunbulak_appendix_a_map.json"
CONSTRAINTS = ROOT / "src/gpc_census/data/constraints.json"
MANUSCRIPT = ROOT / "results/report/main.md"

EXPECTED = {
    "(3,6)": dict(total=4, int_=4, real=0, interference=0),
    "(3,7)": dict(total=10, int_=10, real=0, interference=0),
    "(3,8)": dict(total=38, int_=26, real=1, interference=11),
    "(4,8)": dict(total=22, int_=22, real=0, interference=0),
    "(3,9)": dict(total=58, int_=37, real=1, interference=20),
    "(4,9)": dict(total=103, int_=85, real=2, interference=16),
    "(3,10)": dict(total=113, int_=69, real=2, interference=42),
    "(4,10)": dict(total=159, int_=132, real=2, interference=25),
    "(5,10)": dict(total=292, int_=246, real=4, interference=42),
}
EXPECTED_TOTALS = dict(total=799, int_=631, real=12, interference=156)
EXPECTED_TABLE = {
    "(3,6)": dict(total=4, positive=4, proved_negative=0, solver_negative=0),
    "(3,7)": dict(total=10, positive=10, proved_negative=0, solver_negative=0),
    "(3,8)": dict(total=38, positive=27, proved_negative=11, solver_negative=0),
    "(4,8)": dict(total=22, positive=22, proved_negative=0, solver_negative=0),
    "(3,9)": dict(total=58, positive=38, proved_negative=0, solver_negative=20),
    "(4,9)": dict(total=103, positive=87, proved_negative=1, solver_negative=15),
    "(3,10)": dict(total=113, positive=71, proved_negative=0, solver_negative=42),
    "(4,10)": dict(total=159, positive=134, proved_negative=0, solver_negative=25),
    "(5,10)": dict(total=292, positive=250, proved_negative=0, solver_negative=42),
}
EXPECTED_RANK10 = 564
EXPECTED_DESCENT = dict(padding=235, core=161, primitive=403)
EXPECTED_CERTS = {
    "record_count": 12,
    "farkas_vectors": 192,
    "expanded_clauses": 5141,
    "proof_nodes": 226,
    "maximum_proof_depth": 8,
}


def vertex_path(system: str) -> pathlib.Path:
    n, d = system.strip("()").split(",")
    return VERTICES / f"vertices_{n}_{d}.json"


def parse_table(text: str) -> dict[str, dict[str, int]]:
    """Parse the five numeric columns of the census Markdown table."""
    found: dict[str, dict[str, int]] = {}
    pattern = re.compile(
        r"^\| \$P\^\{\\mathrm\{tab\}\}_\{(\d+),(\d+)\}\$ \| (\d+) \| [^|]+ "
        r"\| (\d+) \| (\d+) \| (\d+) \|$"
    )
    for line in text.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        n, d, total, positive, proved_negative, solver_negative = map(
            int, match.groups()
        )
        found[f"({n},{d})"] = dict(
            total=total,
            positive=positive,
            proved_negative=proved_negative,
            solver_negative=solver_negative,
        )
    return found


def require(text: str, needle: str, errors: list[str]) -> None:
    if needle not in text:
        errors.append(f"manuscript is missing required scoped claim: {needle!r}")


def main() -> int:
    manuscript = MANUSCRIPT.read_text()
    records = [json.loads(line) for line in STATES.read_text().splitlines()]
    errors: list[str] = []

    by_system: dict[str, Counter[str]] = defaultdict(Counter)
    for record in records:
        by_system[record["system"]][record["classified"]] += 1

    data_counts: dict[str, dict[str, int]] = {}
    for system in EXPECTED:
        counts = by_system[system]
        data_counts[system] = dict(
            total=sum(counts.values()),
            int_=counts["DESIGN-INT"],
            real=counts["DESIGN-REAL"],
            interference=counts["INTERFERENCE"],
        )
    if data_counts != EXPECTED:
        errors.append(f"state-atlas system counts drifted: {data_counts}")

    table_counts = parse_table(manuscript)
    if table_counts != EXPECTED_TABLE:
        errors.append(f"manuscript census table disagrees: {table_counts}")

    totals = dict(
        total=len(records),
        int_=sum(row["int_"] for row in data_counts.values()),
        real=sum(row["real"] for row in data_counts.values()),
        interference=sum(row["interference"] for row in data_counts.values()),
    )
    if totals != EXPECTED_TOTALS:
        errors.append(f"state-atlas totals drifted: {totals}")

    # Each record must be exact and match the same indexed vertex exactly.
    seen: set[tuple[str, int]] = set()
    for record in records:
        key = (record["system"], record["index"])
        if key in seen:
            errors.append(f"duplicate state key: {key}")
        seen.add(key)
        if record.get("status") != "OK" or not record.get("closed_form"):
            errors.append(f"state lacks an exact closed form: {key}")

    vertex_keys: set[tuple[str, int]] = set()
    for system in EXPECTED:
        vertices = json.loads(vertex_path(system).read_text())
        if len(vertices) != EXPECTED[system]["total"]:
            errors.append(f"vertex file count drifted for {system}")
        state_rows = {r["index"]: r for r in records if r["system"] == system}
        for index, vertex in enumerate(vertices):
            key = (system, index)
            vertex_keys.add(key)
            record = state_rows.get(index)
            if record is None:
                errors.append(f"missing state record for vertex {key}")
                continue
            if record["integer_form"] != vertex["integer_form"]:
                errors.append(f"integer form disagrees at {key}")
            if record["denominator"] != vertex["denominator"]:
                errors.append(f"denominator disagrees at {key}")
    if seen != vertex_keys:
        errors.append("state records and vertex files are not in exact bijection")

    rank10 = sum(
        1 for record in records if record["system"] in {"(3,10)", "(4,10)", "(5,10)"}
    )
    if rank10 != EXPECTED_RANK10:
        errors.append(f"ten-orbital count drifted: {rank10}")

    # The novelty paragraph uses a disjoint, reproducible descent convention:
    # test empty-orbital padding first, then a frozen occupied orbital.
    vertex_spectra: dict[tuple[int, int], set[tuple[Fraction, ...]]] = {}
    all_spectra: list[tuple[int, int, tuple[Fraction, ...]]] = []
    for system in EXPECTED:
        n, d = map(int, system.strip("()").split(","))
        rows = json.loads(vertex_path(system).read_text())
        spectra = {
            tuple(Fraction(value) for value in row["spectrum"]) for row in rows
        }
        vertex_spectra[n, d] = spectra
        all_spectra.extend((n, d, spectrum) for spectrum in spectra)
    descent = Counter()
    for n, d, spectrum in all_spectra:
        if (
            spectrum[-1] == 0
            and (n, d - 1) in vertex_spectra
            and spectrum[:-1] in vertex_spectra[n, d - 1]
        ):
            descent["padding"] += 1
        elif (
            spectrum[0] == 1
            and (n - 1, d - 1) in vertex_spectra
            and spectrum[1:] in vertex_spectra[n - 1, d - 1]
        ):
            descent["core"] += 1
        else:
            descent["primitive"] += 1
    if dict(descent) != EXPECTED_DESCENT:
        errors.append(f"descent counts drifted: {dict(descent)}")

    certs = json.loads(CERTS.read_text())
    if certs.get("record_count") != EXPECTED_CERTS["record_count"]:
        errors.append("global no-design record count drifted")
    for key in (
        "farkas_vectors",
        "expanded_clauses",
        "proof_nodes",
        "maximum_proof_depth",
    ):
        if certs.get("summary", {}).get(key) != EXPECTED_CERTS[key]:
            errors.append(f"global no-design certificate {key} drifted")

    cert_targets = {
        (f"({row['system'][0]},{row['system'][1]})", row["vertex_index"])
        for row in certs.get("records", [])
    }
    expected_targets = {
        (record["system"], record["index"])
        for record in records
        if record["classified"] == "INTERFERENCE"
        and (
            record["system"] == "(3,8)"
            or (record["system"] == "(4,9)" and record["index"] == 65)
        )
    }
    if cert_targets != expected_targets:
        errors.append("global no-design certificate target set drifted")

    source_map = json.loads(SOURCE_MAP.read_text())
    constraints = json.loads(CONSTRAINTS.read_text())
    if source_map.get("source", {}).get("pdf_sha256") != (
        "42b29e4b27cbf3d38d66d08661be2dd1fef301d97f19a628cda27d31eea0e31d"
    ):
        errors.append("constraint source identity drifted")
    for system in ("(3,9)", "(4,9)", "(3,10)", "(4,10)", "(5,10)"):
        mapped = source_map["systems"][system]
        rows = constraints[system.strip("()")]["inequalities"]
        if len(mapped["rows"]) != len(rows):
            errors.append(f"constraint source-map count drifted for {system}")
            continue
        for mapping, row in zip(mapped["rows"], rows):
            if mapping["coeffs"] != row["coeffs"] or mapping["rhs"] != row["rhs"]:
                errors.append(f"constraint source-map row drifted for {system}")
                break
    if sum(
        source_map["systems"][system]["printed_rows"]
        for system in ("(3,9)", "(4,9)", "(3,10)", "(4,10)", "(5,10)")
    ) != 490:
        errors.append("printed source-row total drifted")

    equation_tags = re.findall(r"\\tag\{([^}]+)\}", manuscript)
    expected_equation_tags = [
        *(str(index) for index in range(1, 18)),
        "A1",
        "A2",
        "B1",
        "C1",
    ]
    if equation_tags != expected_equation_tags:
        errors.append(f"displayed-equation numbering drifted: {equation_tags}")
    if manuscript.count("<!-- Alt text:") != 2:
        errors.append("both manuscript tables must carry source alt text")

    # These are the claim boundaries a referee must see in the prose.
    for claim in (
        "Symbolic 1-RDM reconstruction verifies the preimage in all 799 records",
        "The 144 solver classifications record exploratory CP-SAT and CBC outcomes",
        "**Proposition 2** (Existence of a finite certificate)",
        "**Proposition 3** (Vertex exhaustion of the tabulated H-polytopes)",
        "**Corollary 2** (Four fermions in nine orbitals)",
        "*If every inequality",
        "1\\ge\\lambda_1\\ge\\cdots\\ge\\lambda_{10}\\ge0",
        "235 records descend by padding, 161 by a frozen-core lift, and 403",
        "all 60 rows at $(4,9)$ and all 379 rows at $d=10$",
        "All 490 printed rows match",
        "does not establish that it was the source's intended omitted row",
        "This notation does not assume that a tabulated description equals $\\Pi_{N,d}$",
        "z=\\sqrt{14}\\eta+2",
    ):
        require(manuscript, claim, errors)

    for unbecoming in (
        "claims retained in Paper 1",
        "production pipeline",
        "shipped weights",
        "does not trust the discovery solver",
        "PASS logs",
        "locked environment",
        "current release supplies",
        "The scope is deliberately narrow",
        "Every retained claim",
        "coherent-cancellation",
        "two exchange channels cancel",
        "normalized block",
        "GPC rows",
        "non-Pauli rows",
        "generalized Pauli half-spaces",
        "| $\\Pi_",
    ):
        if unbecoming in manuscript:
            errors.append(f"internal or misleading wording returned: {unbecoming}")

    # Paper 2 and generator material must not leak back into Paper 1.
    for forbidden in (
        "TD-2RDM",
        "rank-11",
        "rank-12",
        "Hamiltonian control",
        "holonomy Galois",
    ):
        if forbidden in manuscript:
            errors.append(f"out-of-scope Paper 2 material returned: {forbidden}")

    if errors:
        print("manuscript consistency FAILED")
        for error in errors:
            print(f"  - {error}")
        return 1

    print(
        "manuscript consistency PASS: 799 state records, nine vertex "
        "bijections, scoped 12/144 no-design evidence, 235/161/403 descent "
        "counts, source concordance, and 564 ten-orbital records"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
