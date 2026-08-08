#!/usr/bin/env python3
"""Emit documentation tables from their artifacts, so no count is typed twice.

The sync invariant: every number that appears in more than one place has
exactly one generator. A table in a doc that restates artifact data is
machine-emitted between markers

    <!-- sync:NAME:start -->  ...  <!-- sync:NAME:end -->

and this script owns what goes between them. ``--write`` regenerates in place,
``--check`` diffs and exits nonzero on drift, which is what
``tests/test_doc_sync.py`` runs.

Adding a block is two steps: put the markers in the doc, and add a BLOCKS entry
whose builder returns the body. A marker pair with no builder, or a builder
whose markers are missing, is an error rather than a silent skip, so the two
halves cannot drift apart either.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
DOCS = ROOT / "docs"

START = "<!-- sync:{name}:start -->"
END = "<!-- sync:{name}:end -->"


def read_json(path: Path):
    return json.loads(path.read_text())


def system_key(system: str) -> tuple[int, int]:
    n, d = (int(part) for part in system.strip("()").split(","))
    return n, d


# ------------------------------------------------------------------ builders


def census_coverage() -> str:
    """The README coverage table, from the census ledger."""
    rows = list(csv.DictReader((DATA / "census_master.csv").read_text().splitlines()))
    by_system: dict[str, Counter] = {}
    for row in rows:
        by_system.setdefault(row["system"], Counter())[row["verdict"]] += 1
    lines = [
        "| System (N, d) | Corners | Design | Interference | Certified closed forms |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    totals = Counter()
    for system in sorted(by_system, key=system_key):
        counts = by_system[system]
        corners = sum(counts.values())
        design = counts["DESIGN-INT"] + counts["DESIGN-REAL"]
        interference = counts["INTERFERENCE"]
        totals.update(
            {"corners": corners, "design": design, "interference": interference}
        )
        n, d = system_key(system)
        lines.append(
            f"| ({n}, {d}) | {corners} | {design} | {interference} | {corners} |"
        )
    lines.append(
        f"| **Total** | **{totals['corners']}** | **{totals['design']}** | "
        f"**{totals['interference']}** | **{totals['corners']}** |"
    )
    return "\n".join(lines)


def bracket_3_11_tally() -> str:
    """The (3,11) bracket tally, from the settlement and the level-five closure."""
    settlement = read_json(DOCS / "bracket_3_11_settlement.json")
    level5 = read_json(DATA / "level5_ressayre_3_11.json")
    tally = settlement["tally"]
    by_certificate = Counter(
        candidate.get("certificate") for candidate in settlement["candidates"]
    )
    certified = level5["rank_eleven"]["certified_rows"]
    lines = [
        "| verdict | candidates |",
        "| --- | ---: |",
        f"| TRUE-VERTEX | {tally['TRUE-VERTEX']} |",
        f"| REFUTED | {tally['REFUTED']} |",
        f"| OPEN | {tally['OPEN']} |",
        f"| **Total** | **{sum(tally.values())}** |",
        "",
        "By certificate:",
        "",
        "| certificate | candidates |",
        "| --- | ---: |",
    ]
    for name in sorted(by_certificate, key=lambda k: (k is None, str(k))):
        lines.append(f"| {name} | {by_certificate[name]} |")
    lines.append("")
    lines.append(
        f"The {certified} certified level-five rows refute "
        f"{by_certificate.get('LEVEL5-RESSAYRE', 0)} of these; validity is proved, "
        f"facetness and completeness of the (3,11) system are not."
    )
    return "\n".join(lines)


def _machine_footer(payload) -> list[str]:
    """State the machine, or state that the artifact does not record one.

    Multi-machine results are reported per machine and never averaged, so an
    artifact with no machine tag gets that said out loud rather than a tidy
    table that invites forgetting it.
    """
    machine = payload.get("machine") or payload.get("context", {}).get("machine")
    lines = [""]
    if machine:
        lines.append(f"_Machine: {machine}. Results are reported per machine, never averaged._")
    else:
        lines.append(
            "_This artifact records no machine tag, so these are one unrecorded "
            "machine's numbers. Runtime factors move across machines while method "
            "ordering holds, so treat the ordering as the finding._"
        )
    lines.append("")
    lines.append(
        "_Caveat, still unmet: there is no independently optimized dense baseline. "
        "Speedups are against this repository's own control arm._"
    )
    return lines


def _fmt(value, digits="%.6g"):
    if isinstance(value, float):
        return digits % value
    return str(value)


def compiled_sparse_benchmark() -> str:
    payload = read_json(DATA / "compiled_sparse_benchmark.json")
    context = payload.get("context", {})
    lines = [
        "| Method | Median runtime (s) | Maximum carrier-matrix error | ODE variables |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in payload["summary"]:
        lines.append(
            f"| {row['method']} | {_fmt(row['median_seconds'])} | "
            f"{_fmt(row['max_error'], '%.3g')} | {row['variables']} |"
        )
    for label, key in (
        ("dense 2-RDM, previous benchmark", "dense2rdm_seconds"),
        ("full FCI, previous benchmark", "full_fci_seconds"),
        ("carrier wavefunction, previous benchmark", "carrier_wavefunction_seconds"),
    ):
        if key in context:
            lines.append(f"| {label} | {_fmt(context[key])} | | |")
    ratios = payload.get("ratios", {})
    if ratios:
        lines.append("")
        lines.append("| Ratio | Value |")
        lines.append("| --- | ---: |")
        for name in sorted(ratios):
            lines.append(f"| {name} | {_fmt(ratios[name], '%.4g')} |")
    return "\n".join(lines + _machine_footer(payload))


def rdm_method_benchmark() -> str:
    payload = read_json(DATA / "rdm_method_benchmark.json")
    lines = [
        "| Method | Real ODE variables | Median runtime (s) | "
        "Median peak RSS (MiB) | Maximum reference error |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in payload["summary"]:
        lines.append(
            f"| {row['method']} | {row['real_ode_variables']} | "
            f"{_fmt(row['median_elapsed_seconds'])} | "
            f"{_fmt(row['median_peak_rss_mib'], '%.2f')} | "
            f"{_fmt(row['max_reference_error'], '%.3g')} |"
        )
    return "\n".join(lines + _machine_footer(payload))


def exact_observable_benchmark() -> str:
    payload = read_json(DATA / "exact_observable_benchmark.json")
    baselines = payload["baseline_medians"]
    lines = [
        "| Method | Median runtime (s) |",
        "| --- | ---: |",
        f"| compiled sparse | {_fmt(payload['compiled_sparse_median'])} |",
    ]
    for name in sorted(baselines):
        lines.append(f"| {name} | {_fmt(baselines[name])} |")
    speedups = payload.get("speedups", {})
    if speedups:
        lines.append("")
        lines.append("| Speedup | Value |")
        lines.append("| --- | ---: |")
        for name in sorted(speedups):
            lines.append(f"| {name} | {_fmt(speedups[name], '%.4g')} |")
    return "\n".join(lines + _machine_footer(payload))


def generalized_observable_benchmark() -> str:
    payload = read_json(DATA / "generalized_observable_benchmark.json")
    lines = [
        "| Cells q | Dense carrier amplitudes | Observable variables | "
        "Observable time (s) | Dense carrier time (s) | Speedup |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in payload["rows"]:
        speedup = row.get("observable_vs_dense_carrier_speedup")
        if speedup is None and row.get("observable_median_seconds"):
            speedup = row["dense_carrier_median_seconds"] / row["observable_median_seconds"]
        lines.append(
            f"| {row['q']} | {row['carrier_complex_amplitudes']} | "
            f"{row['observable_real_variables']} | "
            f"{_fmt(row['observable_median_seconds'])} | "
            f"{_fmt(row['dense_carrier_median_seconds'])} | "
            f"{_fmt(speedup, '%.4g')} |"
        )
    return "\n".join(lines + _machine_footer(payload))


def khovanskii_rank_growth() -> str:
    """Rank growth of the highest-weight semigroup's minimal generators."""
    payload = read_json(DATA / "equivariant_khovanskii_sagbi.json")
    rows = sorted(payload["rank_growth"], key=lambda row: system_key(row["system"]))
    window = rows[0]["window"] if rows else 0
    lines = [
        f"Common degree window `m <= {window}`, all values exact.",
        "",
        "| System | Raw values | Minimal generators | Primitive | Primitive / raw |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| `{row['system']}` | {row['raw_values']} | "
            f"{row['minimal_generators']} | {row['primitive_generators']} | "
            f"{row['primitive_ratio']} |"
        )
    return "\n".join(lines)


def khovanskii_valuations() -> str:
    """The three valuations on one table, compared in a fixed scope."""
    payload = read_json(DATA / "equivariant_khovanskii_sagbi.json")
    lines = [
        "| Valuation | System | Degrees | Status | Raw values | "
        "Minimal generators | Max generator degree | Compression |",
        "| --- | --- | ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    rows = [row for row in payload["v2_by_system"]
            if system_key(row["system"])[0] == 3
            and system_key(row["system"])[1] <= 8]
    rows += payload["v1_v3_by_system"]
    for row in sorted(rows, key=lambda r: (r["valuation"], system_key(r["system"]))):
        lines.append(
            f"| `{row['valuation']}` | `{row['system']}` | "
            f"`m <= {row['reached_degree']}` | {row['status']} | "
            f"{row['raw_values']} | {row['minimal_generators']} | "
            f"{row['max_generator_degree']} | {row['compression_ratio']} |"
        )
    return "\n".join(lines)


def khovanskii_held_out() -> str:
    """The preregistered held-out degree test, scored per system."""
    payload = read_json(DATA / "equivariant_khovanskii_sagbi.json")
    lines = [
        "| System | Fit `m <=` | Test degrees | Fitted generators | "
        "Misses | Verdict |",
        "| --- | ---: | --- | ---: | ---: | --- |",
    ]
    for row in sorted(payload["held_out_degree"],
                      key=lambda r: system_key(r["system"])):
        misses = sum(entry["misses"] for entry in row["per_degree"])
        degrees = ", ".join(str(d) for d in row["test_degrees"]) or "none"
        lines.append(
            f"| `{row['system']}` | {row['fit_degree']} | {degrees} | "
            f"{row['fitted_generators']} | {misses} | {row['verdict']} |"
        )
    return "\n".join(lines)


BLOCKS = {
    "khovanskii-rank-growth": (
        DOCS / "equivariant_khovanskii_sagbi.md", khovanskii_rank_growth),
    "khovanskii-valuations": (
        DOCS / "equivariant_khovanskii_sagbi.md", khovanskii_valuations),
    "khovanskii-held-out": (
        DOCS / "equivariant_khovanskii_sagbi.md", khovanskii_held_out),
    "census-coverage": (ROOT / "README.md", census_coverage),
    "bracket-3-11": (DOCS / "level5_ressayre_3_11.md", bracket_3_11_tally),
    "compiled-sparse-benchmark": (
        DOCS / "compiled_sparse_benchmark.md", compiled_sparse_benchmark),
    "rdm-method-benchmark": (
        DOCS / "rdm_method_benchmark.md", rdm_method_benchmark),
    "exact-observable-benchmark": (
        DOCS / "exact_observable_benchmark.md", exact_observable_benchmark),
    "generalized-observable-benchmark": (
        DOCS / "generalized_observable_benchmark.md", generalized_observable_benchmark),
}


def block_pattern(name: str) -> re.Pattern:
    return re.compile(
        re.escape(START.format(name=name)) + r".*?" + re.escape(END.format(name=name)),
        re.DOTALL,
    )


def rendered(name: str, body: str) -> str:
    return f"{START.format(name=name)}\n{body}\n{END.format(name=name)}"


def apply_blocks(write: bool) -> tuple[list[str], list[str]]:
    drift, missing = [], []
    by_path: dict[Path, str] = {}
    for name, (path, build) in BLOCKS.items():
        text = by_path.get(path)
        if text is None:
            if not path.exists():
                missing.append(f"{name}: {path} does not exist")
                continue
            text = path.read_text()
        pattern = block_pattern(name)
        if not pattern.search(text):
            missing.append(
                f"{name}: markers not found in {path.relative_to(ROOT)}; add "
                f"{START.format(name=name)} and {END.format(name=name)}"
            )
            continue
        updated = pattern.sub(lambda _m: rendered(name, build()), text, count=1)
        if updated != text:
            drift.append(f"{name} in {path.relative_to(ROOT)}")
        by_path[path] = updated
    if write:
        for path, text in by_path.items():
            path.write_text(text)
    return drift, missing


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true", help="regenerate blocks in place")
    group.add_argument("--check", action="store_true", help="fail on drift")
    args = parser.parse_args()

    drift, missing = apply_blocks(write=args.write)
    for problem in missing:
        print(f"MISSING {problem}", file=sys.stderr)
    if args.write:
        if drift:
            for item in drift:
                print(f"rewrote {item}")
        else:
            print(f"all {len(BLOCKS)} emitted blocks were already current")
        return 1 if missing else 0
    if drift:
        for item in drift:
            print(f"DRIFT {item} differs from its artifact; run --write", file=sys.stderr)
    if not drift and not missing:
        print(f"all {len(BLOCKS)} emitted blocks are current")
    return 1 if (drift or missing) else 0


if __name__ == "__main__":
    raise SystemExit(main())
