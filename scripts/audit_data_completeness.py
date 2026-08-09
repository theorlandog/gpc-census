#!/usr/bin/env python3
"""Audit that every published data artifact is accounted for four ways.

The sync invariant: every file under ``results/data/`` must have a
``SHA256SUMS`` entry, a ``PROVENANCE.md`` row, a named generating script, and
at least one guard test. An artifact missing any of those is an orphan, and an
orphan is a number nobody can regenerate or check.

Detection is literal by default: the artifact's file name must appear in the
relevant place. That is deliberate, because inferring coverage is exactly the
drift this audit exists to catch.

Literal matching cannot see two legitimate patterns, so both are handled by an
explicit rule table rather than by loosening the default:

- artifacts whose path is built by f-string, ``vertices_{n}_{d}.json`` and
  friends, whose name never appears as a literal anywhere;
- artifact families documented by one prose block covering the whole family
  rather than a row per file.

Every relaxation is therefore declared in ``FAMILY_RULES`` with the exact
marker that must be present, so renaming a generator or deleting a provenance
block still fails the audit. A family rule that stops matching any file is
itself an error, so the table cannot rot silently.

The audit also runs in reverse. The forward direction catches an artifact that
nobody can regenerate or check. The reverse direction catches the opposite
failure, which is just as bad and was invisible until it happened: a guard test
that names an artifact which does not exist. Such a test does not fail, it
SKIPS, so a study can ship its script, its documentation and its guards while
the numbers those guards would check are absent from the repository. That is
how nine guards for the chemistry study sat inert on main. A skipping guard
proves nothing, so an artifact referenced by a test and missing from disk is an
orphaned guard and fails this audit.

Exit status is 0 when there are no orphans and 1 otherwise, so the same code
backs ``tests/test_data_completeness.py``.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SUMS = DATA / "SHA256SUMS"
PROVENANCE = DATA / "PROVENANCE.md"
SCRIPTS = ROOT / "scripts"
TESTS = ROOT / "tests"

# The bookkeeping itself, not artifacts it describes.
SELF = {"SHA256SUMS", "PROVENANCE.md"}

# Artifact-shaped names a test may reference without the file existing, each
# with the reason. Anything not listed here that a test names must be on disk.
# Empty on purpose: nothing is currently exempt. Add an entry only with a
# reason, and only for a test that deliberately asserts an artifact's absence.
GUARD_REFERENCE_EXEMPTIONS: dict[str, str] = {}

# Suffixes that identify a published artifact rather than an incidental path.
ARTIFACT_SUFFIXES = (".json", ".jsonl", ".csv", ".txt")

# How tests name an artifact under results/data. Interpolated names (f-strings)
# are deliberately not matched: those are the documented families, and the
# forward direction already covers them through FAMILY_RULES.
GUARD_REFERENCE_PATTERNS = (
    # DATA / "name.json"  and  DATA / "sub" / "name.json"
    re.compile(r'DATA\s*/\s*"([^"{}]+?)"(?:\s*/\s*"([^"{}]+?)")?'),
    # ROOT / "results" / "data" / "name.json"
    re.compile(
        r'"results"\s*/\s*"data"\s*/\s*"([^"{}]+?)"(?:\s*/\s*"([^"{}]+?)")?'
    ),
    # a literal path fragment, "results/data/name.json"
    re.compile(r'results/data/([\w./-]+)'),
)

# glob -> the marker that must appear in each of provenance, scripts, tests.
# A None marker means the default literal file-name check still applies.
FAMILY_RULES = (
    {
        "glob": "vertices/vertices_*.json",
        "why": "path built as vertices_{n}_{d}.json; one provenance block covers all nine",
        "provenance": "Vertices: exact enumeration",
        "script": 'f"vertices_{n}_{d}.json"',
        "test": None,
    },
    {
        "glob": "census/census_*_results.txt",
        "why": "path built as census_{tag}_results.txt; raw solver logs behind census_master",
        "provenance": "Census verdicts:",
        "script": 'f"census_{tag}_results.txt"',
        "test": None,
    },
    # Known gap, declared rather than hidden. Both artifacts predate the
    # current scripts and no generator for them survives anywhere in the
    # repository. They keep the checksum, provenance and guard-test
    # requirements; only the generating-script requirement is waived, and the
    # provenance block they point at says so in as many words. Deleting or
    # regenerating them is an open decision, not this audit's to take.
    {
        "glob": "cyclic_window_benchmark.json",
        "why": "no surviving generator; see the UNREGENERABLE block in PROVENANCE.md",
        "provenance": "UNREGENERABLE",
        "script": "UNREGENERABLE_NO_GENERATOR",
        "test": None,
    },
    {
        "glob": "entangling_window_hamiltonian.json",
        "why": "no surviving generator; see the UNREGENERABLE block in PROVENANCE.md",
        "provenance": "UNREGENERABLE",
        "script": "UNREGENERABLE_NO_GENERATOR",
        "test": None,
    },
)

# The literal string a waived family looks for, so the waiver is greppable.
UNREGENERABLE_NO_GENERATOR = "UNREGENERABLE_NO_GENERATOR"


def artifacts() -> list[Path]:
    return sorted(
        path for path in DATA.rglob("*") if path.is_file() and path.name not in SELF
    )


def checksum_names() -> set[str]:
    if not SUMS.exists():
        return set()
    names = set()
    for line in SUMS.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        # sha256sum writes "<digest>  <path>", the path carrying a "*" for
        # binary mode and a "./" prefix from the directory it was run in.
        name = parts[1].strip().lstrip("*")
        if name.startswith("./"):
            name = name[2:]
        names.add(name)
    return names


def text_of(paths) -> str:
    return "\n".join(p.read_text(errors="ignore") for p in paths if p.is_file())


def rule_for(relative: str):
    for rule in FAMILY_RULES:
        if fnmatch.fnmatch(relative, rule["glob"]):
            return rule
    return None


def guard_references() -> dict[str, list[str]]:
    """Artifact names each guard test claims to check, mapped to the tests."""
    references: dict[str, list[str]] = {}
    for test in sorted(TESTS.rglob("*.py")):
        text = test.read_text(errors="ignore")
        for pattern in GUARD_REFERENCE_PATTERNS:
            for match in pattern.findall(text):
                parts = [part for part in _as_tuple(match) if part]
                if not parts:
                    continue
                relative = "/".join(parts)
                if not relative.endswith(ARTIFACT_SUFFIXES):
                    continue
                references.setdefault(relative, []).append(test.name)
    return references


def _as_tuple(match) -> tuple:
    return match if isinstance(match, tuple) else (match,)


def orphaned_guards() -> list[dict]:
    """Guards naming an artifact that is not on disk.

    Such a test skips rather than fails, so without this check the repository
    can carry a study whose numbers were never committed and still report a
    clean audit.
    """
    findings = []
    for relative, tests in sorted(guard_references().items()):
        name = relative.split("/")[-1]
        if name in SELF or name in GUARD_REFERENCE_EXEMPTIONS:
            continue
        if (DATA / relative).is_file():
            continue
        if any(path.name == name for path in artifacts()):
            continue
        findings.append(
            {
                "artifact": relative,
                "referenced_by": sorted(set(tests)),
                "why": (
                    "a guard test names this artifact and it is not under "
                    "results/data, so the guard skips instead of checking"
                ),
            }
        )
    return findings


def audit() -> dict:
    sums = checksum_names()
    corpus = {
        "PROVENANCE.md": PROVENANCE.read_text() if PROVENANCE.exists() else "",
        "generating script": text_of(sorted(SCRIPTS.rglob("*.py"))),
        "guard test": text_of(sorted(TESTS.rglob("*.py"))),
    }
    rule_keys = {
        "PROVENANCE.md": "provenance",
        "generating script": "script",
        "guard test": "test",
    }

    found = sorted(artifacts())
    findings, used_rules = [], set()
    for path in found:
        relative = path.relative_to(DATA).as_posix()
        rule = rule_for(relative)
        if rule is not None:
            used_rules.add(rule["glob"])
        missing = []
        if relative not in sums and path.name not in sums:
            missing.append("SHA256SUMS")
        for label, text in corpus.items():
            marker = rule.get(rule_keys[label]) if rule else None
            if marker is not None:
                if marker not in text:
                    missing.append(f"{label} (family marker {marker!r})")
                continue
            if path.name not in text and relative not in text:
                missing.append(label)
        if missing:
            findings.append({"artifact": relative, "missing": missing})

    stale = [rule["glob"] for rule in FAMILY_RULES if rule["glob"] not in used_rules]
    guards = orphaned_guards()
    return {
        "artifacts_audited": len(found),
        "checksum_entries": len(sums),
        "family_rules": len(FAMILY_RULES),
        "stale_family_rules": stale,
        "orphans": findings,
        "guard_references_checked": len(guard_references()),
        "orphaned_guards": guards,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit the report as JSON")
    args = parser.parse_args()
    report = audit()
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return (
            0
            if not report["orphans"]
            and not report["stale_family_rules"]
            and not report["orphaned_guards"]
            else 1
        )
    print(f"audited {report['artifacts_audited']} artifacts under results/data")
    print(f"SHA256SUMS entries: {report['checksum_entries']}")
    print(f"family rules: {report['family_rules']}")
    for glob in report["stale_family_rules"]:
        print(f"  STALE family rule matches nothing: {glob}")
    print(f"guard references checked: {report['guard_references_checked']}")
    if not report["orphans"]:
        print("no orphans: every artifact has checksum, provenance, generator, guard test")
    else:
        print(f"ORPHANS: {len(report['orphans'])}")
        for finding in report["orphans"]:
            print(f"  {finding['artifact']}: missing {', '.join(finding['missing'])}")
    if not report["orphaned_guards"]:
        print("no orphaned guards: every artifact a test names is on disk")
    else:
        print(f"ORPHANED GUARDS: {len(report['orphaned_guards'])}")
        for finding in report["orphaned_guards"]:
            tests = ", ".join(finding["referenced_by"])
            print(f"  {finding['artifact']}: named by {tests}, not on disk")
    return (
        0
        if not report["orphans"]
        and not report["stale_family_rules"]
        and not report["orphaned_guards"]
        else 1
    )


if __name__ == "__main__":
    sys.exit(main())
