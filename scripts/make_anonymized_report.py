#!/usr/bin/env python3
"""Generate the double-anonymous review copy of the paper.

Reads results/report/main.md (the master document) and writes a
self-contained results/report/anonymized/ folder with author-identifying
material removed per the IOP double-anonymous peer-review guidelines:
the author list, the affiliation/email/ORCID footnote, the preprint DOI
(it resolves to the author's record), and the personal acknowledgment in
Theorem 2. The bibliography and CSL are copied alongside so the folder
renders standalone. The script fails if any identifying string survives,
and it is drift-proof: rerunning it always reproduces the anonymized copy
from the master, so edit main.md only.

Usage: python3 scripts/make_anonymized_report.py [--out DIR]
"""
import argparse
import pathlib
import shutil
import sys

REPORT = pathlib.Path("results/report")
MASTER = REPORT / "main.md"
BIB = REPORT / "references.bib"
CSL = REPORT / "institute-of-physics-numeric.csl"

# Exact-match replacements; each must occur exactly once in the master so
# a silent drift in main.md fails loudly here instead of shipping a
# half-anonymized manuscript.
REPLACEMENTS = [
    # author list (also removes the footnote reference)
    ("  - James Orlando[^1]",
     "  - Anonymized for double-anonymous review"),
    # the preprint DOI resolves to the author's record
    ("  Preprint: [doi:10.5281/zenodo.21313736]"
     "(https://doi.org/10.5281/zenodo.21313736)",
     "  Preprint DOI withheld for double-anonymous review"),
    # personal acknowledgment inside Theorem 2 (restore on acceptance)
    ("(We thank T. Maciazek for prompting this basis-free formulation: a",
     "(A"),
    # self-contained asset paths
    ("bibliography: results/report/references.bib",
     "bibliography: results/report/anonymized/references.bib"),
    ("csl: results/report/institute-of-physics-numeric.csl",
     "csl: results/report/anonymized/institute-of-physics-numeric.csl"),
]

# the affiliation footnote definition is dropped entirely
FOOTNOTE_PREFIX = "[^1]:"

# nothing on this list may survive in the anonymized output
FORBIDDEN = ["orlando", "jamie", "orcid", "zenodo", "derry",
             "we thank", "[^1]"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=pathlib.Path,
                        default=REPORT / "anonymized")
    args = parser.parse_args()

    text = MASTER.read_text()
    for old, new in REPLACEMENTS:
        if text.count(old) != 1:
            print(f"anonymize: expected exactly one occurrence of "
                  f"{old[:60]!r}, found {text.count(old)}; master drifted, "
                  f"update scripts/make_anonymized_report.py")
            return 1
        text = text.replace(old, new)

    lines = [ln for ln in text.splitlines()
             if not ln.startswith(FOOTNOTE_PREFIX)]
    text = "\n".join(lines) + "\n"

    lowered = text.lower()
    leaks = [s for s in FORBIDDEN if s in lowered]
    if leaks:
        print(f"anonymize: identifying strings survived: {leaks}")
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "main.md").write_text(text)
    shutil.copy2(BIB, args.out / BIB.name)
    shutil.copy2(CSL, args.out / CSL.name)
    print(f"anonymized manuscript written to {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
