#!/usr/bin/env python3
"""Generate the double-anonymous review copy of the paper.

Reads results/report/main.md (the master document) and writes a
self-contained results/report/anonymized/ folder with author-identifying
material removed for double-anonymous peer review: the author list, the
affiliation/email/ORCID footnote, author contribution name, repository URL,
and archival DOI. The bibliography and CSL are copied alongside so the folder
renders standalone. The script fails if any identifying string survives,
and it is drift-proof: rerunning it always reproduces the anonymized copy
from the master, so edit main.md only.

Usage: python3 scripts/make_anonymized_report.py [--out DIR]
"""
import argparse
import pathlib
import re
import sys

REPORT = pathlib.Path("results/report")
MASTER = REPORT / "main.md"
BIB = REPORT / "references.bib"
CSL = REPORT / "physics-numeric.csl"

# Exact-match replacements; each must occur exactly once in the master so
# a silent drift in main.md fails loudly here instead of shipping a
# half-anonymized manuscript.
REPLACEMENTS = [
    # author and affiliation block
    ("  - |\n"
     "    James Orlando\\\n"
     "    Independent researcher, Derry, New Hampshire, USA\\\n"
     "    `jamie@orlandonh.com`",
     "  - Anonymized for double-anonymous review"),
    # named CRediT statement
    ("James Orlando: Conceptualization; Methodology; Software; Validation; "
     "Formal analysis; Investigation; Data curation; Writing - original "
     "draft; Writing - review and editing; Project administration.",
     "The author: Conceptualization; Methodology; Software; Validation; "
     "Formal analysis; Investigation; Data curation; Writing - original "
     "draft; Writing - review and editing; Project administration."),
    # public artifacts resolve to the author
    ("The Supplementary Material and machine-readable research-data archive "
     "are submitted with this manuscript. Versioned project materials are "
     "also available from Zenodo [@OrlandoZenodo2026] and the [public source "
     "repository](https://github.com/theorlandog/gpc-census).",
     "The Supplementary Material and machine-readable research-data archive "
     "are submitted with this manuscript. Identifying links to the versioned "
     "archival record and source repository are withheld for "
     "double-anonymous review."),
    # self-contained asset paths
    ("bibliography: results/report/references.bib",
     "bibliography: results/report/anonymized/references.bib"),
    ("csl: results/report/physics-numeric.csl",
     "csl: results/report/anonymized/physics-numeric.csl"),
]

# nothing on this list may survive in the anonymized output
FORBIDDEN = ["orlando", "jamie", "orcid", "zenodo", "derry",
             "we thank"]


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

    lowered = text.lower()
    leaks = [s for s in FORBIDDEN if s in lowered]
    if leaks:
        print(f"anonymize: identifying strings survived: {leaks}")
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "main.md").write_text(text)
    bibliography = BIB.read_text()
    bibliography, count = re.subn(
        r"\n?@\w+\{OrlandoZenodo2026,.*?\n\}\n?",
        "\n",
        bibliography,
        flags=re.DOTALL,
    )
    if count != 1:
        print(
            "anonymize: expected exactly one OrlandoZenodo2026 bibliography "
            f"entry, found {count}"
        )
        return 1
    bibliography = bibliography.strip() + "\n"
    (args.out / BIB.name).write_text(bibliography)
    csl = CSL.read_text().replace(
        "https://github.com/theorlandog/gpc-census/styles/physics-numeric",
        "https://example.org/styles/physics-numeric",
    )
    (args.out / CSL.name).write_text(csl)
    print(f"anonymized manuscript written to {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
