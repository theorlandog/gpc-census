#!/usr/bin/env python3
"""Generate pandoc-crossref markdown from the LaTeX report.

pandoc's LaTeX reader resolves amsthm numbering itself and
pandoc-crossref has no theorem type, so theorem, proposition,
definition and remark references keep their baked numbers as plain
anchor links. Section, equation and table references become live
pandoc-crossref citations ({#eq:x} attributes, [-@sec:x] refs) that
renumber whenever the markdown is rendered with the filter.
"""

import os
import pathlib
import re
import shlex
import subprocess
import sys

# The Makefile points this at the pinned pandoc container image.
PANDOC = shlex.split(os.environ.get("PANDOC", "pandoc"))


def attach_stray_section_labels(md: str) -> str:
    """Move floating []{#sec:x label="sec:x"} spans onto the preceding header.

    pandoc only attaches a \\label to the section header when it directly
    follows \\section{}; labels placed later in the section body come out
    as bare anchor spans that carry no number.
    """
    lines = md.split("\n")
    for i, line in enumerate(lines):
        m = re.fullmatch(r'\[\]\{(#sec:[^ }]+) label="[^"]*"\}', line.strip())
        if not m:
            continue
        lines[i] = ""
        for j in range(i - 1, -1, -1):
            if lines[j].startswith("#") and not lines[j].endswith("}"):
                lines[j] += " {" + m.group(1) + "}"
                break
    return "\n".join(lines)


def tag_equation_labels(md: str) -> str:
    """Turn $$...\\label{eq:x}...$$ into $$...$$ {#eq:x}."""

    def block(m: re.Match) -> str:
        body = m.group(1)
        lab = re.search(r"\\label\{(eq:[^}]+)\}", body)
        if not lab:
            return m.group(0)
        body = body.replace(lab.group(0), "", 1)
        return f"$${body}$$ {{#{lab.group(1)}}}"

    return re.sub(r"\$\$(.*?)\$\$", block, md, flags=re.S)


def table_div_to_caption_attr(md: str) -> str:
    """Unwrap ::: {#tab:x} divs and tag the caption with {#tbl:x}."""

    def div(m: re.Match) -> str:
        body = "\n".join(
            ln[2:] if ln.startswith("  ") else ln for ln in m.group(2).split("\n")
        )
        body, n = re.subn(r"(?m)^: (.*)$", rf": \1 {{#tbl:{m.group(1)}}}", body)
        if n != 1:
            raise SystemExit("tex2md: expected exactly one caption per table div")
        return body

    return re.sub(r"::: \{#tab:([^}]+)\}\n(.*?)\n:::", div, md, flags=re.S)


def rewrite_refs(md: str) -> str:
    def ref(m: re.Match) -> str:
        text, target, rtype = m.groups()
        if target.startswith("tab:"):
            target = "tbl:" + target[4:]
        if target.split(":", 1)[0] in ("sec", "eq", "tbl"):
            # The surrounding prose already says Section/Table, so the
            # crossref prefix is suppressed; \eqref keeps the "eq." prefix.
            return f"[@{target}]" if rtype == "eqref" else f"[-@{target}]"
        return f"[{text}](#{target})"

    return re.sub(
        r"\[((?:\\.|[^\\\]])*)\]"
        r"\(#([a-z]+:[^)]+)\)"
        r'\{reference-type="(ref|eqref)" reference="[^"]*"\}',
        ref,
        md,
    )


def convert(tex: pathlib.Path) -> str:
    md = subprocess.run(
        [*PANDOC, str(tex), "-f", "latex", "-t", "markdown", "-s", "--wrap=none"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    md = attach_stray_section_labels(md)
    md = tag_equation_labels(md)
    md = table_div_to_caption_attr(md)
    md = rewrite_refs(md)

    for leftover in ("reference-type=", "\\label{", "{#tab:", ' label="'):
        if leftover in md:
            raise SystemExit(f"tex2md: unhandled construct left in output: {leftover}")
    return md


def main() -> None:
    tex, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    out.write_text(convert(tex))


if __name__ == "__main__":
    main()
