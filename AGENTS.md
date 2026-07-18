# gpc-census agent instructions

gpc-census is an algorithmic pipeline that constructs exact extremal states for
fermionic natural-occupation-number (moment) polytopes. The accompanying paper
lives in `results/report/main.md`; computed data results live under
`results/data/`.

## House rules

1. Never use em dashes in writing (docs, comments, commit messages).
2. GitHub Actions: only GitHub-owned actions (`actions/*`) are allowed, always
   pinned to a full commit hash with a version comment
   (`uses: actions/checkout@<sha> # v6.0.3`). Install other tooling (uv, etc.)
   with plain shell steps, version pinned.
3. New library dependencies must have been released at least 14 days before
   adoption. `make upgrade` enforces this via `uv lock --upgrade
   --exclude-newer`.
4. Always use a lock file. `uv.lock` is committed and CI installs with
   `uv sync --locked`.
5. Keep third-party dependencies to a minimum. Prefer mature, widely adopted
   libraries; no fly-by-night packages.

## Project requirements (keep these true)

- Python package managed with `pyproject.toml` and uv (hatchling build backend).
- Usable both as a library (`import gpc_census`) and as a CLI (`gpc-census`).
- Ships as a pip package (wheel + sdist via `uv build`) and as an RPM
  (`gpc-census.spec`, Fedora `pyproject-rpm-macros`).
- Runtime dependencies: numpy, scipy, pulp, ortools (pinned). classify.py
  detects ortools at import and falls back to CBC when absent or unpinned,
  so the RPM can omit the ortools dependency; verdicts carry a solver field
  and gpc_census.certify upgrades them to exact certificates.

## Layout

- `src/gpc_census/`: library code (`core.py`) and CLI (`cli.py`, argparse,
  entry point `gpc-census = gpc_census.cli:main`).
- `tests/`: pytest suite. `tests/data/` is the fixture copy of the dataset
  that the tests (and the RPM %check) read; a test asserts it stays
  byte-identical to `results/data/` whenever the latter is present.
- `gpc-census.spec`: RPM spec. `Makefile`: build entry points.
- `results/report/main.md`: the paper. The markdown is the master document,
  in pandoc-crossref syntax; `make report` renders `main.pdf` with pandoc in
  the pinned `pandoc/extra` container image, which bundles pandoc, a matched
  pandoc-crossref, and a TeX engine (podman by default; `CONTAINER=docker`
  to override).
  Section, equation, and table refs are live crossref citations
  (`[-@sec:x]`); theorem-family numbering is literal text since crossref has
  no theorem type, so renumber by hand when inserting theorems. Citations
  resolve via citeproc against `results/report/references.bib` (APS numeric
  style, vendored CSL); `nocite: "@*"` prints uncited entries. The build
  fails if any reference or citation does not resolve. `results/report` is
  excluded from the sdist and RPM; the PDF ships in the release
  `data-output.zip`.
- `results/data/`: computed data results, shipped in the release
  `data-output.zip`.
- `.github/pull_request_template.md`: PR template. The `pr-metadata` workflow
  fills its Metadata section (branch, head SHA, CI build version, diff stats)
  when a PR is opened; keep the `pr-metadata:start/end` markers intact.

## Commands

```sh
make sync     # uv sync (create .venv, install dev deps)
make test     # uv run pytest
make lint     # uv run ruff check
make build    # uv build -> dist/ (wheel + sdist)
make rpm      # sdist + rpmbuild into build/rpm/ (needs rpm-build etc.)
make report   # render results/report/main.pdf with pandoc (containerized)
make upgrade  # upgrade locked deps, excluding releases newer than 14 days
```

## Versioning (CI-owned)

- Single source of truth: `project.version` in `pyproject.toml`. Never hardcode
  a version anywhere else; `gpc_census.__version__` reads `importlib.metadata`.
- GitHub Actions (`.github/workflows/build.yml`) stamps builds with
  `uv version`: tag `vX.Y.Z` gives `X.Y.Z`; `main` gives
  `<base>+main.<short-sha>`; other refs give `<base>+git.<short-sha>`.
- `make rpm` regenerates `build/gpc-census.spec` with the pyproject version;
  the committed spec's `Version:` line is not authoritative.
- Releases: a clean `vX.Y.Z` tag publishes a GitHub release marked latest; a
  suffixed tag (e.g. `vX.Y.Zrc1`) starts as a pre-release for manual
  promotion, and its suffix must be valid PEP 440 (rc1, b1, a1, .post1,
  .dev1) or the version stamp fails the build. A push to main publishes a
  rolling `snapshot-<short-sha>` pre-release and prunes older snapshots.
  Tag names are single use: immutable releases reserve a published tag name
  forever, even after deletion (the bare name `snapshot` is already burned
  this way), so never delete-and-recreate a release tag. Never name a
  release tag after a branch either: a release tag named `main` makes the
  refname ambiguous with the branch. The release job runs directly on the runner
  (no job container) and rebuilds everything fresh (it does not reuse CI
  artifacts): wheel and sdist with uv, RPMs inside a Fedora container and
  the paper inside the pandoc container, both driven by the runner's
  podman, then publishes with the `gh` CLI. It also attaches
  `data-output.zip`: the paper PDF plus `results/data/`, with a signed
  provenance attestation (`actions/attest-build-provenance`) and SHA256SUMS
  of those files inside the zip.

## Gotchas

- `.gitignore` ignores `*.spec` (PyInstaller template pattern) with a
  `!gpc-census.spec` negation. Keep the RPM spec un-ignored.
- Sdist contents are pinned explicitly in `[tool.hatch.build.targets.sdist]`.
  Files that must ship in the sdist (e.g. for the RPM build) go in that list.
- `CLAUDE.md` is a symlink to this file. Edit `AGENTS.md` only.

## Research context

Read docs/RESEARCH.md before working on the science. It encodes the
classification trichotomy, the validation law, campaign state, and the
location of every source document. The short version: never ship a result
that has not passed the structural invariants in gpc_census.validate, and
never change state-solving code without the v_B preflight.
