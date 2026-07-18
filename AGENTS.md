# gpc-census agent instructions

gpc-census is an algorithmic pipeline that constructs exact extremal states for
fermionic natural-occupation-number (moment) polytopes. The accompanying paper
lives in `report/main.tex`.

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
- Zero runtime dependencies so far. Add new ones consciously; each must be
  available as an RPM (`python3dist(...)`) or the RPM build breaks.

## Layout

- `src/gpc_census/`: library code (`core.py`) and CLI (`cli.py`, argparse,
  entry point `gpc-census = gpc_census.cli:main`).
- `tests/`: pytest suite.
- `gpc-census.spec`: RPM spec. `Makefile`: build entry points.
- `report/main.tex`: LaTeX paper, built with `make report` (latexmk).
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
make report   # build report/main.pdf with latexmk
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
- Releases: a tag push publishes a GitHub release marked latest; a push to
  main refreshes a rolling `snapshot` pre-release (never name it after a
  branch: a release tag named `main` makes the refname ambiguous with the
  branch). The release job rebuilds the
  wheel, sdist, and RPMs fresh in a Fedora container (it does not reuse CI
  artifacts) and publishes with the `gh` CLI.

## Gotchas

- `.gitignore` ignores `*.spec` (PyInstaller template pattern) with a
  `!gpc-census.spec` negation. Keep the RPM spec un-ignored.
- Sdist contents are pinned explicitly in `[tool.hatch.build.targets.sdist]`.
  Files that must ship in the sdist (e.g. for the RPM build) go in that list.
- `CLAUDE.md` is a symlink to this file. Edit `AGENTS.md` only.
