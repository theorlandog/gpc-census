# gpc-census

gpc-census implements an algorithmic pipeline that constructs exact extremal states for fermionic natural-occupation-number (moment) polytopes.

It is usable both as a Python library (`import gpc_census`) and as a command-line tool (`gpc-census`).

## Installation

From a source checkout:

```sh
pip install .
```

Or with [uv](https://docs.astral.sh/uv/):

```sh
uv pip install .
```

## Usage

### As a library

```python
from gpc_census import slater_vertices

# Vertices of the Pauli polytope Delta(d=6, n=3): occupation-number
# vectors of Slater determinants for 3 fermions in 6 orbitals.
for vec in slater_vertices(6, 3):
    print(vec)
```

### As a CLI

```sh
gpc-census --version
gpc-census vertices -d 6 -n 3
```

## Development

Requires [uv](https://docs.astral.sh/uv/). Common tasks are wrapped in the Makefile:

```sh
make sync     # create .venv and install dev dependencies
make test     # run pytest
make lint     # run ruff
make upgrade  # upgrade locked dependencies (14-day release cooldown)
```

## Versioning

The canonical version lives in `pyproject.toml` (`project.version`). At runtime,
`gpc_census.__version__` reports whatever version was actually installed (via
`importlib.metadata`).

CI stamps each build before packaging (`uv version <computed>`):

| Ref              | Package version           |
| ---------------- | ------------------------- |
| tag `vX.Y.Z`     | `X.Y.Z`                   |
| `main` branch    | `<base>+main.<short-sha>` |
| other refs / PRs | `<base>+git.<short-sha>`  |

The GitHub Actions workflow (`.github/workflows/build.yml`) runs tests, then
builds and uploads the wheel/sdist and the RPMs as artifacts. Tag pushes
(`vX.Y.Z`) also publish a GitHub release marked latest with those packages
attached; each push to main refreshes a rolling `snapshot` pre-release holding
the packages built from the branch tip.

## Building packages

### Wheel and sdist (pip package)

```sh
make build          # uv build -> dist/*.whl and dist/*.tar.gz
```

The resulting wheel installs with `pip install dist/gpc_census-<version>-py3-none-any.whl`.

### RPM

The spec file is `gpc-census.spec` and builds from the sdist using Fedora's
`pyproject-rpm-macros`. Build dependencies (Fedora/RHEL):

```sh
sudo dnf install rpm-build python3-devel pyproject-rpm-macros \
    python3-hatchling python3-pytest
```

Then:

```sh
make rpm            # builds sdist, then rpmbuild into build/rpm/
```

`make rpm` regenerates `build/gpc-census.spec` with the version from
`pyproject.toml`, so CI-stamped versions flow into the RPM automatically.

Binary RPMs land in `build/rpm/RPMS/noarch/` and the source RPM in
`build/rpm/SRPMS/`. `make srpm` builds only the source RPM.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
