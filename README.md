# gpc-census

Which one-body occupation spectra can an N-fermion pure state have? The
answer is a convex polytope in the ordered simplex, cut out by the
generalized Pauli constraints (GPCs): finitely many linear inequalities
beyond 0 <= lambda <= 1 that Klyachko's solution of the pure-state
N-representability problem attaches to each pair (N, d) of particle and
orbital numbers. Altunbulak and Klyachko computed these systems through
rank d = 10 (Commun. Math. Phys. 282, 287 (2008)) and left two vertices
of the wedge^4 H_9 polytope resolved only numerically. This project
closed both, then classified every vertex of every known system, and
packages the machinery so the results are reproducible with one command.

## The classification

Each polytope vertex is an extremal spectrum, and the census asks how it
can be attained. Exactly one of three things is true:

- DESIGN-INT: a superposition of Slater determinants with nonnegative
  integer weights at the spectrum's natural denominator attains it, with
  a one-hop-free support (no two determinants share N-1 orbitals), so no
  phase cancellation is needed. A combinatorial design, hence the name.
- DESIGN-REAL: real nonnegative weights suffice but integers do not.
- INTERFERENCE: no nonnegative weighting exists at all; complex relative
  phases are forced. The canonical example is the vertex
  v_B = (20:14:14:14:14:4:4:4:4)/23 of wedge^4 H_9, whose extremal state
  requires cos(gamma) = 3/(4*sqrt(14)).

Verdicts are solver certificates (integer stage plus real stage, each a
feasibility proof), and the repository ships the complete census: every
vertex of every determinate system, ranks 6 through 10, with duality
transporting the classification to the complement systems. Highlights:
interference is absent in the N=4 series at rank 8 and first appears at
rank 9; interference fractions stabilize across ranks within each
particle-number series; padding and frozen-core lifts generate most of
each rank's interference vertices from lower-rank originals.

## The pipeline

Four stages, mirroring the mathematics:

1. Constraints: the rank 6-10 systems as validated lookup data
   (`gpc_census.constraints`), including the repaired (3,9) table (one
   inequality of the published version was lost at a page break; see
   `results/data/PROVENANCE.md`). Constraint generation for d >= 11 via
   Klyachko's extremal-edge algorithm is scoped in
   `docs/stage1_klyachko_spec.md` and is the project's research frontier.
2. Vertices: exact enumeration by lrslib in rational arithmetic
   (`gpc_census.polytope`).
3. Classification: the two-stage certificate solver
   (`gpc_census.classify`), CP-SAT reference backend with CBC fallback,
   every verdict labeled with its provenance.
4. States: attain a vertex spectrum with a complex pure state by
   alternating spectral projection, sparsify to minimal support, and
   exactify the numerics into certified closed forms (moduli snap to the
   natural denominator, phases to the p*sqrt(q)/r lattice, verification
   by exact characteristic polynomial identity in sympy)
   (`gpc_census.states`, `gpc_census.exactify`).

Every output passes structural validation (`gpc_census.validate`):
embedding coherence between ranks, frozen-core lifts, particle-hole
self-duality, physicality. These invariants caught every real error made
while building the census; nothing ships without them. The research
narrative and continuity notes live in `docs/RESEARCH.md`; the paper in
`results/report/`; the certified dataset with provenance and checksums
in `results/data/`.

## Quick start

```sh
uv sync --extra cpsat                       # reference solver backend
gpc-census constraints -n 3 -d 9            # the 52-inequality system
gpc-census polytope   -n 3 -d 9             # its 58 vertices, exactly
gpc-census classify   -n 4 -d 9             # verdicts incl. v_A and v_B
uv run python scripts/solve_all.py --preflight   # v_B end to end, gate for campaigns
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
attached; suffixed tags (e.g. `vX.Y.Zrc1`, `vX.Y.Z-beta`) start as pre-releases and are
promoted manually. Each push to main refreshes a rolling `snapshot`
pre-release holding the packages built from the branch tip. Every release
also carries `data-output.zip`: the compiled paper (`results/report/main.pdf`)
and the data results (`results/data/`), together with a signed provenance
attestation and SHA256SUMS of those files.

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

### Moment polytope pipeline

```sh
gpc-census constraints -n 3 -d 9      # tabulated constraint system (ranks 6-10)
gpc-census polytope   -n 3 -d 9      # exact vertex enumeration (requires lrslib)
gpc-census classify   -n 4 -d 9      # design/interference verdicts
gpc-census solve -n 4 -d 9 --den 23 --spectrum 20,14,14,14,14,4,4,4,4

```

Constraint generation beyond rank 10 is future work; version 1 ships the
known systems as a validated lookup table with particle-hole duality
transport. The test suite reproduces the historical census in
`results/data` end to end.
