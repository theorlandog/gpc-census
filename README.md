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
   (`gpc_census.constraints`). Constraint generation for d >= 11 via
   Klyachko's extremal-edge algorithm is scoped in
   `docs/stage1_klyachko_spec.md` and is the project's research frontier.
2. Vertices: exact enumeration by lrslib in rational arithmetic
   (`gpc_census.polytope`).
3. Classification: the two-stage certificate solver
   (`gpc_census.classify`), CP-SAT reference backend with CBC fallback,
   every verdict labeled with its provenance.
4. States: construct and certify the extremal state of each vertex,
   routed by the verdict from stage 3 so no work is wasted. The state
   algorithm is a pipeline of its own:
   1. Route. A design vertex is built directly from its classification
      witness: the design's support is one-hop free, so
      psi = sum sqrt(k_t/den) |t> has a diagonal 1-RDM equal to the
      spectrum, exact by construction, no iterative solve
      (`solve_design_vertex`). Interference (and DESIGN-REAL) vertices go
      to the state solver.
   2. Block-budget preflight. `min_block_count` computes, as pure
      feasibility with no phase solve, how many 2x2 natural-orbital blocks
      the sparse realization needs: 0 for a design, k >= 1 for
      interference, or None when the vertex lies outside the current
      ansatz family (the extended-ansatz frontier), which fails fast
      instead of sweeping.
   3. Weights-first solve. `solve_vertex_exact_first` sweeps block ansatze
      at that budget, fixing exact moduli sqrt(k/den) and solving only the
      phases as a smooth quartic with analytic gradients (no
      eigendecomposition, immune to the second-order degeneracy flatness
      that stalls gradient and moment-matching methods).
   4. Block-size generalization (k x k). A 2x2 block mixes two degenerate
      classes; some vertices need a natural-orbital rotation mixing k >= 3
      classes at once. The general block is a k-mode clique whose canonical
      diagonal (integer occupations) is any vector majorized by the block's
      k eigenvalues: the Schur-Horn theorem (I. Schur, Sitzungsber. Berl.
      Math. Ges. 22, 9 (1923); A. Horn, Amer. J. Math. 76, 620 (1954))
      characterizes exactly these, and `_schur_horn_diagonals` enumerates
      them (the 2x2 split is the k=2 case). Because a clique mixes distinct
      eigenvalues the block is non-degenerate, so `phase_solve_clique`
      matches it by its characteristic-polynomial coefficients (power sums
      and Newton's identities, analytic gradients, no eigendecomposition),
      and `min_clique_count` is the block-size preflight. This is what the
      (4,9) vertices outside the 2x2 family need: idx 24,
      (9:6:5:5:5:2:2:1:1)/9, reconstructs from the spectrum alone with one
      3x3 clique. `max_clique` enables it. On the (4,9) interference
      vertices, k=3 cliques raise certified closed forms from 4 to 10 of 16.
   5. Exactify. Moduli snap to the natural denominator; after gauge-fixing
      the single-particle U(1)^d phase freedom, the residual interference
      phases are recognized (rational multiples of pi, or cosines on the
      p*sqrt(q)/r lattice) and the symbolic state is verified by exact
      characteristic polynomial identity in sympy
      (`gpc_census.exactify`). v_B certifies as a single 14/4 block with a
      pi/8 interference phase. Larger blocks are often easier here: a
      k-clique has a (k-1)(k-2)/2-dimensional Schur-Horn fiber, so for
      k >= 3 a real realization (phases 0/pi) is frequently available and
      certifies with no algebraic recognition at all. idx 24's closed form
      is real, (|0125>+|0134>+|0237>+|0245>)/3 + sqrt(2)/3 |0268> +
      sqrt(3)/3 |0348>. The tight k=2 fiber is the hard case (v_B genuinely
      needs pi/8); the residual across all systems is the small tail whose
      phases are higher-degree algebraic (no real, no p*sqrt(q)/r).
   The historical alternating-projection solver (`attain`) remains as a
   general Tier-A fallback: `scripts/solve_all.py` cascades to it when the
   block solver fails, so every vertex records at least a numeric state.

Every output passes structural validation (`gpc_census.validate`):
embedding coherence between ranks, frozen-core lifts, particle-hole
self-duality, physicality. These invariants caught every real error made
while building the census; nothing ships without them. The research
narrative and continuity notes live in `docs/RESEARCH.md`; the paper in
`results/report/`; the certified dataset with provenance and checksums
in `results/data/`.

## Quick start

```sh
uv sync                                     # install deps (CP-SAT backend included)
gpc-census constraints -n 3 -d 9            # the 52-inequality system
gpc-census polytope   -n 3 -d 9             # its 58 vertices, exactly
gpc-census classify   -n 4 -d 9             # verdicts incl. v_A and v_B
uv run python scripts/solve_all.py --preflight   # reconstruct + certify v_B, gate for campaigns
uv run python scripts/solve_all.py --all         # full state census, routed by verdict
```

The state campaign routes each vertex by verdict (design vertices built
from their witness, interference vertices through the block solver, with
a cascade to `attain` on the frontier), is restartable and checkpointed,
and writes `results/data/states.jsonl`. `--solver`, `--max-blocks`, and
`--workers` tune coverage and parallelism.

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
builds the wheel and sdist. Tag pushes
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

RPM packaging (`gpc-census.spec`, `make rpm`) is currently disabled and not
built in CI: ortools is now a required runtime dependency and is not packaged
for Fedora. The spec is retained, dormant, for future revival.

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
