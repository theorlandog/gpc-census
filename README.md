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

A narrative walkthrough is on the blog:
https://nutfieldsecurity.com/posts/ai-vacation-physics/

## The classification

Each polytope vertex is an extremal spectrum, and the census asks how it
can be attained. Exactly one of three things is true:

- DESIGN-INT: a superposition of Slater determinants with nonnegative
  integer weights at the spectrum's natural denominator attains it, with
  a one-hop-free support (no two determinants share N-1 orbitals), so no
  phase cancellation is needed. A combinatorial design, hence the name.
- DESIGN-REAL: real nonnegative weights suffice but integers do not.
- INTERFERENCE: no nonnegative weighting exists at all; phase
  cancellation is forced (signed real amplitudes can suffice, or complex
  ones are needed). The canonical example is the vertex
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
      needs pi/8); the small tail whose phases are higher-degree algebraic
      (no real, no p*sqrt(q)/r) is closed by a constructive solver that stops
      guessing phases and solves the pinned variables instead. The moduli are
      exactly rational, so wherever two orbitals share an occupation the
      spectrum must split, the off-diagonal 1-RDM magnitude is forced by
      Schur-Horn to an exact algebraic number; each off-diagonal is a closed
      polygon whose relative phase is then an exact arccos, propagated across
      edges and gated by the same characteristic-polynomial certificate. This
      is v_B's cos(gamma) = 3/(4 sqrt(14)) as a general rule, and it certifies
      every interference corner whose support the solver finds
      (`exactify_interference`, which certifies v_B itself). Of the 799 census
      vertices, all 799 now carry a certified closed form; the last two rank-10
      corners were closed off the rational search grid after the attainability
      audit self-corrected (see docs/RESEARCH.md).
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

The CLI has two modes: serve the precomputed results (fast, no solve) and run
the engine to recompute or extend them.

```sh
uv sync                                     # install deps (CP-SAT backend included)

# precomputed results, machine-readable (usable as a data source / library)
gpc-census export   -n 4 -d 9               # constraints + vertices + verdicts + states, JSON
gpc-census export   -n 4 -d 9 --kind states # just the certified closed-form states
gpc-census states   -n 4 -d 9 --index 65    # the closed form for one vertex (v_B)

# run the engine
gpc-census constraints -n 3 -d 9            # the 52-inequality system
gpc-census polytope    -n 3 -d 9            # its 58 vertices, exactly
gpc-census classify    -n 4 -d 9            # verdicts incl. v_A and v_B
gpc-census solve       -n 4 -d 9 --den 23 --spectrum 20,14,14,14,14,4,4,4,4  # certify a state
gpc-census states      -n 4 -d 9 --source hybrid   # lookup where available, solve the rest
gpc-census states      -n 4 -d 9 --source solve --max-cliques 0   # solve everything, push bounds
```

As a library, the precomputed dataset is available without recomputing:

```python
from gpc_census import dataset
dataset.vertices(4, 9)        # extremal spectra (exact)
dataset.classification(4, 9)  # design/interference verdict per vertex
dataset.states(4, 9)          # certified closed-form states from the lookup
dataset.export(4, 9)          # all of the above, one object

# the mode is the provenance: precompute serves the lookup, solve recomputes
# independently (ignoring the lookup), hybrid serves the lookup and solves gaps
dataset.resolve_states(4, 9, mode="hybrid")
```

### What's precomputed

`results/data/states.jsonl` ships a certified closed-form state for every
vertex the campaign has closed so far, so `gpc-census states --source
precompute` (the default) is a lookup, not a solve. Current coverage:

| System (N, d) | Corners | Design | Interference | Certified closed forms |
| --- | ---: | ---: | ---: | ---: |
| (3, 6) | 4 | 4 | 0 | 4 |
| (3, 7) | 10 | 10 | 0 | 10 |
| (3, 8) | 38 | 27 | 11 | 38 |
| (4, 8) | 22 | 22 | 0 | 22 |
| (3, 9) | 58 | 38 | 20 | 58 |
| (4, 9) | 103 | 87 | 16 | 103 |
| (3, 10) | 113 | 71 | 42 | 113 |
| (4, 10) | 159 | 134 | 25 | 159 |
| (5, 10) | 292 | 250 | 42 | 292 |
| **Total** | **799** | **643** | **156** | **799** |

Every design corner (integer and real) is certified by construction. All
156 interference corners are now certified as well; the final two rank-10
corners were closed off the rational search grid by an exact reconstruction
after the attainability audit self-corrected (`exactify_interference`, the
off-diagonal-target constructive solver, certifies the rest).
`scripts/build_states.py --retry-uncertified` closes that gap with more
compute, and `scripts/transport_states.py` closes more of it for free by
transporting a certified sibling's state along the padding, frozen-core lift,
and particle-hole maps (14 vertices so far); the engine call is below.

To recompute or push further, the engine is a single call
(`gpc_census.states.certify_state(n, d, spectrum, max_clique=..., max_cliques=...)`),
and the full campaign (`scripts/solve_all.py --all`) routes each vertex by
verdict, is restartable and checkpointed, and writes `results/data/states.jsonl`.
The compute knobs (`--max-card`, `--max-blocks`, `--max-clique`, `--max-cliques`)
bound the ansatz search, so more compute reaches further vertices.

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
