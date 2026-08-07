# Screening is not a Weyl-orbit invariant

Status: **REFUTED**, with a structural decomposition that survives and is
useful. The hypothesis was that Ressayre screening could be run on one
representative per `S_d` orbit and the verdict transported. It cannot. What
survives is that the trace test splits into an orbit invariant and a
permutation statistic, which removes the dominant cost without removing the
work.

Certifying artifacts: `scripts/screening_orbit_structure.py`,
`results/data/screening_orbit_structure.json`,
`tests/test_screening_orbit_structure.py`.

## The claim that failed, and it was mine

`docs/orbit_canonical_flats.md` reduced candidate ENUMERATION from hyperplanes
to orbits, losslessly, and then said the screening cost "scales with the ORBIT
count only if screening is orbit invariant, which is true for validity but has
not been re-verified". The parenthetical was wrong. Screening is **not** orbit
invariant, and the doc has been corrected.

Measured at rank 7: **25 of the 35 oriented-tau orbits carry mixed verdicts**,
several containing both `certified` and `trace_rejected` in the same orbit. A
representative's verdict says nothing about its orbit.

## Why, exactly

The Ressayre trace test compares two counts:

    (a) negative roots (i < j) with h_j - h_i < 0
    (b) exterior weights strictly below the hyperplane

**Count (b) is `S_d` invariant.** Permuting modes permutes the weights and
"below" is equivariant, so the count belongs to the orbit. Verified on every
row at ranks 6, 7 and 8, with zero failures.

**Count (a) is exactly the INVERSION NUMBER of the sequence `h`**, verified
against `_negative_root_set` directly. Inversions are the archetypal
non-invariant statistic: on one rank-6 orbit the count ranges over all of
0 through 13 while `B` sits fixed at 10.

The underlying reason is not an implementation detail. The dominant chamber
`lambda_0 >= ... >= lambda_{d-1}` IS a fundamental domain, so it is precisely
the thing the Weyl group does not preserve. The invariance of the moment cone
is a statement about the ambient cone; the GPC system is a list of
chamber-facing inequalities, and a chamber-referenced condition cannot be an
orbit invariant. Orbit reduction is lossless for enumerating candidates and
meaningless for deciding them.

## What survives, and it is worth having

    trace test passes  <=>  inv(sigma h) == B,      B invariant per orbit.

So `B` is computed ONCE per orbit rather than once per row, and each row then
costs an inversion count rather than a sweep over all `C(d,3)` weights.

| system | oriented rows | orbits | `B` evaluations saved | trace survivors | survivor fraction |
|---|---|---|---|---|---|
| (3,6) | 724 | 15 | 48x | 64 | 8.84% |
| (3,7) | 10,682 | 35 | 305x | 549 | 5.14% |
| (3,8) | 332,840 | 109 | 3,053x | 6,607 | **1.98%** |

The 549 at rank 7 is an exact cross-check: `docs/fixed_n3_generator_methodology.md`
records 10,133 trace failures out of 10,682 rows, and 10,682 - 10,133 = 549.

**The survivor fraction is falling fast**, 8.8% to 5.1% to 2.0%, so the
expensive determinant stage shrinks as a share of the work exactly where it
matters.

## The consequence for the generator

The pipeline now stands as:

| stage | orbit reduction | status |
|---|---|---|
| enumerate admissible hyperplanes | **applies, lossless** | done to rank 10, 1,337 orbits |
| compute `B` per candidate | **applies**, `B` is an orbit invariant | done |
| trace filter | does not apply, but costs only `inv(h)` | done |
| determinant certification | does not apply | survivors only, ~2% and falling |
| facet reduction | not needed, systems are tens of rows | existing |

So enumeration is solved and screening is not reduced but is cheap per row and
sparse in survivors. The remaining engineering step is to GENERATE the
survivors directly instead of filtering them out of the full orbit: the
permutations of a multiset with a prescribed inversion number are a classical
object, counted by a Gaussian multinomial and generable in output-linear time.
Without that step every row must still be visited, which is what sets the
reachable rank.

## Attacking the determinant: most zeros are structural, not algebraic

The determinant stage is where the remaining cost sits, and it turns out most
of it is not algebra at all.

The tangent matrix has rows indexed by weights below and columns by roots, with
each cell either empty or a signed variable. Its determinant is the sum over
permutation terms, so **if the support bipartite graph has no perfect matching,
every term contains an empty cell and the determinant is identically zero for
structural reasons.** By Hall's theorem the obstruction is a set of rows whose
joint neighbourhood is strictly smaller, which is a compact certificate a
reader can check against the support alone.

| system | reach the determinant | structural zeros | algebraic zeros | nonzero |
|---|---|---|---|---|
| (3,7) | 549 | **474 (86.3%)** | 25 (4.6%) | 50 |
| (3,8) | 6,607 | **6,414 (97.1%)** | 55 (0.8%) | 138 |

At rank 7 the split 474 + 25 = 499 reproduces the methodology doc's
symbolic-zero count exactly, and 50 matches its 49 nonstructural Ressayre rows
plus 1 structural Pauli row.

**The structural share RISES with rank**, 86.3% to 97.1%, and the cost gap
widens sharply: at rank 8 the matching test takes **0.08 seconds for all 6,607
survivors** while the dynamic programming takes **224 seconds for the 193 that
genuinely need it**. The DP is exponential in the matrix size in the worst
case; matching is `O(E sqrt(V))`.

`tangent_determinant_is_structurally_zero` returns the verdict with its Hall
violator and `verify_hall_violator` replays it. Every one of the 474 structural
verdicts at rank 7 is cross-checked against the symbolic route, and every
violator replays.

**The refined pipeline**, in increasing cost per candidate:

1. trace test, `inv(sigma h) == B` with `B` per orbit, kills about 98%;
2. Hall matching, kills 97% of what remains, with a certificate;
3. one evaluation at an integer point: NONZERO is already a proof of
   nonvanishing, so certified rows never need the symbolic route;
4. symbolic expansion, needed only where a matching exists and every
   evaluation vanishes, which is the genuine algebraic-cancellation case:
   **55 of 6,607 at rank 8**.

## Symmetry still on the table, and not yet used

The stabilizer of a candidate in `S_d` permutes the weights below and the roots
simultaneously, so the tangent matrix is equivariant and its determinant
block-diagonalises over isotypic components. The stabilizers are large: over
the 56 rank-8 orbits their orders run 5040, 5040, 1440, 1440, 720, 720, 576,
240, ... down to 2, against `8! = 40320`. That is a genuine second attack on
the residual algebraic cases, and it is NOT implemented or verified here.

## What this does and does not license

**Does.** `B` as an orbit invariant, proved by equivariance and verified on
every row at three ranks. The identification of count (a) with the inversion
number, verified directly against the implementation. A measured survivor
fraction falling through 2%.

**Does not.** It does not rescue verdict transport, which is dead. It does not
reduce the number of CANDIDATES that must be visited: every oriented row still
has to be touched, which is what sets the reachable rank. The Hall filter
reduces the cost of deciding a candidate, not the number of them. The
stabilizer block-diagonalisation is named as available and is unimplemented.
And it does not license the earlier sentence in
`docs/orbit_canonical_flats.md` claiming orbit invariance "is true for
validity"; that sentence was wrong and is corrected there.

## Reproducing

```sh
python scripts/screening_orbit_structure.py --ranks 6 7 8
uv run pytest tests/test_screening_orbit_structure.py
```
