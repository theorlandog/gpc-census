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

## What this does and does not license

**Does.** `B` as an orbit invariant, proved by equivariance and verified on
every row at three ranks. The identification of count (a) with the inversion
number, verified directly against the implementation. A measured survivor
fraction falling through 2%.

**Does not.** It does not rescue verdict transport, which is dead. It does not
reduce the number of determinant certifications. And it does not license the
earlier sentence in `docs/orbit_canonical_flats.md` claiming orbit invariance
"is true for validity"; that sentence was wrong and is corrected there.

## Reproducing

```sh
python scripts/screening_orbit_structure.py --ranks 6 7 8
uv run pytest tests/test_screening_orbit_structure.py
```
