# The root-budget boundary layer: one half verified, one half refuted

**Status:** external claim, tested against the census before adoption. The
positive-layer theorem is confirmed and is a real reduction. The zero-layer
extension is refuted for non-regular normals, which is where essentially all
known mechanisms live.

**Artifact:** `results/data/root_budget_boundary.json`, written by
`scripts/root_budget_boundary.py`.

## The claim

For a dominant one-parameter subgroup `tau` satisfying the
Bulois-Denis-Ressayre grade condition, the positive weight family

```text
Omega_+(tau) = { I in binom([d],N) : sum_{i in I} tau_i > 0 }
```

is an ideal of the Grassmannian dominance order, and summing Condition (B) over
positive grades gives `|Omega_+(tau)| = |Phi(w)| <= binom(d,2)`. An ideal
contains the principal ideal of each member, so every positive weight obeys
`|down I| <= binom(d,2)`. Under the bijection between `N`-subsets and Young
diagrams in the `N x (d-N)` rectangle, `|down I|` is the number of subdiagrams,
so the layer is enumerable without materialising `binom(d,N)` weights.

## The positive half: verified, and it is a real reduction

The enumeration reproduces exactly, in about 7 seconds of standard library:

| `N` | ambient `binom(40,N)` | positive layer | zero-or-positive layer | frontier | max area |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 9,880 | 2,587 | 2,590 | 485 | 74 |
| 4 | 91,390 | 4,634 | 4,638 | 1,537 | 72 |
| 5 | 658,008 | 6,525 | 6,530 | 2,576 | 70 |
| 10 | 847,660,528 | 12,372 | 12,382 | 5,724 | 60 |
| 20 | 137,846,528,820 | **15,359** | 15,371 | 7,446 | 42 |

At `(20,40)` that is a reduction of the positive-side universe by a factor of
about nine million.

Tested against every surviving mechanism at ranks 7, 8 and 9, where the
repository knows the exact population, the positive half holds without
exception:

| System | Mechanisms | `Omega_+` is an ideal | `|Omega_+| <= binom(d,2)` | every positive `|down I| <= binom(d,2)` |
| --- | ---: | ---: | ---: | ---: |
| `(3,7)` | 25 | 25 | 25 | 25 |
| `(3,8)` | 64 | 64 | 64 | 64 |
| `(3,9)` | 191 | 191 | 191 | 191 |

**Regularity is not needed for this half**, which is worth stating because the
source memo assumes it throughout. Below a positive weight every weight is
positive: `J <= I` gives `tau.J >= tau.I > 0` whether or not `tau` has repeated
entries. The largest positive family observed is 36 at `(3,9)`, exactly
`binom(9,2)`, so the budget is attained and the bound is tight.

This also identifies the claim with something the repository already enforces.
`|Omega_+(tau)|` is the orbit invariant `B` of `screening_orbit_structure.md`,
and `B <= binom(d,2)` is exactly the condition that an orbit is not
trace-dead. The new content is not the bound; it is the Young-diagram
consequence, which converts the bound into an enumeration strategy.

## The zero half: refuted for non-regular normals

The memo continues: if `I` may be zero or positive then everything strictly
below it is positive, so `|down I| <= binom(d,2) + 1`.

That step needs `tau.J > tau.I` strictly, which requires `tau` **regular**, all
entries distinct. As soon as two entries coincide, a weight strictly below a
zero weight can itself be zero, and the bound has nothing to stand on.

On the census it fails:

| System | Mechanisms | Uniform bound holds | Worst `|down I|` at a zero weight | Bound |
| --- | ---: | ---: | ---: | ---: |
| `(3,7)` | 25 | 21 | 30 | 22 |
| `(3,8)` | 64 | 54 | 40 | 29 |
| `(3,9)` | 191 | 183 | 56 | 37 |

22 of 280 mechanisms violate it, with the worst case exceeding the bound by
more than half again.

### And regularity is nearly vacuous here

| System | Mechanisms | Regular | Max distinct `tau` levels (out of `d`) |
| --- | ---: | ---: | ---: |
| `(3,7)` | 25 | **0** | 5 |
| `(3,8)` | 64 | **0** | 7 |
| `(3,9)` | 191 | **3** | 9 |

Three of 280. The memo calls the regular case "the largest full-resolution
instance"; in the census it is a thin stratum, and the Levi-block degeneracy
that breaks the zero bound is the normal situation rather than the exception.
So the `(20,40)` zero-layer count of 15,371 describes a regime that the census
suggests is not where the facets are.

## A corrected bound, verified

For any weight with `tau.I >= 0`, everything below it is nonnegative, so
`down I` is contained in `Omega_+ union Omega_0` and

```text
|down I| <= binom(d,2) + |Omega_0|
```

This holds on **every** mechanism at all three ranks, 280 of 280, with tightest
slack 6, 9 and 9.

`|Omega_0|` is the zero-grade dimension the Levi-fusion audit already computes,
and it is available in closed form from the memo's own generating function,
`[y^N x^0] prod_i (1 + y x^{l_i})^{m_i}`, so it can be evaluated per Levi type
before any weight is touched.

The cost is that the corrected bound is **Levi-type dependent** rather than
diagram-only. The memo's version is powerful precisely because `|down I|` is a
function of the diagram alone and so supports one universal truncation; the
corrected version needs a budget computed per Levi type. That is still
implementable, and it degrades exactly where the zero-grade module is large.

## What this means for the proposed architecture

- The positive-layer truncation stands and is the strongest part of the memo.
  A candidate's positive side at `(20,40)` genuinely lives in 15,359 weights.
- Sign control does **not** follow at the stated strength. The zero weights
  carry the tangent-matrix variables, and the uniform bound that was supposed
  to fence them is false off the regular stratum. Any implementation should use
  the per-Levi-type budget, or work with occupancy profiles in the Levi
  quotient, which the memo gestures at but does not make precise.
- The memo's step 2, "reproduce `(3,8)`, `(3,9)`, `(3,10)`, `(4,8)`, `(5,10)`",
  remains the right gate, and this document supplies part of it early: the
  positive-side prediction is already checked against ranks 7 to 9.
- Claims resting on the BDR implementation stay unverified here.
  `docs/bdr_constructor_comparison.md` records that the pinned backend cannot
  be fetched or run in this environment and labels that benchmark INCOMPLETE;
  nothing in this document changes that.

## Provenance note on the uploaded bundle

The shipped `rank40_root_budget_boundary.json` does not byte-match its own
generator's output: the script emits `rectangle` and `area_distribution` fields
that the shipped artifact lacks. Every shared key agrees exactly, so this is a
trimmed copy rather than a different computation, but the bundle's SHA256SUMS
therefore pins a file the script does not reproduce.

## Verification

```sh
uv run scripts/root_budget_boundary.py --skip-census
uv run scripts/root_budget_boundary.py
uv run pytest tests/test_root_budget_boundary.py
```
