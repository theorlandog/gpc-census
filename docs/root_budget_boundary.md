# The root-budget boundary layer: refuted per weight, repaired per profile

**Status:** external claim, tested against the census before adoption. The
positive-layer theorem is confirmed and is a real reduction. The zero-layer
extension is refuted **as a per-weight bound**, and then repaired by the
weighted occupancy-profile formulation of the third memo, which holds 280 of
280. The gap this document opened is therefore closed, by a statement about
profiles rather than about individual determinants.

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

## The follow-up memo: per-Levi budget, Durfee rank, transpose orbits

A second memo (`frontier_orbital_stress.md`) proposes four refinements. All are
reproduced here; three are correct and one has a consequence its author does not
draw.

**Transpose orbits, confirmed exactly.** At half filling the sign-control
universe is closed under partition transpose, since the subdiagram count is
transpose invariant. Recomputed from scratch:

| System | Sign-control universe | Self-conjugate | Transpose orbits | Claimed |
| --- | ---: | ---: | ---: | ---: |
| `(10,20)` | 1,940 | 46 | 993 | 993 |
| `(15,30)` | 8,409 | 89 | 4,249 | 4,249 |
| `(20,40)` | 22,817 | 139 | 11,478 | 11,478 |

**The Durfee bound is correct, and it survives the zero-layer refutation.** A
diagram with Durfee square `r` contains the `r x r` square, so its subdiagram
count is at least `binom(2r,r)`, giving `binom(2r,r) <= budget`. That yields
`r <= 4, 4, 5, 5` at `(13,17)`, `(10,20)`, `(15,30)`, `(20,40)`. Crucially the
answer is the same whether the budget is `binom(d,2)` or `binom(d,2)+1`, so
this claim does not depend on the zero-layer bound refuted above. On the census
the measured maximum Durfee rank of a positive weight is 2 at all three ranks,
against a bound of 3, so the fixed-width Frobenius encoding is sound with room.

**The per-Levi budget is already in the repository.** The memo proposes
replacing `binom(d,2)` by `R_L = sum_{i<j} m_i m_j = (d^2 - sum m_i^2)/2`. That
is exactly `gpc_census.generation.orbit_canonical._max_inversions`, verified
identical on 400 random Levi types, and it is what `orbit_is_trace_dead`
already applies. So this is not a refinement to implement; it is the sharp form
the survivor generator has been using, which is independent evidence that the
budget is the right quantity.

**And the per-Levi budget cuts both ways, which the memo does not say.** It
strictly sharpens the positive half and strictly aggravates the zero half,
because a smaller budget is a stronger claim on both sides:

| System | `R_L` range | Positive family `<= R_L` | Positive weights `<= R_L` | Zero weights `<= R_L + 1` | (global budget was) |
| --- | --- | ---: | ---: | ---: | ---: |
| `(3,7)` | 6 to 19 | 25/25 | 25/25 | **13/25** | 21/25 |
| `(3,8)` | 7 to 27 | 64/64 | 64/64 | **43/64** | 54/64 |
| `(3,9)` | 8 to 36 | 191/191 | 191/191 | **154/191** | 183/191 |

The positive half survives the tightening everywhere. The zero half degrades
from 22 failures out of 280 to 70 out of 280. So adopting the per-Levi budget
is right, and it makes the sign-control claim worse rather than better.

## The profile repair: the zero half, done correctly

A third memo replaces the per-weight statement with one about Levi occupancy
profiles, and that is the right object. For levels `l_1 > ... > l_s` with
multiplicities `m_i`, a determinant is represented only by its profile `k`,
with

```text
M(k) = prod_i binom(m_i, k_i)        multiplicity
g(k) = sum_i l_i k_i                 grade
Q(k) = sum over l >= k of M(l)       weighted upper-ideal mass
```

where `>=` is prefix-sum dominance. Strict dominance strictly raises the grade,
so everything strictly above a nonnegative profile is strictly positive, giving

```text
g(k) >  0   =>   Q(k)        <= R_L
g(k) == 0   =>   Q(k) - M(k) <= R_L
```

**Verified on the census, 280 of 280, on both branches**, with no excess at any
mechanism. The profile-mass sum also reproduces the direct weight count of
`Omega_+` at every mechanism, so the compression is lossless.

The second law is exactly the repair. The per-weight version subtracted one
weight (the "+1"); the profile version subtracts the profile's whole
multiplicity `M(k)`. That is precisely the difference between a regular normal,
where a zero weight is alone in its class, and a degenerate one, where a single
zero profile can stand for an astronomical number of determinants. The zero law
holds at every mechanism where the per-weight bound failed, including all 70
that failed under the tightened per-Levi budget.

So the correct sequence is: the per-weight zero bound is false off the regular
stratum, the per-Levi tightening makes it fail more, and the profile
formulation makes it true. Only the last belongs in an implementation.

### The laws are conditional, and that is useful

They assume the root budget, so they are not unconditional facts about a Levi
type. Applied to `tau = (2,2,0,0,-1,-1,-1,-1)` at `N = 3` they fail, with
`|Omega_+| = 24` against `R_L = 20`. That is the correct verdict: no
arrangement of that multiset can pass the Ressayre trace test, so the normal is
trace-dead. A violation therefore certifies trace-death, which makes the laws
usable as a pruning test and not only as a bound. This is worth stating because
a test written against an arbitrary invented normal will fail for a legitimate
reason and look like a defect in the theorem.

## What this means for the proposed architecture

- The positive-layer truncation stands and is the strongest part of the memo.
  A candidate's positive side at `(20,40)` genuinely lives in 15,359 weights.
- Sign control **does** follow, once stated per profile rather than per weight.
  The per-weight route is false off the regular stratum and the per-Levi
  tightening widens that gap; the weighted-profile laws close it, verified
  280 of 280. The zero weights
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

## The Levi excitation bound: proved, and vacuous where we can test it

A fourth memo proposes bounding the excitation defect

```text
delta(k) = sum_i min(k_i, m_i - k_i)
```

of every positive occupancy profile by `floor(log2 R_L)`. The proof is correct.
Condition (B) forces `M(k) <= R_L` for each positive profile, and

```text
binom(m_i, k_i) = binom(m_i, t_i) >= binom(2 t_i, t_i) >= 2^{t_i},
    t_i = min(k_i, m_i - k_i)
```

so `M(k) >= 2^delta(k)` and hence `delta(k) <= floor(log2 R_L)`. Both
inequalities are verified here, the second on 20,000 random profiles.

**But the bound cannot be tested on anything this repository has.** Since
`min(k_i, m_i - k_i) <= k_i` and `sum_i k_i = N`,

```text
delta(k) <= N        unconditionally, for every profile
```

So the theorem says nothing whenever `N <= floor(log2 R_L)`. On the census:

| System | `N` | `floor(log2 R_L)` range | Mechanisms where the bound can bind | Violations | Max measured `delta` |
| --- | ---: | --- | ---: | ---: | ---: |
| `(3,7)` | 3 | 2 to 4 | **2 of 25** | 0 | 3 |
| `(3,8)` | 3 | 2 to 4 | **2 of 64** | 0 | 3 |
| `(3,9)` | 3 | 3 to 5 | **0 of 191** | 0 | 3 |

276 of 280 mechanisms are vacuous, and rank 9 is entirely vacuous. The four
non-vacuous cases pass, which is the total evidence available.

The measured maxima of 3 are the structural ceiling `delta <= N = 3`, not a
finding about facets. The same reading applies to the memo's own post-hoc check
on 83 mechanisms: those systems have `N <= 5`, so "maximum defect 3 against a
permitted 4 or 5" is consistent with the bound never having applied.

**Consequence for the proposed gate.** The memo's first gate is "audit the
defect bound on the complete rank-7 to rank-10 candidate populations". That
audit cannot falsify the theorem: it is vacuous by construction at those
particle numbers, and a pass would carry no information. At half filling the
bound first bites at `d = 14`, since `N = 7 > 6 = floor(log2 binom(14,2))`, so
the first testable systems are exactly the ones for which no data exists.

None of this disputes the theorem. It is correct, and at `(20,40)` it is a
genuine constraint, cutting a possible defect of 20 down to 9. It simply cannot
be validated below the frontier it was designed for, and the audit that was
proposed to validate it should be replaced or dropped rather than run and
reported as a pass.

## The gate ladder, pinned in the artifact

The proposed validation ladder's entry numbers are recomputed here rather than
quoted, so the ladder's premises are guarded in-repo:

| System | Effective | Ambient `binom(d,N)` | Sign-control universe | Compression | Durfee | Transpose orbits |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `(13,17)` | `(4,17)` | 2,380 | 609 | 3.9x | 4 | n/a |
| `(10,20)` | `(10,20)` | 184,756 | 1,940 | 95.2x | 4 | 993 |
| `(15,30)` | `(15,30)` | 155,117,520 | 8,409 | 18,446.6x | 5 | 4,249 |
| `(20,40)` | `(20,40)` | 137,846,528,820 | 22,817 | 6,041,395.8x | 5 | 11,478 |

`(13,17)` is computed as `(4,17)` by particle-hole normalization, and only the
half-filled rungs admit the internal transpose quotient. Every figure agrees
with the source memo.

What the ladder tests is the WEIGHT ORACLE, not the architecture. Gate A
exercises rectangular geometry, particle-hole transport and weighted nonregular
profiles, all of which are verified above. The Kostant inversion-diagram and
birationality stages, which carry the minimal-facet claim, are exercised by no
rung of this ladder as long as the pinned BDR backend cannot be run here.

## Levi-type and root-budget class counts

The proposal to bucket Levi types by root budget rests on counts that
reproduce exactly:

| `d` | Levi types, `p(d)` | Distinct `R_L` values |
| ---: | ---: | ---: |
| 17 | 297 | 83 |
| 20 | 627 | 118 |
| 30 | 5,604 | 284 |
| 40 | **37,338** | **538** |

So the 37,338 Levi types at `d = 40` really do collapse to 538 budget classes,
and one boundary view can be shared across every type in a bucket.

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
