# Direct cocircuit discovery and the constructive signed-Chow decoder

**Status:** one reformulation and two theorems, all proved, plus an independent
low-rank enumeration backend and a decoder that closes the one operational gap
`docs/signed_chow_projection.md` left open. The reformulation is **not** adopted
as the production enumerator, and the measurement that says so is in this
document rather than left for someone to discover later.

Pre-registration: `docs/prereg_cocircuit_chow_decoder.md`. Certifying
artifacts: `src/gpc_census/generation/cocircuit.py`,
`src/gpc_census/generation/chow_decoder.py`,
`scripts/cocircuit_chow_decoder.py`,
`scripts/verify_cocircuit_chow_decoder_standalone.py`,
`results/data/cocircuit_chow_decoder.json`,
`tests/test_cocircuit_chow_decoder.py`.

## Where this came from, and what was checked before anything was written

An external bundle, `gpc_cocircuit_chow_experiment`, arrived with a
dependency-free C++ cocircuit prototype, a Python profile decoder, and CSVs of
trace survivors at ranks 5 to 8. It was replayed here in full before any
repository code was written: the binary was rebuilt from its source, all four
cocircuit runs and both decoder verifications were rerun, and the regenerated
CSVs were byte-compared against the shipped ones. Every number reproduced.

Two cross-checks tied the bundle to this repository rather than to itself.
**Both were run against a scratch copy of the bundle, which is not vendored
here, so they are PROVENANCE NOTES and not repository artifacts.** They are
recorded because they are why the rest of this work was undertaken, and every
claim the document actually rests on is re-established below by in-repo code.

- The bundle's `(3,7)` and `(3,8)` nonstructural trace-survivor sets came out
  set-identical to the ones this repository generates through
  `scripts/signed_chow_projection.py`, 547 and 6,605 normals, with zero shadow
  mismatches. Not the same counts: the same normals.
- Its `S_d` orbit gates 8, 19, 56, 231 and 1337 are exactly the orbit counts
  already recorded in `results/data/orbit_canonical_flats.json`. They are
  therefore regression targets, not new information, and this document treats
  them as such throughout.

Everything below was then reimplemented from the mathematics inside
`gpc_census`, sharing no code with the bundle, and the in-repo chain replaces
the provenance notes: the cocircuit backend agrees with the flat lattice on all
166,420 `(h, z)` pairs at rank 8, and the trace test is cross-computed in both
conventions, so the survivor sets agree with the repository's own route by
in-repo measurement rather than by an external comparison.

The reverse-search node counts are the sharpest evidence that the
reimplementation is the same algorithm as the prototype: they agree exactly,
and node counts depend on both prunings firing at exactly the same nodes.

## C1. Admissible-hyperplane enumeration is cocircuit enumeration

Let `Omega = {chi_I : I in binom([d], N)}` be the exterior weights and, for a
homogeneous normal `tau`, let `B_0(tau) = {omega in Omega : tau . omega = 0}`.

**Theorem C1.** `B_0(tau)` is admissible if and only if it is a maximal
coplanar subset of `Omega` of linear rank `d-1`, that is the zero set of a
cocircuit of the exterior-weight configuration.

*Proof.* Every weight satisfies `1 . chi_I = N`, which is nonzero, so the
affine hull of any subset of `Omega` misses the origin and
`linear_rank = affine_rank + 1`. The admissibility condition the generator
enforces is `affine_rank(B_0) = d-2`, hence `linear_rank(B_0) = d-1`, hence
`span(B_0)` is a linear hyperplane, and being contained in `tau^perp` it equals
it. No omitted weight lies in that span, since it would then annihilate `tau`
and be in `B_0`, so `B_0` is maximal. Conversely a maximal coplanar subset of
linear rank `d-1` has affine rank `d-2` and its primitive annihilator is
admissible. QED

The content is algorithmic, not geometric: the existing enumerators reach
codimension one THROUGH every lower rank, and C1 says a search is entitled to
go straight there.

`src/gpc_census/generation/cocircuit.py` does that, by lexicographic subset
reverse search with two prunings, a node whose omitted lexicographic
predecessor is already in the current span, and a node whose remaining suffix
cannot lift the span to rank `d-1`. Ranks are taken over the rank-preserving
prime `admissible_flats` already uses, and every emitted zero set is re-derived
from its exact integer normal before the enumeration returns, so the modular
arithmetic is a search device and not part of any claim.

## The populations, reproduced by a second algorithm

<!-- sync:cocircuit-populations:start -->
| System | Weights | Admissible hyperplanes | `S_d` orbits | Reverse-search nodes | Trace survivors | Nonstructural | Pairwise match | Trace test match |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `(3,5)` | 10 | 30 | 4 | 154 | 10 | 8 | yes | yes |
| `(3,6)` | 20 | 362 | 8 | 2775 | 64 | 62 | yes | yes |
| `(3,7)` | 35 | 5341 | 19 | 63136 | 549 | 547 | yes | yes |
| `(3,8)` | 56 | 166420 | 56 | 2042570 | 6607 | 6605 | yes | yes |
<!-- sync:cocircuit-populations:end -->

The `d = 6, 7, 8` hyperplane and orbit counts are sealed repository
populations, so this table is a REGRESSION with an independent sample size of
zero. What is new is the agreement being pair for pair: every cocircuit normal,
carried through the `h_i = S - d tau_i`, `z = N S` bridge, matches a
`(h, z)` pair of `enumerate_admissible_hyperplanes_by_flats`, and the two sets
of pairs are equal at every rank in the table. That is the check
`docs/orbit_canonical_flats.md` records as the one that catches a convention
error, the failure that once reported 15 orbits at rank 6 instead of 8. The
trace test is also cross-computed both ways, in the homogeneous coordinates the
cocircuit search produces and in the repository's centered convention, and the
survivor sets agree at every rank.

The orbit counts here are computed without the repository's canonicalizer at
all: an unoriented `S_d` orbit of a codimension-one flat is just the coordinate
multiset of its primitive normal, up to global sign.

## What C1 does NOT buy, measured rather than assumed

<!-- sync:cocircuit-cost:start -->
At `(3,8)`, one machine, all three reaching the same codimension-one objects:

| Enumerator | Objects it must hold | Seconds |
| --- | ---: | ---: |
| direct cocircuit reverse search | 2042570 search nodes, 166420 hyperplanes | 278.7 |
| materialising closed-flat lattice | 1369357 flats | 140.3 |
| orbit-canonical augmentation | 314 orbit representatives | 10.1 |

Flat lattice by rank: `[1, 56, 1540, 21420, 147630, 467082, 565208, 166420]`. Orbits by rank: `[1, 1, 3, 11, 37, 87, 118, 56]`, expanding to 166420 hyperplanes.

_Machine: Linux-6.18.5-fc-v20-x86_64-with-glibc2.39, Python 3.14.0rc2. Reported per machine, never averaged._

_Caveat, stated because it cuts against the conclusion: the cocircuit backend is a first implementation and the other two have been optimized repeatedly, so the seconds column flatters neither side fairly. The object counts do not depend on optimization effort, and they are the finding._
<!-- sync:cocircuit-cost:end -->

**The direct route is the slowest of the three, and it is the one that scales
worst.** The reason is not the constant factor, it is that the search is
unsymmetric: it emits all 166,420 rank-8 hyperplanes individually, where the
augmentation enumerator holds 314 orbit representatives across the whole
lattice, the empty flat included. At rank 9 that gap is 10,004,154 against 494
at the widest level, and at rank 10 it is 889,205,792 against 2,705. So:

> Direct cocircuit enumeration is an INDEPENDENT CROSS-CHECK at low rank. It is
> not a replacement for `orbit_canonical`, and the recommendation to stop
> improving the orbit route is declined.

That is a correction to the handoff this work came from, which recommended
treating the direct route as the preferred pipeline. The recommendation was
made against the materialising flat lattice, which the direct route does beat
in stored objects; it does not survive comparison with the orbit-native
enumerator the repository already had.

The timings are Python against Python on one machine, recorded in the artifact.
The external C++ prototype does rank 8 in about 13 seconds here, which is a
language difference and not an algorithmic one, and it does not change the
scaling argument: 889 million objects is not a memory budget in any language.

## T4. Equal degrees force genuine symmetry

Recall T1 from `docs/signed_chow_projection.md`: any family of fixed-weight
vectors summing to `c_+(tau)` equals `B_+(tau)`.

**Theorem T4.** If `c_+(tau)_i = c_+(tau)_j` then the transposition `pi` of
orbitals `i` and `j` satisfies `pi(B_+(tau)) = B_+(tau)`.

*Proof.* `pi` permutes the slice, so `pi(B_+)` is a family of fixed-weight
vectors, and its shadow is `pi(c_+)`. Since `c_+` is unchanged by swapping two
equal coordinates, `pi(c_+) = c_+`. T1 applied to `pi(B_+)` gives
`pi(B_+) = B_+`. QED

**Corollary.** Let the equal-degree blocks of `c_+` have sizes `m_1, ..., m_r`.
Then `B_+` is invariant under the full product `S_{m_1} x ... x S_{m_r}`, so
membership of a determinant depends only on its occupancy profile
`k = (k_1, ..., k_r)` with `sum k_a = N` and `0 <= k_a <= m_a`. The unknown is a
subset of the PROFILES, not a subset of the `binom(d,N)` determinants.

This is where the size collapses. At `(3,8)` the slice has 56 determinants and
the profile count never exceeded 41.

## T5. A block-constant separator exists, and it orders the blocks

**Theorem T5.** Let `G = S_{m_1} x ... x S_{m_r}` act as above and let
`bar_tau = (1/|G|) sum_{g in G} g tau`. Then `bar_tau` is constant on each
block, `{omega : bar_tau . omega > 0} = B_+(tau)`, and the block levels are
strictly decreasing in the order of decreasing `c_+`.

*Proof.* Constancy on blocks is immediate from averaging over `G`. For
`omega in B_+`, every `g^{-1} omega` lies in `B_+` by T4, so every term of
`bar_tau . omega = avg_g tau . (g^{-1} omega)` is strictly positive and the
average is positive. For `omega` outside `B_+`, the complement is `G`-invariant
too, so every term is nonpositive and the average is nonpositive. Hence
`bar_tau` separates exactly `B_+`. T2 applies to `bar_tau`, since T2's exchange
injection uses only the family it defines, so `c_+` order implies `bar_tau`
order; distinct blocks have distinct `c_+` by construction, so the levels are
strict. QED

**Corollary (upper ideal).** Write the blocks in decreasing `c_+` order with
levels `t_1 > ... > t_r`. Abel summation gives

```text
sum_a k_a t_a = sum_{j<r} (t_j - t_{j+1}) * (k_1 + ... + k_j) + t_r * N
```

with every `t_j - t_{j+1} > 0`, so the separator value is monotone in
prefix-sum (Gale) dominance. The selected profiles therefore form an UPPER
IDEAL in that order: including a profile includes every profile dominating it,
excluding one excludes every profile it dominates.

## The decoder

For a profile `k`, the number of determinants with that profile containing one
fixed orbital of block `a` is

```text
q_a(k) = C(m_a - 1, k_a - 1) * prod_{b != a} C(m_b, k_b).
```

With `x_k in {0,1}` recording whether the whole profile orbit lies in `B_+`,
the shadow gives one exact equation per block:

```text
sum_k q_a(k) x_k = c_a.
```

`gpc_census.generation.chow_decoder` solves that system in integers, with
upper-ideal propagation from T5 and componentwise attained/reachable bounds
that force or forbid a profile before branching.

**It is an exact decoder, not a heuristic.** Any profile set satisfying the
block equations induces a family whose shadow is `c_+`, because within a block
all vertices have equal degree by T4, and T1 admits exactly one family with a
given shadow. So the search cannot return a wrong answer, and a second solution
would contradict T1; the implementation raises rather than choosing, which
makes T1 a live assertion on every row rather than a comment.

Applying the same decoder to `c_-`, taking the complement for `B_0`, and
computing the exact rational annihilator recovers the primitive normal whenever
`span(B_0)` is a hyperplane, which is exactly the hypothesis of T3. The decoded
positive family fixes the orientation, and both sign families are re-derived
from the recovered normal before it is returned.

<!-- sync:cocircuit-decoder:start -->
| System | Signed rows | Decode failures | Primitive-`tau` failures | Max profile variables | Max states per side | Mean states per side |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `(3,6)` | 62 | 0 | 0 | 10 | 3 | 1.565 |
| `(3,7)` | 547 | 0 | 0 | 18 | 5 | 2.166 |
| `(3,8)` | 6605 | 0 | 0 | 41 | 37 | 6.536 |
<!-- sync:cocircuit-decoder:end -->

This closes the caveat `docs/signed_chow_projection.md` records verbatim as
"recovering `B_+` from a bare `c_+` is the Chow-parameters inversion problem
and is not solved here". It is now solved, constructively, on every
nonstructural trace survivor at ranks 6, 7 and 8.

## No complexity bound, and the measurement that shows why not

<!-- sync:cocircuit-stress:start -->
| System | Trials | Failures | Max search states | Max profile variables |
| --- | ---: | ---: | ---: | ---: |
| `(3,8)` | 500 | 0 | 53 | 56 |
| `(3,10)` | 300 | 0 | 137 | 120 |
| `(3,12)` | 200 | 0 | 253 | 175 |
| `(3,20)` | 50 | 0 | 1717 | 442 |
<!-- sync:cocircuit-stress:end -->

These are random threshold shadows, NOT admissible normals, so nothing in this
table says anything about `(3,20)` geometry; it measures only how the profile
search behaves as the block structure widens. The maximum state count grows
with rank, and it grows faster than the profile count, so
**no polynomial worst-case bound** may be claimed from this evidence. The
state-to-profile ratio moves from 0.95 at `d = 8` to 3.89 at `d = 20`, and four
points do not fix an asymptotic either way. The implementation also
precomputes the profile-dominance relation quadratically in the number of
profiles, which is fine at 442 profiles and would not be at high rank.

## What is not claimed

- **Facetness, redundancy, birationality, completeness.** Cocircuit output is
  discovery data. Every downstream exact gate is unchanged and still mandatory,
  and nothing here produces a facet of anything.
- **A production enumerator.** See the cost table. Rank 9 and above remain
  `orbit_canonical`'s job.
- **Symmetric cocircuit enumeration.** The search here is unsymmetric. Cocircuit
  reverse search modulo `S_d` is not implemented.
- **Decoder complexity.** Measured, not bounded. See above.
- **New information about the polytopes.** The populations reproduced are
  sealed. The decoder rows are the new part.

## The next gate, and an honest assessment of it

The handoff proposes symmetry-native cocircuit reverse search through
TOPCOM, citing Rambau, *Symmetric lexicographic symmetric-subset reverse search
for the enumeration of circuits, cocircuits, and triangulations up to
symmetry*, arXiv:2607.05967 (July 2026). That preprint is real and was
confirmed to exist; TOPCOM's download page carries `TOPCOM-1.2.0b`, a beta, so
"implemented in TOPCOM 1.2" should be read as implemented in a beta rather than
in a stable release.

The proposed blind gate is the orbit ladder 8, 19, 56, 231, 1337 at
`(3,6)..(3,10)`. **Those five numbers are already in
`results/data/orbit_canonical_flats.json`.** An external backend reproducing
them is a regression test of the external backend, not a new result, and the
first rank at which such a run would produce anything the repository does not
have is `(3,11)`, where `orbit_canonical` currently stops.

So the gate is worth running only as a means to `(3,11)`, and it carries costs
the repository's rules make explicit: TOPCOM is a third-party non-Python
dependency that CI cannot install under the current house rules, so it can only
ever be optional discovery infrastructure whose every output must be re-derived
locally. The cheaper experiment, and the one that does not add a dependency, is
to ask whether the canonical-augmentation enumerator's known cost
concentration, 78 percent of `(3,9)` time in the top level, can be attacked
directly. That is a decision for `AGENT_PLAN.md` T2 and it is recorded there.

The bundle's TOPCOM input emitter is deliberately NOT vendored. It would be a
script with no artifact, no guard test, and no way to exercise it in CI, which
is exactly the shape `scripts/audit_data_completeness.py` exists to keep out of
`results/data`. Emitting the exterior weights plus the `d-1` adjacent orbital
transpositions is a short exercise if the `(3,11)` run is ever authorised, and
the equivariance check it needs is the one already written in
`tests/test_orbit_canonical.py`.

## Verification commands

```sh
uv run scripts/cocircuit_chow_decoder.py
uv run scripts/cocircuit_chow_decoder.py --ranks 5 6 --skip-cost --skip-stress
uv run python scripts/verify_cocircuit_chow_decoder_standalone.py
uv run pytest tests/test_cocircuit_chow_decoder.py
```
