# Ordered representation stability for GPC facets

Status: **the finite-generation route is CLOSED in the category the question
proposes, for a proved structural reason, and the two covariant transports it
was meant to organize are EXACT.** The negative half is a theorem about the
category, not a failed computation. The positive half is a 612-of-612 exact
inheritance law that the constructor can use immediately.

Pre-registration: `docs/prereg_ordered_representation_stability.md`, committed
with the measuring script at `1ff8259` before the run. Certifying artifacts:
`scripts/ordered_representation_stability.py`,
`results/data/ordered_representation_stability.json`,
`tests/test_ordered_representation_stability.py`.

## The question

Can primitive chamber-facing GPC facets and their proof objects be organized as
a finitely generated functor on a category of ordered mode embeddings, so that
high-rank facets are generated from finitely many bounded-rank templates?

## Answer

No, and the reason is not that the data is too small. In the fixed-`N` ordered
layer the category of mode embeddings **degenerates to a poset** (Theorem 1),
so every representation-stability theorem one could invoke becomes vacuous
there: finite generation stops being a statement about finitely many templates
plus symmetry and becomes the bare statement that the facet count is eventually
constant in `d`. The census refutes that for any generation degree up to 9.
The natural repair, restoring the symmetric-group action and working over `FI`,
is blocked by a second theorem: chamber-facing GPC rows are not permutation
stable (Theorem 2), and there is no `S_d` invariant convex body whose facets
they could be, because the Weyl saturation of the moment polytope is not convex.
Zero of the 901 nonstructural census rows survive a permutation.

What survives, and is worth keeping, is the inheritance law. The tight padding
extension and the tight frozen-core extension are total, are exact sections of
the two restrictions, compose, and carry published facets to published facets
in **612 of 612** nonstructural instances across all 18 map instances at the
nine non-degenerate consecutive pairs. Every one of the 514 inherited rows in
the census replays exactly from its primitive core.

## Part 1. The category

### Objects and rows

A row of the system at `(N,d)` is an affine inequality `a . lambda <= b` read on
the trace slice `sum(lambda) = N` inside the dominant chamber
`1 >= lambda_1 >= ... >= lambda_d >= 0`. The affine form carries the trace
gauge `(a,b) -> (a + t*1, b + t*N)`, so the canonical representative is the
primitive integral homogeneous vector `tau = prim(N*a - b*1)`, which is
`gpc_census.generation.bdr.tau_from_constraint`. Write `S(N,d)` for the ordered
slice, `V(N,d)` for the valid chamber-facing rows and `F(N,d)` for the published
irredundant system.

### Theorem 1 (chamber rigidity). PROVED.

Call `(N,d)` **non-degenerate** if `S(N,d)` is not contained in
`{lambda_d = 0}`, that is, if some `N`-fermion state in `d` modes has a
full-rank one-body density matrix. The hypothesis is necessary rather than
technical: if every state left mode `d` empty then `(N,d)` would be `(N,d-1)`
wearing an extra coordinate. It is verified exactly at all 15 census systems in
the artifact under `nondegeneracy`, by exhibiting an interior point of the
slice with every coordinate strictly positive.

*Let `1 <= N < d <= e` with `(N,d)` non-degenerate, and let `f : [d] -> [e]` be
an order-preserving injection such that inserting empty modes at the positions
outside `im(f)` carries `S(N,d)` into `S(N,e)`. Then `im(f) = {1, ..., d}`.*

**Proof.** Let `lambda` be in `S(N,d)` and let `mu` be its image, so
`mu_{f(i)} = lambda_i` and `mu_k = 0` for `k` outside `im(f)`. Suppose some `j`
outside `im(f)` has some `k > j` inside `im(f)`. Dominance of `mu` forces
`mu_k <= mu_j = 0`, hence `lambda_{f^{-1}(k)} = 0` for every `lambda` in
`S(N,d)`. Since `lambda` is dominant, that forces `lambda_d = 0` throughout
`S(N,d)`, contradicting non-degeneracy. So `im(f)` is an initial segment. QED

The proof deliberately does not average over the Weyl group. An earlier draft
argued that `S(N,d)` contains the uniform point because the moment polytope is
`S_d` stable and so contains its own barycenter. That is wrong, and Theorem 2
below is exactly why: the ordered slice is a fundamental domain, not an
invariant body, and the Weyl saturation is not convex, so there is no barycentre
argument to make. A positive interior point has to be exhibited, and it is.

**Consequence.** The fixed-`N` ordered category `C_N` has exactly one morphism
`[d] -> [e]` for `e >= d`. It is the poset of ranks, not `OI` and not `FI`.
Empty-orbital insertion is terminal only, and the frozen-core lift is the other
generator but it changes `N`, so it lives in the bigraded layer and is not an
endomorphism of `C_N`.

### What the variance is

`S(N,-)` is covariant: zero padding embeds `S(N,d)` onto the face
`{lambda_{d+1} = 0}` of `S(N,d+1)`, and the frozen-core lift
`lambda -> (1, lambda)` embeds `S(N,d)` onto the face `{lambda_1 = 1}` of
`S(N+1,d+1)`. Rows are functions on the slice, so they **pull back**: the row
functor is naturally contravariant. That is the first thing the brief asks to
be settled, and it is settled by identifying the two faces, which is why the
terminal restriction is total rather than merely observed:

- `rho : V(N,d+1) -> V(N,d)`, `tau -> prim(tau_1, ..., tau_d)`;
- `phi : V(N+1,d+1) -> V(N,d)`, `tau'_i = N*tau_i + tau_1` for `i >= 2`.

Measured on the census: 726 restrictions along `rho`, of which 724 are
nontrivial and every one valid, the other 2 collapsing to a constant row
(P1); 726 along `phi`, with 4 trivial. This is a theorem, and the measurement
is a regression test on the implementation rather than evidence.

### The covariant maps

A covariant map is what finite generation needs, and there is a canonical one
because a valid row admits a strongest valid extension:

- `pi : V(N,d) -> V(N,d+1)`, `pi(tau) = prim(q*tau, p)` where `p/q = t*` and
  `t*` is the largest last entry keeping the row valid on `S(N,d+1)`;
- `psi : V(N,d) -> V(N+1,d+1)`, prepend the smallest coefficient `c*` making
  `c*lambda_1 + a.lambda_{2..} <= b + c` valid on `S(N+1,d+1)`.

Both are defined at every census pair, with no undefined case. The naive
extension `a -> (a,0)`, which is what a coefficient-level reading of "padding"
suggests, is **gauge dependent and not total**: it is valid for only 10 of 31
rows at `(3,8) -> (3,9)` and 0 of 31 at `(5,8) -> (5,9)`. That control is in
the artifact under `naive_padding_control`, and it is the reason tightness is
the right definition.

## Part 2. The inheritance law. CONFIRMED, exact.

Across all nine census systems and their exact particle-hole duals:

| property | scope | result |
|---|---|---|
| `rho` and `phi` total on valid rows | 726 + 726 | every image valid |
| `rho(pi(tau)) = tau` and `phi(psi(tau)) = tau` | 620 | 620 exact |
| `pi` and `psi` carry facets to facets, nonstructural | 612 | **612 exact** |
| `pi . pi` valid and doubly restricting back | 86 | 86 exact |
| particle-hole involutive | 907 | 907 exact |
| particle-hole facet exact | 907 | 904, the 3 misses self-dual, below |
| inherited rows replayed from their core | 514 | **514 exact** |
| published rows that are genuine facets | 907 | **907 exact** |

The five apparent misses in the raw `P4` count are not counterexamples and are
listed individually in the artifact. Three are the structural chamber row
(`lambda_1 <= 1` and its particle-hole dual `lambda_d >= 0`), which the
published tables carry by a per-rank convention rather than by transport. Two
are the Borland-Dennis system `(3,6)`, which ships three equalities and is
therefore not full dimensional, so "facet" does not mean at `(3,6)` what it
means elsewhere. With those two families named and excluded in advance rather
than after inspection, preservation is exact.

The three particle-hole misses are the same two families again. `delta` is
involutive on all 907 rows without exception and lands on a published facet in
904 of them; the three that do not are `(3,6)`, which is degenerate, and the
self-dual systems `(4,8)` and `(5,10)`, where the dual of the structural row
`lambda_1 <= 1` is `lambda_d >= 0` and the published table carries only one of
the pair.

So `F(N,-)` really is a subfunctor of `V(N,-)` over the poset `C_N`, and the
inheritance is provable row by row, not a coefficient pattern. This is what the
brief asked for in Task 2, and the replay of all 514 inherited rows from their
cores is the acceptance criterion "reproduce every known inherited row from its
core".

## Part 3. What the counts say

Descent labels, decided in the fixed order STRUCTURAL, PADDING, FROZEN-CORE,
PRIMITIVE, each inherited label carrying a round-trip witness:

| system | rows | padding | frozen core | structural | primitive |
|---|---|---|---|---|---|
| (3,6) | 1 | 0 | 0 | 0 | 1 |
| (3,7) | 4 | 0 | 0 | 0 | 4 |
| (3,8) | 31 | 4 | 0 | 0 | 27 |
| (3,9) | 52 | 31 | 0 | 0 | 21 |
| (3,10) | 93 | 52 | 0 | 0 | 41 |
| (4,8) | 15 | 4 | 4 | 1 | 6 |
| (4,9) | 60 | 14 | 21 | 1 | 24 |
| (4,10) | 125 | 59 | 20 | 1 | 45 |
| (5,10) | 161 | 59 | 45 | 1 | 56 |

Primitive counts do not fall to zero anywhere, and at `N=3` they are
27, 21, 41 at ranks 8, 9, 10. The held-out test, generators from ranks 6 to 9
closed under `pi`, `psi` and `delta` with nothing added, covers 52 of 93 at
`(3,10)`, 80 of 125 at `(4,10)` and 104 of 161 at `(5,10)`.

That is the whole finite-generation question in the poset category: the
transition maps `pi` are injective, so `F(N,d) = pi(F(N,d-1))` would force
`|F(N,d)| = |F(N,d-1)|`. The census gives 4, 31, 52, 93 at `N=3`. Finite
generation with generation degree at most 9 is therefore **DISPROVED**. For
degree 10 or higher it is **OPEN**, and it is exactly equivalent to eventual
constancy of the facet count.

### The counts are facet counts, certified

That argument counts facets, so a redundant published row would inflate it and
the disproof would be worthless. The repository proves `P = O` at all nine
systems (`docs/census_inner_sandwich.md`), which says the published rows cut the
correct polytope, but it does not say they are irredundant, and no in-repo
result did. So it is certified here directly: a row is a facet exactly when the
vertices at which it is tight span an affine subspace one dimension below the
slice, which is an exact rational rank computation against the census vertex
lists.

**All 907 published rows, at all 15 systems, are facets.** No row is redundant
and none cuts a lower face. The slice dimensions come out at `d - 1` everywhere
except `(3,6)`, which is 3 because Borland-Dennis ships three equalities of
which two are independent of the trace. So `1, 4, 31, 52, 93` at `N=3`,
`15, 60, 125` at `N=4` and `161` at `(5,10)` are genuine facet counts, and the
41 rows at `(3,10)` outside `pi(F(3,9))` are 41 genuine new facets.

### Scorecard

| prediction | verdict | measured |
|---|---|---|
| P1 terminal restriction total | CONFIRMATORY-PASS | 726 of 726 |
| P2 non-terminal restriction fails | PASS | 3048 of 6258 pairs invalid |
| P3 tight extension is a section | PASS | 620 of 620 |
| P4 extension does not preserve facetness | PASS | 5 misses, all explained above |
| P5 primitive counts do not vanish | **FAIL** | 27, 21, 41: nonzero but not monotone |
| P6 width unbounded | PASS | 4, 6, 7, 8, 9 |
| P7 height unbounded | **CONFIRMATORY-FAIL** | spread 2, 3, 12, 12, 12 |
| P8a holdout (3,10) in 30 to 80 percent | PASS | 55.91 percent |
| P8b holdout (4,10) below 100 percent | PASS | 64.0 percent |
| P8c holdout (5,10) below 30 percent | **FAIL** | 64.6 percent |
| P9 composability | PASS | 86 of 86 |

Independent scope, per standing rule R2, is the primitive count in the table
above; the padding and frozen-core rows are functions of their donors and
confirm only the invariance of the derivation.

Three predictions failed and each failure is informative.

**P5** predicted the primitive count would be non-decreasing. It is not: 27 at
rank 8, then 21 at rank 9, then 41 at rank 10. The substantive clause, that the
count never reaches zero, holds. The monotonicity clause was a guess with no
mechanism behind it and it was wrong. The dip matters, because a monotone
sequence would have licensed an extrapolation and the real sequence does not.

**P7** predicted strictly growing coefficient height. Spread stabilizes at 12
across `(3,8)`, `(3,9)`, `(3,10)` while affine height goes 9, 9, 15, so neither
is monotone. This was flagged CONFIRMATORY in advance because the coefficient
ranges were already visible in the source tables, so its failure costs nothing
but it does remove a route: height is not the obstruction.

**P8c** predicted below 30 percent coverage at `(5,10)` and got 64.6 percent.
The prediction was written for a corpus without the derived dual systems; the
prereg was then extended to include them, which opened the `psi` route from
`(4,9)` and the `pi` route from the derived `(5,9)`. The prediction was not
updated to match, so it is scored as written and failed. This is a pre-registration
bookkeeping error on my part, not a surprise about the mathematics.

**P6 passed but is worthless as evidence**, and saying so is more useful than
banking it. Maximum width is 9 at rank 10, but width 9 is reached by inherited
rows just as much as by primitive ones: the artifact records
`max_width_primitive` and `max_width_inherited` as equal at every system. `pi`
raises width by construction, so unbounded width is a property of the transport
maps and not an obstruction to generation. The pre-registered rationale for P6
was wrong even though its prediction came true.

## Part 4. Why the symmetric-group repair fails

The only reason `FI` and `OI` module theory buys anything is the symmetric
group action: finitely many orbits of generators cover infinitely many objects.
Theorem 1 says the ordered category has no such action. The obvious repair is
to drop the ordering and work with the `S_d` stable object. That fails, twice.

### Theorem 2 (no permutation stability). PROVED, and CONFIRMED at the census.

*For `1 <= N < d`, the Weyl saturation of the moment polytope, that is
`{x : sort(x) in S(N,d)}`, is not convex.*

**Proof.** Let `x = (1^N, 0^(d-N))`, the spectrum of a Slater determinant, and
let `y` be `x` with coordinates `N` and `N+1` transposed, which is a permutation
of `x` and so is in the saturation. Their midpoint is
`m = (1^(N-1), 1/2, 1/2, 0^(d-N-1))`, already dominant. If `m` were attainable,
`lambda_1 = ... = lambda_{N-1} = 1` would force those `N-1` orbitals to be
occupied in every determinant of the state, so the state is those `N-1`
orbitals wedged with a one-particle state in the remaining `d-N+1` modes, whose
one-body spectrum is `(1, 0, ..., 0)`. The residual spectrum of `m` is
`(1/2, 1/2, 0, ...)`, so `m` is not attainable. QED

The artifact carries this certificate explicitly at all nine primary systems
under `weyl_saturation_nonconvexity`.

**Corollary.** There is no `S_d` invariant convex body whose facets are the GPC
rows, so "the orbit of a facet" is not a well-formed notion here.

Independently, and this is the sharper measurement, permutations do not
preserve validity at all. By the rearrangement inequality the worst permutation
of a row against a vertex pairs both sorted, so testing permutation stability
costs nothing. The result across all 15 systems: **0 of 901 nonstructural rows
are permutation stable.** The only six survivors in the whole census are the
single structural Pauli row `lambda_1 <= 1` at `(4,8)`, `(4,9)`, `(4,10)`,
`(5,9)`, `(5,10)` and `(6,10)`.

This is the same fact the repository already established from the screening
side, where `docs/orbit_canonical_flats.md` records that orbit reduction is
lossless for ENUMERATING candidate hyperplanes and that screening verdicts are
nevertheless not constant on orbits, and `docs/screening_orbit_structure.md`
measures it. The present result strengthens it: the orbit members of a
chamber-facing facet are not merely differently screened, they are not valid
inequalities.

### The consequence for coefficient-multiset statistics

Grouping census rows by the multiset of their `tau` entries compresses hard:
93 rows at `(3,10)` fall into 16 classes, 52 at `(3,9)` into 8, 161 at `(5,10)`
into 21. That looks like the orbit compression representation stability would
predict, and it is recorded in the artifact as
`coefficient_multiset_classes`. **It must not be read as an orbit count.** By
Theorem 2 the other members of those classes are invalid inequalities, so the
classes are a coefficient-pattern statistic with no group acting behind them.
Treating that 93-to-16 ratio as evidence of representation stability would be
the exact mistake this campaign was run to avoid.

## Part 5. The noetherianity criterion, hypotheses checked

The theorem worth invoking is the Gröbner-category noetherianity theorem of Sam
and Snowden ("Gröbner methods for representations of combinatorial categories"):
for a directed category `C` admitting a Gröbner order on each morphism poset,
representations of `C` over a noetherian ring are noetherian; `OI` is Gröbner
and `FI` is quasi-Gröbner, which is what makes `FI` modules finitely generated
in the useful sense. Checking its hypotheses against the three candidate
categories:

| category | directed | Gröbner | functor exists | verdict |
|---|---|---|---|---|
| `C_N` (ordered, fixed `N`) | yes | yes, trivially | yes, covariantly via `pi` | **hypotheses hold, conclusion vacuous** |
| `OI` | yes | yes | **no** | blocked by Theorem 1 |
| `FI` | yes | quasi | **no** | blocked by Theorem 2 |

The first row is the important one and it is the reason this route closes
rather than merely stalls. `C_N` is a poset, so `Rep(C_N)` is noetherian for
free, and a finitely generated object is precisely one whose transition maps
are eventually surjective. The theorem is true and tells us nothing we did not
already have to prove by hand, because there are no automorphisms to quotient
by and therefore no compression from "finitely many orbits". The whole content
of finite generation collapses onto a single unproved statement about facet
counts.

The other two rows fail at functoriality, not at noetherianity. There is no
`OI` module because non-terminal injections do not act, measured as 3048
invalid restriction pairs out of 6258. There is no `FI` module because `S_d`
does not act, measured as 0 of 901.

**Citation caveat.** The Sam-Snowden statement above is quoted from memory and
must be checked against the published paper before it appears in any manuscript.
Nothing in Parts 1 to 4 depends on it: Theorems 1 and 2 are proved here, and the
disproof of finite generation at generation degree at most 9 is a counting
argument on the census.

## Part 6. The strongest justified statements

**Theorem 1** (proved above). `C_N` is the poset of ranks.

**Theorem 2** (proved above). The Weyl saturation is not convex, and
chamber-facing GPC rows are not permutation stable.

**Inheritance law** (CONFIRMED, exact, ranks 6 to 10, `N = 3, 4, 5` and duals).
`pi(F(N,d-1))` and `psi(F(N-1,d-1))` are subsets of `F(N,d)` for every
non-degenerate census pair, `rho . pi = id`, `phi . psi = id`, and all 514
inherited rows replay from their cores. 612 of 612 nonstructural instances.

**Corollary** (DISPROVED at degree at most 9, OPEN above).
`M_N` is not generated in degree at most 9, since `|F(3,d)| = 4, 31, 52, 93`
strictly increases. Finite generation at any degree is equivalent to eventual
constancy of `|F(N,d)|`.

**Conjecture** (SPECULATIVE, stated so it can be attacked).
`|F(N,d)|` is unbounded in `d` for every fixed `N >= 3`, hence `M_N` is not
finitely generated over `C_N` at any degree. Evidence: 4, 31, 52, 93 at `N=3`
and 15, 60, 125 at `N=4`, with primitive counts 27, 21, 41 and 6, 24, 45. This
is a conjecture, not a finding; four and three data points are not a growth
theorem, and the sequence is not even monotone in its primitive part.

What is deliberately NOT claimed: nothing here says facets stabilize; nothing
here proves they do not, at ranks above 10; the rank-10 rows are conditional on
AK2008 completeness per `results/data/PROVENANCE.md`, though the disproof of
generation degree at most 8 survives without them, since `|F(3,9)| = 52 > 31`.

## Part 7. What this buys the constructor

One measured thing, worth taking. At every census rank the transports supply a
certified block of the next system for free, with a round-trip witness per row
and no candidate screening at all:

| target | free from transport | published | fraction |
|---|---|---|---|
| (3,9) | 31 | 52 | 60 percent |
| (3,10) | 52 | 93 | 56 percent |
| (4,9) | 35 | 60 | 58 percent |
| (4,10) | 79 | 125 | 63 percent |
| (5,10) | 104 | 161 | 65 percent |

That is a seeding layer, not a compression: it removes a constant fraction of
the work at each rank and leaves the asymptotics untouched. It is worth wiring
into the constructor precisely because it is exact and already certified, and
it is worth being clear that it does not touch the rank barrier.

## MERGE NOTES

Not yet applied to the shared ledgers beyond what this branch needs. On
integration:

- `docs/RESEARCH.md`: add under the census program, "Ordered representation
  stability: DISPROVED as a finite-generation route (Theorem 1, the ordered
  category is a poset; Theorem 2, no permutation stability); inheritance law
  CONFIRMED exact, 612 of 612
  [docs/ordered_representation_stability.md;
  results/data/ordered_representation_stability.json]".
- `AGENT_PLAN.md`: the seeding layer of Part 7 belongs in the T2 ladder beside
  the orbit-native constructor, as a candidate-free block, not as a new gate.
- `results/data/PROVENANCE.md` and `SHA256SUMS`: rows added on this branch.

## Reproducing

```sh
uv run scripts/ordered_representation_stability.py
uv run pytest tests/test_ordered_representation_stability.py
```

## Final answer

```text
QUESTION ANSWER:
  No. A precise category and functor exist and are constructed here, but the
  finite-generation theorem is false in the ordered category the question
  proposes, and the failure is structural. Fixing N and imposing the dominant
  chamber forces the category of ordered mode embeddings to be the poset of
  ranks (Theorem 1), which strips representation stability of everything that
  makes it useful: finite generation degenerates to eventual constancy of the
  facet count. Restoring the symmetric group is blocked by Theorem 2, since
  chamber-facing GPC rows are not permutation stable, 0 of 901 in the census,
  and the Weyl saturation of the moment polytope is not convex, so there is no
  invariant body whose facets they could be.

BEST THEOREM OR COUNTEREXAMPLE:
  Theorem 1 (chamber rigidity): the only order-preserving injection inducing a
  map of ordered slices by empty-mode insertion is terminal, so Hom([d],[e])
  is a single morphism and C_N is a poset. Theorem 2 (no permutation
  stability): for 1 <= N < d the point (1^(N-1), 1/2, 1/2, 0^(d-N-1)) is the
  midpoint of two points of the Weyl saturation and is not attainable, so the
  saturation is not convex. Both proved. The positive counterpart is the
  inheritance law: 612 of 612 nonstructural tight extensions of published
  facets are published facets, with 514 of 514 inherited rows replayed from
  their cores.

FINITE DATA RESULT:
  Nine census systems plus six exact duals: 542 primary rows and 365 derived,
  907 in all. Primitive counts at N=3
  are 1, 4, 27, 21, 41 for d = 6..10; at N=4, 4, 6, 24, 45; at N=5, 27, 24,
  56. Held-out rank 10 from rank-9 generators: 52/93, 80/125, 104/161. Nine
  pre-registered predictions scored, six PASS, three FAIL (P5 monotonicity,
  P7 height growth, P8c coverage window), all three failures recorded above
  with mechanism.

MISSING LEMMA:
  A growth or boundedness theorem for |F(N,d)| in d at fixed N. In the only
  category that survives, finite generation is exactly equivalent to eventual
  constancy of that count, so this single statement decides the question and
  nothing in the representation-stability toolkit supplies it.

RANK-20/RANK-40 IMPLICATION:
  None, in the compression sense. This route yields no rank-independent
  presentation. What it does yield at every rank is a certified transported
  block covering 56 to 65 percent of the next system with no candidate
  screening, which is a constant-factor saving and leaves the barrier where it
  was.

NEXT GO/NO-GO DECISION:
  NO-GO on ordered representation stability as posed, and on any variant that
  needs a symmetric-group action on chamber-facing rows. GO, small, on wiring
  pi and psi into the constructor as a certified seeding layer, since both are
  already exact and tested. The open mathematical question moves to the facet
  count |F(N,d)|, which is a different attack and should be scoped separately.
```
