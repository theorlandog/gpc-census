# Orbit-closure moment polytopes of Lambda^3 C^7

Status: 10 of 13 class polytopes CERTIFIED EXACT (4 of 4 at d = 6, 6 of 9 at
d = 7); 3 shipped as certified inner/outer brackets. Emitter
`scripts/lambda3_c7_polytopes.py`, artifact
`results/data/lambda3_c7_polytopes.json`, library `gpc_census.orbit`, guards
`tests/test_lambda3_c7_polytopes.py`.

## What this computes, and how it differs from the census

The census layer answers: which occupation spectra arise from SOME N-fermion
pure state. That is the ambient moment polytope, and for (3,7) it has the ten
vertices in `results/data/vertices/vertices_3_7.json`.

This layer answers the finer question, one polytope per SLOCC class:

    Delta(psi) = { sorted spec(rho_phi) : phi in closure(GL_d . psi) }

Delta(psi) is a convex polytope (Kirwan, Brion, Ness) depending only on the
orbit. It is the fermionic entanglement polytope of Walter, Doran, Gross and
Christandl (Science 340, 1205 (2013)), who computed the four classes of
Lambda^3 C^6. The point of the object: a measured one-body spectrum outside
Delta(psi) certifies that the state is NOT in the closure of psi's class, so
these polytopes are one-body witnesses of genuinely multipartite fermionic
entanglement.

## The classification, recomputed rather than assumed

Lambda^3 C^7 has finitely many GL_7 orbits (Schouten 1931). The nine nonzero
orbits are recovered here computationally, and they have nine DISTINCT orbit
dimensions, so the single integer `orbit_dimension` is a complete orbit
invariant at d = 7. That is what `gpc_census.orbit.classify` relies on.

| class | orbit dim | rank | representative | polytope dim | verdict |
|-------|-----------|------|----------------|--------------|---------|
| I    | 13 | 3 | e123                      | 0 | EXACT |
| II   | 20 | 5 | e123+e145                 | 1 | EXACT |
| III  | 21 | 7 | e123+e145+e167            | 2 | EXACT |
| IV   | 25 | 6 | e123+e145+e246            | 3 | EXACT |
| V    | 26 | 6 | e123+e456                 | 3 | EXACT |
| VI   | 28 | 7 | e123+e145+e167+e246       | 6 | BRACKET |
| VII  | 31 | 7 | e123+e456+e147            | 6 | BRACKET |
| VIII | 34 | 7 | e123+e456+e147+e257       | 6 | BRACKET |
| IX   | 35 | 7 | e123+e456+e147+e257+e367  | 6 | EXACT |

Four classes have rank <= 6 (I, II, IV, V): these are the inherited
Lambda^3 C^6 classes, and their polytopes are the d = 6 answers padded with a
zero, because a rank-r state and its whole orbit live in Lambda^3 W with
dim W = r. Five classes are genuinely full-rank-7 (III, VI, VII, VIII, IX).

Note that the full-rank classes are NOT an upper interval in the dimension
order: class III sits at orbit dimension 21, below the rank-6 classes IV and
V at 25 and 26. Class III is `e_1 ^ omega` with omega a rank-6 two-form, so
it is full rank while being a small orbit.

## Certificates

Both bounds are exact over the rationals. Nothing is sampled, and no verdict
is upgraded past what its certificate proves.

### Inner: attainability

If a state's support is ONE-HOP FREE (no two support triples share two
modes), its 1-RDM is diagonal, and stays diagonal across the whole positive
torus orbit. The moment image of the closure of a torus orbit is exactly the
Newton polytope conv{chi_S}, so every sorted point of that Newton polytope is
an attained spectrum, realised by the explicit state `sum_S sqrt(w_S) e_S`.

The enumeration of one-hop-free supports is COMPLETE, not sampled, and the
reason is worth recording: every one-hop-free support has vanishing
coefficient moduli (`|support| = rank` of its incidence matrix), so the torus
acts transitively on nonzero coefficient vectors over it. The all-ones
trivector is therefore the only state supported there up to the torus, and
the support alone pins down the orbit. There are 270 such supports at d = 6
and 5595 at d = 7, the largest being the seven-triple Fano plane; up to
relabelling of the modes they collapse to 13 canonical supports, and sorted
spectra do not see a relabelling.

Certified degenerations contribute too: a class in the orbit closure of
another contributes its entire polytope. Degenerations are produced
constructively, each an explicit one-parameter limit (Newton faces, which are
torus limits, and `pi(h.psi)`, the limit of `diag(1,...,1,t).h.psi`).

### Outer: validity

Let `c` be a decreasing weight and `phi` any state of the orbit written in a
basis where it has a given support. Then

    <c, diag(rho_phi)> = sum_S |phi_S|^2 <c, chi_S> >= min_S <c, chi_S>,

while Schur-Horn majorisation plus Schur-convexity of `<c, .>` for decreasing
`c` gives `<c, diag(rho_phi)> <= <c, lambda>`. Chaining the two yields

    <c, lambda> >= min_S <c, chi_S>

valid on the orbit, and on the closure by continuity. Taking
`c = (1^m, 0^(d-m))` recovers the Ky Fan / coordinate-subspace form
`lambda_1 + ... + lambda_m >= j` with `j = min_S |S and A|`. General
decreasing `c` is strictly stronger, because the right-hand side is concave
in `c` and so is not implied by its values on the extreme rays. Relabelling
the modes changes the right-hand side, so all supports of the class are
scanned and the strongest bound per `c` is kept.

Inequalities are read only off supports of the class ITSELF. A support
belonging to a degeneration bounds the smaller polytope, not the larger one;
mixing the two would be unsound.

### The sandwich

The verdict is EXACT when every vertex of the outer H-polytope is a convex
combination of certified attainable points, checked by exact rational LP.
Then inner = Delta = outer and the polytope is pinned. Otherwise the result
ships as a BRACKET with the unreached outer vertices listed explicitly.

## Validation gates (all green)

- The exact polytope tools reproduce the published ambient (3,7) vertex set
  from its H-representation.
- All four Lambda^3 C^6 classes come out EXACT, and the OPEN ORBIT reproduces
  the ambient (3,6) polytope exactly. This is the decisive check: the generic
  class must have the whole ambient polytope, and it does.
- At d = 7 the open orbit reproduces the ambient (3,7) polytope exactly, all
  ten vertices.
- The Slater spectrum (1,1,1,0,0,0,0) lies in every class polytope, as it
  must: the Grassmannian is the unique closed orbit in P(V) and so sits in
  every orbit closure.
- Every certified degeneration implies the corresponding polytope inclusion.
- A rank-r class occupies at most r orbitals.

Two bugs were caught by these gates and are worth recording, since both would
have shipped silently:

1. `ambient_rows` initially dropped the constraint system's EQUALITIES. The
   (3,6) system has three (the Borland-Dennis pairings), so every d = 6 orbit
   polytope came out inflated, with vertices that are not ambient-feasible.
2. The seed vertex enumeration counted equality ROWS instead of their RANK.
   The (3,6) trace equality is implied by the three pairing equalities, so the
   tight-subset size was one too small, the vertex set came out EMPTY, and
   every class reported a FALSE EXACT (an empty outer polytope has no
   unreached vertices). This is exactly the failure mode the validation law
   exists to catch.

## Results

Class III (orbit dimension 21) has affine hull `l1 = 1`, `l2 = l3`,
`l4 = l5`, `l2 + l4 + l6 = 1`, `l2 + l4 + l7 = 1`: a frozen mode plus the
even-pairing degeneracy of a two-fermion six-orbital system, which is exactly
what `e_1 ^ omega` should give. Its polytope is the triangle on
(1,1,1,0,0,0,0), (1,1/2,1/2,1/2,1/2,0,0) and (1,1/3,1/3,1/3,1/3,1/3,1/3),
i.e. the whole `l1 = 1` face of the ambient polytope.

Classes IV and V, both rank 6, share an affine hull (the Borland-Dennis
equalities plus `l7 = 0`) but differ by one facet: the tangent class IV
carries `3l1 + 3l2 + 3l3 + 2l4 + 2l5 + 2l6 >= 8`, which the generic rank-6
class V does not. So one-body data separates the two rank-6 classes.

No two classes share a polytope: the containment order computed from the
vertex sets is strict throughout, so at d = 7 the one-body spectrum
distinguishes every SLOCC class from every other, at least to the resolution
of these bounds (the three bracketed classes are reported by their outer
polytopes in that comparison).

Class VI's outer description contains `2l1 + l2 + l3 + l4 >= 3`, a genuinely
new entanglement-polytope inequality of the general dominant-weight form that
the coordinate-subspace family cannot express.

## What is open, and which side the gap is on

Classes VI, VII and VIII ship as brackets. Sweeping the dominant-weight bound
from 3 to 7 changes nothing, so the Schur-Horn/Newton inequality family is
SATURATED: it is exactly tight for six of the nine classes and provably
cannot close the other three. Closing them needs a different tool: Ressayre
inequalities FOR THE ORBIT CLOSURE, or covariant nonvanishing.

RECONCILIATION (2026-08), read this before citing the sentence above. The
repository now also carries `docs/stage1_ressayre_3_7.md` and
`results/data/stage1_ressayre_3_7.json`, a first-principles Ressayre
derivation at (3,7). That manifest does NOT close the three bracketed
classes, and the two results must not be run together:

- The manifest computes the AMBIENT cone. Its four retained rows agree, as an
  oriented set, with the Altunbulak-Klyachko table this layer already consumes
  through `constraints(3,7)`, and feeding those rows alone through
  `gpc_census.exact_polytope` reproduces the ten published (3,7) vertices.
  That is a genuine independent cross-check of the ambient input, and it is
  pinned by `test_ressayre_manifest_agrees_with_ambient_rows`.
- The ambient cone is the polytope of the OPEN orbit, which is class IX, and
  class IX is already EXACT here. So the manifest confirms the one class that
  needed no help and says nothing about the other eight.
- What classes VI, VII and VIII need is Ressayre's theorem applied to a fixed
  orbit closure (the Berenstein-Sjamaar / Ressayre inequalities for a
  subvariety of P(V)), which is a strictly harder object than the ambient
  cone: the ambient case is the special instance where the subvariety is all
  of P(V). None of that machinery is in the repository.

Numerical probing of the unreached outer vertices (Powell descent on
`dist(sorted spec(rho_{g.psi}), target)^2` over g in GL_7, 30 restarts) shows
the gap is on BOTH sides depending on the vertex, and this is evidence, not
certificate:

- Positive floors, so the OUTER bound is loose there: class VI at
  (3/4,3/4,1/2,1/4,1/4,1/4,1/4) floor 3.8e-2, at
  (9/13,9/13,9/13,3/13,3/13,3/13,3/13) floor 3.3e-2, at
  (2/3,2/3,1/2,1/2,1/3,1/3,0) floor 1.7e-2; class VII at
  (7/11,...,3/11) floor 1.9e-2 and (2/3,1/2,1/2,1/2,1/2,1/6,1/6) floor
  1.6e-2.
- Floors at 1e-21, so the point IS attained and the INNER generator is
  incomplete there: class VII at (3/5,3/5,3/5,2/5,2/5,1/5,1/5) and class VIII
  at (5/9,5/9,5/9,1/3,1/3,1/3,1/3). Both are attained by orbit points whose
  1-RDM is not diagonal in any coordinate basis, which is precisely what the
  one-hop-free generator cannot reach. `gpc_census.orbit.rational_spectrum`
  reads exact spectra off the characteristic polynomial and does not need a
  diagonal 1-RDM, so it is the hook for closing these two; the search for an
  exact attaining state has not been done.

Class VIII is one vertex away from EXACT, and that vertex is numerically
attained, so class VIII is the cheapest remaining target.

One structural question stayed open and is recorded rather than guessed:
whether class V (orbit dimension 26) lies in the closure of class VI (orbit
dimension 28). The semicontinuous rank invariants do not obstruct it, and no
explicit degeneration was found: the generic rank-6 degeneration
`pi(h.psi_28)` lands in class IV (dimension 25) for every `h` tried, while
the same construction reaches class V from classes VII, VIII and IX. Since
orbit dimension is lower semicontinuous in `h`, the generic value is the
maximum, which is evidence that class V is NOT in the closure of class VI,
but that argument covers only degenerations of the form
`lim diag(1,...,1,t) h psi` and is not a proof.

## Reproducing

```sh
python scripts/lambda3_c7_polytopes.py     # ~30 s, writes the artifact
uv run pytest tests/test_lambda3_c7_polytopes.py
```
