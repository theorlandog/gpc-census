# Cross-realization invariants of the census (2026-08-07)

Pre-registration: `docs/prereg_symplectic_invariants.md`, written and committed
before either measuring script was run. Artifacts:
`results/data/symplectic_invariants.json` (799 rows, exact),
`results/data/quantization_multiplicities.json` (14 vertices, exact).
Generators: `scripts/symplectic_invariant_census.py`,
`scripts/quantization_multiplicity.py`. Guard test:
`tests/test_symplectic_invariants.py`. Standalone replay with no project
import: `scripts/verify_symplectic_invariants_standalone.py`.

## What this campaign is

The census already describes each extremal object once, through its occupation
spectrum, its support combinatorics, and the facets it saturates. This campaign
computes the same 799 objects in three OTHER presentations and asks which
quantities agree across them.

- LATTICE. The tangent cone `T_v P` at each vertex, inside the weight lattice
  restricted to the affine hull: primitive edge rays, Smith invariant factors,
  lattice index, unimodularity.
- SYMPLECTIC/LIE. The exact projective stabilizer `Stab_{U(d)}([psi])` of the
  certified state, its orbit dimension, and the commutant of the state's own
  1-RDM with the dimension of its orbit through the state.
- REPRESENTATION. The multiplicity of highest weight `k*lambda` in
  `Sym^k(wedge^n C^d)`, which is the quantization of the same moment polytope.

Everything below is exact. Ranks are computed over the algebraic number field
the amplitudes generate; lattice data is integral; multiplicities are integer
counts.

## Scorecard

| prediction | verdict | scope | independent scope |
| --- | --- | --- | --- |
| P1 exactly 645 diagonal 1-RDMs | PASS | 799 | 411 classes |
| P2 Slater vertices are Grassmannian | PASS | 9 | 9 |
| P3 stabilizer never trivial | PASS | 799 | 411 classes |
| P4 toral stabilizer is combinatorial | PASS, 0 violations | 799 | 411 classes |
| P5 cone types | PASS | 799 | 411 classes |
| P6 vertices are multiplicity free | **FAIL**, 4 of 14 | 14 | 14 |
| P7 stabilizer constant on support classes | **FAIL**, 67 of 411 | 799 | 411 classes |

Independent scope is the number of distinct values of the support invariant,
which can merge non-isomorphic supports and therefore understates independence.
The nine systems are not independent of each other: `(4,8)`, `(4,9)`, `(4,10)`
and `(5,10)` reach the census partly through particle-hole duality and lifts,
and both `dim Stab` and the cone type are duality-covariant. The prereg says so
before the run.

## Finding 1. The support determines the toral symmetry, exactly

An orbital phase `phi` acts by `c_S -> exp(i sum_{o in S} phi_o) c_S`, so it
fixes `[psi]` exactly when `A phi` is constant on the support, where `A` is the
0/1 support-incidence matrix. That predicts

    dim (Stab_{U(d)}([psi]) intersect u(1)^d) = d - rank([A_S - A_{S_0}]_S),

a number computed from the discrete support alone, with no amplitude input.

Measured against the exact toral stabilizer of the certified state: **799 of
799, zero violations**. In the determination matrix,
`support_invariant -> toral_stab_dim` is exact over 411 classes and is not
vacuous.

This is the cross-realization statement the campaign was built to test, in the
one place it holds cleanly: a combinatorial invariant equals a symplectic one.

## Finding 2. It does NOT determine the full symmetry, and this is where it breaks

`support_invariant -> stab_dim` fails on **67 of 411 classes**. Supports with
identical invariants carry states whose full projective stabilizers differ, for
example `[6, 10]` within `(3,7)` and `[7, 9]`, `[4, 6]`, `[3, 5]` within
`(3,8)`.

So the honest boundary is: the support fixes the torus, and the amplitudes
carry symmetry beyond it. Across the census `dim Stab` exceeds its toral part
at 445 of 799 vertices, with the excess reaching 16.

Two caveats bound this negative result rather than dissolving it. The support
invariant is not a certified canonical form, so a split class could in
principle hold two genuinely non-isomorphic supports; and 67 splits out of 411
is a lower bound on the failure rate, since collisions can only hide splits,
never manufacture them.

## Finding 3. The GPC polytopes are toric-singular at every vertex above rank 6

- Simplicial tangent cones: **183 of 799**, exactly the simple vertices, which
  independently reproduces the `simple_vertices` column of
  `docs/polytope_invariants.json` system by system.
- Unimodular tangent cones: **3**, and all three are in `(3,6)`, namely the
  Slater vertex, `(2,1,1,1,1,0)/2`, and `(1,1,1,1,1,1)/2`.

From rank 7 up, **no** vertex of any census moment polytope has a unimodular
tangent cone. Lattice indices are large and grow fast: 16807 at `(3,7)`, and a
maximum of `37^8 = 3512479453921` at `(5,10)#144`, whose spectrum
`(34,34,34,17,17,17,8,8,8,8)/37` has denominator 37 and whose invariant factors
are `[1, 37, 37, 37, 37, 37, 37, 37, 37]`.

That last coincidence does NOT generalize. `index = q^(dim-1)` holds at only 29
of the 183 simplicial vertices, so it is a property of the extremes and not a
law; it is recorded here because it looks like one.

Verdict classes split unevenly across cone type: all **12 of 12** DESIGN-REAL
vertices are simplicial, against 120 of 631 DESIGN-INT and 51 of 156
INTERFERENCE. The interference figure agrees with the archive's independently
derived count of 51 interference vertices with the relevant normal-cone
structure. Twelve rows is a small sample and the DESIGN-REAL observation is
recorded as a correlation, not a law.

## Finding 4. The two hardest census vertices are the least symmetric ones

The census minimum `dim Stab_{U(d)}([psi]) = 2` is attained at exactly five
vertices: `(3,8)#32`, `(3,10)#89`, `(3,10)#103`, `(5,10)#233`, `(5,10)#280`.

`(3,10)#89` and `(3,10)#103` are v89 and v103, the two rank-10
cancellation-regime closures that took the most work in the census and that
carry the fiber program's hardest open questions. They are picked out here by
an invariant that knows nothing about how they were found. The maximum, 58, is
`(3,10)#0`, the Slater vertex.

## Finding 5. Quantization has a period, and it is not the spectrum denominator

P6 predicted multiplicity 1 at every admissible `k`. It FAILS at 4 of 14
vertices, and the failure has structure.

| vertex | spectrum | denominator q | period |
| --- | --- | --- | --- |
| `(3,6)#3` | `(1,1,1,1,1,1)/2` | 2 | **4** |
| `(3,7)#5` | `(2,2,1,1,1,1,1)/3` | 3 | **6** |
| `(3,7)#7` | `(1,1,1,1,1,1,0)/2` | 2 | **4** |
| `(3,7)#8` | `(2,2,2,2,2,1,1)/4` | 4 | **8** |

The other ten vertices have period exactly `q`. Three facts make this a
finding rather than noise:

1. Every nonzero multiplicity in the whole run is exactly **1**. The reduced
   space over a census vertex is a point as far as multiplicity growth can see.
   What fails is not the size of the multiplicity but the SUPPORT of the
   sequence.
2. The vanishing is exactly periodic. The model "multiplicity is 1 on multiples
   of `p`, 0 elsewhere" fits all 14 vertices with zero exceptions
   (`period_model_inconsistent: []`).
3. Where the period beats the denominator it is **exactly `2q`**, never
   anything else.

So a census vertex carries an integer invariant, its quantization period, that
is invisible in the polytope: two vertices with the same denominator can have
different periods, and the period is what says at which `k` the vertex is
actually realized by a representation.

`(3,6)#3` is explained by classical invariant theory: `lambda = (1/2)^6` makes
`k*lambda` a power of the determinant, so the multiplicity is the dimension of
the degree-`k` `SL(6)` invariants of `wedge^3 C^6`, and that ring is generated
by a single quartic. Period 4 is the degree of that generator. `(3,7)#7` is
`(3,6)#3` padded with an empty orbital and inherits period 4, which is
consistent with the padding functoriality the census already records. The
mechanism at `(3,7)#5` and `(3,7)#8` is NOT established here.

Coverage is `(3,6)` and `(3,7)` only: all 14 vertices, `k` up to 16, with no
vertex truncated by the budget (`vertices_untested_budget: []`). The four
period-`2q` vertices are each observed over four to eight consecutive
multiples, so the alternation is not a two-point extrapolation. Coverage stops
at rank 7 because the weight count grows quickly in `k` and `d`; a single
`(3,7)` weight at `k = 14` took over fifteen minutes on its own.

## Finding 6. A convention trap, caught by the artifact's own internal check

Differentiating `rho_pq = <psi|a_p^* a_q|psi>` along the flow of `X in u(d)`
gives `d(rho) = [rho, X^T]`. So the subalgebra preserving `rho`, which must
contain the entire projective stabilizer because `U psi = e^{i t} psi` forces
`U rho U^* = rho`, is the commutant of `rho^T` and NOT of `rho`.

The two agree identically for the 645 states whose 1-RDM is real diagonal, and
differ as subspaces of the same dimension for the other 154. The first run used
`rho`, and 504 of the 505 rows it had produced looked fine. Row `(3,10)#73`
did not: it reported `stab_dim = 6` against `commutant_stab_dim = 4`, which is
impossible. One row out of 505 exposed an error that was silently wrong
wherever the two subspaces happened to give equal orbit ranks.

Both forced identities are now checked on every row and asserted by the guard
test, and both are clean across all 799:

- `stab_dim == commutant_stab_dim`: 0 violations.
- `commutant_dim == sum of squared eigenvalue multiplicities of the vertex
  spectrum`: 0 violations. This one is an independent confirmation that each
  certified state's 1-RDM really carries the census spectrum, computed without
  ever diagonalizing `rho`.

## The determination matrix

Ten source fields against six target fields, scored within each system because
`dim Stab` lives in `u(d)` and cannot be compared across ranks. A source whose
classes are all singletons determines everything for free; those rows are
flagged `vacuous` and carry no information. All six `spectrum -> X` rows are
vacuous, because all 799 vertices have distinct spectra.

Non-vacuous exact determinations:

| relation | direction | classes |
| --- | --- | --- |
| `support_invariant -> toral_stab_dim` | combinatorial -> symplectic | 411 |
| `spectrum_multiplicities -> commutant_dim` | spectral -> symplectic | 450 |
| `active_facets -> cone_simplicial` | lattice -> lattice | 141 |
| `cone_lattice_index -> cone_simplicial` | lattice -> lattice | 146 |

Only the first is both non-vacuous and not forced. The second is the identity
`dim = sum m_i^2` and is a validation, not a discovery; the last two are close
to definitional. Every other pairing fails, including every attempt to predict
`stab_dim` or `cone_lattice_index` from combinatorial data.

The user-facing summary: of the four presentations, exactly one pair is tied
together by a clean law on this census, and it is support to torus.

## Independent replay

`scripts/verify_symplectic_invariants_standalone.py` imports nothing from
`gpc_census` and recomputes each quantity by a deliberately different route:
the stabilizer rebuilt from the `pretty` amplitudes from scratch; the toral
dimension as the nullspace of `(phi, c) -> A phi - c*1` instead of the rank of
the row differences of `A`; the cone index from the H-representation by
inverting the active facet normals instead of reading neighbours off the
V-representation; and multiplicities by a dictionary generating-function DP
instead of the memoized suffix search.

Result: **140 checks matched, 0 failures, 2 skipped.** The skips are
`(3,6)#0` and `(3,6)#1`, where more valid inequalities are active at the vertex
than there are facets, so the H-rep route declines to answer rather than guess.
Run separately over every simplicial vertex, the H-rep cone route agrees with
the generator at **181 of 183**, with the same 2 declined and 0 mismatches.

The replay also checks the Weyl dimension identity
`sum_mu m(mu) dim V_mu = dim Sym^k(wedge^n C^d)` at `(3,6)` for `k = 1, 2`
(20 and 210) and `(2,4)` for `k = 2, 3` (21 and 56), which validates the
multiplicity machinery against arithmetic that does not depend on this campaign
at all.

The first version of the H-rep cone route was itself wrong: it replaced the
facet normals with an echelon basis of their span, which changes the dual basis
and hence the index, and it disagreed with the generator at every simple
vertex. The telltale was that its answers were small factorial-looking numbers.
That is recorded because a verifier that agrees for the wrong reason is worse
than no verifier.

## What this does NOT establish

- No claim about `dim M_lambda` from symplectic reduction. First-order tangent
  rank is vacuous at these ranks: `ker d(mu)` has dimension at least
  `2*binom(d,n) - d^2` for trivial reasons, and a pilot at the `(3,6)` Slater
  vertex, whose fiber is a single point, returns a first-order projective
  tangent dimension of 20. The prereg says this before the run. The only
  handle on `M_lambda` used here is multiplicity growth, and it is evidence at
  finitely many `k`, not a proof for all `k`.
- Every stabilizer number is UNREDUCED and FIXED-RHO. Nothing here measures a
  fixed-spectrum fiber.
- P4 holds on this census. It says nothing about rank 11 and above, and it does
  not show the toral part is the whole stabilizer, which Finding 2 refutes.
- Multiplicity 1 at the tested `k` is evidence that the reduced space is a
  point, not a proof.
- Rank-10 rows inherit the standing conditionality of the rank-10 constraint
  lists.

## Next, in the order I would do it

1. Extend the quantization run to `(3,8)` and `(4,8)`. The period is the most
   promising new invariant here and it currently rests on 14 vertices. A
   smarter weight count, or a Littlewood-Richardson route instead of the
   alternating sum, is the enabling work.
2. Explain periods 6 and 8. If the `(3,6)` mechanism generalizes, the period is
   the degree of a generator of a semi-invariant ring, which would tie the
   census to invariant theory the way `(3,6)` already is.
3. Characterize the 67 split classes. They are the exact locus where amplitudes
   carry symmetry the support does not see, which is a sharper target than
   "compute more invariants".
4. Certify the support invariant as a canonical form, so the class counts stop
   being lower bounds and Finding 2 gets an exact failure rate.
