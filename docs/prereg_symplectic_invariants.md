# PRE-REGISTRATION: cross-realization invariants of the census (2026-08-07)

Written BEFORE the measuring scripts were run. Standing rules:
`docs/prereg_template.md`. Certifying artifacts:
`scripts/symplectic_invariant_census.py` -> `results/data/symplectic_invariants.json`,
`scripts/quantization_multiplicity.py` -> `results/data/quantization_multiplicities.json`,
standalone replay `scripts/verify_symplectic_invariants_standalone.py`.

## What is being tested

The census currently describes each vertex through one presentation: the
occupation spectrum, its support combinatorics, and the facets it saturates.
This campaign computes the same 799 objects in three OTHER presentations and
asks whether the presentations agree.

1. LATTICE presentation. Each vertex `v` of the moment polytope `P(n,d)` has a
   tangent cone `T_v P` inside the weight lattice `L_0 = {x in Z^d : sum x = 0}`.
   Recorded per vertex: number of active facets, simpliciality, primitive edge
   ray generators, the Smith normal form of the ray matrix against a basis of
   `L_0`, the lattice index of the cone, and unimodularity.

2. SYMPLECTIC/LIE presentation. Each certified state `psi` has a projective
   stabilizer `Stab_{U(d)}([psi])`. Recorded per state: `dim Stab`, the orbit
   dimension `d^2 - dim Stab`, the commutant `G_rho = {U : U rho U^* = rho}`
   with its dimension, and `dim G_rho . [psi]`. All exact.

3. REPRESENTATION presentation. For a vertex `lambda` with denominator `q`, the
   multiplicity `m_k(lambda)` of the irreducible `U(d)` representation of
   highest weight `k*lambda` inside `Sym^k(wedge^n C^d)`, for `k` a multiple of
   `q`. This is the quantization of the same moment polytope.

The single law under test, in the form that can fail:

> the discrete support combinatorics of a census state determines its
> continuous symmetry, and the vertex's lattice-cone type is not independent
> of that symmetry.

## What is deliberately NOT claimed

`dim M_lambda` for the symplectic reduction `M_lambda = mu^{-1}(diag lambda)/G_lambda`
is NOT computed by first-order tangent rank. A pilot on `(3,6)` index 0 (the
Slater vertex, whose fixed-rho fiber is a single point) returns a first-order
projective tangent dimension of 20, because `ker d(mu)` has dimension at least
`2*binom(d,n) - d^2` for trivial reasons at every state where the linearization
degenerates. Every tangent number of that kind is a vacuous upper bound at the
ranks in this census, and none is reported. `dim M_lambda` is approached ONLY
through prediction P6 below (multiplicity growth), which is exact.

Fiber qualification, per the standing rule: `G_rho` below is the commutant of
the state's OWN 1-RDM, so every stabilizer number here is UNREDUCED and
FIXED-RHO. Nothing in this campaign measures a fixed-spectrum fiber.

## Holdout and its provenance

Not a holdout. This is a census-wide recomputation of all 799 rows from
`results/data/states.jsonl` and `results/data/vertices/`, so R2's derivation
factor applies in full and is stated up front rather than discovered later:

- 799 rows total, 9 systems.
- The `(4,8)`, `(4,9)`, `(4,10)` and `(5,10)` systems reach the census through
  particle-hole duality and lifts from lower rank in parts of the pipeline, so
  their rows are NOT independent of the `(3,d)` rows for any prediction whose
  statement is duality-covariant. `dim Stab` and the cone type ARE
  duality-covariant, so the honest independent scope for P2 to P5 is the set of
  distinct support-isomorphism classes, not 799. The artifact records
  `support_class` per row and every scorecard reports
  `independent_distinct_states` as the number of distinct classes.
- Predictions P1 and P6 are properties of individual rows and systems and are
  scored at full scope, with the derivation factor stated.

## Power caveat, stated before scoring

A PASS on P4 would establish a correspondence on the census as it stands. It
would NOT establish it for any rank above 10, and it would not prove that the
toral stabilizer is the whole stabilizer. A PASS on P6 at the finitely many `k`
tested does not prove multiplicity-freeness for all `k`; it is evidence for the
reduced space being a point, not a proof of it.

## Predictions

P1. Exactly 645 of the 799 certified states have an exactly diagonal 1-RDM
    `rho = diag(lambda)`, and 154 do not. Expectation: PASS, because
    `results/data/uniqueness_census.json` records "the 154 records whose 1-RDM
    is not diagonal", and 631 DESIGN-INT plus 12 DESIGN-REAL is 643, so two of
    the 156 INTERFERENCE states must already be diagonal. Scored exactly. A
    disagreement falsifies either this arithmetic or that artifact.

P2. For the Slater vertex of each of the nine systems (integer form with
    entries in {0, den}), `dim Stab_{U(d)}([psi]) = n^2 + (d-n)^2` and the
    orbit dimension is `2n(d-n)`. Expectation: PASS, this is the Grassmannian
    `Gr(n,d)` and is forced. Pilot confirms it at `(3,6)` (18) and `(5,10)`
    (50). Scored on 9 rows; independent scope 9 systems, but the statement is
    a theorem, so a FAIL means the pipeline is wrong, not the mathematics.

P3. `dim Stab_{U(d)}([psi]) >= 1` for all 799 rows, and the census minimum is
    strictly greater than 1. Expectation: PASS on the first clause (the
    global phase `i*I` always stabilizes projectively) and PASS on the second.
    Pilot minimum over the states sampled so far is 6.

P4. THE CROSS-REALIZATION PREDICTION. Let `A` be the 0/1 support-incidence
    matrix of the state (rows the support determinants, columns the orbitals),
    and let `t_phi(psi) = d - rank_R([A_S - A_{S_0}]_S)` be the dimension of the
    orbital-phase subgroup fixing `[psi]`. Then for every one of the 799 rows,
    the toral part of the stabilizer has dimension exactly `t_phi`, i.e.
    `dim (Stab_{U(d)}([psi]) intersect u(1)^d) = t_phi`, a quantity computed
    from the support combinatorics alone with no amplitude input.
    Expectation: PASS. Rationale: an orbital phase acts on `c_S` by
    `exp(i sum_{o in S} phi_o)`, so it fixes `[psi]` iff `A phi` is constant on
    the support, which is exactly the stated rank condition, independent of the
    amplitudes as long as none vanishes. The content of the prediction is that
    no census state has a hidden amplitude coincidence that enlarges the toral
    stabilizer beyond the combinatorial count.

P5. Cone types. (a) A vertex has a simplicial tangent cone iff it is simple in
    the sense of `docs/polytope_invariants.json`, so the simplicial count is
    exactly 183 = 4+2+7+9+12+25+26+32+66. (b) Unimodular cones are a strict
    minority of the simplicial ones, i.e. fewer than 92 of the 799 vertices
    have a unimodular tangent cone. Expectation: (a) PASS by definition of
    simple, recorded to catch pipeline error; (b) PASS, because
    `docs/RESEARCH_ARCHIVE.md` already records an edge-ray lattice index of
    `23^7` at `v_B`, so large indices exist. (b) is the falsifiable half.

P6. For every vertex of `(3,6)` and `(3,7)` and every admissible `k` tested
    (`k` a positive multiple of the vertex denominator, up to the compute
    budget), the multiplicity of highest weight `k*lambda` in
    `Sym^k(wedge^3 C^d)` equals 1. Expectation: PASS. Rationale: a vertex of
    the moment polytope should have a point as its reduced space, and a point
    quantizes to a one-dimensional space. A FAIL at any vertex is the more
    interesting outcome: it would exhibit a census vertex whose reduced space
    is positive dimensional, and would identify it by name.

P7. Cross-tabulation, exploratory and scored as a single yes/no. Within each
    system, `dim Stab` is constant on support-isomorphism classes. Expectation:
    PASS, forced if P4's mechanism extends to the non-toral part; a FAIL
    localizes exactly where amplitudes carry symmetry that the support does not
    see.

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-FALSIFIED per prediction, decided by the
artifacts alone. Each verdict reports scope, independent scope, and independent
distinct states. No prediction is edited after the run. Numbers are exact
throughout; any row that cannot be settled exactly is reported in its own
`evidence` bucket and never totaled with the exact rows.

## SCORED (2026-08-07, after running; predictions above are unedited)

The authority is the `scorecard` block of
`results/data/symplectic_invariants.json` and of
`results/data/quantization_multiplicities.json`. Discussion in
`docs/symplectic_invariants.md`.

| P | verdict | scope | indep scope | measured |
| --- | --- | --- | --- | --- |
| P1 | PASS | 799 | 411 | 645 diagonal, 154 not, exactly as predicted |
| P2 | PASS | 9 | 9 | all nine Slater vertices Grassmannian, 0 violations |
| P3 | PASS | 799 | 411 | min `dim Stab` = 2 at 5 rows, max 58 |
| P4 | PASS | 799 | 411 | 0 violations of the combinatorial toral prediction |
| P5 | PASS | 799 | 411 | 183 simplicial = 183 simple as predicted; 3 unimodular, bound was 92 |
| P6 | **FAIL** | 14 | 14 | 4 of 14 vertices have a zero multiplicity at an admissible `k` |
| P7 | **FAIL** | 799 | 411 | 67 of 411 support classes carry more than one `dim Stab` |

Independent scope is the count of distinct support-invariant values. The
invariant can merge non-isomorphic supports, so it understates independence,
which is the conservative direction.

### P6, the informative failure

P6 was misformulated, not merely wrong, and the correction is the finding.
It assumed that once `k*lambda` is an integral weight the vertex is realized.
The multiplicity is instead 1 on multiples of a PERIOD `p` and 0 elsewhere,
with `p = q` at ten of the fourteen vertices and `p = 2q` at the other four
(`(3,6)#3`, `(3,7)#5`, `(3,7)#7`, `(3,7)#8`). The periodic model fits all 14
vertices with zero exceptions, and every nonzero multiplicity in the run is
exactly 1, so the reduced-space-is-a-point reading that motivated P6 survives
intact. Only the claim about which `k` realize the vertex was wrong.

Under standing rule R4 this is checked against the corpus: no artifact in the
repository records a quantization period, so this is new information rather
than a prediction already refuted by data in hand.

### P7, the boundary of P4

P4 and P7 together locate exactly where the combinatorics stops. The support
fixes the toral part of the stabilizer with no exceptions in 799 rows, and
fails to fix the full stabilizer on 67 classes. `dim Stab` exceeds its toral
part at 445 of 799 vertices, so the non-toral part is not a rare correction.

### Corrections during the run, per R5

The first run of the census scored P1 to P5 on 505 rows before a forced
identity failed at `(3,10)#73`, exposing a real error: the commutant was taken
of `rho` where the derivative `d(rho) = [rho, X^T]` requires `rho^T`. The 645
real-diagonal rows are provably unaffected, since `rho^T` is literally `rho`
there; the 154 others were recomputed. No verdict from the first run is quoted
anywhere, and none was scored before the fix.
