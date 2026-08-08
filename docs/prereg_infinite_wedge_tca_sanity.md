# PRE-REGISTRATION: the stable highest-weight semigroup and cross-rank facet templates (2026-08-08)

Written BEFORE the measuring script was run, with two disclosed pilots stated
below. Standing rules: `docs/prereg_template.md`. Certifying artifacts:
`scripts/infinite_wedge_tca_sanity.py` ->
`results/data/infinite_wedge_tca_sanity.json`, guard test
`tests/test_infinite_wedge_tca_sanity.py`.

## What is being tested

The stable-rank compression program asks whether all finite-rank fermionic
moment polytopes for a fixed particle number `N` are recoverable from one
`GL_infinity` object, and whether that would compress the constraint
constructor across the rank variable `d`. This campaign does not attempt the
compression theorem. It binds the abstract definitions to certified repository
data at the four places where a wrong convention or a wrong functoriality
assumption would silently invalidate every later argument:

1. the convention tying plethysm multiplicity to occupation spectrum;
2. the recovery of the moment polytope from occurring highest weights;
3. whether the weight monoid is saturated, which decides whether the polytope
   determines the semigroup;
4. whether valid rows transport UP the rank ladder by zero padding, DOWN it by
   restriction, or both, and whether the published systems show any sign of
   template saturation inside the census window.

Definitions used throughout, fixed here before the run.

    Sigma_{N,d} = { (lambda, k) : lambda a partition, ell(lambda) <= d, k >= 1,
                    mult(V_lambda^{GL_d}, Sym^k(wedge^N C^d)) > 0 }

    Sigma_N^stab = { (lambda, k) : <s_lambda, h_k[e_N]> > 0 }   (no d)

    Delta(N,d)  = the ordered fermionic moment polytope at (N,d)

    padding template of a row sum_i c_i lambda_i <= b at rank d:
        the coefficient vector with trailing zeros stripped, divided by the
        positive gcd, paired with the correspondingly scaled b. Two rows at
        different ranks share a padding template exactly when one is the
        zero padding of the other. Its ARITY is the index of the last nonzero
        coefficient. This is the ordered (OI) notion.

    symmetric template of the same row:
        the multiset of nonzero coefficients, divided by the positive gcd,
        paired with the scaled b. Its ARITY is the number of nonzero
        coefficients. This is the unordered (FI, ambient Weyl-invariant cone)
        notion.

## Holdout and its provenance

Every source is already in the corpus. Nothing here is an out-of-sample
measurement and the campaign does not claim to be one. Under R2 the honest
accounting is:

- Multiplicity rows: the 14 `(3,6)` and `(3,7)` census vertices with a measured
  multiplicity sequence in `results/data/quantization_multiplicities.json`,
  93 (vertex, k) pairs. `(3,7)#7` is `(3,6)#3` padded with an empty orbital
  and is NOT independent of it; `docs/quantization_character.md` already
  records that the two carry the same obstruction. Independent distinct
  vertices: 13.
- Vertex sets: the certified census vertex files for the nine systems. The
  `(4,8)`, `(4,9)`, `(4,10)` and `(5,10)` systems reach the census partly
  through particle-hole duality and lifts, so the nine systems are not
  independent of each other.
- Constraint rows: the published Altunbulak-Klyachko tables in
  `src/gpc_census/data/constraints.json`, each already re-certified as an
  exact Ressayre row by `results/data/census_inner_sandwich.json`.

Independent scope for the cross-rank predictions is the number of consecutive
rank pairs, which is 4 at `N = 3` (6 to 7, 7 to 8, 8 to 9, 9 to 10) and 2 at
`N = 4` (8 to 9, 9 to 10). Six pairs, of which the `N = 4` pair 9 to 10 is the
least independent because of the duality route.

## Disclosed pilots

Both were derived by hand before the script existed and are therefore NOT
scored as discoveries.

**Pilot A.** The Borland-Dennis inequality `lambda_4 <= lambda_5 + lambda_6` at
`(3,6)`, zero padded to rank 7, is violated at the certified `(3,7)` vertex
`(3,3,3,3,1,1,1)/5`, where the left side is `3/5` and the right side is `2/5`.
So zero padding does NOT preserve validity in general, and the `(3,6)` system
in particular does not pad up. The three `(3,6)` equalities cannot pad up
either, since `(3,7)` has no equalities.

**Pilot B.** Reading `src/gpc_census/data/constraints.json` while designing
this campaign, the first four rows of the `(3,8)` block are visibly the four
`(3,7)` rows with a zero appended. So the 7 to 8 pair is known to nest before
the run. The 8 to 9, 9 to 10 and both `N = 4` pairs are blind.

## Power caveat, stated before scoring

A PASS anywhere here is a statement about the census window, ranks 6 to 10 and
`N` in {3,4,5}. Nothing measured at four consecutive ranks constrains the
asymptotics in `d`, and no prediction below may be read as a theorem about all
`d`. In particular a template count that grows slowly over four points is not
evidence of polynomial growth, and a template count that does not saturate over
four points is not a proof that it never saturates.

The multiplicity side cannot pass rank 7. That is a hard budget limit recorded
in `docs/symplectic_invariants.md`, not a choice, so predictions P1 to P3 are
permanently capped at 14 vertices unless a better multiplicity algorithm
arrives.

## Predictions

**P1 (forced correctness gate, uninformative as evidence).** For all 93
(vertex, k) pairs, three independently implemented routes agree exactly:
the Murnaghan-Nakayama plethysm pairing `<s_{k lambda}, h_k[e_N]>` in
`gpc_census.plethysm`, the Weyl alternating sum over exact weight counts in
`scripts/quantization_multiplicity.py`, and the stored integers in
`results/data/quantization_multiplicities.json`. Expectation: PASS. Recorded
only as a gate on the code and on the convention. A FAIL means the partition,
dual or transpose convention in one of the two routes is wrong.

**P2 (load-bearing control on P1).** The convention is not vacuous: replacing
`k lambda` by its conjugate partition, or by the reversed complement
`(k - lambda_d, ..., k - lambda_1)`, changes at least one of the 93 answers.
Expectation: PASS. Rationale: a gate that every convention passes proves
nothing, and this project has been bitten by exactly that before
(`docs/census_inner_sandwich.md`, "plumbing that always says yes proves
nothing"). A FAIL means P1 cannot distinguish conventions and must be
redesigned.

**P3 (restatement, already implied by artifacts in hand; not a discovery).**
`Sigma_{N,d}` is NOT saturated. Concretely, for a census vertex with integer
form `lambda_int` and denominator `q`, the pair `(lambda_int, q)` is a lattice
point of `cone(Sigma_{N,d})` by construction, and it lies in `Sigma_{N,d}`
exactly when the quantization period equals `q`. Prediction: exactly 4 of the
14 measured vertices give a lattice point of the cone with multiplicity 0,
and they are exactly the four period-`2q` vertices `(3,6)#3`, `(3,7)#5`,
`(3,7)#7`, `(3,7)#8`. Expectation: PASS, and flagged ALREADY-KNOWN under R4
because `docs/symplectic_invariants.md` Finding 5 records the periods. What is
new is only the semigroup reading, not the numbers.

**P4 (blind above rank 7).** Zero padding preserves validity across every
consecutive census pair at `N = 3` above rank 6 and at `N = 4`: for
`d` in {7,8,9} the published `(3,d)` rows padded with zeros are satisfied by
every certified `(3,d+1)` vertex, and likewise for `(4,8)` into `(4,9)` and
`(4,9)` into `(4,10)`. Expectation: PASS above rank 6, FAIL at the 6 to 7 pair
by Pilot A. Rationale: a Klyachko row is indexed by combinatorial data whose
support is bounded, and Pilot B shows the 7 to 8 pair nests literally; the
`(3,6)` system is expected to be the exception because it is the degenerate
system carrying equalities. A FAIL above rank 6 would be the stronger result:
it would mean the published tables are not monotone and no naive template
family can exist.

**P5 (forced correctness gate on the published tables).** Restriction is
complete downward: for every consecutive pair, the published rank-`d+1`
inequality system restricted to `lambda_{d+1} = 0`, together with the rank-`d`
ordered chamber and trace rows, has exactly the certified rank-`d` vertex set.
Expectation: PASS. This is forced by the theorem
`Delta(N,d) = Delta(N,d+1) intersect {lambda_{d+1} = 0}` plus `P = O` at both
ranks, so a FAIL indicts the tables or this code, not the mathematics. It is
registered because passing it upgrades a citation into a replayed proof that
the `(3,10)` system implies the `(3,6)` through `(3,9)` systems.

**P6 (blind).** No template saturation inside the census window. Define
`new(d)` as the number of published `(N,d)` rows whose padding template does
not already occur at rank `d-1`. Prediction: `new(d) >= 1` at every
consecutive census pair, so the template stock strictly grows at every rank
measured, and the symmetric template count also strictly grows at every rank.
Expectation: PASS. Rationale: the published row counts grow, and no mechanism
is known that would stop new templates appearing. A FAIL, meaning some rank
introduces no new template, would be the first positive evidence for finite
templating and would immediately raise the priority of the whole program.

**P7 (blind, and the one that matters for compression).** The symmetric
template arity is unbounded inside the window: the maximum number of nonzero
coefficients in a published `(3,d)` row strictly increases over
`d = 6,7,8,9,10`. Expectation: PASS. Rationale: this is the quantity that
enters the polynomial bound. If the number of templates were finite with
maximum symmetric arity `a`, the ambient Weyl-invariant moment cone at rank `d`
would have `O(d^a)` facets; unbounded arity inside the window is therefore
the first observable symptom that no finite template list exists. A FAIL,
meaning arity saturates, would be genuine evidence FOR compression and would
identify the saturating arity as the template bound to try to prove.

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-KNOWN per prediction, decided by the
artifact alone. Each verdict reports scope, independent scope and independent
distinct states. No prediction is edited after the run. P1 and P5 are forced
gates and are recorded as gates on the code, never quoted as evidence for the
science. P3 is scored but carries the ALREADY-KNOWN flag under R4.

## Nonclaims fixed before the run

- No prediction here is a statement about all `d`.
- Passing P4 does not prove that zero padding preserves validity at any rank
  above 10, and does not prove it for rows that are not in the published
  tables.
- Passing P5 does not prove facetness or irredundancy of any published row.
  The repository records both as not established, and this campaign does not
  change that.
- Nothing here computes a facet of a rank-11 or higher system, and nothing
  here bears on candidate exhaustiveness.

## SCORED (2026-08-08, after running; predictions above are unedited)

See `docs/infinite_wedge_tca_bridge.md`, section 7, and the `scorecard` block
of `results/data/infinite_wedge_tca_sanity.json`.
