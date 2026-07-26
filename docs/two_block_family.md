# The two-block interference family, and what governs the two-fiber gap (2026-07)

Status: the four-vertex family is described exhaustively below (EXACT
combinatorics, NUMERICAL tangent dimensions). The classification result,
"a cross-block channel is necessary for the two fibers to disagree", is
CONFIRMED census-wide with zero exceptions on 799 states. The converse is
FALSE, with a mechanism. The handoff's proposed counting rule is FALSIFIED.

Certifying artifacts: `scripts/two_block_family.py`,
`results/data/two_block_family.json`. Pre-registration:
`docs/prereg_two_block_inversion.md`, committed before the cross-block counts
were computed.

## Two corrections to the motivating handoff, up front

**The inversion is not sample size one.** The handoff states that
fixed-rho tangent 0 with fixed-spectrum tangent positive occurs only at v103.
`results/data/cascade.jsonl` has **14** such states, in the patterns
(fr, fs) = (0,1) x 11, (0,2) x 2, (0,6) x 1, spanning 7 transport classes.
v103 is the unique inversion of MAGNITUDE 6; it is not the unique inversion.
The sample for the rule was 14, not 1, which is why the rule could be scored
at all.

**Two blocks is not the operative invariant.** The family is real and worth
describing, but block count is not what governs fiber behaviour. The
two-fiber gap population spans block counts 2 through 6 (2:2, 3:21, 4:31,
5:18, 6:1), and two-block vertices are UNDER-represented: 3 of 64 carry a
cross-block channel against 155 of 799 census-wide. See "Does the family have
a closure" below.

## The four members

Exactly four census vertices are INTERFERENCE with a two-valued spectrum.

| vertex | spectrum | supp | dim K | within | cross | fs | fr | gap | cancellation | cascade drop |
|---|---|---|---|---|---|---|---|---|---|---|
| (3,10) v103 | (18^4, 5^6)/34 | 14 | 5 | 2 | 3 | 6 | 0 | 6 | yes | 4 |
| (3,10) v89 | (15^2, 6^8)/26 | 10 | 1 | 1 | 0 | 0 | 0 | 0 | yes | 0 |
| (3,10) v108 | (5^5, 1^5)/10 | 8 | 1 | 0 | 1 | 2 | 1 | 1 | no | 1 |
| (3,8) v32 | (15^4, 6^4)/28 | 7 | 0 | 0 | 1 | 0 | 0 | 0 | no | 0 |

"within" and "cross" are the one-hop channel classes split by whether the two
modes lie in the same eigenvalue block. `fs` and `fr` are the fixed-spectrum
and fixed-rho first-order tangent dimensions modulo gauge.

Arithmetic and gauge data, from the holonomy and gauge-minimization artifacts:

| vertex | state den | ratio | loops | holonomy min polys | gauge lower bound | gauge-minimal certified |
|---|---|---|---|---|---|---|
| v103 | 195364 | 5746 | 5 | x-1, x+1, x+1, x-1, x-1 | 9 | no |
| v89 | 130 | 5 | 1 | x+1 | 4 | no |
| v108 | 10 | 1 | 1 | x | 7 | no |
| v32 | 28 | 1 | 0 | none | 6 | no |

Every holonomy in the family is a sign except v108's, whose cosine is 0: v89
and v103 are real up to gauge, v108 carries a genuine quarter-turn phase, and
v32 is loop-free. The two states with a nontrivial denominator ratio are
exactly the two cancellation-regime members, consistent with the minting
mechanism (`docs/prereg_denominator_arithmetic.md`).

## The mechanism, and what the family shows

The silence-rank lemma (`docs/silence_rank_lemma.md`) says first-order
fixed-SPECTRUM weight deformations see only the WITHIN-block channel
conditions, while fixed-RHO deformations see all of them. Cross-block
channels are therefore the entire first-order difference between the two
fibers. The family reads directly off that:

- **v89** has one channel and it is WITHIN-block. No cross-block condition
  exists, so the two fibers cannot differ, and both are rigid. This is why
  v89 is rigid in BOTH senses while v103 is not.
- **v103** has 3 cross-block channels out of 5, and the fixed-rho fiber sees
  all five while the fixed-spectrum fiber sees only two. Hence the inversion.
- **v108** has its single channel CROSS-block, giving gap 1.
- **v32** has a cross-block channel but incidence kernel 0: there is no
  first-order weight freedom for the condition to cut, so the gap is zero for
  want of anything to differ about.

## What the six v103 dimensions are

The handoff asks what the six fixed-spectrum directions at v103 actually are,
given that the state is fixed-rho rigid and that the sparsification which took
it to support 10 needed a non-block rotation (`results/data/orbit_check.json`),
so they are not stabilizer directions.

Decomposing the tangent by sector (the ambient real coordinates split into
perturbations of the real parts of the amplitudes, the WEIGHT sector, and of
the imaginary parts, the PHASE sector):

    v103   tangent 6 = 3 weight + 3 phase, sectors decouple exactly
    v89    tangent 0 = 0 weight + 0 phase, sectors decouple exactly
    v32    tangent 0 = 0 weight + 0 phase, sectors decouple exactly
    v108   tangent 2, sectors do NOT decouple

This reproduces the lemma's dimension formula term by term at v103:
(dim K - rank J_within^Re) + (|S| - g - rank J_within^Im)
= (5 - 2) + (14 - 9 - 2) = 3 + 3. So the six directions are three
weight deformations inside the five-dimensional incidence kernel that survive
the two within-block real conditions, plus three phase directions beyond the
nine gauge directions that survive the two within-block imaginary conditions.

The v108 line is a SCOPE NOTE on the lemma, not a defect: its proof assumes a
real state, where the weight and phase parts of a deformation are the real and
imaginary sectors. v108 is genuinely complex (holonomy cosine 0), the sectors
mix, and the weight/phase split is not defined there. The lemma's dimension
formula is stated for real states and this is what its hypothesis is doing.

## The scored predictions

Full detail in `results/data/two_block_family.json`, scorecard block. Scope is
799 states; independent scope is counted over transport classes per standing
rule R2 of `docs/prereg_template.md`.

| prediction | registered expectation | verdict | exceptions |
|---|---|---|---|
| P-INV-1 GAP == X | FAIL | **FAIL** | 94 |
| P-INV-2 GAP == 2X | uncertain | **FAIL** | 154 |
| P-INV-3 GAP > 0 iff X > 0 | PASS | **FAIL** | 82 |
| P-INV-3 forward, GAP > 0 implies X > 0 | PASS | **PASS** | 0 |
| P-INV-3 reverse, X > 0 implies GAP > 0 | PASS | **FAIL** | 82 |
| P-INV-4 every inversion has X > 0 | PASS | **PASS** | 0 |
| P-INV-5 cancellation regime only at two blocks | uninformative | **UNINFORMATIVE** | 0 |
| P-INV-6 family rigidity iff dim K == 0 or X == 0 | uncertain | **PASS** | 0 |

**The counting rules are dead.** The handoff proposed that the inversion
magnitude equals the number of cross-block channels. It does not: v103 has 3
cross-block channels and gap 6, while v108 has 1 and gap 1. The factor-two
variant fixes v103 and breaks v108. A bare class count cannot equal a
dimension, because a channel condition is complex, hence two real conditions,
and because those conditions are not independent once several channels share
support.

**One half of P-INV-3 survives, with zero exceptions.**

> **GAP > 0 implies X > 0.** If the fixed-rho and fixed-spectrum first-order
> tangent dimensions differ at a certified state, the support has at least one
> cross-block one-hop channel. Verified on all 799 certified states, 403
> transport classes, zero exceptions. All 73 states with a gap have X > 0.

That is the silence-rank lemma's cross-block clause promoted from a dimension
formula verified at two states to a census-wide NECESSARY CONDITION. It is
the result of this study.

**The converse fails, and the failure has a mechanism.** 82 states have a
cross-block channel and no gap. Every one of them has incidence kernel
dimension 0 and is rigid in both fibers: the cross-block condition is there,
but there is no first-order weight freedom for it to act on. Restricting to
states that deform at all in the fixed-spectrum fiber, the biconditional does
hold with zero exceptions on all 73 -- but that restriction was chosen AFTER
seeing which direction failed, is recorded as POST-HOC in the artifact, and
has had no out-of-sample test.

## Does the family have a closure

No, and the question resolves the framing. Task 3 of the handoff asked whether
the two-block family has a natural closure that grows the sample past four.

- The 60 two-block DESIGN vertices have ZERO one-hop channels, all of them,
  because a design support is one-hop free by definition. So X = 0 and no gap
  is possible. The design end of the family is trivial for a structural
  reason, not an interesting one, and it cannot be a source of new specimens.
- The gap population spans block counts 2 through 6, peaking at 3 and 4.
- Two-block vertices are the WORST source of cross-block channels: 3 of 64,
  against 155 of 799 census-wide.

So the four members span the known fiber behaviours because they happen to
span (X, dim K), not because two blocks causes anything. The natural closure
of the family is the X-stratification, which already has 155 members in 72
transport classes and is the population any future prediction should be
scored on.

The one two-block fact that survives is that both known cancellation-regime
states are two-block. That sample is size two, P-INV-5 is UNINFORMATIVE, and
it is the only reason two-block spectra remain worth targeting in Stage 1
(`docs/GENERATOR_PLAN.md`, targeting rule 1).

## Open

1. Is "GAP > 0 implies X > 0" provable? The lemma gives the mechanism at
   first order for real states in the cancellation regime; the census
   statement holds across regimes and includes complex states, so the proof
   would need the complex bookkeeping the lemma defers.
2. The restricted biconditional is post-hoc. It needs an out-of-sample test,
   which means Stage 1.
3. The v103 germ beyond first order: the three weight and three phase
   directions are first-order data. Which integrate, and whether the germ
   meets supports smaller than 14 in the NATURAL-ORBITAL sense, is untouched
   here. The support-10 cascade endpoint does NOT answer it: that endpoint is
   the library state in rotated orbitals and the rotation is not a gauge
   rotation, so it bounds the orbit-minimal support, not s_Q^NO
   (`results/data/orbit_check.json`, and the support taxonomy in
   RESEARCH.md).
4. Silence-rank genericity remains at sample size two.
