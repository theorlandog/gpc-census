# PRE-REGISTRATION: what governs the fixed-rho / fixed-spectrum gap (2026-07)

Written BEFORE the cross-block channel counts were computed. Standing rules:
`docs/prereg_template.md`. Certifying artifact: `scripts/two_block_family.py`,
`results/data/two_block_family.json`.

## What is being tested

The silence-rank lemma (`docs/silence_rank_lemma.md`) says that first-order
fixed-SPECTRUM weight deformations see only the WITHIN-block channel
conditions, while fixed-RHO deformations see ALL of them. Cross-block channels
are therefore the entire first-order difference between the two fibers.

That is a statement about a mechanism. This pre-registration asks whether it
is also a CLASSIFICATION: does the cross-block channel structure predict, per
state, when and by how much the two fiber dimensions disagree?

## Definitions, fixed before scoring

For a certified state with support S whose 1-RDM is diagonal and equal to
lambda:

- BLOCKS are the eigenvalue blocks of lambda (modes grouped by equal value).
- A one-hop class of mode pair (A, B) is WITHIN-BLOCK when A and B lie in the
  same block, CROSS-BLOCK otherwise.
- `X` = the number of CROSS-BLOCK one-hop classes of S.
- `GAP` = `fixed_spectrum_tangent_dim - fixed_rho_tangent_dim`, both as
  measured by `gpc_census.fiber` and shipped in `results/data/cascade.jsonl`.
- An INVERSION state is one with `fixed_rho_tangent_dim == 0` and
  `fixed_spectrum_tangent_dim > 0`: rigid in one fiber, positive-dimensional
  in the other.

## Sample, and its independence

Scored on all 799 certified states of the census. This is IN-SAMPLE with
respect to the corpus the silence-rank lemma was written from (v89 and v103),
but the lemma was derived from those two only, and the other 797 states were
never consulted when writing it. Per standing rule R2 the counts below are row
counts over independent states: each of the 799 is its own certified state,
not a transport of another within this dataset. Transport relationships
between systems (padding, frozen-core, particle-hole) DO exist inside the
census, so a supplementary independent count over transport-classes is
reported alongside the raw count and is the number to quote.

## Correction to the motivating claim, recorded before scoring

The handoff that motivated this work states that the inversion is "sample size
one, census-wide" at v103. That is wrong, and the correction is recorded here
BEFORE the predictions so it cannot be mistaken for a result of them.
`results/data/cascade.jsonl` has 14 inversion states, in the patterns
(fr, fs) = (0,1) x 11, (0,2) x 2, (0,6) x 1. v103 is the unique one of
magnitude 6; it is not the unique inversion. The sample for the rule below is
therefore 14, not 1.

## Predictions

**P-INV-1.** `GAP == X` on every state. Expectation: **FAIL**. A channel
condition is complex, hence two real conditions, and the gauge dimensions of
the two fibers also differ, so a bare count cannot match a dimension.
Registered because it is the rule as originally proposed and it should be
scored as stated rather than quietly replaced.

**P-INV-2.** `GAP == 2 * X` on every state. Expectation: uncertain, roughly
even odds. Rationale: each cross-block channel contributes a real and an
imaginary condition to the fixed-rho problem and neither to the
fixed-spectrum problem, so a factor of two is what the lemma's bookkeeping
suggests if every cross-block condition is independent. It will fail wherever
those conditions are dependent, and dependence is expected once a state has
several channels on overlapping supports.

**P-INV-3 (the load-bearing one).** `GAP > 0` if and only if `X > 0`. That is,
the two fibers disagree exactly when the support has at least one cross-block
one-hop class. Expectation: **PASS**. This is the direct qualitative content
of the lemma's cross-block clause. If it holds on all 799 it promotes the
lemma from a dimension formula verified at two states to a classification of
which states have a two-fiber gap at all.

**P-INV-4.** Every inversion state has `X > 0`. Expectation: **PASS**, and
implied by P-INV-3, registered separately so the inversion subpopulation is
scored on its own.

**P-INV-5.** The cancellation regime occurs only at two-block spectra.
Expectation: **PASS in sample, UNINFORMATIVE as evidence**. Both known members
are two-block, and they are the only two known, so this prediction cannot fail
on the current corpus and must not be quoted as confirmation. It is registered
to fix the statement that Stage 1 will actually test
(`docs/GENERATOR_PLAN.md`).

**P-INV-6.** Among the four two-block INTERFERENCE vertices, the fiber
signature is determined by (`X`, `dim K`): the two with `dim K = 0` or no
cross-block channel are rigid in both fibers, and the rest are not.
Expectation: uncertain. Registered because the four-vertex family is the
handoff's target and a per-family rule is what would make it generalizable.

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-FALSIFIED per prediction, decided by
`results/data/two_block_family.json` alone. Every verdict reports scope,
independent scope over transport classes, and the full exception list. No
prediction is edited after the run. A FAIL is reported with its mechanism if
one is visible and as a bare FAIL if not.
