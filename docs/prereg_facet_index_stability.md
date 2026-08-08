# PRE-REGISTRATION: the facet-index stability probe (2026-08-08)

Written BEFORE the measuring script was run. Standing rules:
`docs/prereg_template.md`. Certifying artifacts:
`scripts/facet_index_stability.py`,
`results/data/facet_index_stability.json`,
`tests/test_facet_index_stability.py`.

## What is being tested

Two questions, one campaign, because they share a generator and a target: what
object about the facet system is stable in the rank, and can an intrinsic
criterion replace the exact reduction stage.

The motivation is measured, not aesthetic. The BDR comparison
(`docs/bdr_constructor_comparison.md`) found that the oracle ceiling FALLS from
21.41x at rank 7 to 16.38x at rank 8, because the floor no external tool can
remove is composition, which scales with the retained facet count. So the stage
with the worst long-run scaling is redundancy elimination, and the question
worth asking is whether an intrinsic facetness criterion can replace it.

## Part A: the padding and shape census. LANDING, NOT PREDICTED.

Rule R4 applies and is stated here rather than discovered later. The padding
validity table and the normal-entry table below were computed in a scratch
session during analysis, before this pre-registration existed. They are
therefore **not** discoveries and are **not** scored. Part A lands them with a
generator, an artifact, and a guard test, which is what the repository requires
before any of them can be cited.

What scratch found, restated here so the landing run can be checked against it
rather than silently replacing it:

1. Trailing-zero padding is NOT a valid map on facets. Padded lower-rank rows
   are violated at the next rank in 1 of 1, 0 of 4, 21 of 31 and 8 of 52 cases
   for `(3,6)` through `(3,9)` padded into `(3,7)` through `(3,10)`.
2. The primitive `tau` entries of published `N=3` facets are bounded by 7 at
   ranks 8, 9 and 10, while the facet count grows 31, 52, 93.
3. The number of distinct sorted-`tau` shapes grows 1, 1, 9, 8, 16 across ranks
   6 to 10 while the facet count grows 1, 4, 31, 52, 93.

The genuinely predictive part of Part A is the extension to the systems with
`N != 3`, which scratch never examined. That is P6.

Standing caveat on both tables: the rank-9 and rank-10 constraint lists are
conjecturally complete per AK2008 (`results/data/PROVENANCE.md`), so every
statement resting on them is conditional on that completeness.

## Part B: the coefficient-one facetness converse. PREDICTIVE.

AK2008 uses Schubert coefficient one as a sharpness criterion for the complete
irredundant list (`docs/stage1_klyachko_spec.md`, implementation plan step 4).
This repository deliberately refuses to use it as a general facetness or
irredundancy gate (`docs/fixed_n3_generator_methodology.md`), and preserves
coefficient one, positive nonunit, and zero as separate diagnostics.

The refusal has never been tested against data, because exact cyclic-Schubert
evaluation is DEFERRED until after redundancy elimination precisely to avoid
paying for it on every certified candidate. So the coefficients of the
certified rows that reduction REMOVES have never been computed at any rank.
That gap is exactly the converse direction of the criterion.

State of knowledge before this run, checked rather than assumed:

- `(3,7)` retained rows all have exact coefficient 1. **ALREADY ESTABLISHED**
  in `results/data/stage1_ressayre_3_7.json`. Not a prediction here.
- `(3,8)` retained-row coefficients are not stored in any artifact.
  `results/data/stage1_reference_series.json` carries `null` for them because
  cyclic evaluation was disabled on that path.
- Redundant certified rows have never been evaluated at any rank.

The experiment computes the exact cyclic-Schubert coefficient for EVERY
certified row at `(3,7)` and `(3,8)`, 49 and 137 rows, and cross-tabulates
against the retained and redundant partition that exact reduction produces:
4 retained and 45 redundant at rank 7, 31 and 106 at rank 8.

## What a result would and would not establish

A clean separation, coefficient one exactly on the retained set at both ranks,
would be **evidence and not proof** that the criterion could replace exact
reduction. It would be a finite-scope observation at two ranks for `N = 3`, and
promoting it to a generation or pruning rule needs both forward validity and
converse coverage proved, per the standing rule in the methodology doc. It
would justify the next theorem, not a code change.

A single redundant row with coefficient one refutes the converse outright and
vindicates the existing refusal with an exact counterexample. That is the
cheaper and more likely outcome, and it is a first-class result.

Either way this cannot become a facetness gate in code during this sprint.

## Power caveat, stated before scoring

Two ranks, one particle number, one selected cyclic coefficient rather than a
general deformed product. The cyclic test is one coefficient of one recurrence
(`docs/fixed_n3_generator_methodology.md`), so a negative result refutes THAT
criterion and says nothing about the general Belkale-Kumar or Ressayre
deformed-product criterion. A positive result is two data points.

`independent_scope` for Part B is 2 systems. The 137 rank-8 rows are not 137
independent observations: they are the certified rows of one system.

## Predictions

P1. Every RETAINED row at `(3,8)` has exact cyclic-Schubert coefficient 1.
Expectation: PASS, by analogy with the already-established rank-7 result and
the AK2008 sharpness claim. This is predictive at rank 8 because no artifact
records it.

P2. At least one REDUNDANT certified row at `(3,7)` or `(3,8)` has exact
coefficient 1, refuting "coefficient one implies facet". Expectation: PASS,
meaning the criterion fails. Rationale: the local test is one selected cyclic
coefficient, not the general deformed product, so it should not separate
facets from redundant valid Ressayre rows. A FAIL here is the more interesting
outcome and would leave the criterion standing at these two ranks.

P3. Among redundant certified rows at `(3,8)`, both coefficient-one and
non-coefficient-one values occur. Expectation: PASS. This is the sharper form
of P2: it fails if the redundant set is homogeneous in the coefficient.

P4. No certified row at either rank has exact coefficient 0. Expectation: PASS,
because a certified row has a nonzero Ressayre determinant. The cyclic
coefficient is nevertheless a different object, so a zero would be an
informative surprise rather than a contradiction.

P5. Every certified row at both ranks yields a well-formed cyclic problem, so
no row returns `degree_mismatch` or `rhs_not_exterior_weight`. Expectation:
uncertain, leaning PASS. A failure is a scope fact about the cyclic test, not
about the rows.

P6. Out of sample, for the four published systems with `N != 3`, namely
`(4,8)`, `(4,9)`, `(4,10)` and `(5,10)`, the number of distinct sorted-`tau`
shapes is at most one quarter of the facet count in each system.
Expectation: PASS, because the `N = 3` ratio at rank 10 is 16 of 93, which is
0.17. This is the only genuinely out-of-sample prediction in Part A. Note that
a bound on the `tau` ENTRIES is deliberately NOT predicted across `N`, since
`tau_i = N a_i - b` scales with the particle number.

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-ESTABLISHED per prediction, decided by
the artifact alone. Each verdict reports scope and independent scope. No
prediction is edited after the run. Part A's three scratch tables are recorded
with the verdict LANDED, which is not a score.

Resource handling is fail closed: a cyclic evaluation that exceeds its per-row
budget is recorded as `resource_limit_unresolved` and is never converted into a
coefficient claim or into evidence for either direction of P2.

## SCORED (2026-08-08, after running; predictions above are unedited)

Filled in by the run. See `docs/facet_index_stability.md`.
