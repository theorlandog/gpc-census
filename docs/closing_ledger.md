# Closing ledger for paper-arc one (2026-07)

Every pre-registered prediction of the fiber and arithmetic programs with its
verdict, every correction with its cause, and the mortality tally. This is
both the paper's methods appendix and the honest record that the conclusions
survived their own error rate.

Nothing here is softened. The corrections are the argument that the rest can
be believed.

## 1. Scored predictions

Verdict vocabulary: PASS, FAIL, UNINFORMATIVE (antecedent empty by
construction, standing rule R3), ALREADY FALSIFIED (refuted by data in hand
before scoring, R4), NOT SCORED (no positives to score against).

### Census era, at the rank-10 closure vertices

Archive: RESEARCH_ARCHIVE "PRE-REGISTRATION SCORED: v103".

| prediction | verdict |
|---|---|
| P1 | PASS |
| P2 universal <= 2-channel interference | **FAIL**, v103 carries 5 |
| P3 state denominator = spectrum denominator | **FAIL**, ratios 5 and 5746 |
| P4 fiber dimension = incidence kernel dimension | **FAIL** at both closures |
| P5 | PASS |
| P6 | never in scope |

### The (3,11)/(3,12) bracket holdout

`docs/prereg_bracket_3_11_3_12.md`, `results/data/oos_bracket.json`. Original
verdicts left unedited; corrected readings added after the provenance audit.

| prediction | original | corrected |
|---|---|---|
| B1 unified dimension formula | PASS | PASS, SCOPE OVERSTATED |
| B2 padding invariance | PASS | PASS |
| B3 the two fibers agree | FAIL | ALREADY FALSIFIED AT RANK 10 |
| B4 no cancellation-regime member | PASS | UNINFORMATIVE |
| B5 no sparsification | FAIL | INHERITED |
| B6 numerics at new rank | PASS | PASS |

### Denominator arithmetic

`docs/prereg_denominator_arithmetic.md`,
`results/data/denominator_arithmetic.json`.

| prediction | expectation | verdict |
|---|---|---|
| P-DEN-2 ratio 1 implies unimodular | PASS | **FAIL**, 18 exceptions |
| P-DEN-3 unimodular implies ratio 1 | FAIL | **FAIL** as predicted, 1 exception |
| P-DEN-4 same among zero-kernel states | PASS | PASS, 736 states, 0 exceptions |

### The two-block family and the cross-block rule

`docs/prereg_two_block_inversion.md`, `results/data/two_block_family.json`.

| prediction | expectation | verdict |
|---|---|---|
| P-INV-1 GAP == X | FAIL | **FAIL**, 94 |
| P-INV-2 GAP == 2X | uncertain | **FAIL**, 154 |
| P-INV-3 GAP > 0 iff X > 0 | PASS | **FAIL** as a biconditional, 82 |
| P-INV-3 forward, GAP > 0 implies X > 0 | PASS | **PASS, 0 exceptions on 799** |
| P-INV-3 reverse | PASS | **FAIL**, 82 |
| P-INV-4 inversions have X > 0 | PASS | PASS |
| P-INV-5 cancellation only at two blocks | uninformative | UNINFORMATIVE |
| P-INV-6 family rigidity rule | uncertain | PASS, 4 of 4 |

### Natural-orbital gauge minimization

`docs/prereg_gauge_minimization.md`, `docs/gauge_minimization.md`.

| prediction | verdict |
|---|---|
| P-G1 fixed-rho tangent does not predict gauge sparsifiability | NOT SCORED, no positives |
| P-G2 same for fixed-spectrum tangent | NOT SCORED, no positives |

### The trichotomy diagnostic

`docs/prereg_no_support_ratio.md`, `results/data/no_support_ratio.json`.
Registered WITH DISCLOSURE that it is not blind: the aggregate counts already
existed when the prediction was written.

| prediction | expectation | verdict |
|---|---|---|
| P-NSR-1 ratio > 1 exactly on generic interference | FAIL | **FAIL**, 13 |
| P-NSR-2 ratio >= 1 everywhere | PASS | PASS |
| P-NSR-3 the flat records share a property | none | QUALITATIVE, post-hoc |

### The v_B Galois orbit

`results/data/vb_galois_orbit.json`. The last open experiment of the
arithmetic arc.

| prediction | verdict |
|---|---|
| the attaining Galois orbit has size 2 | **CONFIRMED**, both roots attain exactly |

## 2. Corrections, with cause

Ten this program, five this month. Each is a case where a claim in the
repository was wrong and structural evidence, not inspection, caught it.

| # | superseded claim | cause | artifact |
|---|---|---|---|
| 1 | v103 is first-order rigid | fixed-rho reading quoted as fixed-spectrum; the two fibers differ by the cross-block conditions | `docs/silence_rank_lemma.md` |
| 2 | v103's certified support 14 is minimal | conflated s_Q^NO with the orbit-minimal support; the support-10 endpoint is the library state in rotated orbitals and the rotation is not a gauge rotation | `results/data/orbit_check.json` |
| 3 | the census support column bounds s_Q^NO | true only for the 645 gate-passing records; the other 154 bound s_Q^free | `results/data/attainer_gate_audit.json` |
| 4 | exactly one record violates the diagonality gate | 154 do, every interference record but v89 and v103; v_B is ordinary, not a singleton | `results/data/attainer_gate_audit.json` |
| 5 | the bracket holdout scored 36 and 38 independent measurements | 34 of 38 are transports or face embeddings of 17 rank-10 donors; independent sample is 2 states | `oos_bracket.json` provenance block |
| 6 | trisected states are a cubic phase-complexity tier | a printing artifact: (pi + 3 atan t)/3 reduces to pi/3 + atan t and no cubic subextension exists anywhere in the class | `holonomy_fields.json` trisection audit |
| 7 | holonomy radicands are square roots of weight monomials | false at 15 of 82 loops; the mechanism predicts exactly where, confirmed 63/63 | `holonomy_fields.json` |
| 8 | 12 quartic holonomies | 13; 37 + 32 + 12 does not equal 82 | `holonomy_fields.json` |
| 9 | the fixed-rho inversion has sample size one | 14 states in 7 transport classes; v103 is unique in magnitude, not in kind | `two_block_family.json` |
| 10 | v103 has natural-orbital support 79 | 14. My own ledger rotated a record that already passed the gate, and with degenerate lambda eigh returns an arbitrary eigenspace basis rather than the identity | `natural_orbital_summary.json` |

Corrections 4, 5, 9 and 10 all have the same shape: a claim about the library
was carried in memory rather than read off an artifact. That is why the house
rule now reads: any convention claim about the library must cite the audit
artifact that verifies it.

Correction 10 is mine, made this month, and it was caught by scoring a
prediction I expected to fail for a different reason. That is the process
working.

## 3. Mortality: mined laws versus mechanism-first

A MINED law is a pattern read off the corpus and then promoted. A
MECHANISM-FIRST prediction is derived from structure and then tested.

**Mined laws** (P2, P3, P4, B3, B5, P-DEN-2, P-INV-1, P-INV-2, P-INV-3,
P-NSR-1, plus the archive's earlier rounds):

    10 scored, 10 failed or were downgraded.

Every headline law mined from the corpus in these programs died. The archive's
earlier estimate of "about one in two" was optimistic; over the fiber and
arithmetic arcs the mined-law survival rate is zero.

**Mechanism-first predictions**, each derived from a written mechanism before
being tested:

| prediction | mechanism | result |
|---|---|---|
| non-monomial radicands only at class size >= 3 | the quadratic formula adjoins a discriminant, not a monomial | CONFIRMED, 63 of 63 states |
| cross-block channel necessary for a two-fiber gap | the silence-rank lemma's cross-block clause | CONFIRMED, 799 of 799, 0 exceptions |
| the v_B Galois orbit has size 2 | Galois acts on the wall parameter; both roots solve the same quadratic | CONFIRMED, both attain exactly |
| ratio > 1 implies generic interference | sparsity and diagonality cannot share a basis outside design and cancellation | CONFIRMED, 141 of 141 |

    4 scored, 4 confirmed.

This is the single most useful methodological finding of the arc, and it is
what the ordering rule "mechanism before mining" is now based on. It is also
a small sample, and the honest reading is that mechanism-first predictions
were made only where a mechanism was available, which is a selection effect
and not a law about the universe.

## 4. What is proved, what is confirmed, what is open

**PROVED.**

- The silence-rank lemma, first-order, both fibers (`silence_rank_lemma.md`).
- Holonomy Theorem A, channel classes of size at most 2; the adjoined radicand
  is a weight monomial (`holonomy_rationality.md`).
- Proposition B, the reduction step at larger classes, yielding a 2-tower.
- The Rigidity-Rationality design case, covering 643 of 799 states by linear
  algebra alone (`rigidity_rationality.md`).
- Q(exp(i Phi)) is 2-elementary on all 82 loops, certified structurally at
  degree 8 where a group computation does not reach.

**CONFIRMED BY CENSUS**, with scope stated.

- 82 of 82 loop holonomies 2-elementary, degrees 1, 2, 4.
- 740 of 740 fixed-rho rigid records have rational weights, with the power
  warning: all 645 natural-orbital representatives are rigid, so the
  discriminating power lies entirely in the non-natural-orbital population.
- Cross-block necessity, 799 of 799, 403 transport classes.
- ratio > 1 implies generic interference, 141 of 141.

**OPEN.**

- Holonomy Rationality at class sizes 3 and 4 (17 classes census-wide);
  scoped rather than extended, deliberately.
- Rigidity-Rationality gap (a), isolation certificates; gap (b), uniqueness
  modulo gauge. v_B instantiates (b) with degree 2 and orbit size 2.
- Whether the interference case admits any analogue of the design-case proof.
- Every post-hoc pattern listed above, each needing an out-of-sample test that
  only Stage 1 can supply.

## 5. Independence discipline

Per standing rule R2, tallies quoted from this arc report transport classes
alongside raw records:

| tally | records | transport classes |
|---|---|---|
| census | 799 | 403 |
| two-fiber gap | 73 | 41 |
| fixed-rho inversions | 14 | 7 |
| cross-block channel present | 155 | 72 |
| bracket holdout | 38 | 2 independent states |

The bracket row is the cautionary one: 38 measurements, 2 independent states.
