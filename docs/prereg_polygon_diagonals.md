# Pre-registration: the polygon diagonal test (Track B1)

Standing rules in `docs/prereg_template.md`. Committed together with
`scripts/polygon_diagonals.py` BEFORE the script is run (rule R1); the commit
hash is the evidence the predictions are unedited.

Target: `docs/holonomy_two_elementary.md` condition (R1),
`deg_Q(cos^2 Phi) <= 2`, at one-hop class sizes 3 and 4. Proposition B gives a
2-tower there, which is strictly weaker than multiquadratic, and the open
problem is exactly the gap between 2-tower and 2-elementary. This measurement
aims at a MECHANISM that forbids D4, not at more reductions.

## Scope, re-derived before writing

`results/data/holonomy_fields.json` and `docs/holonomy_rationality.md`, both
re-read at `8eff55b`:

- census-wide one-hop class sizes over all 799 records: 123 of size 1, 65 of
  size 2, 16 of size 3, 1 of size 4, so 188 of 205 classes are covered by
  Theorem A and **17 are not**;
- restricted to the 63 loopy states the distribution is 2 / 65 / 16 / 1, so
  every class of size at least 2 lies in a loopy state;
- 82 loops, degrees 37 / 32 / 13, Galois groups C1 37, C2 32, V4 13, all
  2-elementary; (R1) holds 82 of 82 and (R2) holds on 13 of 13 quartics;
- `theorem_A_prime_test.by_class_size` gives, as `size,cos^2 rational`:
  `2,True` 47, `3,False` 15, `3,True` 15, `4,True` 5.

Per rule R2 the 17 large classes are RECORDS, not independent states. Transport
classes are reported alongside the raw count in the scorecard.

## The reduction, proved before measuring

Not a prediction; it is proved in the script's module docstring and pinned by a
test. For a two-element subset `T = {a, b}`,

    |z_a + z_b|^2 = r_a^2 + r_b^2 + 2 t_ab,    t_ab = s_a s_b r_a r_b cos(Phi_ab),

so with `r_a^2` rational under (H1) the diagonal is rational iff `t_ab` is.
Because `(r_a r_b)^2 = p_{u_a} p_{v_a} p_{u_b} p_{v_b}` is rational,

    t_ab^2 = (rational) * cos^2(Phi_ab).

So a rational diagonal FORCES `cos^2(Phi_ab)` rational, which is (R1) for the
within-class holonomy. The polygon question is therefore (R1) on the
within-class difference loops plus the extra demand that the resulting rational
be a square. This is recorded here so the measurement cannot later be presented
as independent evidence for (R1) when it is partly a restatement of it.

## Predictions, committed before the run

| id | prediction | expectation |
|---|---|---|
| P-B1-1 | at least one proper partial sum, in at least one class of size >= 3, has IRRATIONAL squared modulus | IRRATIONAL FOUND |
| P-B1-2 | every MODE-COHERENT proper subset has rational squared modulus | PASS |
| P-B1-3 | every mode-coherent two-element subset's partial sum agrees with the independently computed 2-RDM element `rho^(2)_{AC,BC}` | PASS |
| P-B1-4 | no mode-INCOHERENT proper subset has rational squared modulus | PASS |
| P-B1-5 | every irrational squared modulus has degree exactly 2 over Q | PASS |
| P-B1-6 | `t_ab` rational and `cos^2(Phi_ab)` rational agree pair by pair, except where the rational is a non-square | PASS |

Reasoning behind each, so a hit cannot be read as a lucky guess:

- **P-B1-1.** `holonomy_fields.json` records 15 of 82 loops whose radicand lies
  OUTSIDE the weight group, and the mechanism (Proposition B's corollary)
  predicts that this happens only from class size 3 up, because the quadratic
  formula adjoins a discriminant rather than a monomial. If every diagonal were
  rational the size-3 discriminants would be rational and no non-monomial
  radicand could appear. So irrational diagonals must exist somewhere.
- **P-B1-2 and P-B1-3.** If `T = {a : C in u_a}` then every determinant in the
  subset contains the spectator mode `C`, and the partial sum is, up to one
  overall sign, the 2-RDM entry `rho^(2)_{AC,BC}`. Hypothesis (H2) makes 1-RDM
  channel targets rational for the same structural reason; the prediction is
  that this extends to the 2-RDM entries reached by a spectator mode. This is
  the mechanism the whole track is hunting: if it holds, a size-3 class
  triangulates along a rational diagonal and multiquadratic follows, with
  size 4 by induction on diagonals.
- **P-B1-4.** A subset with no spectator mode is not a matrix element of
  anything, so there is no structural reason for its modulus to be rational.
- **P-B1-5.** By the reduction, an irrational `|z_a + z_b|^2` differs from a
  rational by `2 t_ab` with `t_ab^2` rational, so its degree is 2 unless
  `cos^2` is itself irrational, which the census says happens on 15 loops.
  Degree 4 is therefore possible and the prediction is the sharper one.
- **P-B1-6.** Consistency check on the reduction itself. A disagreement means
  the exact arithmetic is wrong, not that the mathematics is interesting.

**P-B1-2 is the prediction that matters.** A PASS hands over the mechanism and
turns the size-3 case into triangulation. A FAIL kills the cheapest route to
(R1) at large classes and makes the recorded degrees and fields of the
diagonals the induction datum instead, which sharpens the conjecture rather
than killing it. Either outcome is a finished result; the failure mode this
prereg exists to prevent is reporting whichever branch occurred as though it
had been anticipated.

## What would REFUTE the whole track

A mode-coherent subset whose partial sum is irrational AND whose 2-RDM
identification verifies. That would show the spectator-mode decomposition is
real but arithmetically useless, and B1 would return BOUNDED NULL with the
search space being all 17 classes and all their proper subsets.

## Nonclaims

- This measures the SHIPPED representative of each record, not the vertex. The
  154 records that fail the diagonality gate attain their spectrum in another
  basis, and `states.jsonl` is used throughout precisely because the
  natural-orbital ledger is numerical and would destroy the exact closed forms
  the whole arithmetic layer depends on.
- Rationality of a diagonal at every measured class would NOT prove (R1) at
  class sizes 3 and 4 census-wide. It would prove it on the 17 measured
  classes and supply the mechanism for a proof, which then has to be written.
- Nothing here extends Theorem A by bookkeeping, and nothing here licenses
  reading holonomies on a non-integer loop basis.
