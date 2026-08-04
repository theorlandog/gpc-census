# PRE-REGISTRATION: sharp candidate-44 rank-11 boundary

Date: 2026-07.

This document freezes the next rank-11 theorem target before a Schubert, SOS,
or exact-state search is interpreted. The numerical input is
`results/data/rank11_contraction.json`; the Stage-0 source is
`docs/bracket_3_11_settlement.json`.

## Frozen objects

System: `(3,11)`.

Open candidate 44:

    lambda_0 = (1/2,1/2,1/2,1/2,1/2,1/12,1/12,1/12,1/12,1/12,1/12).

Claimed level-5 inequality:

    lambda_2 + lambda_3 + lambda_4 + lambda_5 + lambda_11 <= 2.

The target violates it by exactly `1/12`. Inside the trace-three hyperplane,
the squared norm of the projected normal is `30/11`, so the exact orthogonal
distance squared is

    (1/12)^2 / (30/11) = 11/4320.

The projected spectrum is

    q = (37/72,(29/60)^4,(7/72)^5,1/15).

It has trace three and saturates the claimed inequality exactly. The normalized
real and complex contraction attacks both reach residual `11/4320` and spectrum
`q` numerically. These facts identify a target; they prove neither validity of
the inequality nor exact attainment of `q`.

## Scored fork

### Branch A: validity

PASS only if at least one of the following exact proof objects ships:

1. A finite Schubert-coefficient computation whose independently checkable
   integer output proves the claimed inequality.
2. An exact rational SOS dual identity with a strictly positive lower bound,
   together with exact positive-semidefiniteness checks for every block.

A floating SDP objective, a rounded dual without exact reconstruction, or a
positive multistart floor is not a PASS.

### Branch B: falsification

PASS only if an exact pure three-fermion state is produced at `lambda_0`, with
unit norm and the target 1-RDM characteristic polynomial verified symbolically.
Such a state falsifies the claimed inequality.

A smaller numerical residual, including machine precision without exact
recognition, is not a PASS.

### Boundary witness q

Independently of the fork, `q` becomes an exact attained boundary witness only
when a symbolic state passes unit norm and characteristic-polynomial checks.
The existing numerical contraction state is dense and does not count. Because
the one-hop-free design classifier rejects `q`, any sparse exact witness is
expected to require interference, but that classification is not a proof of
attainment.

## Status and stopping rule

Until Branch A or Branch B passes, candidate 44 remains `OPEN` and the claimed
inequality remains `CLAIMED-NOT-PROVED`. Candidate 44 alone cannot establish
the complete `(3,11)` polytope: adding a valid cut creates new intersections
that must be enumerated and certified.

The first theorem-grade milestone is Branch A plus an exact boundary witness at
`q`, or Branch B. If neither exact object is obtained, publish the normalized
audit and the failed proof/construction attempt without changing the bracket
verdict.

## OUTCOME (2026-08)

Branch A passes, with a recorded deviation. Branch B is dead.

**Result.** The claimed level-5 inequality
`lambda_2 + lambda_3 + lambda_4 + lambda_5 + lambda_11 <= 2`, and the three
sibling rows AK-RMK-4.2.1 and KLY09-EXT-2/3, carry exact Ressayre
certificates. Candidate 44 violates KLY09-EXT-1 by `1/12`, so it is REFUTED,
and with the other three the (3,11) bracket closes at 19 TRUE + 31 REFUTED +
0 OPEN. See `docs/level5_ressayre_3_11.md`,
`results/data/level5_ressayre_3_11.json`.

**Deviation from the frozen fork.** Branch A enumerated two acceptable proof
objects, a Schubert-coefficient computation and an exact rational SOS dual.
The object delivered is neither: it is an exact Ressayre certificate. It meets
what the branch was written to demand, a finite exact computation with
independently checkable integer output, and it is checked by a
standard-library verifier that shares no code with the generator
(`scripts/verify_level5_ressayre_3_11_standalone.py`). But it is a third route
and is recorded as such rather than reported as an unqualified PASS of the
clause as written.

**Branch B.** Dead on both legs. KLY09-EXT-1 is now proved, so no state
attains `lambda_0`; and the design route was independently closed by exact
counting before that (`docs/rank11_open_candidate_designs.md`), which had
already ruled out any one-hop-free witness.

**Boundary witness q.** Still not attained. `q` saturates a now-proved
inequality, so it is a genuine boundary point of the polytope rather than a
target outside it, but no symbolic state has passed the unit-norm and
characteristic-polynomial checks. That part of the preregistration is
unresolved and stays open.

**Stopping rule, honoured.** The document's own caveat stands unchanged:
candidate 44 alone does not establish the complete (3,11) polytope. Adding
these valid cuts creates new intersections that must still be enumerated and
certified, so Stage 1 or a tighter bracket remains required.
