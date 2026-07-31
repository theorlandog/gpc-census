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
