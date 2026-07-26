# The attainment criterion, and a verification record (2026-07)

Companion to `docs/gauge_minimization.md` and `docs/cascade_census.md`. Those
say what the results are. This says what an acceptance test has to check for a
support claim to mean anything, and records which claims from two rounds of
correction survived independent checking.

## The criterion

A state ATTAINS a vertex lambda, in the census convention, when

1. its norm is exactly 1,
2. its 1-RDM is exactly DIAGONAL, and
3. that diagonal equals lambda entrywise.

A characteristic-polynomial identity establishes none of this. det(rho - xI) is
invariant under conjugation, so it holds for every state in the U(d) orbit,
including states whose 1-RDM has been rotated off the diagonal. It is a
consistency check and must never be the acceptance criterion.

Added on top, for anything offered as a NEW state: its 2-RDM spectrum must
differ from the library state's. The 2-RDM transforms by conjugation under a
one-body rotation, so its spectrum is a basis-free fingerprint; a match means
the same state in different orbitals, which is a gauge statement and has to be
reported as one.

Enforced in `scripts/exactify_cascade.py`, pinned by
`tests/test_attainment_criterion.py`, which feeds the checker the support-10
state published for v103 and expects REJECTION.

## Why this was not caught earlier: the library is not uniform

Recomputed from scratch over all 799 certified states
(`gpc_census.cascade.diagonality_error`):

| population | count | 1-RDM |
|---|---|---|
| DESIGN-INT | 631 | exactly diagonal, diagonal = lambda |
| DESIGN-REAL | 12 | exactly diagonal, diagonal = lambda |
| INTERFERENCE, v89 and v103 | 2 | exactly diagonal, diagonal = lambda |
| INTERFERENCE, the rest | 154 | NOT diagonal, worst off-diagonal 0.2222 |

645 of 799 library states attain in the diagonal sense; 154 do not, and their
recorded supports are counted in a basis that is not their natural-orbital
basis. This is a pre-existing property of the library, not something any
cascade introduced, and it is why a spectrum-level check looked adequate for
so long.

## The two loci, in code

`gpc_census.cascade.cascade` now takes a named locus and refuses to guess:

- `natural_orbital` holds rho = diag(lambda) and bounds s_Q^NO. A state whose
  library 1-RDM is not diagonal is not on this locus and the run reports
  `start_off_locus` rather than a number. Census-wide: ZERO of 799 sparsify,
  154 cannot start.
- `fixed_spectrum` holds the spectrum only and bounds s_Q^free. Census-wide:
  71 sparsify, and none of the 71 endpoints is diagonal.

Re-running the exact checker under the corrected criterion:

    certified s_Q^NO attainers                     0 of 71
    exact but spectrum-only (rho not diagonal)    69 of 71
    weights not rational (numerical only)          2 of 71

Cross-referenced with the orbit gate: of those 69, 63 are the library state in
rotated orbitals and 6 are genuinely new states giving vertex-level s_Q^free
bounds.

## Verification record

Two rounds of correction arrived at the right conclusion. Several of the
specific supporting claims did not survive checking, and the difference matters
because it changes what needs fixing.

HELD:

- The acceptance criterion is conjugation invariant and cannot see a rotated
  1-RDM. Confirmed.
- The published v103 support-10 state has rho_{3,9} = -5/68 and
  rho_{1,6} = -sqrt(42)/221, with a diagonal strictly majorized by lambda.
  Confirmed in exact arithmetic.
- s_Q^NO and s_Q^free are different and s_I bounds only the first. Confirmed.
- The natural-orbital gauge prod_k U(m_k) preserves attainment exactly and had
  never been minimized over. Confirmed: random block rotations hold
  rho = diag(lambda) to 1e-16.
- Support-restricted first-order analysis is blind to gauge sparsification.
  Confirmed: a random block rotation at v103 takes its support from 14 to 100.

DID NOT HOLD:

- "Every one of the 799 library states is exactly diagonal." False, 154 are
  not (table above).
- "70 of the 71 improved states END non-diagonal; exactly one ends diagonal,
  v103 at 1.3e-17." That figure is `rdm_max_offdiagonal` in `cascade.jsonl`,
  which is computed on the LIBRARY state BEFORE the cascade runs. Corrected:
  0 of 71 endpoints are diagonal, 70 of 71 were already non-diagonal at the
  start, and v103 is the only one that started diagonal. It ends at
  5/68 = 0.0735.
- "The recognized rationals are not the endpoint's weights; the integer
  relation step landed on a different state." False. The recognized state's
  maximum off-diagonal is exactly 5/68 = 0.073529, matching the numerical
  endpoint's 0.0735. The recognition was correct; accepting a non-diagonal
  state as an attainer was the error.
- "s_Q^NO(v103) <= 10 stands because the endpoint is diagonal." Rests on the
  same misreading. The endpoint is not diagonal, the diagonality-constrained
  cascade finds no improvement at v103 (14 to 14), and gauge minimization finds
  none either. The bound is neither established nor refuted.
- "Not one endpoint is a new state." Holds for 63 of 71; 8 have distinct 2-RDM
  spectra (3.4e-3 to 4.7e-2) and are genuinely new states.

## What stands from the earlier cascade work

The multiplicity gate and the post-drop polish, both of which are real: the
minimal polynomial constrains only the eigenvalue SET, and 33 states are
stopped by a continuation that satisfies it while landing on different
multiplicities.
