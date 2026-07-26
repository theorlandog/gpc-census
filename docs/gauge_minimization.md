# Natural-orbital gauge minimization of the census library (2026-07)

Artifacts: `results/data/gauge_min.jsonl` (one record per vertex),
`results/data/gauge_min_summary.json`. Engine: `src/gpc_census/gauge.py`.
Driver: `scripts/gauge_census.py`. Tests: `tests/test_gauge.py`.
Pre-registration: `docs/prereg_gauge_minimization.md` (committed before the run).

## The question

Natural orbitals are unique only up to unitary rotation WITHIN degenerate
occupation eigenspaces. Every census vertex has a degenerate spectrum, so every
certified state carries a gauge freedom

    G(lambda) = U(m_1) x ... x U(m_k),   m_k the eigenvalue multiplicities,

acting through the N-th compound of the block-diagonal unitary. Every point of
that orbit has 1-RDM exactly diag(lambda) (block-diagonal U commutes with
diag(lambda) identically, so this needs no solver and cannot drift), hence every
point is an equally legitimate natural-orbital representative of the same
vertex. The support the census library records is the support of ONE point of
that family, and it had never been minimized over the family. This run does the
minimization.

The quantity is s_Q^NO, the minimum support over states with rho EXACTLY
diag(lambda). It is the census quantity, the one bounded below by s_I, and it is
NOT the free-basis minimum s_Q^free that the sparsification cascade bounds (see
`docs/cascade_census.md`, and note that the two must never appear in one
inequality).

## Result 1: the library column does not move, at all

    states searched                              799
    states whose support the gauge reduces         0
    amplitudes removed                             0

Not one of the 799 certified states admits a support-reducing rotation in its
own natural-orbital gauge family, under a multi-start exact-closure search. The
library column is therefore not the loose upper bound it was expected to be. It
was never minimized over the gauge, and now that it has been, it is unchanged.

## Result 2: 629 of the 799 are EXACTLY gauge-minimal, certified, with no search

A block-diagonal unitary preserves each determinant's BLOCK PROFILE
(|T ∩ B_1|, ..., |T ∩ B_k|). Its compound therefore preserves the decomposition
into profile classes and the norm each class carries, so for any gauge
representative

    support >= number of occupied profile classes,

refined, since the gauge acts inside a class as a TENSOR PRODUCT of exterior
powers, to

    support >= sum over occupied classes of (best flattening rank of the class).

Both are gauge invariants and both are computed exactly. Where the library
support attains the bound, that state is gauge-minimal with no optimization at
all, and the certificate is a rank computation rather than a search.

    certified gauge-minimal          629 of 799
      DESIGN-INT                     510
      INTERFERENCE                   114
      DESIGN-REAL                      5
    bound not attained (open)        170

For the remaining 170 the search found no reduction but the lower bound is not
attained, so minimality there is UNPROVEN, not proved. Of those 170, 51 have a
class where the bound falls back to the trivial value 1 because the single
active factor sits in Lambda^a with 3 <= a <= m-3, where no cheap invariant is
implemented; those are the weakest rows in the table and the obvious place to
sharpen the bound.

## Why the null result is worth believing (and how far)

A search that finds nothing is worth exactly what its power is worth, so the
power is pinned by a test rather than asserted. Three determinants sharing a
profile class whose flattening has rank 1 form a rank-1 tensor, so a reduction
to a SINGLE determinant provably exists; the search finds it, exactly, at
closure residual 4e-16 (`test_greedy_gauge_drop_recovers_a_known_reduction`).
The matching control is a rank-2 flattening, where a reduction provably does not
exist and the search correctly reports none.

Two earlier methods failed that same test and were removed rather than shipped:
a reweighted-L1 Jacobi sweep, and a Riemannian descent with a smoothed-L1
surrogate whose smoothing made growing a zero amplitude free, so it walked away
from the sparsest point it was handed. The regression for the latter is
`test_descent_does_not_walk_away_from_a_sparse_point`.

What the search actually does is ask, for each support determinant in turn,
whether a block-diagonal U annihilates its amplitude, and answer by a
re-anchored Gauss-Newton whose Jacobian is the EXACT derivative of the compound
action (columns Gamma(U) dGamma(iE_a) psi_0). Success lands on the orbit itself,
not on a truncation of it, so an accepted drop keeps rho = diag(lambda)
identically rather than to a tolerance. The gauge groups involved are not small:
median dimension 27, maximum 100, and 162 of the 170 open cases carry dimension
at least 25.

Limits, stated plainly: multi-start greedy over a curved 27-to-100-dimensional
group is not exhaustive, the drop order is greedy, and a reduction requiring two
simultaneous cancellations that no single-determinant step reaches would be
missed. The 629 certified rows do not depend on the search at all; the 170 open
rows do.

## Scoring the pre-registration

The outcome variable DROP is identically zero across all 799 states, which is
the case `docs/prereg_gauge_minimization.md` named in advance under "What would
make any of this moot". Consequently:

| prediction | status |
|---|---|
| P-G1 fixed-rho tangent does not predict gauge sparsifiability | NOT SCORED (no positives) |
| P-G2 same for fixed-spectrum tangent | NOT SCORED (no positives) |
| P-G3 cascade drop count differs from gauge drop count | PASS, degenerately: 71 states have cascade drops and 0 have gauge drops, so the two counts agree nowhere they could disagree; this is not evidence about the mechanism |
| P-G4 EXCESS > 0 predicts DROP > 0 for a majority | FAIL. 170 states have EXCESS > 0 and none sparsifies |
| P-G5 cancellation-regime membership does not predict | consistent, uninformative: v89 and v103 both DROP 0, like every other state |

P-G4 is a real, if unglamorous, falsification: the slack between the support and
the profile/flattening bound predicts nothing about reducibility, which says the
bound is loose on those 170 rather than that they are reducible. The
pre-registration called P-G4 the weakest of the five before the run, for the
reason that turned out to matter.

## What this settles about v103

    library support                    14
    gauge lower bound                   9
    gauge-minimized support            14   (no reduction found)

so 9 <= (minimum over v103's library gauge family) <= 14, and
s_I(v103) = 6 <= s_Q^NO(v103) <= 14.

The support-10 state of `docs/cascade_census.md` does NOT improve this. Its
1-RDM is not diagonal (exactly, rho_39 = -5/68), so it is not a natural-orbital
representative; it is the certified library state seen in different orbitals,
through a U(10) rotation with off-block entries 0.196 that is not an element of
U(4) x U(6). It bounds s_Q^free and the orbit-minimal support of that one state,
not s_Q^NO.

The first-order picture is consistent with all of this and explains why nobody
saw it earlier: at the v103 certified state all 42 off-diagonal U(4) x U(6)
generators move the state off its support, so the gauge orbit is invisible to
every support-restricted tangent computation in this project. That is a general
fact about degenerate spectra, recorded now in `docs/silence_rank_lemma.md` and
`docs/unified_fiber_dimension.md`.

## Caveats

- The 629 certifications are minimality within the state's OWN gauge family.
  s_Q^NO(lambda) is a vertex-level quantity and could still be smaller at a
  different fiber component; nothing here bounds that below except s_I.
- The 170 open rows are open in both directions.
- Everything is at double precision except the lower bounds, which are integer
  ranks of exactly-built matrices, and the v103 diagonality facts, which are
  exact rationals.
