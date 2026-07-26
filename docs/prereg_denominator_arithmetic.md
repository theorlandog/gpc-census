# Pre-registration: denominator arithmetic (2026-07)

Written BEFORE the census-wide test was run, per the falsification protocol
(RESEARCH_ARCHIVE "SKEPTICAL AUDIT + PRE-REGISTRATION"). The worked examples
at v89 and v103 were computed first and are what motivates the definitions
below; the census-wide claim in P-DEN-2 and P-DEN-3 is the out-of-sample part
and is scored after this file was committed.

Certifying artifact: `scripts/denominator_arithmetic.py`,
`results/data/denominator_arithmetic.json`.

## Definitions (fixed before scoring)

Let a vertex have spectrum lambda with denominator `den_spec` (the smallest
positive integer with `den_spec * lambda` integral), and let the certified
state have support S and squared weights `p_t = w_t / den_state` in lowest
common terms, so `den_state` is the state denominator. Write

    ratio = den_state / den_spec.

**Incidence system.** `A_S` is the (d+1) x |S| integer matrix whose first d
rows are the mode-incidence rows (`A[m,t] = 1` iff mode m is in determinant t)
and whose last row is all ones; the right-hand side is `b = (lambda, 1)`.

**Unimodular incidence solve.** The solve is UNIMODULAR when every nonzero
invariant factor of the Smith normal form of `A_S` equals 1. Equivalently the
lattice `A_S(Z^S)` is saturated in `Z^{d+1}`, so any integral-over-`den_spec`
right-hand side that is solvable at all is solvable over `den_spec`.

**Pinned state.** A state is PINNED when the incidence solution space has
positive dimension (`dim ker A_S > 0`) and the channel conditions cut it down
to isolated points. When `dim ker A_S = 0` the weights are already determined
by incidence alone and there is nothing left for the channels to mint.

## Predictions

**P-DEN-1 (mechanism, in-sample).** At v89 the incidence solution space is a
line over `den_spec = 26`, the single silence condition restricted to that line
loses its quadratic terms, and the surviving linear equation `130 t - 1 = 0`
mints the prime 5, giving `den_state = 130` and `ratio = 5`. SCORED IN-SAMPLE;
this is the observation the definitions were built from, not evidence.

**P-DEN-2 (census-wide, out of sample).** `ratio = 1` implies the incidence
solve is unimodular. Expectation: PASS. Rationale: if the solve is not
unimodular then some invariant factor f > 1 forces a denominator beyond
`den_spec` on the incidence step alone, before any channel condition is
imposed, so the ratio cannot be 1.

**P-DEN-3 (census-wide, out of sample).** The converse: a unimodular incidence
solve implies `ratio = 1`. Expectation: FAIL. Rationale: v89 is already a
counterexample in spirit. Unimodularity governs the incidence step only; a
positive-dimensional kernel leaves the channel conditions free to introduce a
new prime through their own pivot, which is exactly the minting mechanism.
The prediction is registered so the failure is scored rather than narrated
after the fact. Expected shape of the failure: the exceptions are pinned
states, i.e. states with `dim ker A_S > 0`.

**P-DEN-4 (census-wide, out of sample).** Among states with `dim ker A_S = 0`
(nothing left to pin), the incidence solve being unimodular implies
`ratio = 1`. Expectation: PASS. This is P-DEN-3 with the escape route named in
advance and is the statement that should survive.

## Scoring rule

Each prediction is scored PASS or FAIL on the full 799-vertex ledger, with
counts and the full exception list recorded in
`results/data/denominator_arithmetic.json`. A FAIL is a result, not a defect,
and is reported with its mechanism. No prediction is edited after the run.

The two states whose certified squared weights are not rational are excluded
from the ratio predictions and reported separately, since `den_state` is not
defined for them; as of this writing that set is exactly the shipped v_B
record at (4,9) index 65 (see `docs/holonomy_rationality.md`, "The one
non-rational record").
