# PRE-REGISTRATION: the symplectic slice atlas (2026-08-07)

Written BEFORE the atlas script was run, with two disclosed pilot states.
Standing rules: `docs/prereg_template.md`. Certifying artifact:
`scripts/symplectic_slice_atlas.py` -> `results/data/symplectic_slices.json`.

## What is being tested

`docs/symplectic_invariants.md` deliberately reported no fiber tangent
dimension, because at these ranks `ker d(mu)_x` has dimension at least
`2 binom(d,n) - d^2` for trivial reasons: at the `(3,6)` Slater vertex, whose
fixed-rho fiber is provably a single point, it is 20.

The symplectic slice is the object that fixes this. At `x = [psi]` with
`K = U(d)` and the 1-RDM as moment map,

    V_x = T_x(K.x)^omega / T_x(K_mu(x).x),   T_x(K.x)^omega = ker d(mu)_x,

and `T_x(K.x)^omega and T_x(K.x) = T_x(K_mu.x)` is the standard identity that
makes this the Marle-Guillemin-Sternberg local model. `dim V_x = 0` certifies
that the reduced germ is a point outright.

## Forced gates, recorded as code checks and not as evidence

G1. `rank d(mu)_x = dim(K.x)`, because `im d(mu)_x = ann(k_x)`.
G2. `dim ker d(mu)_x = 2(H-1) - dim(K.x)`, the same statement.
G3. `dim T_x(K_mu.x) = dim K_mu - dim K_x`, because `K_x` is contained in
    `K_mu`.

All three are forced. They are scored only as gates. G1 and G2 already failed
once during construction and caught a real bug: pivot COLUMNS of an rref were
used to select independent ROWS, which silently dropped `rank d(mu)` by one at
the Slater vertex.

## Disclosed pilot

Two states were computed before this file was written, to check the
construction against a case with a known answer: `(3,6)#0` (Slater), giving
`ker d(mu) = 20` and `slice_dim = 20`, and `(3,6)#3` (uniform Borland-Dennis),
giving `ker d(mu) = 19` and `slice_dim = 0`. No other state was examined.

## What is deliberately NOT claimed

No reduced-germ dimension where the slice is nonzero. The natural shortcut,
`dim V_x - 2 dim(K_x . v)` at a random `v` in `V_x`, is WRONG: the formula
requires a generic point of `mu_V^{-1}(0)`, and a generic point of `V_x` has a
strictly larger orbit. Applied at the Slater vertex it returns **-10**, which
is impossible, and that is the evidence it is wrong rather than merely
unproven. Getting the reduced germ needs the invariant theory of the slice
representation, which this campaign does not attempt.

Nothing here computes a local moment cone, a Hilbert series, an
orbifold/singular classification, or a Duistermaat-Heckman density.

## Predictions

P1. G1, G2, G3 hold on every state the atlas completes. Expectation: PASS.
    Forced; a FAIL means the code is wrong.

P2. `slice_dim = 0` at `(3,6)#3`, the uniform Borland-Dennis vertex, so its
    reduced germ is certified a point. Expectation: PASS. Already observed in
    the disclosed pilot, so this is recorded, not discovered.

P3. `slice_dim > 0` at `(3,6)#0`, the Slater vertex, even though its reduced
    space IS a point. Expectation: PASS, and it is the point of the campaign:
    a nonzero slice does NOT mean a positive-dimensional reduced space, so the
    slice dimension is a local model input and not itself the answer. Already
    observed in the pilot.

P4. The atlas completes all five pilot states within the compute budget.
    Expectation: uncertain. Exact symbolic linear algebra at `H = 126` for v_B
    and `H = 120` for v89 and v103 may exceed it. Unfinished states are named
    in the artifact's `unfinished` list rather than dropped, and a partial
    atlas is reported as partial.

P5. At v89 and v103, `slice_dim > 0`. Expectation: PASS, weakly held. Both are
    interference states with small stabilizers (`dim Stab = 2`, the census
    minimum), so their `K_mu` orbit is far from filling `ker d(mu)`. A
    `slice_dim = 0` at either would be a strong and surprising rigidity
    statement.

## Scoring rules

PASS / FAIL / UNINFORMATIVE per prediction, decided by the artifact alone.
Exact arithmetic throughout. A prediction about a state the atlas does not
finish is UNINFORMATIVE, never PASS.

## SCORED (2026-08-07, after running; predictions above are unedited)

See `docs/symplectic_slice.md`.
