# The symplectic slice atlas (2026-08-07)

Pre-registration: `docs/prereg_symplectic_slice.md`, written before the atlas
was run, with two disclosed pilot states. Artifact:
`results/data/symplectic_slices.json`. Generator:
`scripts/symplectic_slice_atlas.py`. Library: `gpc_census.symplectic.symplectic_slice`.
Guard test: `tests/test_symplectic_slice.py`.

## Why this exists

`docs/symplectic_invariants.md` reported no fiber tangent dimension at all,
because the raw kernel is vacuous: `ker d(mu)_x` has dimension at least
`2 binom(d,n) - d^2` for trivial reasons, and at the `(3,6)` Slater vertex,
whose fixed-rho fiber is provably a single point, it is **20**.

The symplectic slice is the object that fixes that. At `x = [psi]` with
`K = U(d)` and the 1-RDM as moment map,

    V_x = T_x(K.x)^omega / T_x(K_mu(x).x),      T_x(K.x)^omega = ker d(mu)_x,

which is the Marle-Guillemin-Sternberg local model, using the standard identity
`T_x(K.x)^omega and T_x(K.x) = T_x(K_mu.x)`.

## The result the pilot already gives

| state | `ker d(mu)` | `dim T_x(K_mu.x)` | **`dim V_x`** |
| --- | --- | --- | --- |
| `(3,6)#0` Slater | 20 | 0 | **20** |
| `(3,6)#3` uniform Borland-Dennis `(1/2)^6` | 19 | 19 | **0** |

Two things follow, and they point in opposite directions:

**The Borland-Dennis vertex is certified rigid.** `dim V_x = 0` means the local
model is a point: the reduced germ at that vertex IS a point, and the fixed-rho
fiber near `psi` is exactly the `K_mu`-orbit, nothing more. The raw kernel of
19 said nothing; the slice says everything. This is the first rigidity
statement in the project derived from the honest transverse object rather than
from a support-restricted amplitude count.

**And a nonzero slice does NOT mean a positive-dimensional reduced space.** At
the Slater vertex `dim V_x = 20` while the reduced space is a single point. The
20 dimensions are `K_x = U(3) x U(3)` acting on the symplectic normal to the
Grassmannian orbit, and the quotient collapses them. So the slice dimension is
a local model INPUT, not the answer. Anyone reading `dim V_x` as a fiber
dimension would be wrong at the very first state in the atlas.

## What is NOT claimed, and the shortcut that is wrong

No reduced-germ dimension is reported where the slice is nonzero. The obvious
shortcut is

    dim(V_x // K_x) = dim V_x - 2 dim(K_x . v)   for generic v,

and it is WRONG as applied to a generic `v` in `V_x`: the formula needs a
generic point of `mu_V^{-1}(0)`, and a generic point of `V_x` has a strictly
larger orbit. Measured at the Slater vertex it gives `20 - 2*15 = -10`, a
negative dimension. That is recorded here as a refuted method, not as an open
question: it is not that the number is unproven, it is that the number is
impossible.

Getting the reduced germ needs the invariant theory of the slice
representation. That is the next campaign, not this one.

Also not computed: local moment cones, Hilbert series, orbifold or singular
classification, Duistermaat-Heckman densities, and the finite component group
of `K_x` beyond the monomial elements already recorded in
`quantization_characters.json`.

## The forced gates earned their place again

Three identities are checked on every state and all are forced:
`rank d(mu)_x = dim(K.x)`; `dim ker d(mu)_x = 2(H-1) - dim(K.x)`; and
`dim T_x(K_mu.x) = dim K_mu - dim K_x`.

The first two failed during construction and caught a real bug: pivot COLUMNS
of an `rref` were used to select independent ROWS, which silently dropped
`rank d(mu)` from 18 to 17 at the Slater vertex and inflated the kernel to 21.
Every downstream slice dimension would have been wrong by one. The library
carries a comment at that line.

## Coverage, stated as partial

The pilot brief named five states: the Slater vertex, the uniform
Borland-Dennis vertex, v_B `(4,9)#65`, v89 `(3,10)#89` and v103 `(3,10)#103`.
Exact symbolic linear algebra costs `H = 126` for v_B and `H = 120` for v89 and
v103, against `H = 20` for the two rank-6 states, and the rank-9 and rank-10
states are far more expensive than that ratio suggests because the tangent
basis extraction is an `rref` on a `2H x 2H` symbolic matrix.

Whatever the atlas does not finish is listed in the artifact's `unfinished`
field and is reported as unfinished, never as absent. The predictions about
those states are scored UNINFORMATIVE, not PASS.

## Next

1. The invariant theory of the slice representation, which is what turns
   `dim V_x` into a reduced-germ dimension. The Slater vertex is the test case
   with a known answer.
2. A numerically-seeded exact route for `H` around 120, since the bottleneck is
   symbolic `rref` on the ambient tangent basis rather than anything intrinsic.
3. Only then the local moment cone, which is the input to the slice-cone facet
   oracle. It needs the slice weights under `K_mu`, not just the dimension.
