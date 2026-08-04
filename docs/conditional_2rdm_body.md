# The conditional 2-RDM compatibility body over a pinned fiber

Status: 🟦 order-three theorem implemented and certified; exact primal-dual
certificates for four carriers; order-four boundary established by explicit
counterexample. Library `gpc_census.elliptope`, emitter
`scripts/conditional_2rdm_body.py`, artifact
`results/data/conditional_2rdm_body.json`, guards
`tests/test_conditional_2rdm_body.py`.

## The reframing

`docs/complex_fixed_1rdm_observable_ranges.md` computes an exact interval for
two triangle Hamiltonians via a sextic. That undersells what is happening.
Fix a 1-RDM `gamma` on a pinned face and let

    K(gamma) = conv { Gamma^(2)_Psi : Psi |-> gamma }

be the convex body of 2-RDMs compatible with it. Then

    F_H^+(gamma) = max over K(gamma) of Tr(H Gamma),   F_H^- = -F_{-H}^+,

so the endpoint pair is the SUPPORT FUNCTION of `K(gamma)`, and

    K(gamma) = intersection over H of { Gamma : Tr(H Gamma) <= F_H^+(gamma) }.

The interval solver is therefore an exact support oracle for the conditional
2-RDM body. The claim to lead with is not "we solved two triangle
Hamiltonians"; it is that for a nontrivial GPC-pinned fermionic fiber the
complete convex body of two-body information left undetermined by the 1-RDM is
computable, hence so is the exact range of every two-body observable. That is
an exact conditional N-representability statement.

## The order-three theorem

On this fiber the carrier is three determinants with weights
`w = (3/5, 1/4, 3/20)` fixed by `gamma`, the state is `c` with `|c_i|^2 = w_i`,
and its carrier density matrix `rho = c c^dagger` has `diag(rho) = w`. Because
the collision-free two-body channels expose all three carrier coherences, the
variable part of `K(gamma)` is an affine image of the scaled complex elliptope

    E_w = { rho >= 0 : diag(rho) = w }.

**Theorem (order three).** `E_w = conv { c c^dagger : |c_i|^2 = w_i }`.

Rescaling by `Diag(sqrt(w))` turns `E_w` into the complex correlation-matrix
body. An extreme point of rank `r` requires its `m` blocks `v_i v_i^dagger` to
be linearly independent over R inside the `r^2`-dimensional real space of
Hermitian `r`-by-`r` matrices, so `r^2 <= m` (Li and Tam). At `m = 3` this
forces `r = 1`: every extreme point is a phase matrix, so the spectrahedron is
exactly the convex hull of the fiber's own states.

Hence the support function is an SDP, exactly and not as a relaxation:

    F_H^+ = max { Tr(H rho) : rho >= 0, diag(rho) = w }
          = min { w.y : Diag(y) - H >= 0 }        (dual; Slater holds at Diag(w))
    F_H^- = max { w.y : H - Diag(y) >= 0 }.

Verified operationally: at order three the SDP value and the phase-hull value
agree for all four carriers (`sdp_versus_phase_hull` in the artifact).

## The sextic is the KKT elimination ideal

This is the concrete payoff of the reframing, and it is checked rather than
asserted. At the optimum `M = Diag(y) - H` is singular with kernel `c`, and
stationarity forces `cof_ii(M)` proportional to `w_i`. Eliminating `y` from

    det M = 0,   w_1 cof_00 = w_0 cof_11,   w_2 cof_11 = w_1 cof_22

with objective `t = w.y - E0` yields a univariate polynomial. For the two
complex regressions it reproduces the shipped sextics EXACTLY, up to the
primitive-content scalar (`test_kkt_elimination_reproduces_the_shipped_sextics`).

So the sextic is not a bespoke trigonometric trick; it is the algebraic
solution of a universal 3x3 SDP. The reframing also explains a fact the
original derivation had to discover: only `|h_ij|^2` and `Re(chi)` enter the
KKT system, with `chi = h01 h12 conj(h02)`. That is why the range depends on
the carrier only through its edge magnitudes and the gauge-invariant cycle
flux, and why the unit-edge family depends on `Phi` only through `cos Phi`.

## Dual certificates, and why they are better

The Groebner basis is in shape position, so each `y_i` is a rational
polynomial in `t` and `t` is a root of one univariate polynomial. Optimality
is then witnessed by PSD-ness of `Diag(y) - H`, read off exactly from the
signs of the Hermitian principal minors at the isolated root. Two gains over
the subresultant-and-Sturm route:

1. **Pointwise optimality.** Selecting "the least and greatest real roots"
   needs an argument that every real root was found. A dual PSD witness proves
   optimality at the point, with no exhaustiveness obligation. Empirically the
   selection is clean: of the four real roots in each case, exactly one has all
   minors nonnegative (the maximum), exactly one has the odd-order minors
   flipped (the minimum), and the interior roots certify neither.
2. **Degenerate flux is covered.** Zero and pi flux have `Im(chi) = 0`, which
   the sextic branch engine must exclude as nongeneric and hand to the
   real-symmetric solver. The dual route handles them uniformly. At pi flux the
   optimal dual matrix drops rank and all three 2x2 minors vanish, which the
   sign-at-root machinery reports as sign zero rather than failing.

Certified endpoints, all matching the independently known values:

| carrier | cycle flux | range |
|---|---|---|
| zero flux | `chi = 1` | `[-1, 3(2+sqrt 15)/10]` |
| pi flux | `chi = -1` | `[-3(2+sqrt 15)/10, 1]` |
| quarter flux | `chi = i` | `[-rho, rho]`, `rho ~ 1.545782759984` |
| Pythagorean | `chi = 3/5 + 4i/5` | `~[-1.346174378380, 1.685073545053]` |

The first two reproduce the closed forms stored in
`fixed_1rdm_observable_ranges.json`; the last two reproduce the sextic
endpoints in `complex_fixed_1rdm_observable_ranges.json`.

## Where the theorem stops

At `m = 4` the bound permits `r = 2`, so the order-three argument must not be
extrapolated. `rank_two_extreme_point_4x4` builds four unit vectors in C^2
whose rank-one blocks are the four Pauli directions, hence linearly
independent over R: the resulting correlation matrix has rank two and IS
extreme. Taking `H` to be minus the projector onto the orthogonal complement
of its range gives a witness where

    SDP maximum 0,   true phase-hull maximum -0.0694326889,

a relaxation gap of about `6.94e-2`. So for four or more carrier determinants
the SDP is a genuine OUTER BOUND on the conditional body, not an identity.
This is a real mathematical boundary, recorded rather than papered over.

## The practical census quantity

`hidden_width(H) = F_H^+ - F_H^-` is the two-body information the 1-RDM leaves
undetermined; half of it is the least worst-case error of any estimate of
`<H>` that sees only one-body data. Zero width means the observable is exactly
identifiable from the 1-RDM. The four carriers here have widths 2.7619
(zero and pi flux), 3.0916 (quarter), 3.0312 (Pythagorean), so the quarter-flux
carrier hides the most.

## What this does NOT yet do

Deliberately out of scope here, in the order it should be attacked:

- the flux-dependent energy/pair-current joint body, via the support function
  of `H(Phi) + s J(Phi)` swept over `s`;
- Slater-Condon/FCIDUMP projection onto a certified carrier, so the intervals
  attach to a recognizable electronic Hamiltonian;
- a quasipinning enlargement theorem, without which the exact solver applies
  only to exactly pinned inputs;
- the census-wide `conditional_2rdm_body_atlas` over all certified supports,
  which needs the `m >= 4` regime and therefore an honest exact/relaxed label
  per support.

The `m <= 3` exactness proved here is what licenses the atlas's first tier;
everything above order three needs its own exactness argument or ships as a
bound.

## Reproducing

```sh
python scripts/conditional_2rdm_body.py
uv run pytest tests/test_conditional_2rdm_body.py
```
