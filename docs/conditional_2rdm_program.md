# Conditional 2-RDM program: flux, projection, quasipinning, atlas

Status: steps 4 to 7 of the program in `docs/conditional_2rdm_body.md`.
Library `gpc_census.pinned_physics`, emitter
`scripts/conditional_2rdm_program.py`, artifact
`results/data/conditional_2rdm_program.json`, guards
`tests/test_conditional_2rdm_program.py`.

Read `conditional_2rdm_body.md` first: it proves the order-three theorem that
makes the support function an exact SDP, which is what everything here calls.

## Step 4: the energy and pair-current joint body

The unit-edge family `H(Phi)` carries a gauge-invariant cycle flux, and
`J(Phi) = -dH/dPhi` is its pair-current operator. Since
`(Tr(H rho), Tr(J rho))` is linear in rho and the elliptope is the conditional
body at order three, the joint range is a convex planar body recovered exactly
by sweeping the support function of `cos(theta) H + sin(theta) J`.

🟦 **The hidden current is not small.** Its exact range is

    [-2 sqrt(w1 w2), +2 sqrt(w1 w2)] = [-sqrt(15)/10, +sqrt(15)/10],

INDEPENDENT of the flux. The current couples only to the fluxed edge, so its
expectation is `2 Re(conj(c1) c2 * phase)` and the free relative phase attains
the bound at every flux. So states with exactly the same complete 1-RDM carry
genuinely different two-body pair currents, spanning a fixed interval of width
`sqrt(15)/5`, and no amount of one-body information narrows it.

The envelope currents `J_± = -dF^±/dPhi` are the currents carried by the
extremizing states; they vanish at `Phi = 0` and `Phi = pi` (the conjugation
symmetry `P_{-u}(-z) = P_u(z)` makes those critical points) and grow in
between.

🟦 **Branch exchanges.** The discriminant of the unit-edge sextic factors as

    -(u - 1)(u + 1) * (54675 u^8 + 513324 u^6 - 290250 u^4 + 343800 u^2 - 616225)^3.

Its real roots in `[-1, 1]` are `u = ±1` (the degenerate real-symmetric
fluxes) and two roots of the octic, at `Phi ~ 0.059287` and `Phi ~ 3.082305`.
Those interior values are where the optimizing branch changes, so they are the
geometric non-analyticities of the constrained-search functional, sitting close
to but distinct from the symmetric points.

## Step 5: projection of a physical Hamiltonian

`carrier_projection` performs the Slater-Condon projection: diagonal entries
are ordinary determinant energies, and every off-diagonal entry of this carrier
is a DOUBLE excitation, hence one signed antisymmetrized two-electron integral.
The fermionic signs are computed from occupation-string algebra rather than
assumed, and they reproduce the `+1, -1, +1` stored independently in
`complex_fixed_1rdm_observable_ranges.json`.

🟦 The projection reproduces the full CI submatrix with max absolute error
**0.0**, which is the check that matters: every interval downstream is computed
from that 3x3 matrix. An FCIDUMP round trip is included so the path from an
integral file is exercised.

### Model choice, and two traps that had to be avoided

The benchmark uses a spinless **t-V ring** (6 orbitals, 3 fermions, hopping
`t`, nearest-neighbour `V`, on-site gradient 0.35, optional boundary phase).
Two earlier choices failed for instructive reasons, and both are recorded
because either would have produced confident nonsense:

1. **Hubbard with three electrons.** The ground state is at least a spin
   doublet, so `eigh` returns an arbitrary vector of a degenerate eigenspace
   and its 1-RDM is not well defined. Measured ground-state gaps were `1e-15`.
   `hubbard_integrals` is retained with this caveat in its docstring.
2. **The symmetric ring.** Reflection symmetry forces PAIRS of equal natural
   occupations, which leaves the natural orbitals (hence the carrier)
   ambiguous inside the degenerate subspace: the weight split evenly between
   `(1,3,5)` and its partner `(2,4,5)`. The on-site gradient breaks it.

The genuinely complex benchmark freezes the rational Gaussian boundary phase
`(4+3i)/5`. It breaks time reversal without introducing an inexact model
parameter. At `V = 2` the ground-state gap is 1.435, the minimum occupation
gap is 0.00167, the off-carrier weight is 0.00466, and the gauge-invariant
carrier cycle product has imaginary part -0.00578. Thus this point is neither
a degenerate-basis artifact nor gauge-equivalent to a real carrier.

The complex run caught a basis-convention bug hidden by the real benchmark.
With this module's convention `gamma[p,q] = <a_p^dag a_q>`, a passive CI
coefficient transformation uses `U.T`, not `U.conj().T`; the physical orbital
columns used for integral transformation are `U.conj()`. The old formula
leaves a 0.169 off-diagonal 1-RDM residual at the selected point. The corrected
formula leaves 2.2e-16. A regression asserts both values qualitatively.

## Step 6: quasipinning enlargement

The bound has two independent terms, and conflating them understates the error:

    eps_H = 2 ||H|| (eps + eps^2) + ||w - w0||_1 ||y||_inf

1. The state is not exactly on the carrier. With `eps^2 = ||(I-P) Psi||^2` and
   `Phi = P Psi / ||P Psi||`, expanding `<Psi|H|Psi>` around `<Phi|H|Phi>`
   gives at most `2 ||H|| (eps + eps^2)`.
2. The projected state's fiber weights differ from the pinned ones. Because
   `F^+(w) = min { w.y : Diag(y) - H >= 0 }` is a minimum of LINEAR functions
   of `w` over a feasible set that does not depend on `w`, the interval is
   Lipschitz in `w` with constant `||y||_inf`. The dual certificate from step 3
   is therefore also the exact sensitivity of the interval to the weights, at
   no extra cost.

No quasipinning constant is assumed. The bound is stated in terms of the
measured `eps`, so it is self-contained; any theorem supplying
`eps^2 <= C_gap D(lambda)` plugs into the first term.

Here `||H||` is the centered norm of the FULL fixed-particle Hamiltonian,
equal to half its spectral width. Using only the norm of `P H P` would not
control matrix elements between the carrier and its orthogonal complement.
The complex benchmark exposed that earlier underbound, and the artifact now
records the normalized carrier energy and checks its direct deviation from the
full expectation against the corrected enlargement.

🟦 **Measured, not proved:** across the usable scan the ratio
`D(lambda) / eps^2` stays in `[1.43, 2.00]`, so the off-carrier weight is
controlled linearly by the Borland-Dennis slack. At the well-separated points
it sits at almost exactly 2. This is an empirical regularity on one model
family, not a theorem, and it is labelled that way in the artifact.

🟦 The bound holds at every point where its preconditions hold, and the
prediction is worth having: the ratio of the no-1-RDM-information width to the
pinned width runs from about **11 to 160** over the scan, so pinning plus the
1-RDM narrows the achievable range by one to two orders of magnitude.

### Preconditions are enforced, not assumed

Ten of 12 scan points are usable, including five with genuinely complex
carrier flux. A point counts only when the natural
occupations are nondegenerate AND the carrier is the dominant support. In
particular, both `V = 6` points fail because their three largest determinant
weights are not the pinned carrier. They are reported with their diagnostics
and excluded from the claim.

The selected `V = 2` complex point passes six frozen gates: nondegenerate
ground state, nondegenerate occupations, diagonal reconstructed 1-RDM,
dominant pinned carrier, nonzero gauge-invariant carrier flux, and containment
of the true energy in the quasipinning enlargement. The standalone verifier
imports no project module. It rebuilds the occupation-space Hamiltonian and
the exterior-power basis transformation, and passes 15 of 15 checks. This is
a numerical benchmark with exact model input, not an exact eigensystem.

## Step 7: the census-wide atlas

Every one of the 799 certified census supports is classified:

| carrier order | supports | method |
|---|---|---|
| `m <= 3` | 75 | exact elliptope SDP (proved) |
| `m >= 4` | 724 | elliptope outer bound |

Carrier sizes run from 1 to 14, peaking at 5 and 6. The tiering is forced by
the order-three theorem and its order-four counterexample: above order three
exactness needs a separate proof per support, and it can genuinely fail.

🟦 **Most census coherences are invisible to two-body observables.** Only 72 of
799 supports have a complete coherence graph. A pair of carrier determinants
differing in three or more orbitals cannot be connected by ANY two-body
operator, so those directions of the compatibility body are unprobeable by
two-body measurements. The atlas records, per support, the split of carrier
pairs into single excitations, double excitations, and inaccessible pairs, and
the resulting zero-width dimension `m + 2 * (inaccessible pairs)`.

That is the practically important output: it says which supports have hidden
correlation structure that no two-body observable can ever resolve, as opposed
to structure that is merely hard to measure.

## What remains open

- The `m >= 4` tier is a bound, not an answer. Special larger carriers may
  still be exactly spectrahedral, and proving that per family is the next
  mathematical step; a cycle-constraint or SOS hierarchy is the fallback.
- `D(lambda) / eps^2` bounded is measured on one family. A theorem would let
  the enlargement bound be quoted from the slack alone.
- The complex benchmark is a small flux-threaded lattice model. Extension to
  molecular magnetic fields, complex Bloch orbitals, or spin-orbit terms is
  still open.
- The atlas reports geometry, not widths, for `m >= 4`; per-support widths need
  a canonical observable choice and the outer-bound caveat attached.

## Reproducing

```sh
python scripts/conditional_2rdm_program.py
python scripts/verify_complex_flux_physical_benchmark_standalone.py
uv run pytest tests/test_conditional_2rdm_program.py
```
