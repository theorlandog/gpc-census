# The Silence-Rank Lemma (2026-07)

Status: PROVED below (first-order statements; elementary degenerate
perturbation theory). Dimension formula VERIFIED at both census
cancellation states. Supersedes the P4 scorecard's rigidity reading for
v103; confirms it for v89.

## Setting

Fix a system (N, d), a vertex spectrum λ with distinct values
λ^(1) > ... > λ^(k) on orbital blocks B_1, ..., B_k, and a state
ψ = Σ_{t ∈ S} c_t |t⟩ with support S whose 1-RDM ρ is EXACTLY DIAGONAL
with diag(ρ) = λ (the cancellation regime: every one-hop channel of S is
silent). Assume ψ real (both census instances are; the complex case
changes only bookkeeping). Write:

- A_S : R^S → R^d, the incidence map (A_S p)_o = Σ_{t ∋ o} p_t, and
  K = ker A_S ∩ {Σ δp = 0} the incidence kernel;
- for a channel γ = (p, q) (an orbital pair realized by at least one
  one-hop determinant pair in S), σ_γ(c) = ρ_{pq}(c), the signed channel
  sum. Silence: σ_γ(c) = 0 for all γ. A channel is WITHIN-BLOCK if p, q
  lie in the same B_i, CROSS-BLOCK otherwise;
- J_within, J_all: the Jacobians of the within-block (resp. all) channel
  sums with respect to weight deformations δp at fixed phases
  (δc_t = δp_t / (2 c_t)).

## Scope: every statement here is SUPPORT-RESTRICTED

Added 2026-07. All three parts of the lemma, and every fiber dimension derived
from them anywhere in this project, describe deformations of the amplitudes on
a FIXED support S. They say nothing about the unrestricted fiber, and the gap
is not a technicality:

At a degenerate spectrum the unrestricted fixed-ρ fiber always contains the
whole natural-orbital gauge orbit of ψ, the U(m_1) × ... × U(m_k) orbit of
dimension Σ m_i² − (gauge already quotiented), because a block-diagonal
rotation leaves ρ = diag(λ) identically unchanged. Since every census vertex
has a degenerate spectrum, **the unrestricted fixed-ρ fiber is never a point at
a census vertex**, and no "rigid" verdict below should be read as saying it is.

The two are consistent because that orbit generically LEAVES any fixed support:
at the v103 support-14 state all 42 off-diagonal U(4) × U(6) generators move
the state off its support, which is exactly why a support-restricted
first-order computation sees none of it. The practical consequence is that
support-restricted first-order analysis is structurally blind to gauge
sparsification, and a separate finite-rotation search is required
(`src/gpc_census/gauge.py`, `docs/gauge_minimization.md`).

## Lemma

(i) [Fixed spectrum, weight sector, SUPPORT-RESTRICTED] The first-order
deformations of the weights preserving spec(ρ) = λ are exactly
K ∩ ker J_within^Re.
Cross-block silence imposes no first-order fixed-spectrum condition.

(ii) [Fixed ρ, weight sector, SUPPORT-RESTRICTED] The first-order
deformations preserving ρ itself are exactly K ∩ ker J_all^Re. (Unrestricted,
this fiber also contains the gauge orbit above; see Scope.)

(iii) [Dimension formula, full tangent] The first-order fixed-spectrum
tangent space of the support-restricted fiber, modulo the U(1)^d × U(1)
gauge, has dimension

    dim = (dim K − rank J_within^Re) + (|S| − g − rank J_within^Im),

where g is the rank of the gauge action on S (the number of independent
orbital-phase directions, global phase included).

## Proof

Let δc be a first-order deformation and δρ the induced 1-RDM variation;
δρ is R-linear in δc. Since ρ is diagonal with blocks B_i, its spectral
projectors P_i are coordinate projectors. By first-order degenerate
perturbation theory, the eigenvalues of ρ + ε δρ within the λ^(i)
cluster are λ^(i) + ε · eig(P_i δρ P_i) + O(ε²). Hence spec is
stationary and non-splitting at first order iff P_i δρ P_i = 0 for
every i. The off-block components P_i δρ P_j (i ≠ j) are unconstrained:
they rotate eigenvectors at first order and move eigenvalues only at
O(ε²). This proves the cross-block clause of (i).

The entries of P_i δρ P_i are: the diagonal entries δρ_oo for o ∈ B_i,
and the off-diagonal entries δρ_{pq} for p, q ∈ B_i, which are nonzero
only on channels of S, i.e. exactly the within-block channel variations
δσ_γ. At a real state the deformation splits orthogonally into a weight
part (δc real) and a phase part (δc imaginary), and δρ splits
accordingly into its real and imaginary parts (the diagonal receives
only the real part).

Weight sector: δρ_oo = (A_S δp)_o, so the diagonal conditions are
A_S δp = 0 with the norm condition Σ δp = 0, i.e. δp ∈ K; the
within-block off-diagonal conditions are J_within^Re δp = 0. This is
(i). For (ii), fixing ρ itself additionally requires every off-diagonal
of δρ to vanish, i.e. all channel sums stationary: J_all^Re δp = 0 on K.

Phase sector: δc = i θ ⊙ c contributes nothing to the diagonal (at a
real state, δρ_oo = 2 Re(i · (...)) = 0) and contributes i-parts to the
channel sums; the within-block conditions are J_within^Im θ = 0. The
gauge directions θ_t = Σ_{o ∈ t} φ_o form a rank-g subspace of the
phase sector and lie in the kernel of every within-block condition
(gauge rotates ρ_{pq} by e^{i(φ_p − φ_q)}, preserving zeros). Summing
the two sectors and subtracting gauge gives (iii). ∎

## Verified corollaries (census instances; scripts/reproduce_next_claims.py)

- v89 = (15,15,6^8)/26, |S| = 10, dim K = 1, one channel, WITHIN-block
  (both orbitals at 6/26), rank J_within^Re = 1, rank J_within^Im = 1,
  g = 9. Dimension formula: (1−1) + (10−9−1) = 0. v89 is first-order
  fixed-spectrum rigid, in agreement with the direct tangent
  computation (kernel = gauge = 9) and with the P4 scorecard.
- v103 = (18^4,5^6)/34, |S| = 14, dim K = 5, five channels of which two
  within-block ((1,2) in the 18/34 block, (5,8) in the 5/34 block) and
  three cross-block. rank J_within^Re = 2, rank J_within^Im = 2, g = 9.
  Dimension formula: (5−2) + (14−9−2) = 3 + 3 = 6. Directly measured
  first-order tangent: 6, constant along the fiber at t = 0.05, 0.10
  with eigenvalue multiplicities (4,6) preserved: the germ is a smooth
  6-manifold. v103 is fixed-spectrum DEFORMABLE; the P4 "first-order
  rigid" verdict is the fixed-ρ statement (ii): rank J_all^Re = 5 =
  dim K, which this lemma shows is a statement about a different fiber.

## Remarks

1. The "coincidence of dimensions" at v103 (five silence conditions
   exhausting a five-dimensional kernel) is a fixed-ρ fact. Whether
   rank J_all^Re = dim K is forced for cancellation states or accidental
   is open on a sample of size two.
2. The lemma is first-order. The v103 germ's smoothness is inferred from
   constant tangent rank along continued fiber points (numerical); a
   constant-rank argument in exact arithmetic would upgrade it.
3. The fixed-ρ / fixed-spectrum distinction here is the same one
   recorded in the levi.py bug log (fixed_rho_fiber_dim vs the then
   unqualified fiber_kernel_dim, since renamed incidence_kernel_dim).
   This lemma is the positive result that distinction
   was waiting for: the two fibers differ at first order exactly by the
   cross-block channel conditions. The distinction now lives in the code:
   gpc_census.fiber.fixed_rho_tangent and
   gpc_census.fiber.fixed_spectrum_tangent, with tests/test_fiber_notions.py
   pinning the v103 disagreement and the v89 agreement.
4. Consequence for T1 (Jacobian theorem): in the cancellation regime the
   correct fiber-dimension formula is (iii), not
   dim K − #(touched 1-term classes). The magnitude-target regime's
   formula and (iii) should be unified by the observation that in the
   magnitude regime every active channel behaves as a within-block
   constraint on its modulus with a free phase.
5. Scope correction (2026-07, out of sample): the gap between the two
   fibers is NOT confined to the cancellation regime this lemma is about.
   At (3,11) v43, a magnitude-regime state with one channel and 1-RDM
   off-diagonal 0.2, the fixed-ρ tangent is 1 while the fixed-spectrum
   tangent is 2. The lemma computes the gap in the cancellation regime; it
   does not claim the gap occurs only there.
   [docs/prereg_bracket_3_11_3_12.md, prediction B3, scored FAIL]
6. The cancellation regime is a property of the REPRESENTATIVE, not of the
   vertex. The census-wide cascade walks v103 from its certified
   cancellation-regime state to a support-10 point whose 1-RDM has
   off-diagonal 7.4e-2, i.e. out of the regime entirely. So "the census has
   exactly two cancellation states" counts certified representatives.
   [docs/cascade_census.md]
