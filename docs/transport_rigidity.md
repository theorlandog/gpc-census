# Transport rigidity on glued cyclic carriers

The infinite observability construction admits cross-block density interactions,
which create entanglement by conditional phases.  A stronger physical target is
a non-diagonal two-body interaction that transfers excitation between blocks
while preserving the complete glued carrier.

For two blocks, the most general block-number-preserving cross interaction is

\[
H_{\rm hop}=\sum_{p\ne q,r\ne s}K_{pq,rs}
(a_p^\dagger a_q)(a_r^\dagger a_s)+\mathrm{h.c.}
\]

Carrier preservation is a homogeneous integer linear system: every matrix
element from a carrier determinant to a determinant outside the carrier must
vanish.  The script `scripts/transport_rigidity.py` builds this leakage matrix
with exact fermionic signs and computes its rational rank.

## Exact finite results

For \(S(5,2)\boxtimes S(5,2)\), the leakage matrix has 400 columns and exact
rank 400.  Thus the only cross-block pair-hopping kernel preserving the whole
25-dimensional carrier is zero.  The density interaction \(n_pn_q\) still
preserves and entangles the carrier, but transport is rigidly excluded within
this operator class.

For \(S(4,2)\boxtimes S(4,2)\), the preserving hopping space has dimension 16.
For example, the Hermitian combination

\[
(a_2^\dagger a_0)(a_6^\dagger a_4)+
(a_0^\dagger a_2)(a_4^\dagger a_6)
\]

moves four pairs of product windows without leakage.  However this carrier is
not certified by the present 2-RDM theorem: its pair-incidence rank is 11 rather
than 16, and its simple collision-free coherence graph is empty.

This exposes a **transport-observability tension**:

* the coprime odd cyclic carriers give clean exact 2-RDM observability but can
  be rigid against carrier-wide pair hopping;
* the smallest hopping-flexible cyclic product loses the weight certificate.

The finite rank scan through the listed cases supports a broader rigidity
conjecture for nontrivial cyclic windows with \(d\ge5\), but no general theorem
is claimed yet.

## Direct 2-RDM propagation for diagonal interactions

The density-entangling dynamics can be written directly on the certified
2-RDM chart.  Pair occupations are constant.  Every selected collision-free
tree coherence evolves as

\[
X_{ab}(t)=e^{-i(E_b-E_a)t}X_{ab}(0).
\]

Thus the measured diagonal pair coordinates and \(m-1\) measured tree channels
form a closed amplitude-free propagator.  Wavefunction coefficients need not
be reconstructed during the integration.  Higher coherences and
\(\Gamma^{(3)}\) are evaluated only when requested through the rational path
map.

## Next attack

The obstruction redirects the search toward three possibilities:

1. state-specific scheduled hopping that preserves a trajectory rather than
   the entire carrier;
2. enlarged carriers formed by adjoining leakage images until hopping closes,
   followed by a fresh 2-RDM observability test;
3. non-cyclic factor families whose singleton and pair incidence remain full
   while their one-body normalizers are nontrivial.
