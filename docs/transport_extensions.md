# Transport extensions for exactly observable carriers

This campaign tests three ways around the carrier-wide hopping obstruction for
`S(5,2) x S(5,2)`.  Every rank and invariance statement below is checked by
`scripts/transport_extensions.py` with exact integer or rational linear
algebra.

## 1. State-specific scheduled tangent transport

Carrier-wide preservation requires `(1-P)HP=0` and excludes every nonzero
cross-block hop on this 25-state carrier.  A scheduled trajectory only requires

\[
(1-P)H(t)|\psi(t)\rangle=0.
\]

At the uniform carrier state, the search runs over the 200 Hermitian
combinations `O + O^dagger` of the 400 directed cross hops.  That is half of
the real Hermitian span, which also contains the 200 combinations
`i(O - O^dagger)`; since the conclusion is an existence statement, searching
the smaller family only strengthens it.  The state-leakage map on that family
has rank 75, leaving a 125-dimensional exact kernel.  Its induced tangent action has rank 25.  The certificate records an
8-term Hermitian kernel whose leakage is identically zero on that state and
whose carrier action is nonprojective: only carrier coordinates 22 and 24 move,
with opposite signs.

Thus transport rigidity of the whole carrier does **not** imply tangent
rigidity of every state.  This is an exact first-order local result.  It is not
yet a finite-time scheduled path; continuing the kernel smoothly and avoiding
rank loss along the path remain open.

## 2. Hopping-closure enlargement

Choose the Hermitian hop

\[
(a_1^\dagger a_0)(a_7^\dagger a_5)+
(a_0^\dagger a_1)(a_5^\dagger a_7).
\]

Starting from `S(5,2) x S(5,2)`, repeatedly adjoining every leakage image closes
after four new determinants.  The resulting carrier has dimension 29 and is
invariant under the selected hop.  Its pair-incidence matrix has exact rank 29,
and its collision-free 2-RDM coherence graph is connected, with 75 certified
channels.

Therefore the enlarged carrier supports genuine non-diagonal transport while
retaining exact pure-state reconstruction and exact BBGKY closure from the
fixed-frame 2-RDM.  This is presently the strongest physical construction in
the campaign: transport, interaction, invariance, and observability coexist in
one finite exact model.

## 3. Different families with a one-body normalizer

Let `C` be a fixed core and let the active sector contain either one or two
fermions in `a` active orbitals:

\[
S_{C,k,a}=\{C\cup A: A\subseteq \mathcal A,\ |A|=k\},\qquad k=1,2.
\]

For `k=1`, core-active pair occupations recover every weight.  For `k=2`, the
active-pair occupations are an identity submatrix of pair incidence.  In both
cases the 2-RDM directly supplies all active-sector coherences, so the
collision-free graph is connected.  The full active one-body algebra preserves
the carrier; in particular there are `a(a-1)` nontrivial off-diagonal active
hops.

These fixed-core complete-active-space families are controls rather than a
novel active-space construction.  Their importance here is conceptual: exact
2-RDM observability and transport are compatible.  The cyclic product's
rigidity is a property of that sparse support, not a universal consequence of
observability.

## Physics consequence

The combined results separate three levels of dynamics:

1. **carrier-wide rigidity** on the original sparse cyclic product;
2. **state-specific tangent transport** without carrier-wide invariance;
3. **exact invariant transport** after a minimal four-determinant closure
   enlargement.

The next theorem target is a continuation criterion for the scheduled tangent
kernel.  The next computational target is to propagate a state on the
29-dimensional enlarged carrier directly in certified 2-RDM coordinates and
verify the reconstructed 3-RDM throughout the transport trajectory.
