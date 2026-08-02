# Cyclic-window carriers and exact reduced dynamics

## Scope

This note isolates a constructive infinite family behind the finite 2-RDM
observability atlas.  The result is a carrier-restricted theorem.  It does not
solve generic N-representability, infer an unknown carrier, or claim that every
Hamiltonian preserving a moment-map fiber is two-body realizable.

The new point is compositional: small observable carriers can be glued into
arbitrarily large observable carriers while retaining an exact pure-state
2-RDM coordinate chart and an explicit class of carrier-preserving two-body
Hamiltonians.

## Cyclic-window factors

Fix integers `2 <= N < d`.  On orbitals `Z/dZ`, define

\[
W_k=\{k,k+1,\ldots,k+N-1\}\pmod d,
\qquad k\in\mathbb Z/d\mathbb Z,
\]

and let `S(d,N)={W_k}` be the corresponding determinant carrier.
Let `A(d,N)` be its singleton-incidence matrix,

\[
A_{j,k}=\mathbf 1[j\in W_k].
\]

### Proposition 1: exact singleton rank

\[
\operatorname{rank}A(d,N)=d-\gcd(d,N)+1.
\]

In particular, `A(d,N)` is nonsingular exactly when `gcd(d,N)=1`.

**Proof.**  The matrix is circulant.  At the Fourier character `zeta^r`, its
symbol is

\[
p(\zeta^r)=1+\zeta^r+\cdots+\zeta^{(N-1)r}.
\]

For `r != 0`, this vanishes exactly when `d` divides `Nr` but not `r`.
There are `gcd(d,N)-1` such nonzero characters.  The zero character has
symbol `N`, hence is nonzero.  The stated rank follows.  \(\square\)

### Corollary 2: pair-incidence rank in the coprime family

Let `B(d,N)` be pair incidence.  Summing pair rows containing orbital `j`
gives

\[
\sum_{q\ne j}B_{\{j,q\},k}=(N-1)A_{j,k}.
\]

Therefore `row(A)` is contained in `row(B)`.  If `gcd(d,N)=1`, then
`B(d,N)` has full column rank `d`.

This is only a sufficient argument.  Exact computations show full pair rank
in additional noncoprime examples; no stronger general formula is claimed
here.

### Proposition 3: a collision-free cycle

For each `k`, the adjacent windows `W_k` and `W_{k+1}` differ by dropping
`k` and adding `k+N`.  Choose any common spectator `s` in
`W_k intersection W_{k+1}` and the 2-RDM channel

\[
P=\{k,s\},\qquad Q=\{k+N,s\}.
\]

The ordered dropped/added labels determine `k`; hence no other adjacent
window pair contributes the same carrier coherence to this channel.  After
fixing one spectator rule, these channels give a collision-free spanning
cycle `C_d` in the carrier coherence graph.

For `N=2`, the common spectator is unique.  For `N>2`, a deterministic rule
such as the least common orbital is sufficient.  Fermionic signs affect the
reconstruction coefficient, not collision freedom.

### Theorem 4: cyclic-window observability

If `gcd(d,N)=1`, every pure full-support state on `S(d,N)` is rationally
reconstructed, up to global phase, from its fixed-frame 2-RDM.

**Proof.**  Corollary 2 recovers all coefficient weights from diagonal pair
occupations.  Proposition 3 supplies a connected collision-free coherence
graph.  The support theorem in `docs/rdm_observability_census.md` applies.
\(\square\)

## Disjoint-block gluing

Let `S_r` be determinant carriers on pairwise disjoint orbital blocks `O_r`,
with fixed particle number `N_r` in each block.  Their product carrier is

\[
S=\boxtimes_{r=1}^q S_r
 =\{D_1\cup\cdots\cup D_q:D_r\in S_r\}.
\]

Write `A_r` for singleton incidence in factor `r`.

### Proposition 5: cross-pair incidence factorization

For two factors, restrict product pair incidence to cross-block pairs
`{i,j}` with `i in O_1`, `j in O_2`.  With columns indexed by `(D_a,E_b)`,
this block is exactly

\[
B_{12}^{cross}=A_1\otimes A_2.
\]

Consequently, if both singleton-incidence matrices have full column rank,
then the complete product pair-incidence matrix has full column rank.
For `q` factors, any bipartition of the blocks gives the analogous Kronecker
factorization, and iterated gluing yields full weight recovery when all
factor singleton-incidence matrices are full column rank.

### Labeled coherence graphs

A factor coherence edge is called *label-injective* when its selected
one-hop channel has a dropped/added orbital label not shared by any other
selected factor edge.  The cyclic-window cycle above is label-injective.

For two label-injective factor edges `a--a'` and `b--b'`, the product
coherence between `(a,b)` and `(a',b')` is exposed by the cross-block 2-RDM
channel made from the two dropped orbitals and the two added orbitals.  Its
labels identify both factor edges, so the channel is collision-free.

Thus the product collision-free graph contains a spanning subgraph
isomorphic to the categorical graph product

\[
G_1\times G_2.
\]

The word *contains* is essential: additional collision-free channels may be
present.

### Theorem 6: gluing theorem

Suppose each factor `S_r` has

1. full-column-rank singleton incidence;
2. a connected, label-injective collision-free one-hop graph `G_r`;
3. at least one non-bipartite `G_r`.

Then the disjoint-block product carrier has full-column-rank pair incidence
and a connected collision-free 2-RDM graph.  Hence every pure full-support
state on the product carrier is rationally reconstructible from its
fixed-frame 2-RDM.

**Proof.**  Proposition 5 gives full weight recovery.  The selected product
channels give a spanning categorical product of the `G_r`.  A categorical
product of connected graphs is connected when at least one factor is
non-bipartite.  The full collision-free graph contains this connected
spanning graph and is therefore connected.  Apply the support theorem.
\(\square\)

### Corollary 7: glued odd cyclic windows

Let every factor be `S(d_r,N_r)` with `gcd(d_r,N_r)=1`, and suppose at least
one `d_r` is odd.  Then the product carrier is 2-RDM observable on its pure
full-support chart.

Its carrier dimension is

\[
M=\prod_r d_r.
\]

The containing block-number sector has dimension

\[
\prod_r {d_r\choose N_r},
\]

while the full fixed-total-particle FCI sector has dimension

\[
{\sum_r d_r\choose\sum_r N_r}.
\]

These two ambient dimensions must not be conflated.

## Exact carrier dynamics

Let `H_r(t)` be an at-most-two-body Hamiltonian preserving factor carrier
`S_r`.  On disjoint orbital blocks,

\[
H(t)=\sum_r H_r(t)
\]

is again at most two-body and preserves the product carrier.  On every pure
full-support trajectory in the carrier, Theorem 6 gives an exact rational
map

\[
\Gamma^{(3)}=F_{S,3}(\Gamma^{(2)}),
\]

so the TD-2RDM equation closes exactly after substitution.

This establishes an infinite family of exact carrier-conditional reduced
dynamical systems.  It does **not** yet establish a computational speedup:
that requires an amplitude-free implementation, conditioning bounds, and a
benchmark where product-carrier propagation is practical while propagation
in the relevant ambient sector is not.

## Physical content

The theorem identifies a class of interacting fermionic models for which the
2-RDM is not merely a reduced observable but an exact coordinate system on a
pure-state dynamical manifold.  Within that manifold:

- all higher equal-time correlators are functions of the 2-RDM;
- BBGKY truncation at order two is exact rather than approximate;
- disjoint interacting subsystems can be composed without losing exact
  observability;
- carrier-preserving controls can be designed factorwise and lifted to large
  systems.

The genuinely physical next step is to add nonseparable inter-factor
Hamiltonian terms that still preserve the product carrier and produce
correlation transport between blocks.  The sum-of-factors Hamiltonian proves
closure and scalability of the state space, but by itself describes
independent block dynamics.

## Falsifiable next gates

1. Implement the rational reconstruction DAG and propagate `Gamma^(2)` without
   reading amplitudes.
2. Prove conditioning bounds for cyclic windows and their products.
3. Search for carrier-preserving cross-block two-body couplings and determine
   whether they generate entanglement between factors.
4. Benchmark against full propagation in both the block-number sector and the
   full fixed-particle sector.
5. Apply linear-plus-rank-one determinacy closure to the 28 census supports
   not certified by a collision-free tree.

Only gates 1--4 would justify a claim of a scalable new reduced-dynamics
method.  The present theorem supplies the exact mathematical substrate.
