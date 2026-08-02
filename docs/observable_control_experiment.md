# Observable locality versus controllable locality

## Question

For a determinant carrier \(S=\{D_1,\ldots,D_m\}\), compare

\[
q_{\rm ctrl}(S),
\]

the minimum interaction order whose independently addressable determinant
transitions connect the carrier, with

\[
q_{\rm obs}(S),
\]

the minimum RDM order whose collision-free coherence channels connect it.

The motivating conjecture was equality. The experiment disproves universal
equality and replaces it with a sharper hierarchy.

## Exact graph theorem for control

Define the order-\(q\) excitation graph \(E_q(S)\) by connecting determinants
when

\[
\delta(D_a,D_b)=N-|D_a\cap D_b|\le q.
\]

A number-conserving \(q\)-body operator can have a matrix element between the
two determinants only under this condition. If independent Hermitian edge
controls and independent diagonal controls are available, every connected
component of \(E_q(S)\) is controllable. If the graph is connected, the edge
matrices and diagonal commutators generate \(\mathfrak{su}(m)\).

Thus, under this stated control model,

\[
q_{\rm ctrl}(S)=\min\{q:E_q(S)\text{ is connected}\}.
\]

This is a graph-theoretic control theorem, not a numerical conjecture.

## Observable graph

A \(q\)-RDM can contain the coherence \(X_{ab}=\bar c_a c_b\) only if
\(\delta(D_a,D_b)\le q\). Multiple determinant pairs may contribute to the
same RDM channel. The collision-free observable graph retains only edges that
appear alone in at least one channel.

Therefore its edge set is a subgraph of the control graph, giving the general
inequality

\[
\boxed{q_{\rm obs}(S)\ge q_{\rm ctrl}(S)}.
\]

Equality holds when a connected control subgraph has independently readable
channels. A strict gap measures **readout frustration**: the system can be
coherently controlled at locality \(q\), but its phases cannot be separately
identified by the order-\(q\) RDM because channels collide.

## Experiment

The script sampled or exhaustively enumerated 1,790 carriers in
\(\wedge^3\mathbb C^6\), with carrier sizes two through six. For every
carrier it computed exact control and collision-free observable graphs at
orders one through three. It also propagated a random Hamiltonian at the
minimum control order and recorded configuration entropy and participation.

Results:

| Outcome | Carriers | Fraction |
|---|---:|---:|
| \(q_{\rm obs}=q_{\rm ctrl}\) | 1,620 | 90.50% |
| \(q_{\rm obs}=q_{\rm ctrl}+1\) | 170 | 9.50% |
| \(q_{\rm obs}<q_{\rm ctrl}\) | 0 | 0% |

Pair counts:

| \((q_{\rm ctrl},q_{\rm obs})\) | Count |
|---|---:|
| (1,1) | 893 |
| (1,2) | 163 |
| (2,2) | 717 |
| (2,3) | 7 |
| (3,3) | 10 |

The inequality was never violated, which is a consistency check rather than a
discovery: the containment argument above already proves it. The subgraph
containment itself was verified directly on every carrier of size two through
six in \(\wedge^3\mathbb C^6\), 60,439 carriers, at every order, with no
failure. All strict gaps were one order and arose from channel collisions,
and that gap distribution is the actual finding.

The experiment drops fermionic signs when deciding whether a channel is
collision free, which needs justifying rather than assuming. It is sound here:
a given carrier pair contributes exactly one term of \(\pm1\) to a given
channel, so nothing can cancel. Checked against the exactly signed
implementation in `scripts/census_support_geometry.py` on all 21,679 carriers
of size two through five, the two agree on the collision-free edge sets
themselves at orders two and three, not merely on the resulting order.

## Census implication

The merged census has 28 order-two exceptions and reports no collided
order-two channels on those records. For that collision-free class, the only
reason the order-two graph can be disconnected is missing two-body overlap
edges. Hence their observed reconstruction order is predicted to equal their
control locality under independently addressable controls.

This should be checked record-by-record by the repository campaign script,
but it is already a much narrower and more physical statement than the
original universal conjecture.

## New physical quantities

The experiment identifies three locality scales:

1. **Control locality** \(q_{\rm ctrl}\): minimum interaction order required to
   connect the carrier dynamically.
2. **Observable locality** \(q_{\rm obs}\): minimum RDM order required by the
   collision-free reconstruction certificate.
3. **Readout-frustration gap**
   \[
   \Delta_q=q_{\rm obs}-q_{\rm ctrl}\ge0.
   \]

A fourth quantity, the non-toric phase complexity \(\kappa\), is independent:
it counts intrinsic interference directions, not the locality needed to reach
or read them.

## Practical consequences

- A disconnected \(E_2\) proves that no two-body control family can coherently
  connect the full carrier without leaving it.
- A connected \(E_2\) but disconnected observable graph identifies a
  measurement-design problem rather than a control problem.
- Extra observables, orbital-frame changes, or designed probes should target
  the collided channels in the readout-frustrated class.
- Active-space selection can optimize both control locality and readout gap,
  rather than relying only on energy criteria.
- TD-RDM closure order and Hamiltonian locality can be compared on the same
  support filtration.

## Honest limitation

Random Hamiltonian entropy and participation do not collapse onto
\(q_{\rm ctrl}\) alone; they also depend strongly on carrier size, graph
structure, coupling distribution, and time. The experiment does not establish
that one scalar locality controls all entanglement or transport dynamics.

The robust result is the graph hierarchy and the identification of channel
collisions as the sole source of \(q_{\rm obs}>q_{\rm ctrl}\) in the tested
sample.
